#include "FEM4.h"

// Class constructor for a scalar field
template <int dim>
FEM<dim>::FEM (double Alpha, double delta_tin)
:
fe (FE_Q<dim>(order), 1),
  dof_handler (triangulation),
  quadrature_formula(quadRule)
{
  alpha = Alpha;
  delta_t = delta_tin;

  nodal_solution_names.push_back("D");
  nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
}

//Class destructor
template <int dim>
FEM<dim>::~FEM (){dof_handler.clear ();}


/**************************************************************************************/
/*  Define the problem domain and generate the mesh                                   */
/**************************************************************************************/

template <int dim>
void FEM<dim>::generate_mesh(std::vector<unsigned int> numberOfElements){

  //Define the limits of your domain
    x_min = 0.;
    x_max = 1.;
    y_min = 0.;
    y_max = 1.;
    z_min = 0.;
    z_max = 0.1;

  Point<dim,double> min(x_min,y_min,z_min),
    max(x_max,y_max,z_max);
  GridGenerator::subdivided_hyper_rectangle (triangulation, numberOfElements, min, max);
}


/**************************************************************************************/
/*  This function specifies and stores Dirichilet boundary conditions.                */
/*  After assembling Kglobal and Flocal, theses b.c.s will be automatically applied   */
/*  to resive matrices of the system with techniques supported by deal.II.            */
/*  Note that we are not only define the Drichilet boundary values of a steady state, */
/*  we also defines the boundary changing variable V here, which should be 0 in this  */
/*  assignment since we do not have time-variant b,c,s                                */
/**************************************************************************************/

template <int dim>
void FEM<dim>::define_boundary_conds(){

  //Define the Dirichlet boundary conditions.
  const unsigned int totalNodes = dof_handler.n_dofs(); //Total number of nodes

  for(unsigned int globalNodeIndex = 0; globalNodeIndex < totalNodes; globalNodeIndex++) {
    if(nodeLocation[globalNodeIndex][0] == x_min) {
      boundary_values_of_D[globalNodeIndex] = 300.;
      boundary_values_of_V[globalNodeIndex] = 0.;
    } 
    if(nodeLocation[globalNodeIndex][0] == x_max) {
      boundary_values_of_D[globalNodeIndex] = 310.;
      boundary_values_of_V[globalNodeIndex] = 0.;
    } 
  }

}


/**************************************************************************************/
/*  This function sets up data structures for Finite Element Solving System.          */
/*  We initialize matrix shapes here and assign the degree of freedom to elements.    */
/**************************************************************************************/

template <int dim>
void FEM<dim>::setup_system(){

  //Let deal.II organize degrees of freedom
  dof_handler.distribute_dofs (fe);

  //Fill in the Table "nodeLocations" with the x, y, and z coordinates of each node by its global index
  MappingQ1<dim,dim> mapping;
  std::vector< Point<dim,double> > dof_coords(dof_handler.n_dofs());
  nodeLocation.reinit(dof_handler.n_dofs(),dim);
  DoFTools::map_dofs_to_support_points<dim,dim>(mapping,dof_handler,dof_coords);
  for(unsigned int i=0; i<dof_coords.size(); i++){
    for(unsigned int j=0; j<dim; j++){
      nodeLocation[i][j] = dof_coords[i][j];
    }
  }

  //Specify boundary condtions (call the function)
  define_boundary_conds();

  //Define the size of the global matrices and vectors
  sparsity_pattern.reinit (dof_handler.n_dofs(),
                           dof_handler.n_dofs(),
                           dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress();
  K.reinit (sparsity_pattern);
  M.reinit (sparsity_pattern);
  system_matrix.reinit (sparsity_pattern);
  D_steady.reinit(dof_handler.n_dofs());
  D_trans.reinit(dof_handler.n_dofs());
  V_trans.reinit(dof_handler.n_dofs());
  RHS.reinit(dof_handler.n_dofs());
  F.reinit(dof_handler.n_dofs());

  //Just some notes...
  std::cout << "   Number of active elems:       " << triangulation.n_active_cells() << std::endl;
  std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl;   
}


/**************************************************************************************/
/*  This function assembles K_local, F_local, K_glocal and F_global matrices/vectors. */
/*  It generally applies a quadrature approximation to estimate an itegral over the   */
/*  reference domain, rather than solve it analytically.                              */
/*  Note that we ignore Neumann b.c. here because we assume j = 0 at boundaries.      */
/**************************************************************************************/

template <int dim>
void FEM<dim>::assemble_system(){

  M=0; K=0; F=0;

  FEValues<dim> fe_values(fe,
                          quadrature_formula,
                          update_values | 
                          update_gradients | 
                          update_JxW_values);

  const unsigned int dofs_per_elem = fe.dofs_per_cell;         //This gives you dofs per element
  unsigned int       num_quad_pts = quadrature_formula.size(); //Total number of quad points in the element
  FullMatrix<double> Mlocal (dofs_per_elem, dofs_per_elem);
  FullMatrix<double> Mlocal_lumped (dofs_per_elem, dofs_per_elem);
  FullMatrix<double> Klocal (dofs_per_elem, dofs_per_elem);
  Vector<double>     Flocal (dofs_per_elem);

  std::vector<unsigned int> local_dof_indices (dofs_per_elem); //This relates local dof numbering to global dof numbering
  double                    rho = 3.8151e6;  // specify the specific heat per unit volume

  //loop over elements  
  typename DoFHandler<dim>::active_cell_iterator elem = dof_handler.begin_active (),
    endc = dof_handler.end();
  for (;elem!=endc; ++elem){

    elem->get_dof_indices (local_dof_indices);
    fe_values.reinit(elem); //Retrieve values from current element

    elem->get_dof_indices (local_dof_indices);

    //Loop over local DOFs and quadrature points to populate Mlocal
    Mlocal = 0.;
    for(unsigned int q=0; q<num_quad_pts; q++){
      for(unsigned int A=0; A<fe.dofs_per_cell; A++){
        for(unsigned int B=0; B<fe.dofs_per_cell; B++){
          Mlocal[A][B] += rho * fe_values.shape_value(A,q)
                * fe_values.shape_value(B,q) * fe_values.JxW(q);

        }
      }
    }

    //assemble lumped mass matrix
    Mlocal_lumped = 0.;
    double Mlrow;
    for(unsigned int A=0; A<fe.dofs_per_cell; A++) {
      Mlrow = 0.;
      for(unsigned int B=0; B<fe.dofs_per_cell; B++) {
        Mlocal_lumped[A][A] += Mlocal[A][B];
      }
    }

    FullMatrix<double> kappa(dim,dim);
    kappa[0][0] = 385.;
    kappa[1][1] = 385.;
    kappa[2][2] = 385.;

    //Loop over local DOFs and quadrature points to populate Klocal
    Klocal = 0.;
    for(unsigned int A=0; A<fe.dofs_per_cell; A++){
      for(unsigned int B=0; B<fe.dofs_per_cell; B++){
        for(unsigned int q=0; q<num_quad_pts; q++){
          for(unsigned int i=0; i<dim; i++){
            for(unsigned int j=0; j<dim; j++){
              Klocal[A][B] += fe_values.shape_grad(A,q)[i] * kappa[i][j]
                          * fe_values.shape_grad(B,q)[j] * fe_values.JxW(q);
            }
          }
        }
      }
    }

    for (unsigned int i=0; i<dofs_per_elem; ++i){
      for (unsigned int j=0; j<dofs_per_elem; ++j){
        K.add(local_dof_indices[i], local_dof_indices[j], Klocal[i][j]);
        if(iflumped) // modify it in FEM4.h
          M.add(local_dof_indices[i], local_dof_indices[j], Mlocal_lumped[i][j]);
        else
          M.add(local_dof_indices[i], local_dof_indices[j], Mlocal[i][j]);
      }
    }
  }

}


/**************************************************************************************/
/*  This function directly solves "K_glocal * D = F_global" at steady state by        */
/*  implementing a matrix inverse - which is cost expensive, but of relatively        */
/*  high accuracy.                                                                    */
/**************************************************************************************/

template <int dim>
void FEM<dim>::solve_steady(){

  MatrixTools::apply_boundary_values (boundary_values_of_D, K, D_steady, F, false);

  SparseDirectUMFPACK  A;
  A.initialize(K);
  A.vmult (D_steady, F); //D=K^{-1}*F

  output_steady_results();
}


/**************************************************************************************/
/*  This function apply initial conditions for the transient problem and also solves  */
/*  transient solutions at t = 0.                                                     */
/**************************************************************************************/

template <int dim>
void FEM<dim>::apply_initial_conditions(){

  const unsigned int totalNodes = dof_handler.n_dofs(); 

  for(unsigned int i=0; i<totalNodes; i++){
    if(nodeLocation[i][0] < 0.5){
      D_trans[i] = 300.;
    }
    else{
      D_trans[i] = 300. + 20. * (nodeLocation[i][0] - 0.5); 
    }
  }
  //Find V_0 = M^{-1}*(F_0 - K*D_0)
  system_matrix.copy_from(M);

  // Define the right-hand-side vector (RHS = F_0 - K*D_0)
  // where D_0 = D_trans
  K.vmult(RHS,D_trans); 
  RHS *= -1.;
  RHS.add(1.,F); 

  MatrixTools::apply_boundary_values (boundary_values_of_V, system_matrix, V_trans, RHS, false);

  SparseDirectUMFPACK  A;
  A.initialize(system_matrix);
  A.vmult (V_trans, RHS);

  output_trans_results(0);

  double current_l2norm = l2norm();
  l2norm_results.push_back(current_l2norm);

}


/**************************************************************************************/
/*  This function solves transient problem iteratively from t = 1 to t = 3000 s       */
/*  using a v-method, in Forward Euler, BackwardEuler, and Mid-Point schemes.         */
/*  The solutions are auto-saved per 100 s.                                           */
/**************************************************************************************/

template <int dim>
void FEM<dim>::solve_trans(){

  //Call the function to initialize D_trans and V_trans as D_0 and V_0
  apply_initial_conditions();

  //Define delta_t
  //const double delta_t = 1.;

  const unsigned int totalNodes = dof_handler.n_dofs(); //Total number of nodes
  Vector<double>     D_tilde(totalNodes);

  //Loop over time steps and update D_transient from D_n to D_{n+1} using the V method
  for(unsigned int t_step=1; t_step<3001; t_step++){

    //Update D_tilde. Dtile = delta_t * ( 1 - alpha) * V_trans + D_trans
    //where D_trans is D_n, V_trans = V_n at nth step, and D_tilde is at (n+1)th step
    D_tilde = D_trans;
    D_tilde.add(delta_t * (1.-alpha), V_trans);

    //Update V_trans. V_trans = (M + alpha * delta_t * K)^{-1} * (F - K * D_tilde)
    //Note that after updation, V_trans = V_n+1
    system_matrix.copy_from(M);
    system_matrix.add(alpha * delta_t, K);
    K.vmult(RHS, D_tilde);
    RHS *= -1.;
    RHS.add(1., F);

    //Apply boundary conditions on V_trans before solving the matrix/vector system
    MatrixTools::apply_boundary_values (boundary_values_of_V, system_matrix, V_trans, RHS, false);

    //Solve for V_trans
    SparseDirectUMFPACK  A;
    A.initialize(system_matrix);
    A.vmult (V_trans, RHS);

    //Update D_trans. D_trans = D_tilde + delta_t * alpha * V_trans
    //Ater updation, D_trans = D_n+1
    D_trans = D_tilde;
    D_trans.add(delta_t * alpha, V_trans);

    //Output the results every 100 seconds
    if(t_step%100 == 0){
      output_trans_results(t_step);

      double current_l2norm = l2norm();
      l2norm_results.push_back(current_l2norm);
    }
  }

}


/**************************************************************************************/
/*  This function writes the steady state nodal solutions to "solution_ss.vtk" for    */
/*  assignment submission and results checking.                                       */
/**************************************************************************************/

template <int dim>
void FEM<dim>::output_steady_results (){
  //Write results to VTK file
  std::ofstream output1 ("../sols/vtk/solution_ss.vtk");
  DataOut<dim> data_out; data_out.attach_dof_handler (dof_handler);

  //Add nodal DOF data
  data_out.add_data_vector (D_steady, nodal_solution_names, DataOut<dim>::type_dof_data, nodal_data_component_interpretation);
  data_out.build_patches (); data_out.write_vtk (output1); output1.close();
}


/**************************************************************************************/
/*  This function writes the transient nodal solutions to "solution_time.vtk" for     */
/*  assignment submission and results checking.                                       */
/**************************************************************************************/

template <int dim>
void FEM<dim>::output_trans_results (unsigned int index){
  //This adds an index to your filename so that you can distinguish between time steps

  //Write results to VTK file
  char filename[100];
  snprintf(filename, 100, "../sols/vtk/solution_%d.vtk", index);
  std::ofstream output1 (filename);
  DataOut<dim> data_out; data_out.attach_dof_handler (dof_handler);

  //Add nodal DOF data
  data_out.add_data_vector (D_trans, nodal_solution_names, DataOut<dim>::type_dof_data, nodal_data_component_interpretation);
  data_out.build_patches (); data_out.write_vtk (output1); output1.close();
}


/**************************************************************************************/
/*  This function computes the L2 difference of steady state and transient values at  */
/*  a given time. Note that this L2 error is with respect to a numerical steady state */
/*  solution, but not a exact soluiton. It is used to evaluate the stability of a     */
/*  scheme. We used a quadrature Rule here to interpolate the integral value.         */
/**************************************************************************************/

template <int dim>
double FEM<dim>::l2norm(){
  double l2norm = 0.;

  FEValues<dim> fe_values(fe,
                          quadrature_formula,
                          update_values |
                          update_JxW_values);

  const unsigned int                            dofs_per_elem = fe.dofs_per_cell; //dofs per element
  std::vector<unsigned int> local_dof_indices (dofs_per_elem);
  const unsigned int                            num_quad_pts = quadrature_formula.size();
  double                              u_steady, u_trans;

  //loop over elements  
  typename DoFHandler<dim>::active_cell_iterator elem = dof_handler.begin_active (), 
    endc = dof_handler.end();
  for (;elem!=endc; ++elem){
    elem->get_dof_indices (local_dof_indices);
    fe_values.reinit(elem);

    for(unsigned int q=0; q<num_quad_pts; q++){
      u_steady = 0.; u_trans = 0.;
      for(unsigned int A=0; A<dofs_per_elem; A++){
        u_steady += D_steady[local_dof_indices[A]]*fe_values.shape_value(A, q);
        u_trans += D_trans[local_dof_indices[A]]*fe_values.shape_value(A, q);

      }
      // define the l2norm of the difference between u_steady and u_trans
      l2norm += fe_values.JxW(q) * (u_steady - u_trans) * (u_steady - u_trans);

    }

  }

  return sqrt(l2norm);
}



template class FEM<1>;
template class FEM<2>;
template class FEM<3>;