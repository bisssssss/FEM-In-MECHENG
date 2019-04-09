#include "FEM5.h"

// Class constructor for a vector field
template <int dim>
FEM<dim>::FEM ()
:
fe (FE_Q<dim>(order), dim),
  dof_handler (triangulation),
  quadrature_formula(quadRule),
  face_quadrature_formula(quadRule)
{       
                
  //Nodal Solution names - this is for writing the output file
  for (unsigned int i=0; i<dim; ++i){
    nodal_solution_names.push_back("u");
    nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  }
}

//Class destructor
template <int dim>
FEM<dim>::~FEM (){dof_handler.clear ();}


/**************************************************************************************/
/*  function C returns the components of the 4th order elasticity tensor.             */
/*        C(i,j,k,l) = lambda*1(i,j)*1(k,l) + mu*[1(i,k)*1(j,l)+1(i,l)*1(j,k)]        */
/*    where 1(x, y) is the Kronecker delta                                            */
/**************************************************************************************/

template <int dim>
double FEM<dim>::C(unsigned int i,unsigned int j,unsigned int k,unsigned int l){

  //Define the material parameters of Young's modulus and Poisson's ratio
  double E= 2.0e11,
    nu= 0.3;
  double lambda=(E*nu)/((1.+nu)*(1.-2.*nu)),
    mu=E/(2.*(1.+nu));

  return lambda*(i==j)*(k==l) + mu*((i==k)*(j==l) + (i==l)*(j==k));

}


/**************************************************************************************/
/*  Define the problem domain and generate the mesh                                   */
/**************************************************************************************/

template <int dim>
void FEM<dim>::generate_mesh(std::vector<unsigned int> numberOfElements){

  //Define the limits of your domain
    x_min = 0.,
    x_max = 1.,
    y_min = 0.,
    y_max = 1.,
    z_min = 0.,
    z_max = 1.;

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

  const unsigned int totalDOFs = dof_handler.n_dofs(); //Total number of degrees of freedom
    
    //loop over global degrees of freedom 
    for(unsigned int global_dof=0; global_dof<totalDOFs; global_dof++){

        // Since in both problems, u1 = u2 = u3 = 0 on specified boundary, we do not
        // need to add an additional judging condition of nodalDOF = globalDOF % dim
        if(dofLocation[global_dof][0] == x_min) { // x = 0.
          boundary_values_of_D[global_dof] = 0.;
          //boundary_values_of_V[global_dof] = 0.;
          boundary_values_of_a[global_dof] = 0.;
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

  //Get a vector of global degree-of-freedom x-coordinates
  MappingQ1<dim,dim> mapping;
  std::vector< Point<dim,double> > dof_coords(dof_handler.n_dofs());
  dofLocation.reinit(dof_handler.n_dofs(),dim);
  DoFTools::map_dofs_to_support_points<dim,dim>(mapping,dof_handler,dof_coords);
  for(unsigned int i=0; i<dof_coords.size(); i++){
    for(unsigned int j=0; j<dim; j++){
      dofLocation[i][j] = dof_coords[i][j];
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
  a_trans.reinit(dof_handler.n_dofs());
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
/**************************************************************************************/

template <int dim>
void FEM<dim>::assemble_system(){

  //For volume integration/quadrature points
  FEValues<dim> fe_values (fe,
                             quadrature_formula, 
                             update_values | 
                             update_gradients | 
                             update_JxW_values);

  K=0; F=0; M=0;

  const unsigned int dofs_per_elem = fe.dofs_per_cell;                      //dofs per element
  const unsigned int nodes_per_elem = GeometryInfo<dim>::vertices_per_cell;
  const unsigned int num_quad_pts = quadrature_formula.size();              //Total number of quad points in the element
  FullMatrix<double> Klocal (dofs_per_elem, dofs_per_elem);
  FullMatrix<double> Mlocal (dofs_per_elem, dofs_per_elem);

  std::vector<unsigned int> local_dof_indices (dofs_per_elem);              //relates local dof numbering to global dof numbering
  double rho = 7.6e3;

  //loop over elements  
  typename DoFHandler<dim>::active_cell_iterator elem = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (;elem!=endc; ++elem){

    fe_values.reinit(elem);

    elem->get_dof_indices (local_dof_indices);
    

    Mlocal = 0.;
    for(unsigned int q = 0; q < num_quad_pts; q++) {
      for(unsigned int A = 0; A < dofs_per_elem / dim; A++) {
        for(unsigned int i = 0; i < dim; i++) {
          for(unsigned int B = 0; B < dofs_per_elem / dim; B++) {
            for(unsigned int k = 0; k < dim; k++) {
                  Mlocal[dim*A+i][dim*B+k] += rho * fe_values.shape_value(dim*A+i,q) * (i==k)
                                    * fe_values.shape_value(dim*B+k,q) * fe_values.JxW(q);
            }
          }
        }
      }
    } // for(q)


    Klocal = 0.;
    for(unsigned int q = 0; q < num_quad_pts; q++) {
      for(unsigned int A = 0; A < dofs_per_elem / dim; A++) {
        for(unsigned int i = 0; i < dim; i++) {
          for(unsigned int B = 0; B < dofs_per_elem / dim; B++) {
            for(unsigned int k = 0; k < dim; k++) {
              for(unsigned int j = 0; j < dim; j++) {
                for(unsigned int l = 0; l < dim; l++) {
                  Klocal[dim*A+i][dim*B+k] += fe_values.shape_grad(dim*A+i,q)[j] * C(i,j,k,l)
                                    * fe_values.shape_grad(dim*B+k,q)[l] * fe_values.JxW(q);
                }
              }
            }
          }
        }
      }
    } // for(q)

    //Assemble local K and F into global K and F
    for(unsigned int i=0; i<dofs_per_elem; i++){
      for(unsigned int j=0; j<dofs_per_elem; j++){
        K.add(local_dof_indices[i], local_dof_indices[j], Klocal[i][j]);
        M.add(local_dof_indices[i], local_dof_indices[j], Mlocal[i][j]);
      }
    }
  } // for(elem)

}


/**************************************************************************************/
/*  This function directly solves "K_glocal * D_steady = F_global" at steady state by */
/*  implementing a matrix inverse - which is cost expensive, but of relatively        */
/*  high accuracy.                                                                    */
/**************************************************************************************/

template <int dim>
void FEM<dim>::solve_steady(){

  MatrixTools::apply_boundary_values (boundary_values_of_D, K, D_steady, F, false);

  SparseDirectUMFPACK  A;
  A.initialize(K);
  A.vmult (D_steady, F); //D_steady=K^{-1}*F

  output_steady_results();

}


/**************************************************************************************/
/*  This function apply initial conditions for the transient problem and also solves  */
/*  transient solutions at t = 0.                                                     */
/**************************************************************************************/

template <int dim>
void FEM<dim>::apply_initial_conditions(){

  const unsigned int totalNodes = dof_handler.n_dofs(); 

  D_trans = 0.;
  V_trans = 0.;
  for(unsigned int i=0; i<totalNodes; i++){
    if(i % 3 == 0)
    D_trans[i] = 0.01 * dofLocation[i][0]; 
  }

  system_matrix.copy_from(M);

  // Define the right-hand-side vector
  // where D_0 = D_trans
  K.vmult(RHS,D_trans); 
  RHS *= -1.;
  RHS.add(1.,F); 

  MatrixTools::apply_boundary_values (boundary_values_of_a, system_matrix, a_trans, RHS, false);
  SparseDirectUMFPACK  A;

  A.initialize(system_matrix);
  A.vmult (a_trans, RHS);

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
  const double delta_t = 1.e-6;
  double beta = 0.25, gamma = 0.5;

  const unsigned int totalNodes = dof_handler.n_dofs(); //Total number of nodes
  Vector<double>     D_tilde(totalNodes);
  Vector<double>     V_tilde(totalNodes);

  //Loop over time steps and update D_transient from D_n to D_{n+1} using the V method
  for(unsigned int t_step=1; t_step<1001; t_step++){

    //Update D_tilde
    D_tilde = 0.;
    D_tilde.add(1., D_trans);
    D_tilde.add(delta_t, V_trans);
    double const_a = delta_t * delta_t / 2. * (1. - 2. * beta);
    D_tilde.add(const_a, a_trans);

    //Update V_tilde
    V_tilde = 0;
    V_tilde.add(1., V_trans);
    double const_v = delta_t * ( 1 - gamma);
    V_tilde.add(const_v, a_trans);

    //Update a_trans
    system_matrix.copy_from(M);
    double const_a1 = delta_t * delta_t * beta;
    system_matrix.add(const_a1, K);

    K.vmult(RHS, D_tilde);
    RHS *= -1.;
    RHS.add(1., F);

    //Apply boundary conditions on a_trans before solving the matrix/vector system
    MatrixTools::apply_boundary_values (boundary_values_of_a, system_matrix, a_trans, RHS, false);

    //Solve for a_trans{n+1}
    SparseDirectUMFPACK  A;
    A.initialize(system_matrix);
    A.vmult (a_trans, RHS);

    //Update D_trans.
    D_trans = 0;
    D_trans.add(1., D_tilde);
    D_trans.add(const_a1, a_trans);

    //Update V_trans.
    V_trans = 0.;
    V_trans.add(1., V_tilde);
    double const_v2 = delta_t * gamma;
    V_trans.add(const_v2, a_trans);

    //Output the results every 100 seconds
    if(t_step == 350 || t_step == 700 || t_step == 1000){

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

  std::string filename;
  filename = "../sols/vtk/solution_ss.vtk";

  //Write results to VTK file
  std::ofstream output1 (filename);
  DataOut<dim> data_out; data_out.attach_dof_handler (dof_handler);

  //Add nodal DOF data
  data_out.add_data_vector (D_steady,
                            nodal_solution_names,
                            DataOut<dim>::type_dof_data,
                            nodal_data_component_interpretation);
  data_out.build_patches();
  data_out.write_vtk(output1);
  output1.close();
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
