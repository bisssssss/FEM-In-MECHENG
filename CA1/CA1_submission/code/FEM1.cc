#include "FEM1.h"

//Class constructor for a vector field
template <int dim>
FEM<dim>::FEM(unsigned int order,unsigned int problem)
:
fe (FE_Q<dim>(order), dim), 
  dof_handler (triangulation)
{
  basisFunctionOrder = order;
  if(problem == 1 || problem == 2 || problem == 3 || problem == 4){
    prob = problem;
  }
  else{
    std::cout << "Error: problem number should be 1, 2, 3 or 4.\n";
    exit(0);
  }

  for (unsigned int i=0; i<dim; ++i){
    nodal_solution_names.push_back("u");
    nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  }
}

//Class destructor
template <int dim>
FEM<dim>::~FEM(){
  dof_handler.clear();
}

//Find the value of xi at the given node (using deal.II node numbering)
template <int dim>
double FEM<dim>::xi_at_node(unsigned int dealNode){
  double xi;

  if(dealNode == 0){
    xi = -1.;
  }
  else if(dealNode == 1){
    xi = 1.;
  }
  else if(dealNode <= basisFunctionOrder){
    xi = -1. + 2.*(dealNode-1.)/basisFunctionOrder;
  }
  else{
    std::cout << "Error: you input node number "
        << dealNode << " but there are only " 
        << basisFunctionOrder + 1 << " nodes in an element.\n";
    exit(0);
  }

  return xi;
}


/**************************************************************************************/
/*  basis_function defines the polynomial function to evaluate state values.          */
/*  basis_gradient defines the gradient of basis polynomial functions wrt. xi.        */
/*  node :  specifies which node the basis function corresponds to,                   */
/*  xi   :  point (in the bi-unit domain) where the function is being evaluated.      */
/*  Note that in real computing, these functions are NOT called every time for the    */
/*  same "node" and "xi". Instead, the function and gradient values are stored in     */
/*  basis_func_map and basis_grad_map.                                                */
/**************************************************************************************/

template <int dim>
double FEM<dim>::basis_function(unsigned int node, double xi){

  double value = 1.;
  for(unsigned int Ni = 0; Ni < basisFunctionOrder + 1; Ni++){
    if(Ni != node) {
        value *= (xi - xi_at_node(Ni))/(xi_at_node(node) - xi_at_node(Ni));
    } //if
  } //for(Ni)

  return value;
}

template <int dim>
double FEM<dim>::basis_gradient(unsigned int node, double xi){
  
  double value = 0., LIdiv = 0., LIden = 1.;
  for(unsigned int Ni = 0; Ni < basisFunctionOrder + 1; Ni ++){

    if(Ni != node){
      double LIP = 1.;
      for(unsigned int Nj = 0; Nj < basisFunctionOrder + 1; Nj ++){
        if(Nj != node && Nj != Ni) LIP *= (xi - xi_at_node(Nj));
      } //for(Nj)
      LIdiv = LIdiv + LIP;
      LIden = LIden * (xi_at_node(node) - xi_at_node(Ni));
    } //if
  } //for(Ni)
  value = LIdiv / LIden;

  return value;
}


/**************************************************************************************/
/*  Define the problem domain and generate the mesh                                   */
/**************************************************************************************/

template <int dim>
void FEM<dim>::generate_mesh(unsigned int numberOfElements){

  L = 0.1; //Define the limits of domain
  double x_min = 0.;
  double x_max = L;

  Point<dim,double> min(x_min), max(x_max);
  std::vector<unsigned int> meshDimensions (dim,numberOfElements);
  GridGenerator::subdivided_hyper_rectangle (triangulation, meshDimensions, min, max);
}


/**************************************************************************************/
/*  This function specifies and stores Dirichilet boundary conditions.                */
/*  After assembling Kglobal and Flocal, theses b.c.s will be automatically applied   */
/*  to resive matrices of the system with techniques supported by deal.II.            */
/**************************************************************************************/

template <int dim>
void FEM<dim>::define_boundary_conds(){
  const unsigned int totalNodes = dof_handler.n_dofs();

  //Apply Dirichile b.c.s
  for(unsigned int globalNode=0; globalNode<totalNodes; globalNode++){
    if(nodeLocation[globalNode] == 0) boundary_values[globalNode] = g1;
    if(nodeLocation[globalNode] == L){
      if(prob == 1 || prob == 3) boundary_values[globalNode] = g2;
    }
  }
}


/**************************************************************************************/
/*  This function sets up data structures for Finite Element Solving System.          */
/*  The parameters and quadrature rule are defined here as well.                      */
/**************************************************************************************/

template <int dim>
void FEM<dim>::setup_system(){

  //Define constants for problem (Dirichlet boundary values)
  E = 1.e11; Area = 1.e-4; // 1D-Bar specifications
  g1 = 0; g2 = 0.001; f_bar = 1.e11; f_hat = 1.e12; h_in = 1.e6; // B.C. parameters

  dof_handler.distribute_dofs (fe);   //Let deal.II organize degrees of freedom

  MappingQ1<dim,dim> mapping;
  std::vector< Point<dim,double> > dof_coords(dof_handler.n_dofs());
  nodeLocation.resize(dof_handler.n_dofs());
  DoFTools::map_dofs_to_support_points<dim,dim>(mapping,dof_handler,dof_coords);
  for(unsigned int i=0; i<dof_coords.size(); i++){
    nodeLocation[i] = dof_coords[i][0];
  }

  define_boundary_conds(); //Specify boundary condtions

  sparsity_pattern.reinit (dof_handler.n_dofs(), dof_handler.n_dofs(),
         dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress();
  K.reinit (sparsity_pattern);
  F.reinit (dof_handler.n_dofs());
  D.reinit (dof_handler.n_dofs());

  //Define quadrature rule
  quadRule = 5; //Number of quadrature points
  quad_points.resize(quadRule); quad_weight.resize(quadRule);

  quad_points[0] = -1./3. * sqrt(5. + 2. * sqrt(10. / 7.));
  quad_points[1] = -1./3. * sqrt(5. - 2. * sqrt(10. / 7.));
  quad_points[2] = 0.;
  quad_points[3] = - quad_points[1];
  quad_points[4] = - quad_points[0];

  quad_weight[0] = (322. - 13. * sqrt(70.)) / 900.;
  quad_weight[1] = (322. + 13. * sqrt(70.)) / 900.;
  quad_weight[2] = 128. / 225.;
  quad_weight[3] = quad_weight[1];
  quad_weight[4] = quad_weight[0];

  //store repetitive computed basis functions and gradients in two 2D-arraies
  basis_func_map.resize(basisFunctionOrder + 1);
  basis_grad_map.resize(basisFunctionOrder + 1);
  for(unsigned int A = 0; A < basisFunctionOrder + 1; A++) {
    basis_func_map[A].resize(quadRule);
    basis_grad_map[A].resize(quadRule);
    for(unsigned int qi = 0; qi < quadRule; qi++) {
      basis_func_map[A][qi] = basis_function(A, quad_points[qi]);
      basis_grad_map[A][qi] = basis_gradient(A, quad_points[qi]);
    }
  }

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

  K=0; F=0;

  const unsigned int        dofs_per_elem = fe.dofs_per_cell; //This gives you number of degrees of freedom per element
  FullMatrix<double>        Klocal (dofs_per_elem, dofs_per_elem);
  Vector<double>            Flocal (dofs_per_elem);
  std::vector<unsigned int> local_dof_indices (dofs_per_elem);
  double                    h_e, x, f;

  //loop over elements  
  typename DoFHandler<dim>::active_cell_iterator elem = dof_handler.begin_active(), endc = dof_handler.end();
  for (;elem!=endc; ++elem){

    elem->get_dof_indices (local_dof_indices);
    h_e = nodeLocation[local_dof_indices[1]] - nodeLocation[local_dof_indices[0]]; //element length
    std::vector<unsigned int> local2global  = local_dof_indices;

    //Loop over local DOFs and quadrature points to populate F_local and K_local.
    Flocal = 0.;
    for(unsigned int A=0; A<dofs_per_elem; A++){
      for(unsigned int q=0; q<quadRule; q++){
        x = 0.;
        for(unsigned int B=0; B<dofs_per_elem; B++){
          x += nodeLocation[local_dof_indices[B]] * basis_func_map[B][q];
        }
        if(prob == 1 || prob == 2)
          Flocal[A] += quad_weight[q] * f_bar * x * basis_func_map[A][q] * Area * h_e / 2;
        else if(prob == 3 || prob == 4)
          Flocal[A] += quad_weight[q] * f_hat * x * (L-x) * basis_func_map[A][q] * Area * h_e / 2;
      }
    }
    
    //Add nonzero Neumann condition for prob 2 or 4, if applicable
    if(prob == 2 || prob == 4){
      if(nodeLocation[local_dof_indices[1]] == L)
          Flocal[1] += h_in;
    }

    //Loop over local DOFs and quadrature points to populate Klocal
    Klocal = 0;
    for(unsigned int A=0; A<dofs_per_elem; A++){
      for(unsigned int B=0; B<dofs_per_elem; B++){
        for(unsigned int q=0; q<quadRule; q++){
          Klocal[A][B] += quad_weight[q] * basis_grad_map[A][q]
                                * basis_grad_map[B][q] * E * Area / h_e * 2;
        } // for(q)
      } // for(B)
    } // for(A)

    //Assemble local K and F into global K and F
    for(unsigned int A=0; A<dofs_per_elem; A++){
      F[local_dof_indices[A]] += Flocal[A]; 
      for(unsigned int B=0; B<dofs_per_elem; B++){
        K.add(local_dof_indices[A], local_dof_indices[B], Klocal[A][B]);
      }
    }

  } // for(elem)

  //Apply Dirichlet boundary conditions
  MatrixTools::apply_boundary_values (boundary_values, K, D, F, false);
}


/**************************************************************************************/
/*  This function directly solves "K_glocal * D = F_global" by implementing a matrix  */
/*  inverse - which is cost expensive, but of relatively high accuracy.               */
/**************************************************************************************/

template <int dim>
void FEM<dim>::solve(){

  SparseDirectUMFPACK  A;
  A.initialize(K);
  A.vmult (D, F); //D=K^{-1}*F

}


/**************************************************************************************/
/*  This function writes nodal solutions to "solution_Prob_NumMesh_Order.vtk" for     */
/*  assignment submission and results checking.                                       */
/**************************************************************************************/

template <int dim>
void FEM<dim>::output_results (){

  std::string filename = "../sols/vtk/solution_p" + std::to_string(prob) 
                    + "_" + std::to_string(triangulation.n_active_cells()) + "_";
  switch(basisFunctionOrder) {
    case 1: filename += "linear.vtk"; break;
    case 2: filename += "quadratic.vtk"; break;
    case 3: filename += "cubic.vtk"; break;
    default: filename += "undefinedOrder.vtk"; break;
  }
  std::ofstream output1(filename);
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);

  //Add nodal DOF data
  data_out.add_data_vector(D, nodal_solution_names, DataOut<dim>::type_dof_data,
         nodal_data_component_interpretation);
  data_out.build_patches();
  data_out.write_vtk(output1);
  output1.close();
}


/**************************************************************************************/
/*  The four functions below are related to l2 and h1 error norm implementaion.       */
/*  Function uexact(x) and duexact(x) computes exact values and its derivatives       */
/*  at a global coordinate x, and the two error norm functions sums error over        */
/*  each element, returns a normalized value.                                         */
/**************************************************************************************/

template <int dim>
double FEM<dim>::l2norm_of_error(){
  
  double l2norm = 0.;

  const unsigned int        dofs_per_elem = fe.dofs_per_cell;
  std::vector<unsigned int> local_dof_indices (dofs_per_elem);
  double u_h, x, h_e, u_exact;

  //loop over elements  
  typename DoFHandler<dim>::active_cell_iterator elem = dof_handler.begin_active (), 
    endc = dof_handler.end();
  for (;elem!=endc; ++elem){

    //Retrieve the effective "connectivity matrix" for this element
    elem->get_dof_indices (local_dof_indices);
    h_e = nodeLocation[local_dof_indices[1]] - nodeLocation[local_dof_indices[0]]; //Element length

    for(unsigned int q=0; q<quadRule; q++){
      x = 0.; u_h = 0.;

      for(unsigned int B=0; B<dofs_per_elem; B++){
        x += nodeLocation[local_dof_indices[B]]*basis_func_map[B][q];
        u_h += D[local_dof_indices[B]]*basis_func_map[B][q];
      }
      u_exact =  uexact(x);
      l2norm += quad_weight[q] * (u_exact - u_h) * (u_exact - u_h) * h_e / 2; //(h_e/2) corresponds to a chain rule
    } //for(q)
  } //for(elem)

  return sqrt(l2norm);
}

template <int dim>
double FEM<dim>::h1norm_of_error(){
  
  double h1norm = 0.;

  const unsigned int        dofs_per_elem = fe.dofs_per_cell;
  std::vector<unsigned int> local_dof_indices (dofs_per_elem);
  double u_exact, u_h, x, h_e, u_h_x, du_exact;

  //loop over elements  
  typename DoFHandler<dim>::active_cell_iterator elem = dof_handler.begin_active (), 
    endc = dof_handler.end();
  for (;elem!=endc; ++elem){

    //Retrieve the effective "connectivity matrix" for this element
    elem->get_dof_indices (local_dof_indices);

    h_e = nodeLocation[local_dof_indices[1]] - nodeLocation[local_dof_indices[0]]; //Element length

    for(unsigned int q=0; q<quadRule; q++){
      x = 0.; u_h = 0.; u_h_x = 0.;

      for(unsigned int B=0; B<dofs_per_elem; B++){
        x += nodeLocation[local_dof_indices[B]]*basis_func_map[B][q];
        u_h += D[local_dof_indices[B]]*basis_func_map[B][q];
        u_h_x += D[local_dof_indices[B]]*basis_grad_map[B][q]/(h_e/2); //(h_e/2) corresponds to a chain rule
      }
      u_exact =  uexact(x);
      du_exact = duexact(x);
      h1norm += quad_weight[q] * ((u_exact - u_h) * (u_exact - u_h)
               + (L * L) * (du_exact - u_h_x) * (du_exact - u_h_x)) * h_e / 2;
    }
  }

  return sqrt(h1norm);
}

template <int dim>
double FEM<dim>::uexact(double x_in){

  double u_exact, ux_0;

  switch(prob){
    case 1:
      ux_0 = (g2 - g1 + f_bar / E * (pow(L, 3) / 6.)) / L;
      u_exact = - f_bar / (6. * E) * pow(x_in, 3) + ux_0 * x_in + g1;
      break;
    case 2:
      ux_0 = (h_in + Area * f_bar / 2. * pow(L, 2)) / (E * Area);
      u_exact = - f_bar / (6. * E) * pow(x_in, 3) + ux_0 * x_in + g1;
      break;
    case 3:
      ux_0 = (g2 - g1 + f_hat / E * (pow(L, 4) / 12.)) / L;
      u_exact = - f_hat / E * (L * pow(x_in, 3) / 6. - 1. / 12. * pow(x_in, 4)) + ux_0 * x_in + g1;
      break;
    case 4:
      ux_0 = (h_in + Area * f_hat / 6. * pow(L, 3)) / (E * Area);
      u_exact = - f_hat / E * (L * pow(x_in, 3) / 6. - 1. / 12. * pow(x_in, 4)) + ux_0 * x_in + g1;
      break;
  } // switch(prob)

  return u_exact;

}

template <int dim>
double FEM<dim>::duexact(double x_in){

  double du_exact, ux_0;

  switch(prob){
    case 1:
      ux_0 = (g2 - g1 + f_bar / E * (pow(L, 3) / 6.)) / L;
      du_exact = - f_bar / (2. * E) * pow(x_in, 2) + ux_0;
      break;
    case 2:
      ux_0 = (h_in + Area * f_bar / 2. * pow(L, 2)) / (E * Area);
      du_exact = - f_bar / (2. * E) * pow(x_in, 2) + ux_0;
      break;
    case 3:
      ux_0 = (g2 - g1 + f_hat / E * (pow(L, 4) / 12.)) / L;
      du_exact = - f_hat / E * (L * pow(x_in, 2) / 2. - 1. / 3. * pow(x_in, 3)) + ux_0;
      break;
    case 4:
      ux_0 = (h_in + Area * f_hat / 6. * pow(L, 3)) / (E * Area);
      du_exact = - f_hat / E * (L * pow(x_in, 2) / 2. - 1. / 3. * pow(x_in, 3)) + ux_0;
      break;
  } // switch(prob)

  return du_exact;

}


template class FEM<1>;
template class FEM<2>;
template class FEM<3>;
