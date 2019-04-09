#include "FEM2b.h"

// Class constructor for a scalar field
template <int dim>
FEM<dim>::FEM (unsigned int order)
:
fe (FE_Q<dim>(QIterated<1>(QTrapez<1>(),order)),1),
  dof_handler (triangulation)
{
  basisFunctionOrder = order;

  nodal_solution_names.push_back("D");
  nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
}

//Class destructor
template <int dim>
FEM<dim>::~FEM (){
  dof_handler.clear ();
}


/**************************************************************************************/
/*  function xi_at_node returns xi value with an 3d input node index.                 */
/*  The input dealNode should be within [0, (order +1)^3 -1 ]                         */
/*  The output xi is a vector with three elements: xi_x, xi_y, xi_z.                  */
/**************************************************************************************/

template <int dim>
std::vector<double> FEM<dim>::xi_at_node(unsigned int dealNode){
  std::vector<double> xi (dim, 0.0);

  //Node on vertices
  if(dealNode == 0){
    xi[0] = -1.; xi[1] = -1.; xi[2] = -1.;
  }
  else if(dealNode == 1){
    xi[0] = 1.; xi[1] = -1.; xi[2] = -1.;
  }
  else if(dealNode == 2){
    xi[0] = -1.; xi[1] = 1.; xi[2] = -1.;
  }
  else if(dealNode == 3){
    xi[0] = 1.; xi[1] = 1.; xi[2] = -1.;
  }
  else if(dealNode == 4){
    xi[0] = -1.; xi[1] = -1.; xi[2] = 1.;
  }
  else if(dealNode == 5){
    xi[0] = 1.; xi[1] = -1.; xi[2] = 1.;
  }
  else if(dealNode == 6){
    xi[0] = -1.; xi[1] = 1.; xi[2] = 1.;
  }
  else if(dealNode == 7){
    xi[0] = 1.; xi[1] = 1.; xi[2] = 1.;
  }
  // Node on edges
  else if(dealNode >= 8 && dealNode < 8 + 12 * (basisFunctionOrder-1)){
    int rem = (dealNode-8) / (basisFunctionOrder-1);
    int mod = (dealNode-8) % (basisFunctionOrder-1);
    if(rem == 0){
      xi[0] = -1.; xi[1] = -1. + 2.* (mod+1) / basisFunctionOrder; xi[2] = -1.;
    }
    else if(rem == 1){
      xi[0] = 1.; xi[1] = -1. + 2.* (mod+1) / basisFunctionOrder; xi[2] = -1.;
    }
    else if(rem == 2){
      xi[1] = -1.; xi[0] = -1. + 2.* (mod+1) / basisFunctionOrder; xi[2] = -1.;
    }
    else if(rem == 3){
      xi[1] = 1.; xi[0] = -1. + 2.* (mod+1) / basisFunctionOrder; xi[2] = -1.;
    }
    else if(rem == 4){
      xi[0] = -1.; xi[1] = -1. + 2.* (mod+1) / basisFunctionOrder; xi[2] = 1.;
    }
    else if(rem == 5){
      xi[0] = 1.; xi[1] = -1. + 2.* (mod+1) / basisFunctionOrder; xi[2] = 1.;
    }
    else if(rem == 6){
      xi[1] = -1.; xi[0] = -1. + 2.* (mod+1) / basisFunctionOrder; xi[2] = 1.;
    }
    else if(rem == 7){
      xi[1] = 1.; xi[0] = -1. + 2.* (mod+1) / basisFunctionOrder; xi[2] = 1.;
    }
    else if(rem == 8){
      xi[0] = -1.; xi[1] = -1.; xi[2]= -1. + 2.* (mod+1) / basisFunctionOrder;
    }
    else if(rem == 9){
      xi[0] = 1.; xi[1] = -1.; xi[2]= -1. + 2.* (mod+1) / basisFunctionOrder;
    }
    else if(rem == 10){
      xi[0] = -1.; xi[1] = 1.; xi[2]= -1. + 2.* (mod+1) / basisFunctionOrder;
    }
    else if(rem == 11){
      xi[0] = 1.; xi[1] = 1.; xi[2]= -1. + 2.* (mod+1) / basisFunctionOrder;
    }
  }
  // Node on faces
  else if(dealNode >= 8 + 12 * (basisFunctionOrder-1)
              && dealNode < 8 + 12 * (basisFunctionOrder-1) 
              + 6 * (basisFunctionOrder-1) * (basisFunctionOrder-1)) {
    int prev = 8 + 12 * (basisFunctionOrder-1);
    int rem = (dealNode-prev) / ((basisFunctionOrder-1) * (basisFunctionOrder-1));
    int mod = (dealNode-prev) % ((basisFunctionOrder-1) * (basisFunctionOrder-1));
    int rem1 = mod / (basisFunctionOrder-1);
    int mod1 = mod % (basisFunctionOrder-1);
    if(rem == 0){
      xi[0] = -1.;
      xi[1] = -1. + 2.* (mod1+1) / basisFunctionOrder;
      xi[2] = -1. + 2.* (rem1+1) / basisFunctionOrder;
    }
    else if(rem == 1){
      xi[0] = 1.;
      xi[1] = -1. + 2.* (mod1+1) / basisFunctionOrder;
      xi[2] = -1. + 2.* (rem1+1) / basisFunctionOrder;
    }
    else if(rem == 2){
      xi[1] = -1.;
      xi[0] = -1. + 2.* (mod1+1) / basisFunctionOrder;
      xi[2] = -1. + 2.* (rem1+1) / basisFunctionOrder;
    }
    else if(rem == 3){
      xi[1] = 1.;
      xi[0] = -1. + 2.* (mod1+1) / basisFunctionOrder;
      xi[2] = -1. + 2.* (rem1+1) / basisFunctionOrder;
    }
    else if(rem == 4){
      xi[2] = -1.;
      xi[0] = -1. + 2.* (mod1+1) / basisFunctionOrder;
      xi[1] = -1. + 2.* (rem1+1) / basisFunctionOrder;
    }
    else if(rem == 5){
      xi[2] = 1.;
      xi[0] = -1. + 2.* (mod1+1) / basisFunctionOrder;
      xi[1] = -1. + 2.* (rem1+1) / basisFunctionOrder;
    }
  }
  // Node inside element
  else if(dealNode >= 8 + 12 * (basisFunctionOrder-1) 
              + 6 * (basisFunctionOrder-1) * (basisFunctionOrder-1) &&
       dealNode < (basisFunctionOrder+1)*(basisFunctionOrder+1)*(basisFunctionOrder+1)){
    int prev = 8 + 12 * (basisFunctionOrder-1) 
              + 6 * (basisFunctionOrder-1) * (basisFunctionOrder-1);
    int rem = (dealNode-prev) / ((basisFunctionOrder-1) * (basisFunctionOrder-1));
    int mod = (dealNode-prev) % ((basisFunctionOrder-1) * (basisFunctionOrder-1));
    int rem1 = mod / (basisFunctionOrder-1);
    int mod1 = mod % (basisFunctionOrder-1);
    xi[0] = -1. + 2.* (mod1+1) / basisFunctionOrder;
    xi[1] = -1. + 2.* (rem1+1) / basisFunctionOrder;
    xi[2] = -1. + 2.* (rem+1) / basisFunctionOrder;
  }
  else{
    std::cout << "Error: you input node number "
        << dealNode << " but there are only " 
        << (basisFunctionOrder+1)*(basisFunctionOrder+1)*(basisFunctionOrder+1) - 1 
        << " nodes in an element.\n";
    exit(0);
  }

  return xi;
}


/**************************************************************************************/
/*  Similar to function xi_at_node, this function returns xi value                    */
/*  with an 1d input node index.                                                      */
/*  The input dealNode should be within [0, order]                                    */
/*  The output xi is a scalar value: xi                                               */
/**************************************************************************************/

template <int dim>
double FEM<dim>::xi_at_node_1d(unsigned int dealNode){
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
/*  basis_gradient defines the gradient of basis polynomial functions wrt. xi (3x1).  */
/*  node :  specifies which node the basis function corresponds to,                   */
/*  xi   :  point (in the bi-unit domain) where the function is being evaluated.      */
/*  Note that in real computing, these functions are NOT called every time for the    */
/*  same "node" and "xi". Instead, the function and gradient values are stored in     */
/*  basis_func_map and basis_grad_map.                                                */
/**************************************************************************************/

template <int dim>
double FEM<dim>::basis_function(unsigned int node, double xi_1, double xi_2, double xi_3){
 
  double value = 1.;

  std::vector<double> idx = xi_at_node(node);

  for(unsigned int Ni = 0; Ni < basisFunctionOrder+1; Ni++){
        if(idx[0] != xi_at_node_1d(Ni)) {
          value *= (xi_1 - xi_at_node_1d(Ni)) /(idx[0] - xi_at_node_1d(Ni));
        }
        if(idx[1] != xi_at_node_1d(Ni)) {
          value *= (xi_2 - xi_at_node_1d(Ni)) /(idx[1] - xi_at_node_1d(Ni));
        }
        if(idx[2] != xi_at_node_1d(Ni)) {
          value *= (xi_3 - xi_at_node_1d(Ni)) /(idx[2] - xi_at_node_1d(Ni));
        }
    } //for(Ni)

  return value;
}

//Define basis function gradient
template <int dim>
std::vector<double> FEM<dim>::basis_gradient(unsigned int node, double xi_1, double xi_2, double xi_3){

  std::vector<double> values(dim, 0.0); 
  std::vector<double> idx = xi_at_node(node);
  double LIdiv1 = 0., LIdiv2 = 0., LIdiv3 = 0.;
  double LIden1 = 1., LIden2 = 1., LIden3 = 1.;
  double param1 = 1., param2 = 1., param3 = 1.;

  for(unsigned int Ni = 0; Ni < basisFunctionOrder + 1; Ni++){

      if(idx[0] != xi_at_node_1d(Ni)) {
        double LIP1 = 1.;
        for(unsigned int Nk = 0; Nk < basisFunctionOrder + 1; Nk ++){
          if(xi_at_node_1d(Nk) != idx[0] && Nk != Ni) 
            LIP1 *= (xi_1 - xi_at_node_1d(Nk));
        } //for(Nk)
        param1 *= (xi_1 - xi_at_node_1d(Ni)) /(idx[0] - xi_at_node_1d(Ni));
        LIdiv1 += LIP1;
        LIden1 *= (idx[0] - xi_at_node_1d(Ni));
      }

      if(idx[1] != xi_at_node_1d(Ni)) {
        double LIP2 = 1.;
        for(unsigned int Nk = 0; Nk < basisFunctionOrder + 1; Nk ++){
          if(xi_at_node_1d(Nk) != idx[1] && Nk != Ni) 
            LIP2 *= (xi_2 - xi_at_node_1d(Nk));
        } //for(Nk)
        param2 *= (xi_2 - xi_at_node_1d(Ni)) /(idx[1] - xi_at_node_1d(Ni));
        LIdiv2 += LIP2;
        LIden2 *= (idx[1] - xi_at_node_1d(Ni));
      }

      if(idx[2] != xi_at_node_1d(Ni)) {
        double LIP3 = 1.;
        for(unsigned int Nk = 0; Nk < basisFunctionOrder + 1; Nk ++){
          if(xi_at_node_1d(Nk) != idx[2] && Nk != Ni) 
            LIP3 *= (xi_3 - xi_at_node_1d(Nk));
        } //for(Nk)
        param3 *= (xi_3 - xi_at_node_1d(Ni)) /(idx[2] - xi_at_node_1d(Ni));
        LIdiv3 += LIP3;
        LIden3 *= (idx[2] - xi_at_node_1d(Ni));
      }

  } //for(Ni)

  values[0] = LIdiv1 / LIden1 * param2 * param3;
  values[1] = LIdiv2 / LIden2 * param1 * param3;
  values[2] = LIdiv3 / LIden3 * param2 * param1;

  return values;
}


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
    z_max = 1.;

  Point<dim,double> min(x_min,y_min,z_min),
    max(x_max,y_max,z_max);
  GridGenerator::subdivided_hyper_rectangle (triangulation, numberOfElements, min, max);

}


/**************************************************************************************/
/*  This function specifies and stores Dirichilet boundary conditions.                */
/*  After assembling Kglobal and Flocal, theses b.c.s will be automatically applied   */
/*  to resive matrices of the system with techniques supported by deal.II.            */
/**************************************************************************************/

template <int dim>
void FEM<dim>::define_boundary_conds(){

  //Define the Dirichlet boundary conditions.
  const unsigned int totalNodes = dof_handler.n_dofs(); //Total number of nodes

  double c0 = 1./3.;
    for(unsigned int globalNode=0; globalNode<totalNodes; globalNode++){     
      if(nodeLocation[globalNode][0] == x_min) { // x = 0
        boundary_values[globalNode] = 300. * (1. + c0 
            * (nodeLocation[globalNode][1] + nodeLocation[globalNode][2]));
      }
      if(nodeLocation[globalNode][0] == x_max) // x = 0.04
        boundary_values[globalNode] = 310. * (1. + c0 
                * (nodeLocation[globalNode][1] + nodeLocation[globalNode][2]));
      //std::cout << boundary_values[globalNode] << std::endl;
    }
}


/**************************************************************************************/
/*  This function sets up data structures for Finite Element Solving System.          */
/*  The parameters and quadrature rule are defined here as well.                      */
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
  basis_func_map.resize(fe.dofs_per_cell);
  basis_grad_map.resize(fe.dofs_per_cell);
  for(unsigned int A = 0; A < fe.dofs_per_cell; A++) {
    basis_func_map[A].resize(quadRule);
    basis_grad_map[A].resize(quadRule);
    for(unsigned int q1 = 0; q1 < quadRule; q1++) {
      basis_func_map[A][q1].resize(quadRule);
      basis_grad_map[A][q1].resize(quadRule);
      for(unsigned int q2 = 0; q2 < quadRule; q2++) {
        basis_func_map[A][q1][q2].resize(quadRule);
        basis_grad_map[A][q1][q2].resize(quadRule);
        for(unsigned int q3 = 0; q3 < quadRule; q3++) {
          basis_func_map[A][q1][q2][q3] = 
                  basis_function(A, quad_points[q1], quad_points[q2], quad_points[q3]);
          basis_grad_map[A][q1][q2][q3].resize(dim);
          for (unsigned j = 0; j < dim; j++){
            basis_grad_map[A][q1][q2][q3][j]= basis_gradient(A, quad_points[q1], 
                  quad_points[q2], quad_points[q3])[j];
          }
        }
      }
    }
  }

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

  K=0; F=0;

  const unsigned int        dofs_per_elem = fe.dofs_per_cell; //This gives you number of degrees of freedom per element
  FullMatrix<double>      Klocal (dofs_per_elem, dofs_per_elem);
  Vector<double>            Flocal (dofs_per_elem);
  std::vector<unsigned int> local_dof_indices (dofs_per_elem);

  //loop over elements  
  typename DoFHandler<dim>::active_cell_iterator elem = dof_handler.begin_active(), 
    endc = dof_handler.end();
  for (;elem!=endc; ++elem){

    elem->get_dof_indices (local_dof_indices);

    //Loop over local DOFs and quadrature points to populate Flocal and Klocal.
    FullMatrix<double> Jacobian(dim,dim);
    double detJ, f = 0.;

    //Loop over local DOFs and quadrature points to populate Flocal
    Flocal = 0.;
    for(unsigned int q1 = 0; q1 < quadRule; q1++) {
      for(unsigned int q2 = 0; q2 < quadRule; q2++) {
        for(unsigned int q3 = 0; q3 < quadRule; q3++) {
          Jacobian = 0.;
          for(unsigned int i = 0; i < dim; i++) {
            for(unsigned int j = 0; j < dim; j++) {
              for(unsigned int A = 0; A < dofs_per_elem; A++) {
                Jacobian[i][j] += nodeLocation[local_dof_indices[A]][i] * basis_grad_map[A][q1][q2][q3][j];
              }
            }
          }
          detJ = Jacobian.determinant();
          for(unsigned int A = 0; A < dofs_per_elem; A ++) {
            Flocal[A] += quad_weight[q1] * quad_weight[q2] * quad_weight[q3]
                                       * basis_func_map[A][q1][q2][q3] * f * detJ;
          }
        }
      }
    }

    //Loop over quadrature  points and local DOFs to populate Jacobian and Flocal 
    FullMatrix<double> invJacob(dim,dim), kappa(dim,dim);

    //"kappa" is the conductivity tensor
    kappa = 0.;
    kappa[0][0] = 385.;
    kappa[1][1] = 385.;
    kappa[2][2] = 385.;

    //Loop over local DOFs and quadrature points to populate Klocal
    Klocal = 0.;
    for(unsigned int q1 = 0; q1 < quadRule; q1++) {
      for(unsigned int q2 = 0; q2 < quadRule; q2++) {
        for(unsigned int q3 = 0; q3 < quadRule; q3++) {
          Jacobian = 0.;
          for(unsigned int i = 0; i < dim; i++) {
            for(unsigned int j = 0; j < dim; j++) {
              for(unsigned int A = 0; A < dofs_per_elem; A++) {
                Jacobian[i][j] += nodeLocation[local_dof_indices[A]][i] * basis_grad_map[A][q1][q2][q3][j];
              }
            }
          }
          detJ = Jacobian.determinant(); //determinant of Jacobian
          
          invJacob.invert(Jacobian); //inverse of Jacobian
          for(unsigned int A = 0; A < dofs_per_elem; A ++) {
            for(unsigned int B = 0; B < dofs_per_elem; B++) {
              for(unsigned int i = 0; i < dim; i++) {
                for(unsigned int j = 0; j < dim; j++) {
                  for(unsigned int I = 0; I < dim; I++) {
                    for(unsigned int J = 0; J < dim; J++) {
                        Klocal[A][B] += quad_weight[q1] * quad_weight[q2] * quad_weight[q3] * detJ * kappa[I][J]
                                  * basis_grad_map[A][q1][q2][q3][i] * invJacob[i][I]
                                  * basis_grad_map[B][q1][q2][q3][j] * invJacob[j][J];
                    }
                  }
                }
              }
            }
          } // for(A)
        } // for(q3)
      } // for(q2)
    } // for(q1)

    //assemble local K and F to global K and F
    for(unsigned int A = 0; A < dofs_per_elem; A++) {
      F[local_dof_indices[A]] += Flocal[A]; 
      for(unsigned int B = 0; B < dofs_per_elem; B++) {
        K.add(local_dof_indices[A], local_dof_indices[B], Klocal[A][B]);
      }
    }

  } // for(elem)






















    FullMatrix<double> Jacobian(dim,dim);
    for(unsigned int q1 = 0; q1 < quadRule; q1++) {
      for(unsigned int q2 = 0; q2 < quadRule; q2++) {
        for(unsigned int q3 = 0; q3 < quadRule; q3++) {
          Jacobian = 0.;
          for(unsigned int i = 0; i < dim; i++) {
            for(unsigned int j = 0; j < dim; j++) {
              for(unsigned int A = 0; A < dofs_per_elem; A++) {
                Jacobian[i][j] += nodeLocation[local_dof_indices[A]][i] * basis_grad_map[A][q1][q2][q3][j];
              }
            }
          }
          for(unsigned int A = 0; A < dofs_per_elem; A++) {
            for(unsigned int i = 0; i < dim; i++) {
              //std::cout << q1*quadRule +q2*quadRule+q3 << " " << A << " " << i << " " << basis_grad_map[A][q1][q2][q3][i] << std::endl;
            }
          }
          double detJ = Jacobian.determinant(); //determinant of Jacobian
          std::cout << q1 << " " << q2 << " " << q3 << " " << detJ << " " << quad_weight[q1] * quad_weight[q2] * quad_weight[q3] * detJ << std::endl;
        }
      }
    }
    std::cout << dofs_per_elem << std::endl;






















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
/*  This function writes nodal solutions to "solution_P3_Order.vtk" for               */
/*  assignment submission and results checking.                                       */
/**************************************************************************************/

template <int dim>
void FEM<dim>::output_results (){

  std::string filename = "../sols/vtk/solution_p3_order" 
              + std::to_string(basisFunctionOrder) + ".vtk";

  //Write results to VTK file
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


template class FEM<1>;
template class FEM<2>;
template class FEM<3>;
