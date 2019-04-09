//Include files
//Data structures and solvers
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
//Mesh related classes
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
//Finite element implementation classes
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
//Standard C++ libraries
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>

using namespace dealii;

template <int dim>
class FEM
{
 public:
  //Class functions
  FEM(unsigned int order);  // Class constructor 
  ~FEM(); //Class destructor

  //Functions to find the value of xi at the given node (using deal.II node numbering)
  std::vector<double> xi_at_node(unsigned int dealNode);
  double xi_at_node_1d(unsigned int dealNode);

  //Define your 3D basis functions and derivatives
  double basis_function(unsigned int node, 
			double xi_1,
			double xi_2,
			double xi_3);
  std::vector<double> basis_gradient(unsigned int node, 
				     double xi_1,
				     double xi_2,
				     double xi_3);

  //Solution steps
  void generate_mesh(std::vector<unsigned int> numberOfElements);
  void define_boundary_conds();
  void setup_system();
  void assemble_system();
  void solve();
  void output_results();

  //Class objects
  Triangulation<dim>   triangulation; //mesh
  FESystem<dim>        fe; 	      //FE element
  DoFHandler<dim>      dof_handler;   //Connectivity matrices

  //Gaussian quadrature
  unsigned int	      quadRule;    //quadrature rule, i.e. number of quadrature points
  std::vector<double> quad_points; //vector of Gauss quadrature points
  std::vector<double> quad_weight; //vector of the quadrature point weights
  std::vector<std::vector<std::vector<std::vector<double>>>> basis_func_map;   //store the calculated basis function values
  std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> basis_grad_map;   //store the calculated basis gradient values

  //Data structures
  SparsityPattern    	        sparsity_pattern; //Sparse matrix pattern
  SparseMatrix<double>	        K; 		  //Global stiffness (sparse) matrix
  Vector<double>       	        D, F; 		  //Global vectors - Solution vector (D) and Global force vector (F)
  Table<2,double>	        nodeLocation;     //Table of the coordinates of nodes by global dof number
  std::map<unsigned int,double> boundary_values;  //Map of dirichlet boundary conditions 
  unsigned int                  basisFunctionOrder;  //Case specification
  double                        x_min, x_max, y_min, y_max, z_min, z_max;   //Boundary conditions

  //solution name array
  std::vector<std::string> nodal_solution_names;
  std::vector<DataComponentInterpretation::DataComponentInterpretation> nodal_data_component_interpretation;
};
