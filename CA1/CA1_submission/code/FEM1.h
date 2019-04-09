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
  FEM (unsigned int order,unsigned int problem); // Class constructor 
  ~FEM(); //Class destructor

  //Function to find the value of xi at the given node (using deal.II node numbering)
  double xi_at_node(unsigned int dealNode);

  //Define your 1D basis functions and derivatives
  double basis_function(unsigned int node, double xi);
  double basis_gradient(unsigned int node, double xi);

  //Solution steps
  void generate_mesh(unsigned int numberOfElements);
  void define_boundary_conds();
  void setup_system();
  void assemble_system();
  void solve();
  void output_results();

  //Function to calculate the l2 and h1 norm of the error in the finite element sol'n vs. the exact solution
  double l2norm_of_error();
  double h1norm_of_error();
  double uexact(double x_in);
  double duexact(double x_in);

  //Class objects
  Triangulation<dim>   triangulation; //mesh
  FESystem<dim>        fe;        //FE element
  DoFHandler<dim>      dof_handler;   //Connectivity matrices

  //Gaussian quadrature
  unsigned int          quadRule;    //quadrature rule, i.e. number of quadrature points
  std::vector<double> quad_points; //vector of Gauss quadrature points
  std::vector<double> quad_weight; //vector of the quadrature point weights
  std::vector<std::vector<double>> basis_func_map;   //store the calculated basis function values
  std::vector<std::vector<double>> basis_grad_map;   //store the calculated basis gradient values
    
  //Data structures
  SparsityPattern       sparsity_pattern; //Sparse matrix pattern
  SparseMatrix<double>  K;     //Global stiffness (sparse) matrix
  Vector<double>        D, F;      //Global vectors - Solution vector (D) and Global force vector (F)
  std::vector<double>   nodeLocation;  //Vector of the x-coordinate of nodes by global dof number
  std::map<unsigned int,double> boundary_values;  //Map of dirichlet boundary conditions
  double                L, E, Area;   //1D bar properties
  double                g1, g2, f_bar, f_hat, h_in;   //boundary conditions
  unsigned int          prob, basisFunctionOrder;   //case specifications

  //solution name array
  std::vector<std::string> nodal_solution_names;
  std::vector<DataComponentInterpretation::DataComponentInterpretation> nodal_data_component_interpretation;
};

