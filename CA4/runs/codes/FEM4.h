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
using namespace dealii;

//Define the order of the basis functions (Lagrange polynomials)
//and the order of the quadrature rule globally
const unsigned int order = 1;
const unsigned int quadRule = 2;
const bool iflumped = false;

template <int dim>
class FEM
{
 public:
  //Class functions
  FEM (double Alpha, double delta_tin); // Class constructor 
  ~FEM(); //Class destructor

  //Solution steps
  void generate_mesh(std::vector<unsigned int> numberOfElements);
  void define_boundary_conds();
  void setup_system();
  void assemble_system();

  void solve_steady();
  void apply_initial_conditions();
  void solve_trans();
  void output_steady_results();
  void output_trans_results(unsigned int index);

  //Calculate the l2norm of the difference between the steady state and transient solution
  double l2norm();

  //Class objects
  Triangulation<dim> triangulation; //mesh
  FESystem<dim>      fe;            //FE element
  DoFHandler<dim>    dof_handler;   // Connectivity matrices

  QGauss<dim>  	     quadrature_formula; //Quadrature

  //Data structures
  SparsityPattern      sparsity_pattern;                   //Sparse matrix pattern
  SparseMatrix<double> M, K, system_matrix;                //Global stiffness matrix - Sparse matrix - used in the solver
  Vector<double>       D_steady, D_trans, V_trans, F, RHS; //Global vectors - Solution vector (D) and Global force vector (F)

  Table<2,double>	        nodeLocation;	      //Table of the coordinates of nodes by global dof number
  std::map<unsigned int,double> boundary_values_of_D; //Map of dirichlet boundary conditions for the temperature
  std::map<unsigned int,double> boundary_values_of_V; //Map of dirichlet boundary conditions for the time derivative of temperature

  std::vector<double> l2norm_results; //A vector to store the l2norms calculated in the time loop in solve_trans()
  double	      alpha; 	      //Specifies the Euler method, 0 <= alpha <= 1
  double x_min, x_max, y_min, y_max, z_min, z_max; // boundary domain
  double delta_t;             // time step



  //solution name array
  std::vector<std::string> nodal_solution_names;
  std::vector<DataComponentInterpretation::DataComponentInterpretation> nodal_data_component_interpretation;
};

