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

//Define the order of the basis functions (Lagrange polynomials)
//and the order of the quadrature rule globally
const unsigned int order = 1;
const unsigned int quadRule = 5;
const double PI = std::acos(-1);
const bool mybasis = false;

template <int dim>
class FEM
{
 public:
  //Class functions
  FEM (unsigned int prob); // Class constructor 
  ~FEM(); //Class destructor

  //Function to calculate components of the elasticity tensor
  double C(unsigned int i,unsigned int j,unsigned int k,unsigned int l);

  //Solution steps
  void generate_mesh(std::vector<unsigned int> numberOfElements);
  void define_boundary_conds();
  void setup_system();
  void assemble_system();
        double torque();
  void solve();
  void output_results();

  //Class objects
  Triangulation<dim> triangulation; //mesh
  FESystem<dim>      fe;            //FE element
  DoFHandler<dim>    dof_handler;   //Connectivity matrices

  //NEW - deal.II quadrature
  QGauss<dim>   quadrature_formula;      //Quadrature
  QGauss<dim-1> face_quadrature_formula; //Face Quadrature

  //Data structures
  SparsityPattern      sparsity_pattern; //Sparse matrix pattern
  SparseMatrix<double> K;                //Global stiffness matrix - Sparse matrix - used in the solver
  Vector<double>       D, F;             //Global vectors - Solution vector (D) and Global force vector (F)

  Table<2,double>               dofLocation;     //Table of the coordinates of dofs by global dof number
  std::map<unsigned int,double> boundary_values; //Map of dirichlet boundary conditions
  unsigned int                  problem; //Problem 1 or 2
  double                        x_min, x_max, y_min, y_max, z_min, z_max;
  
  //solution name array
  std::vector<std::string> nodal_solution_names;
  std::vector<DataComponentInterpretation::DataComponentInterpretation> nodal_data_component_interpretation;
};


