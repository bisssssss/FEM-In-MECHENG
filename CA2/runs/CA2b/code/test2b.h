#ifndef TEST2B_H_
#define TEST2B_H_
#include <cassert>
using namespace dealii;

void test_xi_at_node();
void test_basis_function();
void test_basis_gradient();
void test_klocal();
void test_flocal();
void test_all();

void test_all() {

	test_xi_at_node();
	// This function has already been tested in "test2a.h"
	// test_xi_at_node_1d();
	
	test_basis_function();
	test_basis_gradient();

	test_klocal();
	test_flocal();

	std::cout << "All test passed!" << std::endl;
	std::cout << "====================================================" << std::endl;
	std::cout << std::endl;

}


// test xi_at_node function with order = 1, 2, 3
void test_xi_at_node() {

	FEM<3> problemObject(1);
	std::vector<double> basis1_0 = problemObject.xi_at_node(0);
	assert(basis1_0[0] == -1.);
	assert(basis1_0[1] == -1.);
	assert(basis1_0[2] == -1.);
	basis1_0 = problemObject.xi_at_node(2);
	assert(basis1_0[1] = 1.);
	basis1_0 = problemObject.xi_at_node(4);
	assert(basis1_0[2] = 1.);

	FEM<3> problemObject2(2);
	std::vector<double> basis2 = problemObject2.xi_at_node(17);
	assert(basis2[0] == 1.);
	assert(basis2[1] == -1.);
	assert(basis2[2] == 0.);
	basis2 = problemObject2.xi_at_node(21);
	assert(basis2[0] == 1.);
	assert(basis2[1] == 0.);
	assert(basis2[2] == 0.);
	basis2 = problemObject2.xi_at_node(26);
	assert(basis2[0] == 0.);
	assert(basis2[1] == 0.);
	assert(basis2[2] == 0.);

	FEM<3> problemObject3(3);
	std::vector<double> basis3 = problemObject3.xi_at_node(33);
	assert(basis3[0] == -1.);
	assert(std::abs(basis3[1] -(1./3.)) < 1e-4);
	assert(std::abs(basis3[2] -(-1./3.)) < 1e-4);

	basis3 = problemObject3.xi_at_node(60);
	assert(std::abs(basis3[0] -(-1./3.)) < 1e-4);
	assert(std::abs(basis3[1] -(-1./3.)) < 1e-4);
	assert(std::abs(basis3[2] -(1./3.)) < 1e-4);

	std::cout << "xi_at_node function tests PASSED!" << std::endl;
}

// test basis functions with order = 1, 2, 3
void test_basis_function() {
	
	FEM<3> problemObject(1);

	double a = problemObject.basis_function(0, -1., -1., -1.);
	assert(a == 1.);
	a = problemObject.basis_function(6, -1., 1., 1.);
	assert(a == 1.);
	a = problemObject.basis_function(0, -1., -1., 1.);
	assert(a == 0.);

	// when order = 1, test whether the node function returns 0 at (1, 1, 1)
	// if the node is not at (1, 1, 1), and returns 1 when the node is at (node_x, node_y, node_z)
	for(unsigned int i = 0; i < (problemObject.basisFunctionOrder + 1) * (problemObject.basisFunctionOrder + 1)
	                             * (problemObject.basisFunctionOrder + 1); i++) {
		std::vector<double> basis1 = problemObject.xi_at_node(i);
		assert(problemObject.basis_function(i, basis1[0], basis1[1], basis1[2]) == 1.);
		if(basis1[0] != 1. || basis1[1] != 1. || basis1[2] != 1.){
			assert(std::abs(problemObject.basis_function(i, 1., 1., 1.)) < 1.e-4 );
		}
	}

	FEM<3> problemObject2(2);

	// when order = 2, test whether the node function returns 0 at (1, 1, 1) 
	// if the node is not at (1, 1, 1), and returns 1 when the node is at (node_x, node_y, node_z).
	for(unsigned int i = 0; i < (problemObject.basisFunctionOrder + 1) * (problemObject.basisFunctionOrder + 1)
	                             * (problemObject.basisFunctionOrder + 1); i++) {
		std::vector<double> basis1 = problemObject.xi_at_node(i);
		assert(problemObject.basis_function(i, basis1[0], basis1[1], basis1[2]) == 1.);
		if(basis1[0] != 1. || basis1[1] != 1. || basis1[2] != 1.){
			assert(std::abs(problemObject.basis_function(i, 1., 1., 1.)) < 1.e-4 );
		}
	}

	FEM<3> problemObject3(3);

	// when order = 3, test whether the node function returns 0 at (1/3, -1/3, 1) 
	// if the node is not at (1/3, -1/3, 1), and returns 1 when the node is at (node_x, node_y, node_z).
	for(unsigned int i = 0; i < (problemObject3.basisFunctionOrder + 1) * (problemObject3.basisFunctionOrder + 1) ; i++) {
		std::vector<double> basis1 = problemObject3.xi_at_node(i);
		assert(problemObject3.basis_function(i, basis1[0], basis1[1], basis1[2]) == 1.);
		if(std::abs(basis1[0] - 1./3.) > 1.e-4 || std::abs(basis1[1] + 1./3.) > 1.e-4
					|| std::abs(basis1[2] - 1.) > 1.e-4) {
			assert(std::abs(problemObject3.basis_function(i, 1./3., -1./3., 1.)) < 1.e-4 );
		}
	}

	std::cout << "basis_function tests PASSED!" << std::endl;
}

// test basis gradients with order = 1, 2
void test_basis_gradient() {
	
	FEM<3> problemObject(1);

	std::vector<double> a = problemObject.basis_gradient(0, -1., -1., -1.);
	assert(a[0] == -.5);
	assert(a[1] == -.5);
	assert(a[2] == -.5);
	a = problemObject.basis_gradient(6, 0., 0., 0.5);
	assert(a[0] == -0.1875);
	assert(a[1] == 0.1875);
	assert(a[2] == .125);
	
	FEM<3> problemObject2(2);

	a = problemObject2.basis_gradient(26, 0.5, 0., -0.5);
	assert(a[0] == -0.75);
	assert(a[1] == 0.);
	assert(a[2] == 0.75);

	a = problemObject2.basis_gradient(18, 0.5, 0.5, 0.5);
	assert(std::abs(a[0] - 0.) < 1.e-4);
	assert(std::abs(a[1] + 0.09375) < 1.e-4);
	assert(std::abs(a[2] - 0.046875) < 1.e-4);
	
	std::cout << "basis_gradient tests PASSED!" << std::endl;
}

// Test first-order local K (4x4 matrix) with analytically precalculated exact values.
// Here we test local K because the global K is more messy if using the node
// index of deal.II compared to what we learned in class. The values of global
// matrix is not exactly located near diagonals, thus needing reorganization,
// and compare the elements seems to be difficult. However, since the steps from
// local K to global K is pretty simple, we can check the values of local K
// instead, and finally check the solved nodal values.
void test_klocal() {

	// Initialize Coding structure
	std::vector<unsigned int> num_of_elems(3);
    num_of_elems[0] = 10;
    num_of_elems[1] = 10;
    num_of_elems[2] = 10;
    unsigned int dim = 3;

	FEM<3> problemObject(1);
	problemObject.generate_mesh(num_of_elems);
    problemObject.setup_system();

    // Jacobian and its determinant and inverse
    FullMatrix<double> Jacobian(dim,dim), invJacob(dim,dim);
    double detJ;

    // parameters used in problem 1 and 2
    FullMatrix<double> kappa(dim,dim);
    kappa = 0.;
    kappa[0][0] = 385.;
    kappa[1][1] = 385.;
    kappa[2][2] = 385.;

    std::vector<unsigned int> local_dof_indices = {0,1,2,3,4,5,6,7}; // selected element
	unsigned int dofs_per_elem = 8;
    FullMatrix<double>        Klocal (dofs_per_elem, dofs_per_elem);
    Klocal = 0.;
    for(unsigned int q1 = 0; q1 < problemObject.quadRule; q1++) {
      for(unsigned int q2 = 0; q2 < problemObject.quadRule; q2++) {
        for(unsigned int q3 = 0; q3 < problemObject.quadRule; q3++) {
          Jacobian = 0.;
          for(unsigned int i = 0; i < dim; i++) {
            for(unsigned int j = 0; j < dim; j++) {
              for(unsigned int A = 0; A < dofs_per_elem; A++) {
                Jacobian[i][j] += problemObject.nodeLocation[local_dof_indices[A]][i] 
                							* problemObject.basis_grad_map[A][q1][q2][q3][j];
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
                        Klocal[A][B] += problemObject.quad_weight[q1] * problemObject.quad_weight[q2] 
                        			* problemObject.quad_weight[q3] * detJ * kappa[I][J]
                                  * problemObject.basis_grad_map[A][q1][q2][q3][i] * invJacob[i][I]
                                  * problemObject.basis_grad_map[B][q1][q2][q3][j] * invJacob[j][J];
                    }
                  }
                }
              }
            }
          } // for(A)
        } // for(q3)
      } // for(q2)
    } // for(q1)

    // compare the klocal with values that are pre-calculated by human.
    // Note that the matrix is 8x8, thus compare every values is expensive
    // hence I randomly select three nodes for comparision.
    assert(std::abs(Klocal[0][0] - 0.898333) < 1e-4);
    assert(std::abs(Klocal[4][1] + 0.417083) < 1e-4);
    assert(std::abs(Klocal[7][6] - 0.1925) < 1e-4);

	std::cout << "local K tests PASSED!" << std::endl;
	
}

// Similarly, this function test a first-order local F (4x1 vector) wih analytically
// precalcualted exact values. Since local F is a constant, we can choose any 
// element to test. 
void test_flocal() {

	// Initialize Coding structure
	std::vector<unsigned int> num_of_elems(3);
    num_of_elems[0] = 10;
    num_of_elems[1] = 10;
    num_of_elems[2] = 10;
    unsigned int dim = 3;

	FEM<3> problemObject(1);
	problemObject.generate_mesh(num_of_elems);
    problemObject.setup_system();

    // Jacobian and its determinant
    FullMatrix<double> Jacobian(dim,dim);
    double detJ;
    
    std::vector<unsigned int> local_dof_indices = {0,1,2,3,4,5,6,7}; // selected element
	unsigned int dofs_per_elem = 8;
	double f = -333.; // a test f

	Vector<double>            Flocal (dofs_per_elem);
    Flocal = 0.; // initialize Flocal

    for(unsigned int q1 = 0; q1 < problemObject.quadRule; q1++) {
      for(unsigned int q2 = 0; q2 < problemObject.quadRule; q2++) {
        for(unsigned int q3 = 0; q3 < problemObject.quadRule; q3++) {
          Jacobian = 0.;
          for(unsigned int i = 0; i < dim; i++) {
            for(unsigned int j = 0; j < dim; j++) {
              for(unsigned int A = 0; A < dofs_per_elem; A++) {
                Jacobian[i][j] += problemObject.nodeLocation[local_dof_indices[A]][i] 
                						* problemObject.basis_grad_map[A][q1][q2][q3][j];
              }
            }
          }
          detJ = Jacobian.determinant();
          for(unsigned int A = 0; A < dofs_per_elem; A ++) {
            Flocal[A] += problemObject.quad_weight[q1] * problemObject.quad_weight[q2] 
            						* problemObject.quad_weight[q3]
                                       * problemObject.basis_func_map[A][q1][q2][q3] * f * detJ;
          }
        }
      }
    }

    for(unsigned int A = 0; A < dofs_per_elem; A++)
    	assert(std::abs(Flocal[A] + 2.664e-06) < 1e-8);
    
	std::cout << "local F tests PASSED!" << std::endl;

}

#endif