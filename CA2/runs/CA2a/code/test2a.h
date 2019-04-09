#ifndef TEST2A_H_
#define TEST2A_H_
#include <cassert>
using namespace dealii;

void test_xi_at_node();
void test_xi_at_node_1d();
void test_basis_function();
void test_basis_gradient();
void test_klocal();
void test_flocal();
void test_all();
void test_l2norm();
void test_uexact();

// Similar to main(), it runs all test cases
void test_all() {

	test_xi_at_node();
	test_xi_at_node_1d();

	test_basis_function();
	test_basis_gradient();

	test_klocal();
	test_flocal();

	test_l2norm();
	test_uexact();

	std::cout << "All test passed!" << std::endl;
	std::cout << "====================================================" << std::endl;
	std::cout << std::endl;

}

// test xi_at_node function with order = 1, 2, 3
void test_xi_at_node() {

	FEM<2> problemObject(1,1);
	std::vector<double> basis1_0 = problemObject.xi_at_node(0);
	assert(basis1_0[0] == -1.);
	assert(basis1_0[1] == -1.);
	basis1_0 = problemObject.xi_at_node(2);
	assert(basis1_0[1] = 1.);

	FEM<2> problemObject2(2,1);
	std::vector<double> basis2 = problemObject2.xi_at_node(7);
	assert(basis2[0] == 0.);
	assert(basis2[1] == 1.);

	FEM<2> problemObject3(3,1);
	std::vector<double> basis3 = problemObject3.xi_at_node(3);
	assert(basis3[0] == 1.);
	assert(basis3[1] == 1.);

	basis3 = problemObject3.xi_at_node(14);
	assert(std::abs(basis3[0] -(-1./3.)) < 1e-4);
	assert(std::abs(basis3[1] -(1./3.)) < 1e-4);

	std::cout << "xi_at_node function tests PASSED!" << std::endl;
}

// test xi_at_node_1d function with order = 1, 2, 3
void test_xi_at_node_1d() {

	FEM<2> problemObject(1,1);
	double basis1 = problemObject.xi_at_node_1d(0);
	assert(basis1 == -1.);
	basis1 = problemObject.xi_at_node_1d(1);
	assert(basis1 == 1.);

	FEM<2> problemObject3(3,1);
	double basis3 = problemObject3.xi_at_node_1d(2);
	assert(std::abs(basis3 -(-1./3.)) < 1e-4);

	std::cout << "xi_at_node_1d function tests PASSED!" << std::endl;
}

// test basis functions with order = 1, 2, 3
void test_basis_function() {

	FEM<2> problemObject(1,1);

	double a = problemObject.basis_function(0, -1., -1.);
	assert(a == 1.);
	a = problemObject.basis_function(1, 1., -1.);
	assert(a == 1.);
	a = problemObject.basis_function(0, 1., 0.5);
	assert(a == 0.);

	// when order = 1, test whether the node function returns 0 at (1, 1)
	// if the node is not at (1, 1), and returns 1 when the node is at (node_x, node_y)
	for(unsigned int i = 0; i < (problemObject.basisFunctionOrder + 1) * (problemObject.basisFunctionOrder + 1) ; i++) {
		std::vector<double> basis1 = problemObject.xi_at_node(i);
		assert(problemObject.basis_function(i, basis1[0], basis1[1]) == 1.);
		if(basis1[0] != 1. || basis1[1] != 1.){
			assert(std::abs(problemObject.basis_function(i, 1., 1)) < 1.e-4 );
		}
	}

	FEM<2> problemObject2(2,1);

	a = problemObject2.basis_function(6, 0., -1.);
	assert(a == 1.);
	a = problemObject2.basis_function(6, 1., -1.);
	assert(a == 0.);
	a = problemObject2.basis_function(8, 0., 0.);
	assert(a == 1.);

	// when order = 2, test whether the node function returns 0 at (1, 1) 
	// if the node is not at (1, 1), and returns 1 when the node is at (node_x, node_y).
	for(unsigned int i = 0; i < (problemObject2.basisFunctionOrder + 1) * (problemObject2.basisFunctionOrder + 1) ; i++) {
		std::vector<double> basis1 = problemObject2.xi_at_node(i);
		assert(problemObject2.basis_function(i, basis1[0], basis1[1]) == 1.);
		if(basis1[0] != 1. || basis1[1] != 1.){
			assert(std::abs(problemObject2.basis_function(i, 1., 1)) < 1.e-4 );
		}
	}

	FEM<2> problemObject3(3,1);

	// when order = 3, test whether the node function returns 0 at (1/3, -1/3) 
	// if the node is not at (1/3, -1/3), and returns 1 when the node is at (node_x, node_y).
	for(unsigned int i = 0; i < (problemObject3.basisFunctionOrder + 1) * (problemObject3.basisFunctionOrder + 1) ; i++) {
		std::vector<double> basis1 = problemObject3.xi_at_node(i);
		assert(problemObject3.basis_function(i, basis1[0], basis1[1]) == 1.);
		if(std::abs(basis1[0] - 1./3.) > 1.e-4 || std::abs(basis1[1] + 1./3.) > 1.e-4) {
			//std::cout << abs(problemObject3.basis_function(i, 1./3., -1./3.)) << std::endl;
			//std::cout << basis1[0] << " " << basis1[1] << " " << i << std::endl;
			assert(std::abs(problemObject3.basis_function(i, 1./3., -1./3.)) < 1.e-4 );
		}
	}

	std::cout << "basis_function tests PASSED!" << std::endl;
}

// test basis gradients with order = 1, 2
void test_basis_gradient() {

	FEM<2> problemObject(1,1);

	std::vector<double> a = problemObject.basis_gradient(0, -1., -1.);
	assert(a[0] == -.5);
	assert(a[1] == -.5);
	a = problemObject.basis_gradient(1, 1., -1.);
	assert(a[0] == .5);
	assert(a[1] == -.5);

	FEM<2> problemObject2(2, 1);
	a = problemObject2.basis_gradient(8, 0., 1.);
	assert(a[0] == 0.);
	assert(a[1] == -2.);
	a = problemObject2.basis_gradient(5, 0.5, 0.5);
	//std::cout << a[0] << " " << a[1]<< std::endl;
	assert(std::abs(a[0] - 0.75) < 1.e-4);
	assert(std::abs(a[1] + 0.375) < 1.e-4);

	std::cout << "basis_gradient tests PASSED!" << std::endl;
}

// Test first-order  local K (4x4 matrix) with analytically precalculated exact values.
// Here we test local K because the global K is more messy if using the node
// index of deal.II compared to what we learned in class. The values of global
// matrix is not exactly located near diagonals, thus needing reorganization,
// and compare the elements seems to be difficult. However, since the steps from
// local K to global K is pretty simple, we can check the values of local K
// instead, and finally check the solved nodal values.
void test_klocal() {

	// Initialize Coding structure
	std::vector<unsigned int> num_of_elems(2);
    num_of_elems[0] = 10;
    num_of_elems[1] = 10;
    unsigned int dim = 2;

	FEM<2> problemObject(1, 1);
	problemObject.generate_mesh(num_of_elems);
    problemObject.setup_system();

    // Jacobian and its determinant
    FullMatrix<double> Jacobian(dim,dim);
    double detJ;
    FullMatrix<double> invJacob(dim,dim), kappa(dim,dim);
    
    // parameters used in problem 1 and 2
    kappa = 0.;
    kappa[0][0] = 385.;
    kappa[1][1] = 385.;

    std::vector<unsigned int> local_dof_indices = {1, 4, 3, 5}; // selected element
	unsigned int dofs_per_elem = 4;
    FullMatrix<double>        Klocal (dofs_per_elem, dofs_per_elem);
    Klocal = 0.;
    for(unsigned int q1 = 0; q1 < problemObject.quadRule; q1++) {
      for(unsigned int q2 = 0; q2 < problemObject.quadRule; q2++) {
        Jacobian = 0.;
        for(unsigned int i = 0; i < dim; i++) {
          for(unsigned int j = 0; j < dim; j++) {
            for(unsigned int A = 0; A < dofs_per_elem; A++) {
              Jacobian[i][j] += problemObject.nodeLocation[local_dof_indices[A]][i] 
              			* problemObject.basis_grad_map[A][q1][q2][j];
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
                    Klocal[A][B] += problemObject.quad_weight[q1] 
                    			* problemObject.quad_weight[q2] * detJ * kappa[I][J] 
                                * problemObject.basis_grad_map[A][q1][q2][i] * invJacob[i][I]
                                * problemObject.basis_grad_map[B][q1][q2][j] * invJacob[j][J];
                  }
                }
              }
            }
          }
        } // for(A)
      } // for(q2)
    } // for(q1)

    std::vector<std::vector<double>>  Kexact = {
    							{390.347, -318.16, 122.986, -195.174},
    							{-318.16, 390.347, -195.174, 122.986},
    							{122.986, -195.174, 390.347, -318.16},
    							{-195.174, 122.986, -318.16, 390.347}};

    for(unsigned int A=0; A<dofs_per_elem; A++){
      for(unsigned int B=0; B<dofs_per_elem; B++){
        assert(std::abs(Klocal[A][B] - Kexact[A][B]) < 1e-2);
      } // for(B)
    } // for(A)

	std::cout << "local K tests PASSED!" << std::endl;
}

// Similarly, this function test a first-order local F (4x1 vector) wih analytically
// precalcualted exact values. Since local F is a constant, we can choose any 
// element to test. 
void test_flocal() {

	// Initialize Coding structure
	std::vector<unsigned int> num_of_elems(2);
    num_of_elems[0] = 10;
    num_of_elems[1] = 10;
    unsigned int dim = 2;

	FEM<2> problemObject(1, 1);
	problemObject.generate_mesh(num_of_elems);
    problemObject.setup_system();

    // Jacobian and its determinant
    FullMatrix<double> Jacobian(dim,dim);
    double detJ;
    
    std::vector<unsigned int> local_dof_indices = {1, 4, 3, 5}; // selected element
	unsigned int dofs_per_elem = 4;
	double f = 333.; // a test f

	Vector<double>            Flocal (dofs_per_elem);
    Flocal = 0.; // initialize Flocal

    for(unsigned int q1 = 0; q1 < problemObject.quadRule; q1++) {
      for(unsigned int q2 = 0; q2 < problemObject.quadRule; q2++) {
        Jacobian = 0.;
        for(unsigned int i = 0; i < dim; i++) {
          for(unsigned int j = 0; j < dim; j++) {
            for(unsigned int A = 0; A < dofs_per_elem; A++) {
              Jacobian[i][j] += problemObject.nodeLocation[local_dof_indices[A]][i] 
              			* problemObject.basis_grad_map[A][q1][q2][j];
            }
          }
        }
        detJ = Jacobian.determinant();
        for(unsigned int A = 0; A < dofs_per_elem; A ++) {
          Flocal[A] += problemObject.quad_weight[q1] * problemObject.quad_weight[q2] 
          				* problemObject.basis_func_map[A][q1][q2] * f * detJ;
        }
      }
    }

    for(unsigned int A = 0; A < dofs_per_elem; A++)
    	assert(std::abs(Flocal[A] - 0.001998) < 1e-4);

	std::cout << "local F tests PASSED!" << std::endl;
}

// The following function is very naive that only check if the error goes below 1e-5.
// Note that this function only tests problem2.
void test_l2norm() {

	std::vector<unsigned int> num_of_elems(2);
    num_of_elems[0] = 10;
    num_of_elems[1] = 10;

	FEM<2> problemObject(1, 2); // For problem 2
    problemObject.generate_mesh(num_of_elems);
    problemObject.setup_system();
    problemObject.assemble_system();
    problemObject.solve();
    problemObject.output_results();

    assert(problemObject.l2norm_of_error() < 1.e-5);

    std::cout << "Naive l2 error norm tests PASSED!" << std::endl;
}

// This function check whether the uexact is implmented correctly.
void test_uexact(){
	FEM<2> problemObject(1, 1);
	assert(problemObject.uexact(0., 0.) == 100.);
	assert(problemObject.uexact(1., 0.) == 100. - 10000. / 4. / 385.);

	std::cout << "uexact tests PASSED!" << std::endl;
}



#endif