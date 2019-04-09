#ifndef TEST_H_
#define TEST_H_
//#include "FEM1.h"
#include <cassert>
using namespace dealii;

void test_basis_function();
void test_basis_gradient();
void test_klocal();
void test_flocal();
void test_compare_values();
void test_all();
void test_l2norm();
void test_h1norm();

// test basis functions of order 3
void test_basis_function() {

	FEM<1> problemObject(3,3);

	double a = problemObject.basis_function(0, -1.);
	assert(a == 1.);
	a = problemObject.basis_function(1, 1.);
	assert(a == 1.);
	a = problemObject.basis_function(0, 1.);
	assert(a == 0.);
	for(unsigned int i = 2; i < problemObject.basisFunctionOrder + 1; i++) {
		double loc = (i-1) * 2. / (problemObject.basisFunctionOrder) - 1.;
		assert(problemObject.basis_function(i, loc) == 1.);
	}
	std::cout << "basis_function tests PASSED!" << std::endl;
}

// test basis gradients of order 2 and 3
void test_basis_gradient() {

	FEM<1> problemObject(2,2);

	double a = problemObject.basis_gradient(0, -1.);
	assert(a == -1.5);
	a = problemObject.basis_gradient(1, 1.);
	assert(a == 1.5);
	a = problemObject.basis_gradient(0, 1.);
	assert(a == 0.5);

	FEM<1> problemObject2(3,2);

	a = problemObject2.basis_gradient(0, 0.);
	assert(abs(a - (-1./3.*1./3/((-4.)/3.*(-2.)/3.*(-2.)))) < 1e-8);
	a = problemObject2.basis_gradient(1, 0.);
	assert(abs(a - (-1./3.*1./3/((4.)/3.*(2.)/3.*(2.)))) < 1e-8);
	a = problemObject2.basis_gradient(0, 1.);
	assert(abs(a - ((-1./3.*1./3 + 1)/((-4.)/3.*(-2.)/3.*(-2.)))) < 1e-8);

	std::cout << "basis_gradient tests PASSED!" << std::endl;
}

// Test local K (4x4 matrix), with analytically precalculated exact values.
// Here we test local K because the global K is more chaotic if using the node
// index of deal.II compared to what we learned in class. The values of global
// matrix is not exactly located near diagonals, thus needing reorganization,
// and compare the elements seems to be difficult. However, since the steps from
// local K to global K is pretty simple, we can check the values of local K
// instead, and finally check the solved nodal values.
void test_klocal() {

	FEM<1> problemObject(3,1);
	problemObject.generate_mesh(10);
    problemObject.setup_system();
	// Following code is copied from FEM1.h
	unsigned int dofs_per_elem = 4;
    FullMatrix<double>        Klocal (dofs_per_elem, dofs_per_elem);
    for(unsigned int A=0; A<dofs_per_elem; A++){
      for(unsigned int B=0; B<dofs_per_elem; B++){
        for(unsigned int q=0; q<problemObject.quadRule; q++){
          Klocal[A][B] += problemObject.quad_weight[q] * problemObject.basis_gradient(A, problemObject.quad_points[q])
                                *  problemObject.basis_gradient(B, problemObject.quad_points[q]);
        } // for(q)
      } // for(B)
    } // for(A)

    std::vector<std::vector<double>>  Kexact (4);
    Kexact[0] = {1.85, -0.1675, -2.3625, 0.6750};
    Kexact[2] = {-2.3625, 0.6750, 5.4000, -3.7125};
    Kexact[3] = {0.6750, -2.3625, -3.7125, 5.40000};
    Kexact[1] = {-0.1625, 1.8500, 0.6750, -2.3625};
    for(unsigned int A=0; A<dofs_per_elem; A++){
      for(unsigned int B=0; B<dofs_per_elem; B++){
        assert(abs(Klocal[A][B] - Kexact[A][B]) < 1e-6);
      } // for(B)
    } // for(A)

    std::cout << "K_local tests PASSED!" << std::endl;
}

// Similarly, this function test a local F (3x1 vector) wih analytically precalcualted exact values.
// Since local F is not a constant, I choose the first element to test. Note that the implementation
// does not consider it as a special element, thus the test has generality.
void test_flocal() {
	 
	FEM<1> problemObject(2,1);
	problemObject.generate_mesh(10);
    problemObject.setup_system();
	// Following code is copied from FEM1.h
	unsigned int dofs_per_elem = 3;
	double x;
	double h_e = problemObject.L/10.;
	std::vector<unsigned int> local_dof_indices = {0,2,1};
	std::vector<double> nodeLocation (21);
	nodeLocation[0] = 0.; nodeLocation[1] = problemObject.L/10.; nodeLocation[2] = nodeLocation[1] /2.;
	std::vector<double> Flocal (4);
    for(unsigned int A=0; A<dofs_per_elem; A++){
      for(unsigned int q=0; q<problemObject.quadRule; q++){
        x = 0.;
        for(unsigned int B=0; B<dofs_per_elem; B++){
          x += nodeLocation[local_dof_indices[B]] * problemObject.basis_func_map[B][q];
        }
        Flocal[A] += problemObject.quad_weight[q] * problemObject.f_bar * x 
      			* problemObject.basis_func_map[A][q] * problemObject.Area * h_e / 2;
      }
    }

    std::vector<double> fexact = {50, 400/3, 1700/3};
    for(unsigned int A=0; A<dofs_per_elem; A++)
    	assert(abs(Flocal[A] - fexact[A]) < 1e-6);

    std::cout << "F_local tests PASSED!" << std::endl;

}

// Here we typically compare the values of solved nodal with exact ones, to ensure that we have
// correctly implement global K and global F. This is because that with different index system, 
// to directly compare global K or global F becomes expensive. Considering the non-complicated 
// transformation between global and local ones, I think compare the solved values is an efficiet
// way cause if the evaluated values are wrong, then there must be some problems in global K or
// global F.
void test_compare_values() {

	FEM<1> problemObject(2,1);
	problemObject.generate_mesh(10);
    problemObject.setup_system();
    problemObject.assemble_system();
    problemObject.solve();
    std::vector<double> Dexact = {0, 1.165e-4, 0.000232, 0.0003455, 0.000456,  0.0005625,
    							0.000664, 0.0007595, 0.000848, 0.0009285,  0.001};
    assert(abs(Dexact[0] - problemObject.D[0]) < 1e-6);
    for(unsigned int i = 1; i < 11; i++)
    	assert(abs(Dexact[i] - problemObject.D[2*i-1]) < 1e-6);

    std::cout << "solved D tests PASSED!" << std::endl;

}

// The following two are very naive tests that only check if the error goes below 1e-5.
void test_l2norm() {

	FEM<1> problemObject(2,4);
	problemObject.generate_mesh(10);
    problemObject.setup_system();
    problemObject.assemble_system();
    problemObject.solve();
    assert(problemObject.l2norm_of_error() < 1e-5);

	FEM<1> problemObject2(3,3);
	problemObject2.generate_mesh(100);
    problemObject2.setup_system();
    problemObject2.assemble_system();
    problemObject2.solve();
    assert(problemObject2.l2norm_of_error() < 1e-5);

    std::cout << "Naive L2 error tests PASSED!" << std::endl;

}

void test_h1norm() {

	FEM<1> problemObject(1,4);
	problemObject.generate_mesh(10);
    problemObject.setup_system();
    problemObject.assemble_system();
    problemObject.solve();
    assert(problemObject.h1norm_of_error() < 1e-5);

	FEM<1> problemObject2(3,2);
	problemObject2.generate_mesh(100);
    problemObject2.setup_system();
    problemObject2.assemble_system();
    problemObject2.solve();
    assert(problemObject2.h1norm_of_error() < 1e-5);

    std::cout << "Naive H1 error tests PASSED!" << std::endl;

}

// Similar to main()
void test_all() {

	test_basis_function();
	test_basis_gradient();

	test_klocal();
	test_flocal();

	test_compare_values();

	test_l2norm();
	test_h1norm();

	std::cout << "All tests passed!" << std::endl;
	std::cout << std::endl;

}

#endif
