/**
	test_pn_solver.cpp

	Purpose: Test of the implementation of the PNSolver class.

	@author Thomas Caleb

	@version 1.0 17/11/2023
*/

#include "pch.h"
#include "pn_solver.cpp"

using namespace DACE;
using namespace std;

TEST(TestPNSolver, Setters) {
	// Init DACE
	DA::init(2, 3);

	// Init
	vectordb x0(SIZE_VECTOR + 2, 1.0);
	vectordb xg(SIZE_VECTOR + 2, 1.0);
	vectorDA u(3, 1.0);
	SpacecraftParameters spacecraft_p;
	SolverParameters solver_p;
	Dynamics dynamics = get_tbp_SUN_lt_dynamics();
	AULSolver aul_solver(solver_p, spacecraft_p, dynamics);
	PNSolver pn_solver(aul_solver);
	vector<vectordb> list_x(solver_p.N() + 1, vectordb(3));
	vector<vectordb> list_u(solver_p.N(), vectordb(3));
	double ToF = 10.0;
	double homotopy_coefficient = 0.5;
	double huber_loss_coefficient = 1e-2;

	// Set
	pn_solver.set_list_x(list_x);
	pn_solver.set_list_u(list_u);

	// Tests
	for (size_t i = 0; i < solver_p.N() + 1; i++) {
		EXPECT_EQ(pn_solver.list_x()[i], list_x[i]);
	}
	for (size_t i = 0; i < solver_p.N(); i++) {
		EXPECT_EQ(pn_solver.list_u()[i], list_u[i]);
	}
}
