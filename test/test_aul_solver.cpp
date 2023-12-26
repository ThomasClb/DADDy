/**
	test_aul_solver.cpp

	Purpose: Test of the implementation of the AULSolver class.

	@author Thomas Caleb

	@version 1.0 17/11/2023
*/

#include "pch.h"
#include "aul_solver.cpp"

using namespace DACE;
using namespace std;

TEST(TestAULSolver, Setters) {
	// Init DACE
	DA::init(2, 3);

	// Init
	vectordb x0(SIZE_VECTOR + 2, 1.0);
	vectordb xg(SIZE_VECTOR + 2, 1.0);
	vectorDA u(3, 1.0);
	Dynamics dynamics = get_tbp_SUN_low_thrust_dynamics();
	SpacecraftParameters spacecraft_p(dynamics.constants());
	SolverParameters solver_p;
	AULSolver aul_solver(solver_p, spacecraft_p, dynamics);
	vector<vectordb> list_lambda(solver_p.N(), vectordb(3));
	vector<vectordb> list_mu(solver_p.N(), vectordb(3));
	double ToF = 10.0;
	double homotopy_coefficient = 0.5;
	double huber_loss_coefficient = 1e-2;

	// Set
	aul_solver.set_homotopy_coefficient(homotopy_coefficient);
	aul_solver.set_ToF(ToF);
	aul_solver.set_huber_loss_coefficient(huber_loss_coefficient);
	DDPSolver ddp_solver = aul_solver.DDPsolver();
	SolverParameters solver_parameters = ddp_solver.solver_parameters();

	// Tests
	EXPECT_EQ(solver_parameters.homotopy_coefficient(), homotopy_coefficient);
	EXPECT_EQ(solver_parameters.ToF(), ToF);
	EXPECT_EQ(solver_parameters.huber_loss_coefficient(), huber_loss_coefficient);
}