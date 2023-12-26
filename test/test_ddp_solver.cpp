/**
	test_ddp_solver.cpp

	Purpose: Test of the implementation of the DDPSolver class.

	@author Thomas Caleb

	@version 1.0 17/11/2023
*/

#include "pch.h"
#include "ddp_solver.cpp"

using namespace DACE;
using namespace std;


TEST(TestDDPSolver, Setters) {
	// Init DACE
	DA::init(2, 3);

	// Init
	vectordb x0(SIZE_VECTOR + 2, 1.0);
	vectordb xg(SIZE_VECTOR + 2, 1.0);
	vectorDA u(3, 1.0);
	Dynamics dynamics = get_low_trust_2bp_SUN_dynamics();
	SpacecraftParameters spacecraft_p(dynamics.constants());
	SolverParameters solver_p;
	DDPSolver ddp_solver(solver_p, spacecraft_p, dynamics);
	vector<vectordb> list_lambda(solver_p.N(), vectordb(3));
	vector<vectordb> list_mu(solver_p.N(), vectordb(3));
	double ToF = 10.0;
	double homotopy_coefficient = 0.5;
	double huber_loss_coefficient = 1e-2;

	// Set
	ddp_solver.set_list_lambda(list_lambda);
	ddp_solver.set_list_mu(list_mu);
	ddp_solver.set_homotopy_coefficient(homotopy_coefficient);
	ddp_solver.set_ToF(ToF);
	ddp_solver.set_huber_loss_coefficient(huber_loss_coefficient);

	// Tests
	EXPECT_EQ(ddp_solver.solver_parameters().homotopy_coefficient(), homotopy_coefficient);
	EXPECT_EQ(ddp_solver.solver_parameters().ToF(), ToF);
	EXPECT_EQ(ddp_solver.solver_parameters().huber_loss_coefficient(), huber_loss_coefficient);
	for (size_t i = 0; i < solver_p.N(); i++) {
		vectordb lambda_i = ddp_solver.solver_parameters().list_lambda()[i];
		vectordb lambda_target_i = list_lambda[i];
		vectordb mu_i = ddp_solver.solver_parameters().list_mu()[i];
		vectordb mu_target_i = list_mu[i];
		for (size_t j = 0; j < mu_target_i.size(); j++) {
			EXPECT_EQ(lambda_i[j], lambda_target_i[j]);
			EXPECT_EQ(mu_i[j], mu_target_i[j]);
		}
	}
}
