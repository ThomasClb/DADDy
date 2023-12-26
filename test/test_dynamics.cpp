/**
	test_dynamics.cpp

	Purpose: Test of the implementation of the dynamics
	methods.

	@author Thomas Caleb

	@version 1.0 16/11/2023
*/

#include "pch.h"
#include "dynamics.cpp"

using namespace DACE;
using namespace std;

// Acceleration

// Sun two-body problem
TEST(TestDynamics, Acc2BSUN) {
	// Init DACE
	size_t nb_variables = 3;
	DA::init(2, nb_variables);

	// Init
	Constants constants;
	vectorDA x0(SIZE_VECTOR + 2, 1.0 + 2.0 * DA(2));
	vectorDA u(3, 1.0/ constants.thrustu() + DA(1));
	double t = 0; double dt = 1.0;
	SpacecraftParameters sp(constants);
	vectorDA xf = acceleration_tbp_SUN_low_thrust(x0, u, t, sp, constants);

	// Tests
	EXPECT_EQ(xf.size(), x0.size());

	// For doubles
	vectordb xf_db = acceleration_tbp_SUN_low_thrust(x0.cons(), u.cons(), t, sp, constants);

	// Tests
	EXPECT_EQ(xf_db.size(), x0.cons().size());
	EXPECT_EQ(xf_db, xf.cons());
}

// Dynamics class

// Test dyanmics constructor + Getters
TEST(TestDynamics, DynamicsConsGet) {
	// Init DACE
	size_t nb_variables = 3;
	DA::init(2, nb_variables);

	// Init
	Dynamics dynamics = get_tbp_SUN_low_thrust_dynamics();
	Constants constants = dynamics.constants();
	vectordb xg(SIZE_VECTOR + 2, 1.0);
	vectorDA x0(SIZE_VECTOR + 2, 1.0 + 2.0 * DA(2));
	vectorDA u(3, 1.0 / constants.thrustu() + DA(1));
	double t = 0; double dt = 1.0;
	double ToF = 1.0;
	SpacecraftParameters spacecraft_p(constants);
	SolverParameters solver_p;
	solver_p.set_ToF(ToF);
	vectorDA xf = dynamics.dynamic()(x0, u, spacecraft_p, constants, solver_p);
	DA ctg = dynamics.cost_to_go()(x0, u, spacecraft_p, constants, solver_p);
	vectorDA eq = dynamics.equality_constraints()(x0, u, spacecraft_p, constants, solver_p);
	vectorDA ineq = dynamics.inequality_constraints()(x0, u, spacecraft_p, constants, solver_p);
	DA tc = dynamics.terminal_cost()(x0, xg, spacecraft_p, constants, solver_p);
	vectorDA teq = dynamics.terminal_equality_constraints()(
		x0, xg, spacecraft_p, constants, solver_p);
	vectorDA tineq = dynamics.terminal_inequality_constraints()(
		x0, xg, spacecraft_p, constants, solver_p);
	vectordb xf_db = dynamics.dynamic_db()(x0.cons(), u.cons(), spacecraft_p, constants, solver_p);
	double ctg_db = dynamics.cost_to_go_db()(x0.cons(), u.cons(), spacecraft_p, constants, solver_p);
	vectordb eq_db = dynamics.equality_constraints_db()(x0.cons(), u.cons(), spacecraft_p, constants, solver_p);
	vectordb ineq_db = dynamics.inequality_constraints_db()(x0.cons(), u.cons(), spacecraft_p, constants, solver_p);
	double tc_db = dynamics.terminal_cost_db()(x0.cons(), xg, spacecraft_p, constants, solver_p);
	vectordb teq_db = dynamics.terminal_equality_constraints_db()(
		x0.cons(), xg, spacecraft_p, constants, solver_p);
	vectordb tineq_db = dynamics.terminal_inequality_constraints_db()(
		x0.cons(), xg, spacecraft_p, constants, solver_p);

	// Tests
	EXPECT_EQ(xf.size(), x0.size());
	EXPECT_EQ(xf.size(), solver_p.Nx());
	EXPECT_EQ(eq.size(), solver_p.Neq());
	EXPECT_EQ(ineq.size(), solver_p.Nineq());
	EXPECT_EQ(teq.size(), solver_p.Nteq());
	EXPECT_EQ(tineq.size(), solver_p.Ntineq());

	// db vs DA
	EXPECT_EQ(xf_db, xf.cons());
	EXPECT_EQ(ctg_db, ctg.cons());
	EXPECT_EQ(eq_db, eq.cons());
	EXPECT_EQ(ineq_db, ineq.cons());
	EXPECT_EQ(tc_db, tc.cons());
	EXPECT_EQ(teq_db, teq.cons());
	EXPECT_EQ(tineq_db, tineq.cons());
}

// Test dyanmics copy constructor + Getters
TEST(TestDynamics, CopyDynamicsConsGet) {
	// Init DACE
	size_t nb_variables = 3;
	DA::init(2, nb_variables);

	// Init
	Dynamics dynamics = get_tbp_SUN_low_thrust_dynamics();
	Constants constants = dynamics.constants();
	vectordb xg(SIZE_VECTOR + 2, 1.0);
	vectorDA x0(SIZE_VECTOR + 2, 1.0 + 2.0 * DA(2));
	vectorDA u(3, 1.0 / constants.thrustu() + DA(1));
	double t = 0; double dt = 1.0;
	double ToF = 1.0;
	SpacecraftParameters spacecraft_p(constants);
	SolverParameters solver_p;
	solver_p.set_ToF(ToF);
	Dynamics copy_dynamics = dynamics;
	vectorDA xf = copy_dynamics.dynamic()(x0, u, spacecraft_p, constants, solver_p);
	DA ctg = copy_dynamics.cost_to_go()(x0, u, spacecraft_p, constants, solver_p);
	vectorDA eq = copy_dynamics.equality_constraints()(x0, u, spacecraft_p, constants, solver_p);
	vectorDA ineq = copy_dynamics.inequality_constraints()(x0, u, spacecraft_p, constants, solver_p);
	DA tc = copy_dynamics.terminal_cost()(x0, xg, spacecraft_p, constants, solver_p);
	vectorDA teq = copy_dynamics.terminal_equality_constraints()(
		x0, xg, spacecraft_p, constants, solver_p);
	vectorDA tineq = copy_dynamics.terminal_inequality_constraints()(
		x0, xg, spacecraft_p, constants, solver_p);
	vectordb xf_db = copy_dynamics.dynamic_db()(x0.cons(), u.cons(), spacecraft_p, constants, solver_p);
	double ctg_db = copy_dynamics.cost_to_go_db()(x0.cons(), u.cons(), spacecraft_p, constants, solver_p);
	vectordb eq_db = copy_dynamics.equality_constraints_db()(x0.cons(), u.cons(), spacecraft_p, constants, solver_p);
	vectordb ineq_db = copy_dynamics.inequality_constraints_db()(x0.cons(), u.cons(), spacecraft_p, constants, solver_p);
	double tc_db = copy_dynamics.terminal_cost_db()(x0.cons(), xg, spacecraft_p, constants, solver_p);
	vectordb teq_db = dynamics.terminal_equality_constraints_db()(
		x0.cons(), xg, spacecraft_p, constants, solver_p);
	vectordb tineq_db = copy_dynamics.terminal_inequality_constraints_db()(
		x0.cons(), xg, spacecraft_p, constants, solver_p);

	// Tests
	EXPECT_EQ(xf.size(), x0.size());
	EXPECT_EQ(xf.size(), solver_p.Nx());
	EXPECT_EQ(eq.size(), solver_p.Neq());
	EXPECT_EQ(ineq.size(), solver_p.Nineq());
	EXPECT_EQ(teq.size(), solver_p.Nteq());
	EXPECT_EQ(tineq.size(), solver_p.Ntineq());

	// db vs DA
	EXPECT_EQ(xf_db, xf.cons());
	EXPECT_EQ(ctg_db, ctg.cons());
	EXPECT_EQ(eq_db, eq.cons());
	EXPECT_EQ(ineq_db, ineq.cons());
	EXPECT_EQ(tc_db, tc.cons());
	EXPECT_EQ(teq_db, teq.cons());
	EXPECT_EQ(tineq_db, tineq.cons());
}
