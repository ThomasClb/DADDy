/**
	test_cases.cpp

	Purpose: Executes a test case.

	@author Thomas Caleb

	@version 1.0 07/12/2023
*/

#include "test_cases.h"

using namespace DACE;
using namespace std::chrono;
using namespace std;

double real_constraints(
	vectordb const& x_goal, PNSolver const& PNsolver) {

	// Unpack
	DDPSolver ddp_solver = PNsolver.DDPsolver();
	SolverParameters solver_parameters = ddp_solver.solver_parameters();
	SpacecraftParameters spacecraft_parameters = ddp_solver.spacecraft_parameters();
	Dynamics dynamics = ddp_solver.dynamics();
	unsigned int N = solver_parameters.N();
	unsigned int Nx = solver_parameters.Nx();
	unsigned int Nu = solver_parameters.Nu();
	unsigned int Neq = solver_parameters.Neq();
	unsigned int Nineq = solver_parameters.Nineq();
	unsigned int Nteq = solver_parameters.Nteq();
	unsigned int Ntineq = solver_parameters.Ntineq();

	// Init
	double max_norm = 1e-15;

	// Update path constraints

	// Loop on all steps
	for (size_t i = 0; i < N; i++) {

		// Get DA x, u
		vectordb x = PNsolver.list_x()[i];
		vectordb u = PNsolver.list_u()[i];

		// Constraints evaluations
		vectordb eq_eval = dynamics.equality_constraints_db()(
			x, u, spacecraft_parameters, dynamics.constants(), solver_parameters);
		vectordb ineq_eval = dynamics.inequality_constraints_db()(
			x, u, spacecraft_parameters, dynamics.constants(), solver_parameters);

		// Continuity constraints
		vectordb x_kp1_eval = dynamics.dynamic_db()(
			x, u, spacecraft_parameters, dynamics.constants(), solver_parameters) - PNsolver.list_x()[i + 1];
		x_kp1_eval[SIZE_VECTOR] /= spacecraft_parameters.initial_mass(); // Normalize mass
		eq_eval = eq_eval.concat(x_kp1_eval);


		// Assign
		vectordb eq = eq_eval;
		vectordb ineq = ineq_eval;
		for (size_t j = 0; j < Neq + Nx; j++) {
			if (max_norm < abs(eq[j]))
				max_norm = abs(eq[j]);
		}
		for (size_t j = 0; j < Nineq; j++) {
			if (max_norm < ineq[j])
				max_norm = ineq[j];
		}
	}

	// Update terminal constraints

	// Get DA x, u
	vectordb x = PNsolver.list_x()[N];

	// Constraints evaluations
	vectordb teq_eval = dynamics.terminal_equality_constraints_db()(
		x, x_goal, spacecraft_parameters, dynamics.constants(), solver_parameters);
	vectordb tineq_eval = dynamics.terminal_inequality_constraints_db()(
		x, x_goal, spacecraft_parameters, dynamics.constants(), solver_parameters);

	// Assign
	vectordb teq = teq_eval;
	vectordb tineq = tineq_eval;
	for (size_t j = 0; j < Nteq; j++) {
		if (max_norm < abs(teq[j]))
			max_norm = abs(teq[j]);
	}
	for (size_t j = 0; j < Ntineq; j++) {
		if (max_norm < tineq[j])
			max_norm = tineq[j];
	}

	return max_norm;
}

void run_test_cases(int argc, char** argv) {

	// Input check
	if (argc < 2) {
		cout << "Wrong number of arguments." << endl;
		cout << "Requested number : 1" << endl;
		cout << "0 - Test case number from 0 to 8." << endl;
		return;
	}

	// Unpack inputs
	unsigned int test_case = atoi(argv[1]);

	// Excute correct test case

	// Double integrator
	if (test_case == 0) {}

	// TBP
	
	// SUN centered

	// Earth-mars transfer
	else if (test_case == 1)
		tbp_SUN_lt_earth_to_mars(argc, argv);

	// Earth centered

	// GTO to GEO
	else if (test_case == 2) {
		tbp_EARTH_lt_leo_to_geo(argc, argv);
	}

	// LEO to LEO
	else if (test_case == 3) {
		tbp_EARTH_lt_leo_to_leo(argc, argv);
	}

	// MEO to MEO
	else if (test_case == 4) {
		tbp_EARTH_lt_meo_to_meo(argc, argv);
	}

	// CR3BP

	// Earth-Moon

	// Halo L2 to Halo L1
	else if (test_case == 5)
		cr3bp_EARTH_MOON_lt_haloL2_to_haloL1(argc, argv);

	// Halo L2 (Gateway) to DRO
	else if (test_case == 6) {
		cr3bp_EARTH_MOON_lt_nrho_to_dro(argc, argv);
	}

	// Lyapunov L2 to Lyapunov L1
	else if (test_case == 7) {
		cr3bp_EARTH_MOON_lt_lyapunovL1_to_lyapunovL2(argc, argv);
	}

	// Manyrev DRO to DRO
	else if (test_case == 8) {
		cr3bp_EARTH_MOON_lt_dro_to_dro(argc, argv);
	}

	// Manyrev GSO to NRHO
	else if (test_case == 9) {}
}