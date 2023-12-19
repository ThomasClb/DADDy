/**
	low-thrust_cr3bp.cpp

	Purpose: Low-thrust Halo L2 to Halo L1 transfer execution script.

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
			x, u, spacecraft_parameters, solver_parameters);
		vectordb ineq_eval = dynamics.inequality_constraints_db()(
			x, u, spacecraft_parameters, solver_parameters);

		// Continuity constraints
		vectordb x_kp1_eval = dynamics.dynamic_db()(
			x, u, spacecraft_parameters, solver_parameters) - PNsolver.list_x()[i + 1];
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
		x, x_goal, spacecraft_parameters, solver_parameters);
	vectordb tineq_eval = dynamics.terminal_inequality_constraints_db()(
		x, x_goal, spacecraft_parameters, solver_parameters);

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
