/**
	daddy.cpp

	Purpose: Implementation of the DADDy solver class.

	@author Thomas Caleb

	@version 1.0 08/03/2024
*/

#include "daddy.h"

using namespace DACE;
using namespace std;
using namespace std::chrono;

// Constructors
DADDy::DADDy() : 
		solver_parameters_(), dynamics_(), spacecraft_parameters_(),
		AULsolver_(), PNsolver_() {}
DADDy::DADDy(
	SolverParameters const& solver_parameters,
	SpacecraftParameters const& spacecraft_parameters,
	Dynamics const& dynamics) :
		solver_parameters_(solver_parameters), dynamics_(dynamics), spacecraft_parameters_(spacecraft_parameters),
		PNsolver_(), AULsolver_(solver_parameters, spacecraft_parameters, dynamics) {}

// Copy constructor
DADDy::DADDy(
	DADDy const& solver) : 
		solver_parameters_(solver.solver_parameters_), dynamics_(solver.dynamics_),
		spacecraft_parameters_(solver.spacecraft_parameters_), AULsolver_(solver.AULsolver_),
		PNsolver_(solver.PNsolver_) {}

// Destructors
DADDy::~DADDy() {}

// Getters
const AULSolver DADDy::AULsolver() const { return AULsolver_; }
const PNSolver DADDy::PNsolver() const { return PNsolver_; }
const vector<vectordb> DADDy::list_x() const { return PNsolver_.list_x(); }
const vector<vectordb> DADDy::list_u() const { return PNsolver_.list_u(); }
const double DADDy::cost() const { return PNsolver_.cost(); }
const double DADDy::runtime() const { return runtime_; }
const double DADDy::PN_runtime() const { return PN_runtime_; }
const double DADDy::AUL_runtime() const { return AUL_runtime_; }
const size_t DADDy::PN_n_iter() const { return PN_n_iter_; }
const size_t DADDy::AUL_n_iter() const { return AUL_n_iter_; }
const size_t DADDy::DDP_n_iter() const { return DDP_n_iter_; }



// Computes the constraints with propagated dynamics.
double DADDy::real_constraints(
	vectordb const& x_goal) {

	// Unpack
	unsigned int N = solver_parameters_.N();
	unsigned int Nx = solver_parameters_.Nx();
	unsigned int Nu = solver_parameters_.Nu();
	unsigned int Neq = solver_parameters_.Neq();
	unsigned int Nineq = solver_parameters_.Nineq();
	unsigned int Nteq = solver_parameters_.Nteq();
	unsigned int Ntineq = solver_parameters_.Ntineq();

	// Init
	double max_norm = 1e-15;

	// Update path constraints

	// Loop on all steps
	for (size_t i = 0; i < N; i++) {

		// Get DA x, u
		vectordb x = PNsolver_.list_x()[i];
		vectordb u = PNsolver_.list_u()[i];

		// Constraints evaluations
		vectordb eq_eval = dynamics_.equality_constraints_db()(
			x, u, spacecraft_parameters_, dynamics_.constants(), solver_parameters_);
		vectordb ineq_eval = dynamics_.inequality_constraints_db()(
			x, u, spacecraft_parameters_, dynamics_.constants(), solver_parameters_);

		// Continuity constraints
		vectordb x_kp1_eval = dynamics_.dynamic_db()(
			x, u, spacecraft_parameters_, dynamics_.constants(), solver_parameters_) - PNsolver_.list_x()[i + 1];
		x_kp1_eval[SIZE_VECTOR] /= spacecraft_parameters_.initial_mass(); // Normalize mass
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
	vectordb x = PNsolver_.list_x()[N];

	// Constraints evaluations
	vectordb teq_eval = dynamics_.terminal_equality_constraints_db()(
		x, x_goal, spacecraft_parameters_, dynamics_.constants(), solver_parameters_);
	vectordb tineq_eval = dynamics_.terminal_inequality_constraints_db()(
		x, x_goal, spacecraft_parameters_, dynamics_.constants(), solver_parameters_);

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

// Performs solving given a starting point,
// initial controls and a final state.
// Inspired from ALTRO (Julia).
// See: https://github.com/RoboticExplorationLab/Altro.jl
void DADDy::solve(
	vectordb const& x0,
	vector<vectordb> const& list_u_init,
	vectordb const& x_goal,
	bool const& fuel_optimal,
	bool const& pn_solving) {
	// Unpack
	unsigned int verbosity = solver_parameters_.verbosity();
	unsigned int Nx = solver_parameters_.Nx();
	unsigned int Nu = solver_parameters_.Nu();
	double AUL_tol = solver_parameters_.AUL_tol();
	double DDP_tol = solver_parameters_.DDP_tol();

	// Run DDP
	auto start = high_resolution_clock::now();
	vectordb homotopy_sequence = solver_parameters_.homotopy_coefficient_sequence();
	vectordb huber_loss_coefficient_sequence = solver_parameters_.huber_loss_coefficient_sequence();
	vector<vectordb> list_u_init_ = list_u_init;
	AUL_n_iter_ = 0;
	DDP_n_iter_ = 0;
	for (size_t i = 0; i < homotopy_sequence.size(); i++) {
		AULsolver_.set_homotopy_coefficient(homotopy_sequence[i]);
		AULsolver_.set_huber_loss_coefficient(huber_loss_coefficient_sequence[i]);
		if (i != 0) {
			for (size_t j=0; j<list_u_init.size(); j++) {
				list_u_init_[j] = AULsolver_.list_u()[j] + DDP_tol*DDP_tol; // Small perturbation
			}
		}
		AULsolver_.solve(x0, list_u_init_, x_goal);
		DDP_n_iter_ += AULsolver_.DDP_n_iter();
		AUL_n_iter_ += AULsolver_.AUL_n_iter();

		if (!fuel_optimal)
			break;
	}

	// PN test
	auto start_inter = high_resolution_clock::now();
	PNsolver_ = PNSolver(AULsolver_);
	if (pn_solving)
		PNsolver_.solve(x_goal);
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	auto duration_AUL = duration_cast<microseconds>(start_inter - start);
	auto duration_PN = duration_cast<microseconds>(stop - start_inter);
	PN_runtime_ = static_cast<double>(duration_PN.count()) / 1e6;
	AUL_runtime_ = static_cast<double>(duration_AUL.count()) / 1e6;
	runtime_ = static_cast<double>(duration.count()) / 1e6;
	PN_n_iter_ = PNsolver_.n_iter();

	// Output
	if (verbosity <= 1) {
		cout << endl;
		cout << "Optimised" << endl;
		cout << "	Total runtime : " + to_string(runtime_) + "s" << endl;
		cout << "	AUL solver runtime : " + to_string(AUL_runtime_) + "s" << endl;
		cout << "	PN solver runtime : " + to_string(PN_runtime_) + "s" << endl;
		cout << "	FINAL COST [-] : " << PNsolver_.cost() << endl;
		cout << "	FINAL ERROR [-] : " << real_constraints(x_goal) << endl;
	}
}
