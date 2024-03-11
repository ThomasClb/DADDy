/**
	aul_solver.cpp

	Purpose: Implementation of the AULSolver class.

	@author Thomas Caleb

	@version 1.0 17/11/2023
*/

#include "aul_solver.h"

using namespace DACE;
using namespace std;
using namespace std::chrono;

// Constructors
AULSolver::AULSolver() : DDPsolver_(DDPSolver()),
	list_x_(vector<vectordb>(0)), list_u_(vector<vectordb>(0)), cost_(0),
	list_eq_(vector<vectordb>(0)), list_ineq_(vector<vectordb>(0)),
	list_lambda_(vector<vectordb>(0)), list_mu_(vector<vectordb>(0)),
	AUL_n_iter_(0), DDP_n_iter_(0),
	list_x_mem_(vector<vector<vectordb>>(0)), list_u_mem_(vector<vector<vectordb>>(0)) {}

AULSolver::AULSolver(
	SolverParameters const& solver_parameters,
	SpacecraftParameters const& spacecraft_parameters,
	Dynamics const& dynamics) : DDPsolver_(
		DDPSolver(solver_parameters, spacecraft_parameters, dynamics)),
	list_x_(vector<vectordb>(0)), list_u_(vector<vectordb>(0)), cost_(0),
	list_eq_(vector<vectordb>(0)), list_ineq_(vector<vectordb>(0)),
	list_lambda_(solver_parameters.list_lambda()), list_mu_(solver_parameters.list_mu()),
	AUL_n_iter_(0), DDP_n_iter_(0),
	list_x_mem_(vector<vector<vectordb>>(0)), list_u_mem_(vector<vector<vectordb>>(0)) {}

// Copy constructor
AULSolver::AULSolver(
	AULSolver const& solver) : DDPsolver_(solver.DDPsolver_),
	list_x_(solver.list_x_), list_u_(solver.list_u_), cost_(solver.cost_),
	list_eq_(solver.list_eq_), list_ineq_(solver.list_ineq_),
	list_lambda_(solver.list_lambda_), list_mu_(solver.list_mu_),
	AUL_n_iter_(solver.AUL_n_iter_), DDP_n_iter_(solver.DDP_n_iter_),
	list_x_mem_(solver.list_x_mem_), list_u_mem_(solver.list_u_mem_) {}

// Destructors
AULSolver::~AULSolver() {}

// Getters
const DDPSolver AULSolver::DDPsolver() const { return DDPsolver_; }
const vector<vectordb> AULSolver::list_x() const { return list_x_; }
const vector<vectordb> AULSolver::list_u() const { return list_u_; }
const vector<vectordb> AULSolver::list_eq() const { return list_eq_; }
const vector<vectordb> AULSolver::list_ineq() const { return list_ineq_; }
const vectordb AULSolver::teq() const { return teq_; }
const vectordb AULSolver::tineq() const { return tineq_; }
const double AULSolver::cost() const { return cost_; }
const vector<vectordb> AULSolver::list_lambda() const { return list_lambda_; }
const vector<vectordb> AULSolver::list_mu() const { return list_mu_; }
const unsigned int AULSolver::AUL_n_iter() const { return AUL_n_iter_; }
const unsigned int AULSolver::DDP_n_iter() const { return DDP_n_iter_; }
const vector<vector<vectordb>> AULSolver::list_x_mem() const { return list_x_mem_; }
const vector<vector<vectordb>> AULSolver::list_u_mem() const { return list_u_mem_; }


// Setters
void AULSolver::set_ToF(double const& ToF) {
	DDPsolver_.set_ToF(ToF);
}
void AULSolver::set_homotopy_coefficient(double const& homotopy_coefficient) {
	DDPsolver_.set_homotopy_coefficient(homotopy_coefficient);
}
void AULSolver::set_huber_loss_coefficient(double const& huber_loss_coefficient) {
	DDPsolver_.set_huber_loss_coefficient(huber_loss_coefficient);
}

// Compute the value of the max constraint.
double AULSolver::get_max_constraint_() {
	// Unpack parameters
	SolverParameters solver_parameters = DDPsolver_.solver_parameters();
	unsigned int N = solver_parameters.N();
	unsigned int Neq = solver_parameters.Neq();
	unsigned int Nineq = solver_parameters.Nineq();
	unsigned int Nteq = solver_parameters.Nteq();
	unsigned int Ntineq = solver_parameters.Ntineq();

	// Loop on all steps
	double max = -1e15;
	for (size_t i = 0; i < N; i++) {

		// Unpack
		vectordb eq_i = list_eq_[i];
		vectordb ineq_i = list_ineq_[i];

		// Find max
		for (size_t j = 0; j < Neq; j++) {
			double abs_eq_j = abs(eq_i[j]);
			if (abs_eq_j > max)
				max = abs_eq_j;
		}
		for (size_t j = 0; j < Nineq; j++) {
			double ineq_j = ineq_i[j];
			if (ineq_j > max)
				max = ineq_j;
		}
	}

	// Terminal constraints
	for (size_t j = 0; j < Nteq; j++) {
		double abs_teq_j = abs(teq_[j]);
		if (abs_teq_j > max)
			max = abs_teq_j;
	}
	for (size_t j = 0; j < Ntineq; j++) {
		double ineq_j = tineq_[j];
		if (ineq_j > max)
			max = ineq_j;
	}

	return max;
}

// Update dual state in Augmented Lagrangian formulation.
void AULSolver::update_lambda_() {
	// Unpack parameters
	SolverParameters solver_parameters = DDPsolver_.solver_parameters();
	unsigned int N = solver_parameters.N();
	unsigned int Nx = solver_parameters.Nx();
	unsigned int Nu = solver_parameters.Nu();
	unsigned int Neq = solver_parameters.Neq();
	unsigned int Nineq = solver_parameters.Nineq();
	unsigned int Nteq = solver_parameters.Nteq();
	unsigned int Ntineq = solver_parameters.Ntineq();
	double lambda_ub = solver_parameters.lambda_parameters()[1];
	
	// Init output
	vector<vectordb> list_lambda;
	list_lambda.reserve(list_lambda.size());

	// Loop on path constraints
	for (size_t i = 0; i < N; i++) {
		
		// Unpack
		vectordb lambda = list_lambda_[i];
		vectordb mu = list_mu_[i];
		vectordb ineq = list_ineq_[i];
		vectordb eq = list_eq_[i]; 
		
		// Iterate on constraints
		double buff = 0.0;
		for (size_t j = 0; j < Neq; j++) {
			buff = lambda[j] + mu[j] * eq[j];
			if (buff > 0)
				lambda[j] = min(buff, lambda_ub);
			else
				lambda[j] = max(buff, -lambda_ub);
		}
		for (size_t j = 0; j < Nineq; j++) {
			buff = max(0.0, lambda[j + Neq] + mu[j + Neq] * ineq[j]);
			lambda[j + Neq] = min(buff, lambda_ub);
		}

		// Assign
		list_lambda.push_back(lambda);
	}

	// Terminal contraints
	
	// Unpack
	vectordb lambda = list_lambda_[N];
	vectordb mu = list_mu_[N];
	vectordb ineq = tineq_;
	vectordb eq = teq_;

	// Iterate on constraints
	double buff = 0.0;
	for (size_t j = 0; j < Nteq; j++) {
		buff = lambda[j] + mu[j] * eq[j];
		if (buff > 0)
			lambda[j] = min(buff, lambda_ub);
		else
			lambda[j] = max(buff, -lambda_ub);
	}
	for (size_t j = 0; j < Ntineq; j++) {
		buff = max(0.0, lambda[j + Nteq] + mu[j + Nteq] * ineq[j]);
		lambda[j + Nteq] = min(buff, lambda_ub);
	}

	// Assign
	list_lambda.push_back(lambda);
	DDPsolver_.set_list_lambda(list_lambda);
}

// Update penalty in Augmented Lagrangian formulation.
void AULSolver::update_mu_() {
	// Unpack parameters
	SolverParameters solver_parameters = DDPsolver_.solver_parameters();
	unsigned int N = solver_parameters.N();
	unsigned int Nx = solver_parameters.Nx();
	unsigned int Nu = solver_parameters.Nu();
	unsigned int Neq = solver_parameters.Neq();
	unsigned int Nineq = solver_parameters.Nineq();
	unsigned int Nteq = solver_parameters.Nteq();
	unsigned int Ntineq = solver_parameters.Ntineq();
	vectordb mu_parameters = solver_parameters.mu_parameters();
	double mu_init = mu_parameters[0]; double mu_ub = mu_parameters[1];
	double mu_factor = mu_parameters[2];

	vector<vectordb> list_mu;
	list_mu.reserve(list_mu_.size());
	double buff = 0.0;
	for (size_t i = 0; i < N; i++) {
		// Unpack
		vectordb mu = list_mu_[i];

		// Iterate on constraints
		for (size_t j = 0; j < Neq + Nineq; j++) {
			buff = mu_factor * mu[j];
			mu[j] = max(0, min(buff, mu_ub));
		}

		// Assign
		list_mu.push_back(mu);
	}

	// Unpack
	vectordb mu = list_mu_[N];

	// Iterate on constraints
	for (size_t j = 0; j < Nteq + Ntineq; j++) {
		buff = mu_factor * mu[j];
		mu[j] = max(0, min(buff, mu_ub));
	}

	// Assign
	list_mu.push_back(mu);
	DDPsolver_.set_list_mu(list_mu);
}

// Performs AUL solving given a starting point,
// initial controls and a final state.
// Inspired from ALTRO (Julia).
// See: https://github.com/RoboticExplorationLab/Altro.jl
void AULSolver::solve(
	vectordb const& x0,
	vector<vectordb> const& list_u_init,
	vectordb const& x_goal) {
	// Unpack parameters
	SolverParameters solver_parameters = DDPsolver_.solver_parameters();
	unsigned int N = solver_parameters.N();
	unsigned int Neq = solver_parameters.Neq();
	unsigned int Nineq = solver_parameters.Nineq();
	unsigned int Nteq = solver_parameters.Nteq();
	unsigned int Ntineq = solver_parameters.Ntineq();
	double AUL_tol = solver_parameters.AUL_tol();
	int AUL_max_iter = solver_parameters.AUL_max_iter();
	vectordb mu_parameters = solver_parameters.mu_parameters();
	vectordb lambda_parameters = solver_parameters.lambda_parameters();
	unsigned int verbosity = solver_parameters.verbosity();
	unsigned int saving_iterations = solver_parameters.saving_iterations();
	Constants constants = DDPsolver_.dynamics().constants();

	// Init DDPsolver
	DDPsolver_.set_ToF(x_goal[x_goal.size() - 1]);
	DDPsolver_.set_recompute_dynamics(true);

	// Output
	auto start_aul = high_resolution_clock::now();
	if (verbosity < 1) {
		cout << "############################################################################" << endl;
		cout << "#                                                                          #" << endl;
		cout << "#                             START AUL SOLVING                            #" << endl;
		cout << "#                                                                          #" << endl;
		cout << "############################################################################" << endl << endl;
		cout << "Homotopy coefficient [0, 1] : " << solver_parameters.homotopy_coefficient() << endl;
		cout << "Huber-loss coefficient [0, 1] : " << solver_parameters.huber_loss_coefficient();
		cout << endl << endl << endl;
	}
	else if (verbosity < 2) {
		cout << endl;
		cout << "AUL solving - Homotopy coefficient [0, 1] : " << solver_parameters.homotopy_coefficient() << endl;
		cout << "            - Huber-loss coefficient [0, 1] : " << solver_parameters.huber_loss_coefficient() << endl;
		cout << "	ITERATION [-], DDP ITERATIONS [-], RUNTIME [s], FINAL MASS [kg], MAX CONSTRAINT [-]" << endl;
	}

	// Init lists dual state and penalty factors lists
	double lambda_0 = lambda_parameters[0];
	double mu_0 = mu_parameters[0];
	list_lambda_ = vector<vectordb>();
	list_mu_ = vector<vectordb>();
	list_lambda_.reserve(N + 1); list_mu_.reserve(N + 1);
	for (size_t i = 0; i < N; i++) {
		list_lambda_.emplace_back(Neq + Nineq, lambda_0);
		list_mu_.emplace_back(Neq + Nineq, mu_0);
	}
	list_lambda_.emplace_back(Nteq + Ntineq, lambda_0);
	list_mu_.emplace_back(Nteq + Ntineq, mu_0);

	// Set dual state and penalty factors lists
	DDPsolver_.set_list_lambda(list_lambda_);
	DDPsolver_.set_list_mu(list_mu_);
	list_u_ = list_u_init;

	// Init loop variables
	bool loop = true;
	double cost = 1e15;
	cost_ = cost;
	AUL_n_iter_ = 0;
	size_t counter_rejected = 0;
	double max_constraint = cost;
	while (loop && AUL_n_iter_ < AUL_max_iter) {

		// Solve DDP problem
		auto start = high_resolution_clock::now();
		DDPsolver_.solve(x0, list_u_, x_goal);
		auto stop = high_resolution_clock::now();
		auto duration_mapping = duration_cast<microseconds>(stop - start);

		// Check constraints and that the solver is not stuck
		double max_constraint_new = DDPsolver_.get_max_constraint_();
		if (abs((max_constraint_new - max_constraint)/ max_constraint) < AUL_tol) {
			DDPsolver_.set_recompute_dynamics(false);
		}
		else {DDPsolver_.set_recompute_dynamics(true);}
		max_constraint = max_constraint_new;

		// Store results
		list_x_ = DDPsolver_.list_x();
		list_u_ = DDPsolver_.list_u();
		list_eq_ = DDPsolver_.list_eq();
		list_ineq_ = DDPsolver_.list_ineq();
		teq_ = DDPsolver_.teq();
		tineq_ = DDPsolver_.tineq();
		cost_ = DDPsolver_.cost();

		// Output
		if (verbosity < 1) {
			cout << AUL_n_iter_ << " - RUNTIME [s] : "
				<< to_string(static_cast<double>(duration_mapping.count()) / 1e6) << ", "
				<< "FINAL MASS [kg] : " << DDPsolver_.list_x()[N][SIZE_VECTOR] * constants.massu() << ", "
				<< "MAX CONSTRAINT [-] : " << max_constraint << endl << endl;
		}
		else if (verbosity < 2) {
			cout << "	" << AUL_n_iter_ << ", " << DDPsolver_.n_iter()
				<< ", "	<< to_string(static_cast<double>(duration_mapping.count()) / 1e6)
				<< ", " << DDPsolver_.list_x()[N][SIZE_VECTOR] * constants.massu()
				<< ", " << max_constraint << endl;
		}

		// Update dual state and penalities
		update_lambda_(); update_mu_();		
		list_lambda_ = DDPsolver_.solver_parameters().list_lambda();
		list_mu_ = DDPsolver_.solver_parameters().list_mu();

		// Check constraints
		bool force_continue_loop = max_constraint > AUL_tol;

		// Stopping conditions
		bool force_stop_loop = AUL_n_iter_ > AUL_max_iter;
		loop = !force_stop_loop && force_continue_loop;
		cost = DDPsolver_.cost();
		DDP_n_iter_ += DDPsolver_.n_iter();
		AUL_n_iter_++;

		// Save iterations
		if (saving_iterations > 2) { // Save all DDP iterations
			for (size_t k = 0; k < DDPsolver_.list_u_mem().size(); k++) {
				list_x_mem_.push_back(DDPsolver_.list_x_mem()[k]);
				list_u_mem_.push_back(DDPsolver_.list_u_mem()[k]);
			}
		}
		else if (saving_iterations > 2) { // Save all AUL iterations
			list_x_mem_.push_back(list_x_);
			list_u_mem_.push_back(list_u_);
		}
		else if (saving_iterations > 1 && !loop) { // Save last AUL iteration
			list_x_mem_.push_back(list_x_);
			list_u_mem_.push_back(list_u_);
		}
	}

	// Output
	auto stop_aul = high_resolution_clock::now();
	auto duration_aul = duration_cast<microseconds>(stop_aul - start_aul);
	if (verbosity < 1) {
		cout << "Runtime : " + to_string(static_cast<double>(duration_aul.count()) / 1e6) + "s" << endl;
	}
	else if (verbosity < 2) {
		cout << "Runtime : " + to_string(static_cast<double>(duration_aul.count()) / 1e6) + "s" << endl;
	}
}
