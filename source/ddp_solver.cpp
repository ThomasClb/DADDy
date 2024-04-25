/**
	ddp_solver.cpp

	Purpose: Implementation of the DDPSolver class.

	@author Thomas Caleb

	@version 1.0 17/11/2023
*/

#include "ddp_solver.h"

using namespace DACE;
using namespace std;
using namespace std::chrono;

// Empty constructor
DDPSolver::DDPSolver() : solver_parameters_(),
	spacecraft_parameters_(Constants()), dynamics_(),
	list_x_(vector<vectordb>(0)), list_u_(vector<vectordb>(0)), cost_(0),
	list_eq_(vector<vectordb>(0)), list_ineq_(vector<vectordb>(0)),
	rho_(0.0), d_rho_(0.0),
	list_x_mem_(vector<vector<vectordb>>(0)), list_u_mem_(vector<vector<vectordb>>(0)) {}

// Constructor
DDPSolver::DDPSolver(
	SolverParameters const& solver_parameters,
	SpacecraftParameters const& spacecraft_parameters,
	Dynamics const& dynamics) : solver_parameters_(solver_parameters),
	spacecraft_parameters_(spacecraft_parameters), dynamics_(dynamics),
	list_x_(vector<vectordb>(0)), list_u_(vector<vectordb>(0)), cost_(0),
	list_eq_(vector<vectordb>(0)), list_ineq_(vector<vectordb>(0)),
	rho_(0.0), d_rho_(0.0), alpha_(1.0),
	list_x_mem_(vector<vector<vectordb>>(0)), list_u_mem_(vector<vector<vectordb>>(0)) {}

// Copy constructor
DDPSolver::DDPSolver(
	DDPSolver const& solver) : solver_parameters_(solver.solver_parameters_),
	spacecraft_parameters_(solver.spacecraft_parameters()), dynamics_(solver.dynamics()),
	list_x_(solver.list_x_), list_u_(solver.list_u_), cost_(solver.cost_),
	list_eq_(solver.list_eq_), list_ineq_(solver.list_ineq_),
	rho_(0.0), d_rho_(0.0), alpha_(1.0),
	list_x_mem_(solver.list_x_mem_), list_u_mem_(solver.list_u_mem_) {}

// Destructors
DDPSolver::~DDPSolver() {}

// Getters
SolverParameters DDPSolver::solver_parameters() { return solver_parameters_; }
const SpacecraftParameters DDPSolver::spacecraft_parameters() const {
	return spacecraft_parameters_;
}
const Dynamics DDPSolver::dynamics() const { return dynamics_; }
const vector<vectordb> DDPSolver::list_x() const { return list_x_; }
const vector<vectordb> DDPSolver::list_u() const { return list_u_; }
const vector<vectordb> DDPSolver::list_eq() const { return list_eq_; }
const vector<vectordb> DDPSolver::list_ineq() const { return list_ineq_; }
const vectordb DDPSolver::teq() const { return teq_; }
const vectordb DDPSolver::tineq() const { return tineq_; }
const double DDPSolver::cost() const { return cost_; }
const DA DDPSolver::tc_eval() const { return tc_eval_; }
const unsigned int DDPSolver::n_iter() const { return n_iter_; }
const vector<vector<vectordb>> DDPSolver::list_x_mem() const { return list_x_mem_; }
const vector<vector<vectordb>> DDPSolver::list_u_mem() const { return list_u_mem_; }

// Setters
void DDPSolver::set_list_lambda(vector<vectordb> const& list_lambda) {
	solver_parameters_.set_list_lambda(list_lambda);
}
void DDPSolver::set_list_mu(vector<vectordb> const& list_mu) {
	solver_parameters_.set_list_mu(list_mu);
}
void DDPSolver::set_ToF(double const& ToF) {
	solver_parameters_.set_ToF(ToF);
}
void DDPSolver::set_homotopy_coefficient(double const& homotopy_coefficient) {
	solver_parameters_.set_homotopy_coefficient(homotopy_coefficient);
}
void DDPSolver::set_huber_loss_coefficient(double const& huber_loss_coefficient) {
	solver_parameters_.set_huber_loss_coefficient(huber_loss_coefficient);
}
void DDPSolver::set_recompute_dynamics(bool const& recompute_dynamics) {
	recompute_dynamics_ = recompute_dynamics;
}

// Returns the Augmented lagrangian cost-to-go:
// AUL_ctg = ctg + Lambda^T * c(x) + 0.5 * c(x)^T * I_mu * c(x)
vectorDA DDPSolver::get_AUL_cost_to_go(
	vectorDA const& x_star_DA, vectorDA const& u_star_DA, size_t const& index) {
	// Unpack parameters
	unsigned int N = solver_parameters_.N();
	unsigned int Nx = solver_parameters_.Nx();
	unsigned int Nu = solver_parameters_.Nu();
	unsigned int Neq = solver_parameters_.Neq();
	unsigned int Nineq = solver_parameters_.Nineq();
	vectordb lambda = solver_parameters_.list_lambda()[index];
	vectordb mu = solver_parameters_.list_mu()[index];

	// Init output
	vectorDA constraints_eval;
	constraints_eval.reserve(Neq + Nineq + 1);

	// Evaluate ctg
	DA ctg_eval = dynamics_.cost_to_go()(
		x_star_DA, u_star_DA,
		spacecraft_parameters_, dynamics_.constants(), solver_parameters_);
	
	// Constraints evaluations
	vectorDA eq_eval = dynamics_.equality_constraints()(
		x_star_DA, u_star_DA,
		spacecraft_parameters_, dynamics_.constants(), solver_parameters_);
	vectorDA ineq_eval = dynamics_.inequality_constraints()(
		x_star_DA, u_star_DA,
		spacecraft_parameters_, dynamics_.constants(), solver_parameters_);

	// Assign to output
	for (size_t i = 0; i < Neq; i++) {
		constraints_eval.push_back(eq_eval[i]);
	}
	for (size_t i = 0; i < Nineq; i++) {
		constraints_eval.push_back(ineq_eval[i]);
	}

	// Set mu=0 for inactive constraints
	for (size_t i = Neq; i < Neq + Nineq; i++) {
		if (constraints_eval[i].cons() < 0 && lambda[i] <= 0)
			mu[i] *= 0.0;
	}

	// AUL formulation
	// Lambda^T * c(x) + 0.5 * c(x)^T * I_mu * c(x)
	for (size_t i = 0; i < Neq + Nineq; i++) {
		DA constraints_eval_i = constraints_eval[i];
		ctg_eval += constraints_eval_i*(lambda[i] + (0.5 * mu[i]) *  constraints_eval_i);
	}

	// Assign
	constraints_eval.push_back(ctg_eval);

	return constraints_eval;
}

// Returns the Augmented lagrangian terminal cost:
// AUL_tc = tc + Lambda^T * c(x) + 0.5 * c(x)^T * I_mu * c(x)
vectorDA DDPSolver::get_AUL_terminal_cost(
	vectorDA const& x_star_DA, vectordb const& x_goal) {
	// Unpack parameters
	unsigned int N = solver_parameters_.N();
	unsigned int Nx = solver_parameters_.Nx();
	unsigned int Nu = solver_parameters_.Nu();
	unsigned int Nteq = solver_parameters_.Nteq();
	unsigned int Ntineq = solver_parameters_.Ntineq();
	vectordb lambda = solver_parameters_.list_lambda()[N];
	vectordb mu = solver_parameters_.list_mu()[N];

	// Init output
	vectorDA constraints_eval;
	constraints_eval.reserve(Nteq + Ntineq + 1);

	// Evaluate terminal cost
	DA tc_eval = dynamics_.terminal_cost()(
		x_star_DA, x_goal,
		spacecraft_parameters_, dynamics_.constants(), solver_parameters_);

	// Constraints evaluations
	vectorDA teq_eval = dynamics_.terminal_equality_constraints()(
		x_star_DA, x_goal,
		spacecraft_parameters_, dynamics_.constants(), solver_parameters_);
	vectorDA tineq_eval = dynamics_.terminal_inequality_constraints()(
		x_star_DA, x_goal,
		spacecraft_parameters_, dynamics_.constants(), solver_parameters_);
	
	// Assign to output
	for (size_t i = 0; i < Nteq; i++) {
		constraints_eval.push_back(teq_eval[i]);
	}
	for (size_t i = 0; i < Ntineq; i++) {
		constraints_eval.push_back(tineq_eval[i]);
	}

	// Set mu=0 for inactive constraints
	for (size_t i = Nteq; i < Nteq + Ntineq; i++) {
		if (constraints_eval[i].cons() < 0 && lambda[i] <= 0)
			mu[i] *= 0.0;
	}

	// AUL formulation
	// Lambda^T * c(x) + 0.5 * c(x)^T * I_mu * c(x)
	for (size_t i = 0; i < Nteq + Ntineq; i++) {
		DA constraints_eval_i = constraints_eval[i];
		tc_eval += constraints_eval_i * (lambda[i] + (0.5 * mu[i]) * constraints_eval_i);
	}

	// Assign
	constraints_eval.push_back(tc_eval);

	return constraints_eval;
}

// Returns the Augmented lagrangian cost-to-go:
// AUL_ctg = ctg + Lambda^T * c(x) + 0.5 * c(x)^T * I_mu * c(x)
// Double version
vectordb DDPSolver::get_AUL_cost_to_go(
	vectordb const& x_star_DA, vectordb const& u_star_DA, size_t const& index) {
	// Unpack parameters
	unsigned int N = solver_parameters_.N();
	unsigned int Nx = solver_parameters_.Nx();
	unsigned int Nu = solver_parameters_.Nu();
	unsigned int Neq = solver_parameters_.Neq();
	unsigned int Nineq = solver_parameters_.Nineq();
	vectordb lambda = solver_parameters_.list_lambda()[index];
	vectordb mu = solver_parameters_.list_mu()[index];

	// Init output
	vectordb constraints_eval;
	constraints_eval.reserve(Neq + Nineq + 1);

	// Evaluate ctg
	double ctg_eval = dynamics_.cost_to_go_db()(
		x_star_DA, u_star_DA,
		spacecraft_parameters_, dynamics_.constants(), solver_parameters_);

	// Constraints evaluations
	vectordb eq_eval = dynamics_.equality_constraints_db()(
		x_star_DA, u_star_DA,
		spacecraft_parameters_, dynamics_.constants(), solver_parameters_);
	vectordb ineq_eval = dynamics_.inequality_constraints_db()(
		x_star_DA, u_star_DA,
		spacecraft_parameters_, dynamics_.constants(), solver_parameters_);

	// Assign to output
	for (size_t i = 0; i < Neq; i++) {
		constraints_eval.push_back(eq_eval[i]);
	}
	for (size_t i = 0; i < Nineq; i++) {
		constraints_eval.push_back(ineq_eval[i]);
	}

	// Set mu=0 for inactive constraints
	for (size_t i = Neq; i < Neq + Nineq; i++) {
		if (constraints_eval[i] < 0 && lambda[i] <= 0)
			mu[i] *= 0.0;
	}

	// AUL formulation
	// Lambda^T * c(x) + 0.5 * c(x)^T * I_mu * c(x)
	for (size_t i = 0; i < Neq + Nineq; i++) {
		double constraints_eval_i = constraints_eval[i];
		ctg_eval += constraints_eval_i * (lambda[i] + (0.5 * mu[i]) * constraints_eval_i);
	}

	// Assign
	constraints_eval.push_back(ctg_eval);

	return constraints_eval;
}

// Returns the Augmented lagrangian terminal cost:
// AUL_tc = tc + Lambda^T * c(x) + 0.5 * c(x)^T * I_mu * c(x)
// Double version
vectordb DDPSolver::get_AUL_terminal_cost(
	vectordb const& x_star_DA, vectordb const& x_goal) {
	// Unpack parameters
	unsigned int N = solver_parameters_.N();
	unsigned int Nx = solver_parameters_.Nx();
	unsigned int Nu = solver_parameters_.Nu();
	unsigned int Nteq = solver_parameters_.Nteq();
	unsigned int Ntineq = solver_parameters_.Ntineq();
	vectordb lambda = solver_parameters_.list_lambda()[N];
	vectordb mu = solver_parameters_.list_mu()[N];

	// Init output
	vectordb constraints_eval;
	constraints_eval.reserve(Nteq + Ntineq + 1);

	// Evaluate terminal cost
	double tc_eval = dynamics_.terminal_cost_db()(
		x_star_DA, x_goal,
		spacecraft_parameters_, dynamics_.constants(), solver_parameters_);

	// Constraints evaluations
	vectordb teq_eval = dynamics_.terminal_equality_constraints_db()(
		x_star_DA, x_goal,
		spacecraft_parameters_, dynamics_.constants(), solver_parameters_);
	vectordb tineq_eval = dynamics_.terminal_inequality_constraints_db()(
		x_star_DA, x_goal,
		spacecraft_parameters_, dynamics_.constants(), solver_parameters_);

	// Assign to output
	for (size_t i = 0; i < Nteq; i++) {
		constraints_eval.push_back(teq_eval[i]);
	}
	for (size_t i = 0; i < Ntineq; i++) {
		constraints_eval.push_back(tineq_eval[i]);
	}

	// Set mu=0 for inactive constraints
	for (size_t i = Nteq; i < Nteq + Ntineq; i++) {
		if (constraints_eval[i] < 0 && lambda[i] <= 0)
			mu[i] *= 0.0;
	}

	// AUL formulation
	// Lambda^T * c(x) + 0.5 * c(x)^T * I_mu * c(x)
	for (size_t i = 0; i < Nteq + Ntineq; i++) {
		double constraints_eval_i = constraints_eval[i];
		tc_eval += constraints_eval_i * (lambda[i] + (0.5 * mu[i]) * constraints_eval_i);
	}

	// Assign
	constraints_eval.push_back(tc_eval);

	return constraints_eval;
}

// Increases the regulation to ensure Quu + rho*I
// is symeteric positive definite.
// Inspired from ALTRO (Julia).
// See: https://github.com/RoboticExplorationLab/Altro.jl
void DDPSolver::increase_regularisation_() {
	// Unpack
	vectordb rho_parameters = solver_parameters_.backward_sweep_regulation_parameters();
	double rho_init = rho_parameters[0];
	double rho_min = rho_parameters[1];
	double rho_max = rho_parameters[2];
	double rho_factor = rho_parameters[3];

	// Update d_rho_
	d_rho_ = max(d_rho_ * rho_factor, rho_factor);
	rho_ = min(rho_max, max(rho_ * d_rho_, rho_min));
}

// Increases the regulation to ensure Quu + rho*I
// is symeteric positive definite.
// Inspired from ALTRO (Julia).
// With a safe guard.
// See: https://github.com/RoboticExplorationLab/Altro.jl
void DDPSolver::increase_regularisation_(matrixdb const& Quu) {
	// Unpack
	vectordb rho_parameters = solver_parameters_.backward_sweep_regulation_parameters();
	double rho_init = rho_parameters[0];
	double rho_min = rho_parameters[1];
	double rho_max = rho_parameters[2];
	double rho_factor = rho_parameters[3];
	
	// Update d_rho_
	d_rho_ = max(d_rho_ * rho_factor, rho_factor);

	// If the reg factor is already too large
	if (rho_ >= rho_max) {
		// The small (or largest negative in absolute value) eigenvalue
		// is smaller than the norm of the symetric matrix.
		// ie |S|_2 = |D|_2 >= |eig_max|
		rho_ = frobenius_norm_(Quu); 
		
	}
	else {
		// Set rho
		rho_ = min(rho_max, max(rho_ * d_rho_, rho_min));
	}
}

// Decreases the regulation to ensure Quu + rho*I
// is symeteric positive definite.
// Inspired from ALTRO (Julia).
// See: https://github.com/RoboticExplorationLab/Altro.jl
void DDPSolver::decrease_regularisation_() {
	// Unpack
	vectordb rho_parameters = solver_parameters_.backward_sweep_regulation_parameters();
	double rho_init = rho_parameters[0];
	double rho_min = rho_parameters[1];
	double rho_max = rho_parameters[2];
	double rho_factor = rho_parameters[3];

	// Set rho and d_rho
	d_rho_ = min(d_rho_ / rho_factor, 1/rho_factor);
	rho_ = max(min(rho_ * d_rho_, rho_max), rho_min);
}

// Computes the expected cost after backward sweep for linesearch.
// Inspired from ALTRO (Julia).
// See: https://github.com/RoboticExplorationLab/Altro.jl
double DDPSolver::expected_cost_(
	double const& alpha) {
	// Unpack parameters
	unsigned int N = solver_parameters_.N();

	// Init
	double dv_1 = 0.0;
	double dv_2 = 0.0;

	// Add expected costs
	for (size_t i = 0; i < N; i++) {
		matrixdb k_i = list_k_[N - 1 - i];
		matrixdb Qu_i = list_Qu_[N - 1 - i];
		dv_1 += (Qu_i * k_i).at(0, 0);
		dv_2 += 0.5*(k_i.transpose() * Qu_i.transpose()).at(0, 0);
	}
	return -1.0* alpha * (dv_1 + alpha * dv_2);
}

// Evaluates the convergence of DDP optimisation
bool DDPSolver::evaluate_convergence_(double const& d_cost) {
	// Unpack parameters
	double tol = solver_parameters_.DDP_tol();
	unsigned int max_iter = solver_parameters_.DDP_max_iter();

	// Converged
	if ((d_cost < tol && d_cost >= 0.0 && (n_iter_ >= 3 || d_cost==0))
		|| n_iter_ >= max_iter)
		return true;

	return false;
}

// Compute the value of the max constraint.
double DDPSolver::get_max_constraint_() {
	// Unpack parameters
	SolverParameters solver_parameters = solver_parameters_;
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
		double tineq_j = tineq_[j];
		if (tineq_j > max)
			max = tineq_j;
	}

	return max;
}

// Performs the DDP backward sweep, that consists in the computation
// of the gains corrections.
// Inspired from ALTRO (Julia).
// See: https://github.com/RoboticExplorationLab/Altro.jl
void DDPSolver::backward_sweep_() {
	// Unpack parameters
	unsigned int N = solver_parameters_.N();
	unsigned int Nx = solver_parameters_.Nx();
	unsigned int Nu = solver_parameters_.Nu();
	vectordb rho_parameters = solver_parameters_.backward_sweep_regulation_parameters();
	double rho_max = rho_parameters[2];

	// Evaluate terminal cost derivatives
	vector<matrixdb> der_terminal_cost = deriv_x(
		vectorDA{ tc_eval_ }, Nx, true);

	// Init Vx0, Vxx0
	matrixdb Vx0 = der_terminal_cost[0].transpose();
	matrixdb Vxx0 = der_terminal_cost[1];

	// Init vectors
	vectorDA vect_ctg(1);

	// Backward loop
	bool success = false;
	while (!success) {
		// Init Vx, Vxx
		matrixdb Vx = Vx0;
		matrixdb Vxx = Vxx0;

		// Init gains
		list_Qu_ = vector<matrixdb>(0); list_Qu_.reserve(N);
		list_K_ = vector<matrixdb>(0); list_K_.reserve(N);
		list_k_ = vector<matrixdb>(0); list_k_.reserve(N);
		success = true;
		for (int j = N - 1; j >= 0; j--) {

			// Get derivatives
			vect_ctg[0] = list_ctg_eval_[j];
			vector<matrixdb> der_dynamic = deriv_xu(
				list_dynamic_eval_[j], Nx, Nu, false);
			vector<matrixdb> der_cost_to_go = deriv_xu(
				vect_ctg, Nx, Nu, true);

			// Unpack derivatives and transpose if needed
			matrixdb der_dynamic_0(der_dynamic[0]);
			matrixdb der_dynamic_1(der_dynamic[1]);
			matrixdb der_dynamic_1_t(der_dynamic_1.transpose());
			matrixdb Vx_t = Vx.transpose();

			// Get Q derivatives
			matrixdb Qx(der_cost_to_go[0] + Vx_t * der_dynamic_0);
			matrixdb Qu(der_cost_to_go[1] + Vx_t * der_dynamic_1);
			matrixdb Qxx(der_cost_to_go[2] + der_dynamic_0.transpose() * Vxx * der_dynamic_0);
			matrixdb Quu(der_cost_to_go[3] + der_dynamic_1_t * Vxx * der_dynamic_1);
			matrixdb Qux(der_cost_to_go[4] + der_dynamic_1_t * Vxx * der_dynamic_0);

			// Regularisation
			for (size_t k = 0; k < Nu; k++) {
				Quu.at(k, k) += rho_;
			}

			// Check that Quu is definite positive (hence invertible)
			if (!is_def_pos_(Quu)) {
				if (rho_ == rho_max)
					increase_regularisation_(Quu); // safe version
				else
					increase_regularisation_();
				success = false; break;
			}

			// Compute Cholesky Factorisation to ease solving
			matrixdb Luu = cholesky_(Quu);

			// Compute gains
			matrixdb Qux_t = Qux.transpose();
			matrixdb K = -1.0 * solve_cholesky_(Luu, Qux);
			matrixdb k = -1.0 * solve_cholesky_(Luu, Qu.transpose());

			// Store gains and Qu
			list_Qu_.push_back(Qu);
			list_K_.push_back(K);
			list_k_.push_back(k);

			// Get Vx and Vxx for next step
			Vx = Qx.transpose() + Qux_t * k;
			Vxx = Qxx + Qux_t * K;
		}
	}
	decrease_regularisation_();
}

// Performs the DDP backward sweep, that consists in the computation
// of the gains corrections, with hessian.
// Inspired from ALTRO (Julia).
// See: https://github.com/RoboticExplorationLab/Altro.jl
void DDPSolver::backward_sweep_hessian_() {
	// Unpack parameters
	unsigned int N = solver_parameters_.N();
	unsigned int Nx = solver_parameters_.Nx();
	unsigned int Nu = solver_parameters_.Nu();
	vectordb rho_parameters = solver_parameters_.backward_sweep_regulation_parameters();
	double rho_max = rho_parameters[2];

	// Evaluate terminal cost derivatives
	vector<matrixdb> der_terminal_cost = deriv_x(
		vectorDA{ tc_eval_ }, Nx, true);

	// Init Vx0, Vxx0
	matrixdb Vx0 = der_terminal_cost[0].transpose();
	matrixdb Vxx0 = der_terminal_cost[1];

	// Init vectors
	vectorDA vect_ctg(1);

	// Backward loop
	bool success = false;
	while (!success) {
		// Init Vx, Vxx
		matrixdb Vx = Vx0;
		matrixdb Vxx = Vxx0;

		// Init gains
		list_Qu_ = vector<matrixdb>(0); list_Qu_.reserve(N);
		list_K_ = vector<matrixdb>(0); list_K_.reserve(N);
		list_k_ = vector<matrixdb>(0); list_k_.reserve(N);
		success = true;
		for (int j = N - 1; j >= 0; j--) {

			// Get derivatives
			vect_ctg[0] = list_ctg_eval_[j];
			vector<matrixdb> der_dynamic = deriv_xu(
				list_dynamic_eval_[j], Nx, Nu, true);
			vector<matrixdb> der_cost_to_go = deriv_xu(
				vect_ctg, Nx, Nu, true);

			// Unpack derivatives and transpose if needed
			matrixdb der_dynamic_0(der_dynamic[0]);
			matrixdb der_dynamic_1(der_dynamic[1]);
			matrixdb der_dynamic_1_t(der_dynamic_1.transpose());
			matrixdb Vx_t = Vx.transpose();

			// Get Q derivatives
			matrixdb Qx(der_cost_to_go[0] + Vx_t * der_dynamic_0);
			matrixdb Qu(der_cost_to_go[1] + Vx_t * der_dynamic_1);
			matrixdb Qxx(der_cost_to_go[2] + der_dynamic_0.transpose() * Vxx * der_dynamic_0);
			matrixdb Quu(der_cost_to_go[3] + der_dynamic_1_t * Vxx * der_dynamic_1);
			matrixdb Qux(der_cost_to_go[4] + der_dynamic_1_t * Vxx * der_dynamic_0);

			// Add hessians
			/**/
			for (size_t k=0; k < Nx; k++) {
				matrixdb fxx_k = der_dynamic[2 + k];
				matrixdb fuu_k = der_dynamic[2 + k + Nx];
				matrixdb fux_k = der_dynamic[2 + k + Nx*2];
				double Vx_k = Vx.at(k, 0);
				Qxx = Qxx + Vx_k * fxx_k;
				Quu = Quu + Vx_k * fuu_k;
				Qux = Qux + Vx_k * fux_k;
			}
			
			// Regularisation
			for (size_t k = 0; k < Nu; k++) {
				Quu.at(k, k) += rho_;
			}

			// Check that Quu is definite positive (hence invertible)
			if (!is_def_pos_(Quu)) {
				if (rho_ == rho_max)
					increase_regularisation_(Quu); // safe version
				else
					increase_regularisation_();
				success = false; break;
			}

			// Compute Cholesky Factorisation to ease solving
			matrixdb Luu = cholesky_(Quu);

			// Compute gains
			matrixdb Qux_t = Qux.transpose();
			matrixdb K = -1.0 * solve_cholesky_(Luu, Qux);
			matrixdb k = -1.0 * solve_cholesky_(Luu, Qu.transpose());

			// Store gains and Qu
			list_Qu_.push_back(Qu);
			list_K_.push_back(K);
			list_k_.push_back(k);

			// Get Vx and Vxx for next step
			Vx = Qx.transpose() + Qux_t * k;
			Vxx = Qxx + Qux_t * K;
		}
	}
	decrease_regularisation_();
}

// Performs the DDP backward sweep, that consists in the computation
// of the gains corrections, Q is evaluated using DA, thus, with hessian.
// Inspired from ALTRO (Julia).
// See: https://github.com/RoboticExplorationLab/Altro.jl
void DDPSolver::backward_sweep_DA_Q_() {
	// Unpack parameters
	unsigned int N = solver_parameters_.N();
	unsigned int Nx = solver_parameters_.Nx();
	unsigned int Nu = solver_parameters_.Nu();
	vectordb rho_parameters = solver_parameters_.backward_sweep_regulation_parameters();
	double rho_max = rho_parameters[2];

	// Evaluate terminal cost derivatives
	vector<matrixdb> der_terminal_cost = deriv_x(
		vectorDA{ tc_eval_ }, Nx, true);

	// Init useful vector vectors
	vectorDA dx = id_vector(vectordb(Nx, 0.0), 0, 0, Nx);
	vectorDA dxu(Nu + Nx); for (size_t k = 0; k < Nx; k++) { dxu[k] = dx[k]; }
	matrixDA dx_mat(Nx, 1); dx_mat.setcol(0, dx);
	vectorDA df(Nu + Nx, 0.0); vectorDA vect_Qkp1(1);

	// Backward loop
	bool success = false;
	while (!success) {
		// Init V
		DA Vkp1 = tc_eval_;

		// Init gains
		list_Qu_ = vector<matrixdb>(0); list_Qu_.reserve(N);
		list_K_ = vector<matrixdb>(0); list_K_.reserve(N);
		list_k_ = vector<matrixdb>(0); list_k_.reserve(N);
		success = true;
		for (int j = N - 1; j >= 0; j--) {	 
			
			// Get dxkp1 = f(xk, uk) - xkp1
 			vectorDA dxkp1 = list_dynamic_eval_[j] - list_x_[j + 1];
			for (size_t k = 0; k < Nx; k++) { df[k] = dxkp1[k]; }

			// Get Qkp1 = V \cdot dxkp1 + lk
			DA Qkp1 = Vkp1.eval(df) + list_ctg_eval_[j];
			vect_Qkp1[0] = Qkp1;

			// Get Q derivatives
			vector<matrixdb> der_Qk = deriv_xu(
				vect_Qkp1, Nx, Nu, true);
			matrixdb Qu(der_Qk[1]), Quu(der_Qk[3]), Qux(der_Qk[4]);

			// Regularisation
			for (size_t k = 0; k < Nu; k++) {
				Quu.at(k, k) += rho_;
			}

			// Check that Quu is definite positive (hence invertible)
			if (!is_def_pos_(Quu)) {
				if (rho_ == rho_max)
					increase_regularisation_(Quu); // safe version
				else
					increase_regularisation_();
				success = false; break;
			}

			// Compute Cholesky Factorisation to ease solving
			matrixdb Luu = cholesky_(Quu);

			// Compute gains
			matrixdb K = -1.0 * solve_cholesky_(Luu, Qux);
			matrixdb k = -1.0 * solve_cholesky_(Luu, Qu.transpose());

			// Store gains and Qu
			list_Qu_.push_back(Qu);
			list_K_.push_back(K);
			list_k_.push_back(k);

			// Make du
			matrixDA du_mat = K * dx_mat + k;
			for (size_t k = 0; k < Nu; k++) { dxu[k + Nx] = du_mat.at(k, 0); }

			// Update Vkp1
			Vkp1 = Qkp1.eval(dxu);
		}
	}
	decrease_regularisation_();
}

// Performs the DDP forward pass, that consists in the computation
// of the new states and control after correction.
// Inspired from ALTRO (Julia).
// DA only for automatic differentiation.
// See: https://github.com/RoboticExplorationLab/Altro.jl
void DDPSolver::forward_pass_(
	vector<vectordb> const& list_x, vector<vectordb> const& list_u, vectordb const& x_goal) {
	// Unpack parameters
	double tol = solver_parameters_.DDP_tol();
	unsigned int N = solver_parameters_.N();
	unsigned int Nx = solver_parameters_.Nx();
	unsigned int Nu = solver_parameters_.Nu();
	unsigned int Neq = solver_parameters_.Neq();
	unsigned int Nineq = solver_parameters_.Nineq();
	unsigned int Nteq = solver_parameters_.Nteq();
	unsigned int Ntineq = solver_parameters_.Ntineq();
	vectordb alpha_parameters = solver_parameters_.line_search_parameters();
	double z_min = alpha_parameters[0];
	double z_max = alpha_parameters[1];
	double alpha_factor = alpha_parameters[2];
	int max_iter = alpha_parameters[3];
	matrixdb error_mat(Nx, 1);

	// Init loop variables
	bool success = false; size_t counter = 0;
	double alpha = 1.0;
	while (!success) {

		// Init lists
		vector<vectordb> list_x_star(list_x_), list_u_star(list_u_);
		vector<vectordb> list_eq_eval, list_ineq_eval;
		vectordb teq_eval, tineq_eval;

		// Reserve space
		list_x_star.reserve(N);
		list_u_star.reserve(N);
		list_eq_eval.reserve(N);
		list_ineq_eval.reserve(N);

		// Rollout
		double cost = 0; double z = 1e15;
		for (size_t i = 0; i < N; i++) {

			// Get state error
			vectordb error = list_x_star[i] - list_x[i];
			error_mat.setcol(0, error);

			// Retrieve gains
			matrixdb k = list_k_[N - 1 - i];
			matrixdb K = list_K_[N - 1 - i];

			// Get control correction
			vectordb correction = vectordb((alpha * k + K * error_mat).getcol(0));
			list_u_star.push_back(list_u[i] + correction);

			// Identity vector building
			vectordb x_star = list_x_star[i];
			vectordb u_star = list_u_star[i];

			// Evaluate dynamics
			vectordb dynamic_eval = dynamics_.dynamic_db()(
				x_star, u_star,
				spacecraft_parameters_, dynamics_.constants(), solver_parameters_);

			// Evaluate AUL ctg and store constraints
			vectordb constraints = get_AUL_cost_to_go(
				x_star, u_star, i);
			if (Neq == 0)
				list_eq_eval.emplace_back();
			else
				list_eq_eval.emplace_back(constraints.extract(0, Neq - 1));
			if (Nineq == 0)
				list_ineq_eval.emplace_back();
			else
				list_ineq_eval.emplace_back(constraints.extract(Neq, Neq + Nineq - 1));
			double ctg_eval = constraints[Neq + Nineq];

			// Update cost and append list_x_star
			list_x_star.emplace_back(dynamic_eval);
			cost += ctg_eval;
		}

		// AUL terminal cost and store constraints
		vectordb x_star = list_x_star[N];
		vectordb constraints = get_AUL_terminal_cost(
			x_star, x_goal);
		if (Nteq == 0)
			teq_eval = vectordb(0);
		else
			teq_eval = constraints.extract(0, Nteq - 1);
		if (Ntineq == 0)
			tineq_eval = vectordb(0);
		else
			tineq_eval = constraints.extract(Nteq, Nteq + Ntineq - 1);
		double terminal_cost_eval = constraints[Nteq + Ntineq];

		// Update cost 
		cost += terminal_cost_eval;

		// Get expected cost
		double expected_cost = expected_cost_(alpha);

		// If step too small
		if (expected_cost < tol*tol) { // Hard coded: TO DO
			list_x_ = list_x;
			list_u_ = list_u;
			increase_regularisation_();
			break;
		}
		else {
			z = (cost_ - cost) / expected_cost;

			// Check z \in interval
			if (z > z_min && z < z_max) {
				// Save all double lists
				list_x_ = list_x_star;
				list_u_ = list_u_star;
				list_eq_ = list_eq_eval;
				list_ineq_ = list_ineq_eval;
				teq_ = teq_eval;
				tineq_ = tineq_eval;
				cost_ = cost;

				// Recompute DA dynamics, ctg, and tc
				// Make DA lists
				for (size_t i = 0; i < N; i++) {

					// Declare variables
					vectorDA x_star_DA, u_star_DA;

					// Identity vectors building
					x_star_DA = id_vector(list_x_star[i], 0, 0, Nx - 1);
					u_star_DA = id_vector(list_u_star[i], 0, Nx, Nu + Nx);

					// Evaluate dynamics from sractch
					list_dynamic_eval_[i] = dynamics_.dynamic()(
						x_star_DA, u_star_DA,
						spacecraft_parameters_, dynamics_.constants(),
						solver_parameters_);

					// Get constraints
					vectorDA constraints(get_AUL_cost_to_go(
						x_star_DA, u_star_DA, i));
					DA ctg_eval = constraints[Neq + Nineq];
					list_ctg_eval_[i] = ctg_eval;
				}

				// AUL terminal cost and store constraints
				vectorDA x_star_DA = id_vector(list_x_[N], 0, 0, Nx - 1);
				vectorDA constraints = get_AUL_terminal_cost(
					x_star_DA, x_goal);
				tc_eval_ = constraints[Nteq + Ntineq];
				
				break;
			}
			else {
				// Decrease line search parameter
				alpha *= alpha_factor;

				// Check iteration number
				if (counter > max_iter) {
					list_x_ = list_x;
					list_u_ = list_u;
					break;
				}
			}
		}
		counter++;
	}
}

// Performs the DDP forward pass, that consists in the computation
// of the new states and control after correction.
// Inspired from ALTRO (Julia).
// DA only for automatic differentiation.
// // The linesearch is tweaked to implement a memory from one iteration to the other.
// See: https://github.com/RoboticExplorationLab/Altro.jl
void DDPSolver::forward_pass_ls_(
	vector<vectordb> const& list_x, vector<vectordb> const& list_u, vectordb const& x_goal) {
	// Unpack parameters
	double tol = solver_parameters_.DDP_tol();
	unsigned int N = solver_parameters_.N();
	unsigned int Nx = solver_parameters_.Nx();
	unsigned int Nu = solver_parameters_.Nu();
	unsigned int Neq = solver_parameters_.Neq();
	unsigned int Nineq = solver_parameters_.Nineq();
	unsigned int Nteq = solver_parameters_.Nteq();
	unsigned int Ntineq = solver_parameters_.Ntineq();
	vectordb alpha_parameters = solver_parameters_.line_search_parameters();
	double z_min = alpha_parameters[0];
	double z_max = alpha_parameters[1];
	double alpha_factor = alpha_parameters[2];
	int max_iter = alpha_parameters[3];
	matrixdb error_mat(Nx, 1);

	// Init loop variables
	bool success = false; size_t counter = 0;
	while (!success) {

		// Init lists
		vector<vectordb> list_x_star(list_x_), list_u_star(list_u_);
		vector<vectordb> list_eq_eval, list_ineq_eval;
		vectordb teq_eval, tineq_eval;

		// Reserve space
		list_x_star.reserve(N);
		list_u_star.reserve(N);
		list_eq_eval.reserve(N);
		list_ineq_eval.reserve(N);

		// Rollout
		double cost = 0; double z = 1e15;
		for (size_t i = 0; i < N; i++) {

			// Get state error
			vectordb error = list_x_star[i] - list_x[i];
			error_mat.setcol(0, error);

			// Retrieve gains
			matrixdb k = list_k_[N - 1 - i];
			matrixdb K = list_K_[N - 1 - i];

			// Get control correction
			vectordb correction = vectordb((alpha_ * k + K * error_mat).getcol(0));
			list_u_star.push_back(list_u[i] + correction);

			// Identity vector building
			vectordb x_star = list_x_star[i];
			vectordb u_star = list_u_star[i];

			// Evaluate dynamics
			vectordb dynamic_eval = dynamics_.dynamic_db()(
				x_star, u_star,
				spacecraft_parameters_, dynamics_.constants(), solver_parameters_);

			// Evaluate AUL ctg and store constraints
			vectordb constraints = get_AUL_cost_to_go(
				x_star, u_star, i);
			if (Neq == 0)
				list_eq_eval.emplace_back();
			else
				list_eq_eval.emplace_back(constraints.extract(0, Neq - 1));
			if (Nineq == 0)
				list_ineq_eval.emplace_back();
			else
				list_ineq_eval.emplace_back(constraints.extract(Neq, Neq + Nineq - 1));
			double ctg_eval = constraints[Neq + Nineq];

			// Update cost and append list_x_star
			list_x_star.emplace_back(dynamic_eval);
			cost += ctg_eval;
		}

		// AUL terminal cost and store constraints
		vectordb x_star = list_x_star[N];
		vectordb constraints = get_AUL_terminal_cost(
			x_star, x_goal);
		if (Nteq == 0)
			teq_eval = vectordb(0);
		else
			teq_eval = constraints.extract(0, Nteq - 1);
		if (Ntineq == 0)
			tineq_eval = vectordb(0);
		else
			tineq_eval = constraints.extract(Nteq, Nteq + Ntineq - 1);
		double terminal_cost_eval = constraints[Nteq + Ntineq];

		// Update cost 
		cost += terminal_cost_eval;

		// Get expected cost
		double expected_cost = expected_cost_(alpha_);

		// If step too small
		if (expected_cost < tol * tol) {
			list_x_ = list_x;
			list_u_ = list_u;
			increase_regularisation_();
			break;
		}
		else {
			z = (cost_ - cost) / expected_cost;

			// Check z \in interval
			if (z > z_min && z < z_max) {
				// Save all double lists
				list_x_ = list_x_star;
				list_u_ = list_u_star;
				list_eq_ = list_eq_eval;
				list_ineq_ = list_ineq_eval;
				teq_ = teq_eval;
				tineq_ = tineq_eval;
				cost_ = cost;

				// Recompute DA dynamics, ctg, and tc
				// Make DA lists
				for (size_t i = 0; i < N; i++) {

					// Declare variables
					vectorDA x_star_DA, u_star_DA;

					// Identity vectors building
					x_star_DA = id_vector(list_x_star[i], 0, 0, Nx - 1);
					u_star_DA = id_vector(list_u_star[i], 0, Nx, Nu + Nx);

					// Evaluate dynamics from sractch
					list_dynamic_eval_[i] = dynamics_.dynamic()(
						x_star_DA, u_star_DA,
						spacecraft_parameters_, dynamics_.constants(),
						solver_parameters_);

					// Get constraints
					vectorDA constraints(get_AUL_cost_to_go(
						x_star_DA, u_star_DA, i));
					DA ctg_eval = constraints[Neq + Nineq];
					list_ctg_eval_[i] = ctg_eval;
				}

				// AUL terminal cost and store constraints
				vectorDA x_star_DA = id_vector(list_x_[N], 0, 0, Nx - 1);
				vectorDA constraints = get_AUL_terminal_cost(
					x_star_DA, x_goal);
				tc_eval_ = constraints[Nteq + Ntineq];

				// Update linesearch parameter
				alpha_ = min(alpha_ / alpha_factor, 1.0);
				break;
			}
			else {
				// Decrease line search parameter
				alpha_ *= alpha_factor;

				// Check iteration number
				if (counter > max_iter) {
					list_x_ = list_x;
					list_u_ = list_u;
					break;
				}
			}
		}
		counter++;
	}
}

// Performs the DDP forward pass, that consists in the computation
// of the new states and control after correction using the DA mapping
// The linesearch is tweaked to implement a memory from one iteration to the other.
// The linesearch computation are done with floats (in DA vectors), the DA mappings are computed
// only when the linesearch is ended.
// Inspired from ALTRO (Julia).
// See: https://github.com/RoboticExplorationLab/Altro.jl
void DDPSolver::forward_pass_ls_DA_(
	vector<vectordb> const& list_x, vector<vectordb> const& list_u, vectordb const& x_goal) {
	// Unpack parameters
	double tol = solver_parameters_.DDP_tol();
	double tol_DA = solver_parameters_.AUL_tol(); // TEST
	unsigned int N = solver_parameters_.N();
	unsigned int Nx = solver_parameters_.Nx();
	unsigned int Nu = solver_parameters_.Nu();
	unsigned int Neq = solver_parameters_.Neq();
	unsigned int Nineq = solver_parameters_.Nineq();
	unsigned int Nteq = solver_parameters_.Nteq();
	unsigned int Ntineq = solver_parameters_.Ntineq();
	vectordb alpha_parameters = solver_parameters_.line_search_parameters();
	double z_min = alpha_parameters[0];
	double z_max = alpha_parameters[1];
	double alpha_factor = alpha_parameters[2];
	int max_iter = alpha_parameters[3];
	matrixdb error_mat(Nx, 1);

	// Init loop variables
	bool success = false; size_t counter = 0;
	while (!success) {

		// Init lists
		vector<vectordb> list_x_star(list_x_), list_u_star(list_u_);
		vector<vectordb> list_error, list_correction;
		vector<bool> list_condition_radius;
		
		// Reserve space
		list_x_star.reserve(N);
		list_u_star.reserve(N);
		list_condition_radius.reserve(N);
		list_error.reserve(N);
		list_correction.reserve(N);

		// Rollout
		double cost = 0; double z = 1e15;
		for (size_t i = 0; i < N; i++) {

			// Get state error
			vectordb error = list_x_star[i] - list_x[i];
			list_error.push_back(error);
			error_mat.setcol(0, error);
			double error_norm = error.vnorm();

			// Retrieve gains
			matrixdb k = list_k_[N - 1 - i];
			matrixdb K = list_K_[N - 1 - i];

			// Get control correction
			vectordb correction = vectordb((alpha_ * k + K * error_mat).getcol(0));
			list_correction.push_back(correction);
			double correction_norm = correction.vnorm();

			// Get dynamics convergence radius
			vectorDA dynamic_eval_prev = list_dynamic_eval_[i];
			double dynamics_radius = convRadius(dynamic_eval_prev, tol_DA);

			// Declare variables
			vectordb dynamic_eval, x_star, u_star;

			// Check condition
			list_condition_radius.push_back(
				max(correction_norm, error_norm) < dynamics_radius);

			// Checks if [dx, du] belongs to the convergence radius of the dynamics
			// Thus, if they can be approximated
			if (list_condition_radius[i]) {

				// Get DA perturbed vector
				vectordb dx_star = error;
				vectordb du_star = correction;

				// Make identity vectors
				x_star = list_x[i] + error;
				u_star = list_u[i] + correction;
				list_u_star.emplace_back(u_star.cons());

				// Concatenate two vectors to get [dx, du]
				dx_star.reserve(Nu);
				for (size_t k = 0; k < Nu; k++) { dx_star.push_back(du_star[k]); }

				// Evaluate previous evaluation at [dx, du]
				dynamic_eval = dynamic_eval_prev.eval(dx_star);
			}
			else {

				// Identity vectors building
				list_u_star.emplace_back(list_u[i] + correction);
				x_star = list_x_star[i];
				u_star = list_u_star[i];

				// Evaluate dynamics from sractch
				dynamic_eval = dynamics_.dynamic_db()(
					x_star, u_star,
					spacecraft_parameters_, dynamics_.constants(),
					solver_parameters_);
			}

			// Get constraints
			vectordb constraints(get_AUL_cost_to_go(
				x_star, u_star, i));
			double ctg_eval = constraints[Neq + Nineq];

			// Update cost and append list_x_star
			list_x_star.emplace_back(dynamic_eval);
			cost += ctg_eval;
		}

		// AUL terminal cost and store constraints
		vectordb x_star = list_x_star[N];
		vectordb constraints = get_AUL_terminal_cost(
			x_star, x_goal);

		// Update cost 
		cost += constraints[Nteq + Ntineq];

		// Get expected cost
		double expected_cost = expected_cost_(alpha_);

		// If step too small
		if (expected_cost < tol*tol) {
			list_x_ = list_x;
			list_u_ = list_u;
			increase_regularisation_();
			break;
		}
		else {
			z = (cost_ - cost) / expected_cost;

			// Check z \in interval
			if (z > z_min && z < z_max) {

				// Init lists
				list_ctg_eval_ = vectorDA();
				list_eq_ = vector<vectordb>();
				list_ineq_ = vector<vectordb>();

				// Reserve space
				list_ctg_eval_.reserve(N);
				list_eq_.reserve(N);
				list_ineq_.reserve(N);

				// Rollout
				cost_ = 0; 
				double z = 1e15;

				// Make DA lists
				for (size_t i = 0; i < N; i++) {

					// Get control correction
					vectordb correction = list_correction[i];

					// Declare variables
					vectorDA dynamic_eval, x_star_DA, u_star_DA;

					// Checks if [dx, du] belongs to the convergence radius of the dynamics
					// Thus, if they can be approximated
					if (list_condition_radius[i]) {

						// Get dynamics convergence radius
						vectorDA dynamic_eval_prev = list_dynamic_eval_[i];

						// Get DA perturbed vector
						vectorDA dx_star_DA = id_vector(list_error[i], 0, 0, Nx - 1);
						vectorDA du_star_DA = id_vector(correction, 0, Nx, Nu + Nx);

						// Make identity vectors
						x_star_DA = list_x[i] + dx_star_DA;
						u_star_DA = list_u[i] + du_star_DA;
						list_u_.emplace_back(u_star_DA.cons());

						// Concatenate two vectors to get [dx, du]
						dx_star_DA.reserve(Nu);
						for (size_t k = 0; k < Nu; k++) { dx_star_DA.push_back(du_star_DA[k]); }

						// Evaluate previous evaluation at [dx, du]
						dynamic_eval = dynamic_eval_prev.eval(dx_star_DA);
					}
					else {

						// Identity vectors building
						list_u_.emplace_back(list_u[i] + correction);
						x_star_DA = id_vector(list_x_[i], 0, 0, Nx - 1);
						u_star_DA = id_vector(list_u_[i], 0, Nx, Nu + Nx);

						// Evaluate dynamics from sractch
						dynamic_eval = dynamics_.dynamic()(
							x_star_DA, u_star_DA,
							spacecraft_parameters_, dynamics_.constants(),
							solver_parameters_);
					}
					list_dynamic_eval_[i] = dynamic_eval;

					// Get constraints
					vectorDA constraints(get_AUL_cost_to_go(
						x_star_DA, u_star_DA, i));
					if (Neq == 0)
						list_eq_.emplace_back();
					else
						list_eq_.emplace_back(constraints.extract(0, Neq - 1).cons());
					if (Nineq == 0)
						list_ineq_.emplace_back();
					else
						list_ineq_.emplace_back(constraints.extract(Neq, Neq + Nineq - 1).cons());
					DA ctg_eval = constraints[Neq + Nineq];
					list_ctg_eval_.push_back(ctg_eval);

					// Update cost and append list_x_star
					list_x_.emplace_back(dynamic_eval.cons());
					cost_ += ctg_eval.cons();
				}

				// AUL terminal cost and store constraints
				vectorDA x_star_DA = id_vector(list_x_[N], 0, 0, Nx - 1);
				vectorDA constraints = get_AUL_terminal_cost(
					x_star_DA, x_goal);
				if (Nteq == 0)
					teq_ = vectordb(0);
				else
					teq_ = constraints.extract(0, Nteq - 1).cons();
				if (Ntineq == 0)
					tineq_ = vectordb(0);
				else
					tineq_ = constraints.extract(Nteq, Nteq + Ntineq - 1).cons();
				tc_eval_ = constraints[Nteq + Ntineq];

				// Update cost 
				cost_ += tc_eval_.cons();

				// Update linesearch parameter
				alpha_ = min(alpha_ / alpha_factor, 1.0);
				break;
			}
			else {
				// Decrease line search parameter
				alpha_ *= alpha_factor;

				// Check iteration number
				if (counter > max_iter) {
					list_x_ = list_x;
					list_u_ = list_u;
					break;
				}
			}
		}
		counter++;
	}
}

// Performs DDP solving given a starting point,
// initial controls and a final state.
// Inspired from ALTRO (Julia).
// See: https://github.com/RoboticExplorationLab/Altro.jl
void DDPSolver::solve(
	vectordb const& x0,
	vector<vectordb> const& list_u_init,
	vectordb const& x_goal) {
	// Unpack parameters
	double tol = solver_parameters_.DDP_tol();
	unsigned int max_iter = solver_parameters_.DDP_max_iter();
	unsigned int N = solver_parameters_.N();
	unsigned int Nx = solver_parameters_.Nx();
	unsigned int Nu = solver_parameters_.Nu();
	unsigned int Neq = solver_parameters_.Neq();
	unsigned int Nineq = solver_parameters_.Nineq();
	unsigned int Nteq = solver_parameters_.Nteq();
	unsigned int Ntineq = solver_parameters_.Ntineq();
	unsigned int DDP_type = solver_parameters_.DDP_type();
	bool bs_reg = solver_parameters_.backward_sweep_regulation();
	unsigned int verbosity = solver_parameters_.verbosity();
	unsigned int saving_iterations = solver_parameters_.saving_iterations();

	// Init regularisation
	rho_ = solver_parameters_.backward_sweep_regulation_parameters()[0];
	d_rho_ = 0.0;
	
	// Init lists
	list_ctg_eval_ = vector<DA>();
	list_eq_ = vector<vectordb>(); list_ineq_ = vector<vectordb>();

	// Reserve
	list_ctg_eval_.reserve(N);
	list_eq_.reserve(N);
	list_ineq_.reserve(N);

	// Separte dynamic
	auto start = high_resolution_clock::now();
	if (recompute_dynamics_) {
		list_x_ = vector<vectordb>();
		list_dynamic_eval_ = vector<vectorDA>();
		list_dynamic_eval_.reserve(N);
		list_x_.reserve(N + 1);
		list_x_.push_back(x0);
	}

	// Make first trajectory evaluation.

	// Init cost
	cost_ = 0; n_iter_ = 0.0;
	for (size_t i = 0; i < N; i++) {
		// Make DA maps
		vectorDA x_da = id_vector(list_x_[i], 0, 0, Nx - 1);
		vectorDA u_da = id_vector(list_u_init[i], 0, Nx, Nx + Nu);

		// Get x_k+1
		if (recompute_dynamics_) {
			vectorDA x_kp1 = dynamics_.dynamic()(x_da, u_da,
				spacecraft_parameters_, dynamics_.constants(),
				solver_parameters_);
			list_x_.emplace_back(x_kp1.cons());
			list_dynamic_eval_.push_back(x_kp1);
		}

		// Get cost to go

		// Evaluate AUL ctg
		vectorDA ctg_constraint = get_AUL_cost_to_go(
			x_da, u_da, i);
		if (Neq == 0)
			list_eq_.emplace_back(0);
		else
			list_eq_.emplace_back(ctg_constraint.extract(0, Neq - 1).cons());
		if (Nineq == 0)
			list_ineq_.emplace_back(0);
		else
			list_ineq_.emplace_back(ctg_constraint.extract(Neq, Neq + Nineq - 1).cons());
		DA ctg_eval = ctg_constraint[Neq + Nineq];
		list_ctg_eval_.push_back(ctg_eval);

		// Update cost
		cost_ += ctg_eval.cons();
	}
	list_u_ = list_u_init;
	
	// Evaluate AUL tc
	vectorDA tc_constraint = get_AUL_terminal_cost(
		id_vector(list_x_[N], 0, 0, Nx - 1), x_goal);
	if (Nteq == 0)
		teq_ = vectordb(0);
	else
		teq_ = tc_constraint.extract(0, Nteq - 1).cons();
	if (Ntineq == 0)
		tineq_ = vectordb(0);
	else
		tineq_ = tc_constraint.extract(Nteq, Nteq + Ntineq - 1).cons();
	tc_eval_ = tc_constraint[Nteq + Ntineq];

	// Update cost
	cost_ += tc_eval_.cons();

	// Compute costs
	double tc = tc_eval_.cons();
	double ctg = (cost_ - tc_eval_.cons());

	// Output
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	if (verbosity < 1)
		cout << "    " << 0 << " - " << cost_ << ", " << ctg << ", " << tc 
		<< ", "<< to_string(static_cast<double>(duration.count()) / 1e6) << endl;
	
	// DDP solving
	double cost_last = 0.0;
	vector<vectordb> list_x(list_x_), list_u(list_u_init);
	bool loop = true;
	while (loop) {
		// Update looping variables
		n_iter_++;
		cost_last = cost_;

		int nb_method = 3;

		// Get times
		auto start = high_resolution_clock::now();

		// Backward sweep (< nb_method =  classic method from ALTRO)
		if (DDP_type < nb_method)
			backward_sweep_();

		// Backward sweep (< 2*nb_method = classic method from ALTRO, with hessian)
		else if (DDP_type < 2*nb_method)
			backward_sweep_hessian_();

		// Backward sweep (< 3*nb_method = Qk DA direct evaluation)
		else
			backward_sweep_DA_Q_();

		// Forward pass init
		list_x_ = vector<vectordb>(1, x0);
		list_u_ = vector<vectordb>();

		// Forward pass (0 = classic method from ALTRO)
		if (DDP_type % nb_method == 0)
			forward_pass_(list_x, list_u, x_goal);

		// Forward pass (1 = tweaked linesearch)
		else if (DDP_type % nb_method == 1)
			forward_pass_ls_(list_x, list_u, x_goal);

		// Forward pass (2 = use of the DA mappings the dynamics + Linesearch)
		else if (DDP_type % nb_method == 2)
			forward_pass_ls_DA_(list_x, list_u, x_goal);

		// Compute costs
		tc = tc_eval_.cons();
		ctg = (cost_ - tc_eval_.cons());

		// Evaluate convergence 
		double d_cost = (cost_last - cost_)/max(abs(cost_last), abs(cost_));
		loop = !evaluate_convergence_(d_cost);

		// Save iterations
		if (saving_iterations > 2) {
			list_x_mem_.push_back(list_x_);
			list_u_mem_.push_back(list_u_);
		}

		// Update states and control
		list_x = list_x_; list_u = list_u_;
		
		// Output
		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(stop - start);
		if (verbosity < 1)
			cout << "    " << n_iter_ 
				<< " - " << cost_ << ", " << ctg
				<< ", " << tc << ", "<< to_string(static_cast<double>(duration.count()) / 1e6) << endl;
	}
}

// Evaluates the convergence radius of a DA vector.
double convRadius(vectorDA const& x, double const& tol) {
	size_t n = x.size(); double radius = 1e15;
	for (size_t i = 0; i < n; i++) {
		double value = convRadius(x[i], tol);
		if (value < radius)
			radius = value;
	}
	return radius;
}
