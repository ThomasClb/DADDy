/**
	pn_solver.cpp

	Purpose: Implementation of the PNSolver class.

	@author Thomas Caleb

	@version 1.0 17/11/2023
*/

#include "pn_solver.h"

using namespace DACE;
using namespace std;
using namespace std::chrono;

// Empty constructors
PNSolver::PNSolver() : AULsolver_(AULSolver()),
	list_x_(vector<vectordb>(0)), list_u_(vector<vectordb>(0)),
	cost_(0),
	list_der_cost_(vector<vector<matrixdb>>(0)),
	list_eq_(vector<vectordb>(0)), list_ineq_(vector<vectordb>(0)),
	teq_(vectordb(0)), tineq_(vectordb(0)),
	list_der_eq_(vector<vector<matrixdb>>(0)),
	list_der_ineq_(vector<vector<matrixdb>>(0)),
	der_teq_(vector<matrixdb>(0)), der_tineq_(vector<matrixdb>(0)) {}

// Constructors
PNSolver::PNSolver(AULSolver const& AULsolver) : AULsolver_(AULsolver),
	list_x_(AULsolver.list_x()), list_u_(AULsolver.list_u()),
	cost_(AULsolver.cost()),
	list_der_cost_(vector<vector<matrixdb>>(AULsolver.list_eq().size() + 1)),
	list_eq_(AULsolver.list_eq()), list_ineq_(AULsolver.list_ineq()),
	teq_(AULsolver.teq()), tineq_(AULsolver.tineq()),
	list_der_eq_(vector<vector<matrixdb>>(AULsolver.list_eq().size())),
	list_der_ineq_(vector<vector<matrixdb>>(AULsolver.list_ineq().size())),
	der_teq_(vector<matrixdb>(0)), der_tineq_(vector<matrixdb>(0)) {}

// Copy constructor
PNSolver::PNSolver(
	PNSolver const& solver) : AULsolver_(solver.AULsolver_),
	list_x_(solver.list_x_), list_u_(solver.list_u_),
	cost_(solver.cost_),
	list_der_cost_(vector<vector<matrixdb>>(solver.list_eq_.size() + 1)),
	list_eq_(solver.list_eq_), list_ineq_(solver.list_ineq_),
	teq_(solver.teq_), tineq_(solver.tineq_),
	list_der_eq_(vector<vector<matrixdb>>(solver.list_eq_.size())),
	list_der_ineq_(vector<vector<matrixdb>>(solver.list_ineq_.size())),
	der_teq_(vector<matrixdb>(0)), der_tineq_(vector<matrixdb>(0)) {}

// Destructors
PNSolver::~PNSolver() {}

// Getters
const AULSolver PNSolver::AULsolver() const { return AULsolver_; }
const DDPSolver PNSolver::DDPsolver() const { return AULsolver_.DDPsolver(); }
const vector<vectordb> PNSolver::list_x() const { return list_x_; }
const vector<vectordb> PNSolver::list_u() const { return list_u_; }
const vector<vectordb> PNSolver::list_eq() const { return list_eq_; }
const vector<vectordb> PNSolver::list_ineq() const { return list_ineq_; }
const vectordb PNSolver::teq() const { return teq_; }
const vectordb PNSolver::tineq() const { return tineq_; }
const double PNSolver::cost() const { return cost_; }

// Setters
void PNSolver::set_list_x(vector<vectordb> const& list_x) { list_x_ = list_x;}
void PNSolver::set_list_u(vector<vectordb> const& list_u) { list_u_ = list_u;}

// Solves the optimisation problem with a projected Newton method
// Inspired from ALTRO (Julia).
// See: https://github.com/RoboticExplorationLab/Altro.jl
void PNSolver::solve(vectordb const& x_goal) {
	// Unpack
	DDPSolver ddp_solver = DDPsolver();
	SolverParameters solver_parameters = ddp_solver.solver_parameters();
	unsigned int N = solver_parameters.N();
	unsigned int Nx = solver_parameters.Nx();
	unsigned int Nu = solver_parameters.Nu();
	unsigned int Neq = solver_parameters.Neq();
	unsigned int Nineq = solver_parameters.Nineq();
	unsigned int Nteq = solver_parameters.Nteq();
	unsigned int Ntineq = solver_parameters.Ntineq();
	size_t max_iter = solver_parameters.PN_max_iter();
	double constraint_tol = solver_parameters.PN_tol();
	double cv_rate_threshold = solver_parameters.PN_cv_rate_threshold();
	unsigned int verbosity = solver_parameters.verbosity();
	Constants constants = ddp_solver.dynamics().constants();

	// Output
	if (verbosity < 1) {
		cout << "############################################################################" << endl;
		cout << "#                                                                          #" << endl;
		cout << "#                              START PN SOLVING                            #" << endl;
		cout << "#                                                                          #" << endl;
		cout << "############################################################################" << endl << endl << endl;
	}
	else if (verbosity < 2) {
		cout << endl << "PN solving" << endl;
	}

	// Evaluate constraints
	update_constraints_(x_goal);
	double prev_violation = get_max_constraint_();

	// Init loop
	double violation(prev_violation);
	for (size_t i = 0; i < max_iter; i++) {
		// Output
		if (verbosity < 1) {
			if (i % 5 == 0)
				cout << i << " - " << prev_violation << " - " << list_x_[N][SIZE_VECTOR] * constants.massu() << endl;
		}

		// Check termination constraints
		if (violation < constraint_tol)
			break;

		// Build active constraint vector and its gradient
		if (i != 0)
			update_constraints_(x_goal);
		pair<vectordb, pair<vector<matrixdb>, tri_vector_size_t>> d_block_D = get_d_block_D_();
		vectordb d = d_block_D.first;
		vector<matrixdb> block_D = d_block_D.second.first;

		// Make Sigma = D * D^t as a tridiagonal matrix
		sym_tridiag_matrixdb tridiag_sigma = get_block_sigma_sq_(
			block_D, d_block_D.second.second);

		// Compute block tridiagonal Cholesky factorisation of Sigma.
		sym_tridiag_matrixdb tridiag_L_sigma = cholesky_(tridiag_sigma);

		// Line search loop
		double cv_rate = 1e15; double violation_mem = prev_violation;
		for (size_t j = 0; j < 10; j++) {
			// Termination checks
			if (violation < constraint_tol || cv_rate < cv_rate_threshold)
				break;

			// Line search
			violation = line_search_(x_goal, tridiag_L_sigma, block_D, d, violation);

			// Update cv_rate
			cv_rate = log(violation) / log(prev_violation);
			prev_violation = violation;
		}
	}

	return;
}

// Iterative line search for PN.
// Inspired from ALTRO (Julia).
// See: https://github.com/RoboticExplorationLab/Altro.jl
double PNSolver::line_search_(
	vectordb const& x_goal,
	sym_tridiag_matrixdb const& tridiag_L,
	vector<matrixdb> const& block_D,
	vectordb const& d_0,
	double const& violation_0) {
	// Unpack
	DDPSolver ddp_solver = DDPsolver();
	SolverParameters solver_parameters = ddp_solver.solver_parameters();
	double constraint_tol = solver_parameters.PN_tol();
	size_t max_iter = solver_parameters.PN_max_iter();
	double alpha = solver_parameters.PN_alpha();
	double gamma = solver_parameters.PN_gamma();
	unsigned int N = solver_parameters.N();
	unsigned int Nx = solver_parameters.Nx();
	unsigned int Nu = solver_parameters.Nu();
	matrixdb mat_z(d_0.size(), 1);

	// Init loope
	double violation(violation_0); vectordb d(d_0);
	vector<vectordb> list_x = list_x;
	vector<vectordb> list_u = list_u_;
	for (size_t i = 0; i < 10; i++) {
		// Solve Sigma * z = d <=> L * L^t * z = d using
		// block tridiagonal Cholesky factorisation.
		vectordb z = solve_cholesky_(tridiag_L, d);
		mat_z.setcol(0, z);	
		vectordb corr(N * (Nu + Nx));

		// Compute correction
		size_t counter_z = 0; size_t counter_corr = 0;
		for (size_t j = 0; j < block_D.size(); j++) {
			vectordb z_j;
			matrixdb D_j = block_D[j];
			size_t N_j = D_j.nrows();
			
			if (N_j != 0)
				z_j = z.extract(counter_z, N_j + counter_z - 1);
				
			matrixdb mat_z_j(z_j.size(), 1); mat_z_j.setcol(0, z_j);
			vectordb corr_j = (mat_z_j.transpose() * D_j ).getrow(0);

			// Assign
			for (size_t k = 0; k < corr_j.size(); k++) {
				corr[counter_corr + k] += alpha * corr_j[k];
			}
			counter_z += N_j;
			counter_corr += corr_j.size() - Nx;
		}			

		// Update x, u
		pair<vector<vectordb>,
			 vector<vectordb>> list_x_u = update_list_x_u_(corr); // TO DO : modify

		// Evaluate constraints
		pair<
			vector<vectordb>,
			vector<vectordb>> list_eq_ineq = update_constraints_double_( // TO DO : DA method
			x_goal, list_x_u.first, list_x_u.second);

		// Get the max constraint
		violation = get_max_constraint_(
			list_eq_ineq.first, list_eq_ineq.second);

		// Check exit
		if (violation <= violation_0) {
			// Update state
			set_list_x(list_x_u.first); set_list_u(list_x_u.second);;

			// Update constraints
			update_constraints_double_(list_eq_ineq);

			break;
		}
		else {
			// Restore constraints
			violation = violation_0;
			d = d_0;
		}

		// Update
		alpha *= gamma; 
	}
	return violation;
}

// Computes the maximum constraints given eq in ineq constraints
double PNSolver::get_max_constraint_(
	vector<vectordb> const& list_eq,
	vector<vectordb> const& list_ineq) {
	// Unpack parameters
	DDPSolver ddp_solver = DDPsolver();
	SolverParameters solver_parameters = ddp_solver.solver_parameters();
	unsigned int N = solver_parameters.N();
	unsigned int Nx = solver_parameters.Nx();
	unsigned int Neq = solver_parameters.Neq();
	unsigned int Nineq = solver_parameters.Nineq();
	unsigned int Nteq = solver_parameters.Nteq();
	unsigned int Ntineq = solver_parameters.Ntineq();

	// Loop on all steps
	double max = -1e15;
	for (size_t i = 0; i < N; i++) {

		// Unpack
		vectordb eq_i = list_eq[i];
		vectordb ineq_i = list_ineq[i];

		// Find max
		for (size_t j = 0; j < Neq + Nx; j++) {
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
	vectordb teq = list_eq[N];
	vectordb tineq = list_ineq[N];
	for (size_t j = 0; j < Nteq; j++) {
		double abs_teq_j = abs(teq[j]);
		if (abs_teq_j > max)
			max = abs_teq_j;
	}
	for (size_t j = 0; j < Ntineq; j++) {
		double tineq_j = tineq[j];
		if (tineq_j > max)
			max = tineq_j;
	}

	return max;
}

// Computes the maximum constraints using the attributes
double PNSolver::get_max_constraint_() {
	// Unpack parameters
	DDPSolver ddp_solver = DDPsolver();
	SolverParameters solver_parameters = ddp_solver.solver_parameters();
	unsigned int N = solver_parameters.N();
	unsigned int Nx = solver_parameters.Nx();
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
		for (size_t j = 0; j < Neq + Nx; j++) {
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

// Updates the values of controls and states given corrections
// Return the list of equalities, and inequalities
pair<
	vector<vectordb>,
	vector<vectordb>> PNSolver::update_list_x_u_(
	vectordb const& correction) {
	// Unpack
	DDPSolver ddp_solver = DDPsolver();
	SolverParameters solver_parameters = ddp_solver.solver_parameters();
	SpacecraftParameters spacecraft_parameters = ddp_solver.spacecraft_parameters();
	unsigned int N = solver_parameters.N();
	unsigned int Nx = solver_parameters.Nx();
	unsigned int Nu = solver_parameters.Nu();

	// Init
	vector<vectordb> list_x(list_x_), list_u(list_u_);

	// First control
	vectordb du_i = correction.extract(
		0, Nu - 1);
	list_u[0] -= du_i;

	// Last state
	vectordb dx_i = correction.extract(
		Nu + (N - 1) * (Nu + Nx), correction.size() - 1);
	list_x[N] -= dx_i;

	// Loop on all states in [1, N-1]
	for (size_t i = 1; i < N; i++) {
		// Add correction
		vectordb dx_i = correction.extract(
			Nu + (i - 1) * (Nu + Nx), i * (Nu + Nx) - 1);
		vectordb du_i = correction.extract(
			i * (Nu + Nx), i * (Nu + Nx) + Nu - 1);
		list_x[i] -= dx_i;
		list_u[i] -= du_i;
	}

	// Return
	return pair<vector<vectordb>, vector<vectordb>>(list_x, list_u);
}

// Computes the new constraints given states and controls
void PNSolver::update_constraints_(
	vectordb const& x_goal) {
	// Unpack
	DDPSolver ddp_solver = DDPsolver();
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

	// Make id matrix
	matrixdb id_Nx(Neq + Nx, Nx, 0.0), null_Nx(Nineq, Nx, 0.0);
	for (size_t i = 0; i < Nx; i++) { id_Nx.at(Neq + i, i) = -1.0; }

	// Update path constraints

	// Loop on all steps
	list_dynamics_ = vector<vectorDA>(); list_dynamics_.reserve(N);
	for (size_t i = 0; i < N; i++) {

		// Get DA x, u
		vectorDA x_DA = id_vector(list_x_[i], 0, 0, Nx - 1);
		vectorDA u_DA = id_vector(list_u_[i], 0, Nx, Nu + Nx);

		// Constraints evaluations
		vectorDA eq_eval = dynamics.equality_constraints()(
			x_DA, u_DA, spacecraft_parameters, dynamics.constants(), solver_parameters);
		vectorDA ineq_eval = dynamics.inequality_constraints()(
			x_DA, u_DA, spacecraft_parameters, dynamics.constants(), solver_parameters);

		// Continuity constraints
		vectorDA x_kp1_eval = dynamics.dynamic()(
			x_DA, u_DA,
			spacecraft_parameters, dynamics.constants(), solver_parameters); // TO DO : test radius + evaluate 
		list_dynamics_.push_back(x_kp1_eval);
		x_kp1_eval -= list_x_[i + 1];

		// Add continuity constraints
		eq_eval.reserve(Nx);
		for (size_t k = 0; k < Nx; k++) {
			eq_eval.push_back(x_kp1_eval[k]);
		}

		// Get derivatives
		vector<matrixdb> der_eq = deriv_xu(
			eq_eval, Nx, Nu, false);
		vector<matrixdb> der_ineq = deriv_xu(
			ineq_eval, Nx, Nu, false);

		// Dependency with kp1
		der_eq.reserve(1); der_ineq.reserve(1);
		der_eq.push_back(id_Nx);
		der_ineq.push_back(null_Nx);

		// Assign
		list_eq_[i] = eq_eval.cons();
		list_ineq_[i] = ineq_eval.cons();
		list_der_eq_[i] = der_eq;
		list_der_ineq_[i] = der_ineq;
	}

	// Update terminal constraints

	// Get DA x, u
	vectorDA x_DA = id_vector(list_x_[N], 0, 0, Nx - 1);

	// Constraints evaluations
	vectorDA teq_eval = dynamics.terminal_equality_constraints()(
		x_DA, x_goal, spacecraft_parameters, dynamics.constants(), solver_parameters);
	vectorDA tineq_eval = dynamics.terminal_inequality_constraints()(
		x_DA, x_goal, spacecraft_parameters, dynamics.constants(), solver_parameters);

	// Get derivatives
	vector<matrixdb> der_teq = deriv_xu(
		teq_eval, Nx, Nu, false);
	vector<matrixdb> der_tineq = deriv_xu(
		tineq_eval, Nx, Nu, false);

	// Assign
	teq_ = teq_eval.cons();
	tineq_ = tineq_eval.cons();
	der_teq_ = der_teq;
	der_tineq_ = der_tineq;
}

// Computes the new constraints given states and controls without DA
// Return the list of equalities, and inequalities
pair<
	vector<vectordb>,
	vector<vectordb>> PNSolver::update_constraints_double_(
	vectordb const& x_goal,
	vector<vectordb> const& list_x,
	vector<vectordb> const& list_u) {
	// Unpack
	DDPSolver ddp_solver = DDPsolver();
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

	// Init lists
	vector<vectordb> list_eq, list_ineq;
	list_eq.reserve(N + 1);
	list_ineq.reserve(N + 1);

	// Update path constraints

	// Loop on all steps
	for (size_t i = 0; i < N; i++) {

		// Get x, u
		vectordb x = list_x[i];
		vectordb u = list_u[i];

		// Constraints evaluations
		vectordb eq_eval = dynamics.equality_constraints_db()(
			x, u, spacecraft_parameters, dynamics.constants(), solver_parameters);
		vectordb ineq_eval = dynamics.inequality_constraints_db()(
			x, u, spacecraft_parameters, dynamics.constants(), solver_parameters);
		eq_eval.reserve(Nx);

		// Continuity constraints
		vectordb x_kp1_eval = (dynamics.dynamic_db()(
			x, u, spacecraft_parameters, dynamics.constants(), solver_parameters) - list_x[i + 1]); // TO DO : use DA, test

		// Add continuity constraints
		for (size_t k = 0; k < Nx; k++) {
			eq_eval.emplace_back(x_kp1_eval[k]);
		}

		// Assign
		list_eq.emplace_back(eq_eval);
		list_ineq.emplace_back(ineq_eval);
	}

	// Update terminal constraints

	// Get x
	vectordb x = list_x[N];

	// Constraints evaluations
	vectordb teq_eval = dynamics.terminal_equality_constraints_db()(
		x, x_goal, spacecraft_parameters, dynamics.constants(), solver_parameters);
	vectordb tineq_eval = dynamics.terminal_inequality_constraints_db()(
		x, x_goal, spacecraft_parameters, dynamics.constants(), solver_parameters);

	// Assign
	list_eq.emplace_back(teq_eval);
	list_ineq.emplace_back(tineq_eval);

	return pair<
		vector<vectordb>,
		vector<vectordb>>(list_eq, list_ineq);
}

// Updates the constraints without DA
// Given the list of equalities, and inequalities
void PNSolver::update_constraints_double_(
	pair<
		vector<vectordb>,
		vector<vectordb>>  const& list_eq_ineq) {
	// Unpack
	DDPSolver ddp_solver = DDPsolver();
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

	// Update path constraints

	// Loop on all steps
	vector<vectordb> list_eq = list_eq_ineq.first;
	vector<vectordb> list_ineq = list_eq_ineq.second;
	for (size_t i = 0; i < N; i++) {
		// Assign
		list_eq_[i] = list_eq[i];
		list_ineq_[i] = list_ineq[i];
	}

	// Update terminal constraints
	teq_ = list_eq[N];
	tineq_ = list_ineq[N];
}

// Returns the vector of active constraints and their gradients
// first it the active constraints vector
// second is a pair with the list of gradients of constraints first
// second.second is the list of active constraints.
pair<
	vectordb,
	pair<
	vector<matrixdb>,
	tri_vector_size_t>> PNSolver::get_d_block_D_() {
	// Unpack
	DDPSolver ddp_solver = DDPsolver();
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
	double active_constraint_tol = solver_parameters.PN_active_constraint_tol();

	// Init d and D
	vectordb d; vector<matrixdb> block_D;
	tri_vector_size_t list_active_constraint_index;
	d.reserve(N + 1); block_D.reserve(N + 1);
	list_active_constraint_index.reserve(N + 1);

	// Path constraints
	for (size_t i = 0; i < N; i++) {
		// Init lists
		vector<size_t> list_active_eq_index, list_active_ineq_index, list_active_c_index;
		list_active_eq_index.reserve(Neq);
		list_active_ineq_index.reserve(Nineq);
		list_active_c_index.reserve(Nx);

		// Unpack
		vectordb eq_i(list_eq_[i]), ineq_i(list_ineq_[i]);

		// Equality constraints
		for (size_t j = 0; j < Neq; j++) {
			double abs_eq_j = abs(eq_i[j]);
			if (abs_eq_j > active_constraint_tol) {
				// Assign to d
				d.push_back(abs_eq_j);

				// Save index
				list_active_eq_index.push_back(j);
			}
		}

		// Inequality constraints
		for (size_t j = 0; j < Nineq; j++) {
			double ineq_j = ineq_i[j];
			if (ineq_j > active_constraint_tol) {
				// Assign to d
				d.push_back(ineq_j);

				// Save index
				list_active_ineq_index.push_back(j);
			}
		}

		// Continuity constraints
		for (size_t j = Neq; j < Neq + Nx; j++) {
			double eq_j = eq_i[j];
			if (abs(eq_j) > active_constraint_tol) {
				// Assign to d
				d.push_back(eq_j);

				// Save index
				list_active_c_index.push_back(j);
			}
		}

		// Build gradient
		size_t Neq_active(list_active_eq_index.size());
		size_t Nineq_active(list_active_ineq_index.size());
		size_t Nc_active(list_active_c_index.size());
		matrixdb D_i(
			Neq_active + Nineq_active + Nc_active,
			(2 * Nx + Nu));
		vector<matrixdb> der_eq = list_der_eq_[i];
		vector<matrixdb> der_ineq = list_der_ineq_[i];

		// Equality constraints
		for (size_t j = 0; j < list_active_eq_index.size(); j++) {
			// Get derivatives
			size_t active_eq_index_j = list_active_eq_index[j];
			vectordb dx(der_eq[0].getrow(active_eq_index_j));
			vectordb du(der_eq[1].getrow(active_eq_index_j));
			vectordb dxkp1(der_eq[2].getrow(active_eq_index_j));

			// Assign to D_i
			for (size_t k = 0; k < Nx; k++) {
				D_i.at(j, k) = dx[k];
				D_i.at(j, k + Nx + Nu) = dxkp1[k];
			}
			for (size_t k = 0; k < Nu; k++) {
				D_i.at(j, k + Nx) = du[k];
			}
		}

		// Inequality constraints
		for (size_t j = 0; j < list_active_ineq_index.size(); j++) {
			// Get derivatives
			size_t active_ineq_index_j = list_active_ineq_index[j];
			vectordb dx(der_ineq[0].getrow(active_ineq_index_j));
			vectordb du(der_ineq[1].getrow(active_ineq_index_j));
			vectordb dxkp1(der_ineq[2].getrow(active_ineq_index_j));

			// Assign to D_i
			size_t index_j = j + Neq_active;
			for (size_t k = 0; k < Nx; k++) {
				D_i.at(index_j, k) = dx[k];
				D_i.at(index_j, k + Nx + Nu) = dxkp1[k];
			}
			for (size_t k = 0; k < Nu; k++) {
				D_i.at(index_j, k + Nx) = du[k];
			}
		}

		// Continuity constraints
		for (size_t j = 0; j < list_active_c_index.size(); j++) {
			// Get derivatives
			size_t active_c_index_j = list_active_c_index[j];
			vectordb dx(der_eq[0].getrow(active_c_index_j));
			vectordb du(der_eq[1].getrow(active_c_index_j));
			vectordb dxkp1(der_eq[2].getrow(active_c_index_j));

			// Assign to D_i
			size_t index_j = j + Neq_active + Nineq_active;
			for (size_t k = 0; k < Nx; k++) {
				D_i.at(index_j, k) = dx[k];
				D_i.at(index_j, k + Nx + Nu) = dxkp1[k];
			}
			for (size_t k = 0; k < Nu; k++) {
				D_i.at(index_j, k + Nx) = du[k];
			}
		}

		// x0 is not a variable
		if (i == 0) {
			if (D_i.nrows() != 0)
				D_i = D_i.submat(0, Nx, D_i.nrows() - 1, 2 * Nx + Nu - 1);
			else
				D_i = matrixdb(0, Nx + Nu);
		}

		// Assign
		block_D.push_back(D_i);
		list_active_constraint_index.push_back(vector<vector<size_t>>{
			list_active_eq_index, list_active_ineq_index, list_active_c_index});
	}

	// Terminal constraints
	vector<size_t> list_active_teq_index, list_active_tineq_index;
	list_active_teq_index.reserve(Nteq);
	list_active_tineq_index.reserve(Ntineq);

	// Equality constraints
	for (size_t j = 0; j < Nteq; j++) {
		double teq_j = teq_[j];
		if (abs(teq_j) > active_constraint_tol) {
			// Assign to d
			d.push_back(teq_j);

			// Save index
			list_active_teq_index.push_back(j);
		}
	}

	// Inequality constraints
	for (size_t j = 0; j < Ntineq; j++) {
		double tineq_j = tineq_[j];
		if (tineq_j > active_constraint_tol) {
			// Assign to d
			d.push_back(tineq_j);

			// Save index
			list_active_tineq_index.push_back(j);
		}
	}

	// Build gradient
	size_t Nteq_active(list_active_teq_index.size());
	size_t Ntineq_active(list_active_tineq_index.size());
	matrixdb D_i(
		Nteq_active + Ntineq_active,
		Nx);

	// Equality constraints
	for (size_t j = 0; j < list_active_teq_index.size(); j++) {
		vectordb dx(der_teq_[0].getrow(list_active_teq_index[j]));
		D_i.setrow(j, dx);
	}

	// Inequality constraints
	for (size_t j = 0; j < list_active_tineq_index.size(); j++) {
		vectordb dx(der_tineq_[0].getrow(list_active_tineq_index[j]));
		D_i.setrow(j + Nteq_active, dx);
	}

	// Assign
	block_D.push_back(D_i);
	list_active_constraint_index.push_back(vector<vector<size_t>>{
		list_active_teq_index, list_active_tineq_index, vector<size_t>()});

	// Return d & D
	return pair<
		vectordb,
		pair<
			vector<matrixdb>,
			tri_vector_size_t>>(
		d, pair<vector<matrixdb>, tri_vector_size_t>( block_D, list_active_constraint_index));
}

// Return the matrix Sigma = D_a * D_a^t 
// Where is D_a without the active constraints.
// Using tridiagonal symetric block computation.
// TO DO: add reference.
sym_tridiag_matrixdb PNSolver::get_block_sigma_sq_(
	vector<matrixdb> const& block_Delta,
	tri_vector_size_t const& list_active_index) {
	// Unpack
	DDPSolver ddp_solver = DDPsolver();
	SolverParameters solver_parameters = ddp_solver.solver_parameters();
	unsigned int N = solver_parameters.N();
	unsigned int Nx = solver_parameters.Nx();
	unsigned int Nu = solver_parameters.Nu();

	//  Init lists
	vector<matrixdb> list_diag, list_subdiag;
	list_diag.reserve(block_Delta.size());
	list_subdiag.reserve(block_Delta.size() - 1);

	// Path elements
	for (size_t i = 0; i < N; i++) {

		// Unpack
		matrixdb Delta_i = block_Delta[i];
		vector<vector<size_t>> list_active_index_i = list_active_index[i];
		vector<size_t> list_active_index_c_i = list_active_index_i[2];
		size_t Neq_i(list_active_index_i[0].size());
		size_t Nineq_i(list_active_index_i[1].size());
		size_t Nc_i(list_active_index_c_i.size());

		// Unpack contraints gradient blocks
		matrixdb A_i(Neq_i + Nineq_i, Nx), C_i(Nc_i, Nx);
		matrixdb B_i(Neq_i + Nineq_i, Nu), D_i(Nc_i, Nu);

		if (Neq_i + Nineq_i != 0) {
			if (i == 0)
				B_i = Delta_i.submat(
					0, 0,
					Neq_i + Nineq_i - 1, Nu - 1);
			else {
				A_i = Delta_i.submat(
					0, 0,
					Neq_i + Nineq_i - 1, Nx - 1);
				B_i = Delta_i.submat(
					0, Nx,
					Neq_i + Nineq_i - 1, Nx + Nu - 1);
			}
		}
		if (Nc_i != 0) {
			if (i == 0)
				D_i = Delta_i.submat(
					Neq_i + Nineq_i, 0,
					Neq_i + Nineq_i + Nc_i - 1, Nu - 1);
			else {
				C_i = Delta_i.submat(
					Neq_i + Nineq_i, 0,
					Neq_i + Nineq_i + Nc_i - 1, Nx - 1);
				D_i = Delta_i.submat(
					Neq_i + Nineq_i, Nx,
					Neq_i + Nineq_i + Nc_i - 1, Nx + Nu - 1);
			}
		}
		matrixdb B_i_t = B_i.transpose(); matrixdb D_i_t = D_i.transpose();
		matrixdb A_i_t = A_i.transpose(); matrixdb C_i_t = C_i.transpose();

		// Block formulation
		matrixdb R_i, S_i, T_i;
		matrixdb W_i, Y_i;
		if (i == 0) {
			// Diagonal
			R_i = B_i * B_i_t;
			S_i = D_i * D_i_t;
			T_i = D_i * B_i_t;
		}
		else {
			// Diagonal
			R_i = B_i * B_i_t + A_i * A_i_t;
			S_i = D_i * D_i_t + C_i * C_i_t;
			T_i = D_i * B_i_t + C_i * A_i_t;

			// Tridiag
			W_i = -1.0 * (A_i); Y_i = -1.0 * (C_i);
		}
		size_t size_S = S_i.ncols();
		for (size_t j = 0; j < size_S; j++) { S_i.at(j, j) += 1.0; }

		// Make diagonal matrix 
		size_t size_R = R_i.ncols();
		size_t size_diag_im1 = 0;
		matrixdb diag_i(size_R + size_S);
		for (size_t j = 0; j < size_R; j++) {

			// Assign
			for (size_t k = 0; k < size_R; k++) {
				diag_i.at(j, k) = R_i.at(j, k);
			}
			for (size_t k = 0; k < size_S; k++) {
				diag_i.at(j, k + size_R) = T_i.at(k, j);
			}
		}
		for (size_t j = 0; j < size_S; j++) {

			// Assign
			for (size_t k = 0; k < size_R; k++) {
				diag_i.at(size_R + j, k) = T_i.at(j, k);
			}
			for (size_t k = 0; k < size_S; k++) {
				diag_i.at(size_R + j, k + size_R) = S_i.at(j, k);
			}
		}

		// Subdiag
		matrixdb subdiag_i;
		if (i != 0) {
			size_diag_im1 = list_diag[i - 1].ncols();
			size_t size_diag_i = diag_i.ncols();
			vector<size_t> list_active_index_c_im1 = list_active_index[i - 1][2];
			size_t Nc_im1 = list_active_index_c_im1.size();
			subdiag_i = matrixdb(size_diag_i, size_diag_im1);
			if (size_diag_i != 0 && Nc_im1 != 0 && size_diag_im1 != 0) {
				for (size_t j = 0; j < Nc_im1; j++) {
					// Get cols
					size_t active_index_c_im1_j = list_active_index_c_im1[j];
					vectordb col_W_j = W_i.getcol(active_index_c_im1_j);
					col_W_j.reserve(Y_i.nrows());
					for (size_t k = 0; k < Y_i.nrows(); k++) {
						col_W_j.push_back(Y_i.at(k, active_index_c_im1_j));
					}

					// Assign
					subdiag_i.setcol(size_diag_im1 - Nc_im1 + j, col_W_j);
				}
			}
		}
	
		// Assign
		list_diag.push_back(diag_i);
		if (i != 0)
			list_subdiag.push_back(subdiag_i);
	}

	// Last elements
	size_t i = N;
	matrixdb Delta_i = block_Delta[i];
	vector<vector<size_t>> list_active_index_i = list_active_index[i];
	size_t Nteq(list_active_index_i[0].size());
	size_t Ntineq(list_active_index_i[1].size());
	matrixdb A_i(Nteq + Ntineq, Nx, 0.0);
	if (Nteq + Ntineq != 0) {
		A_i = Delta_i.submat(
			0, 0,
			Nteq + Ntineq - 1, Nx - 1);
	}
	matrixdb A_i_t = A_i.transpose();

	// Block formulation

	// Diagonal
	matrixdb R_i = A_i * A_i_t;

	// Make diagonal matrix 
	size_t size_R = R_i.ncols();
	matrixdb diag_i(R_i);

	// Tridiag
	size_t size_diag_im1 = list_diag[i - 1].ncols();
	size_t size_diag_i = diag_i.ncols();
	vector<size_t> list_active_index_c_im1 = list_active_index[i - 1][2];
	size_t Nc_im1 = list_active_index_c_im1.size();
	matrixdb subdiag_i(size_diag_i, size_diag_im1);
	matrixdb W_i = -1.0 * A_i;
	if (size_diag_i != 0 && Nc_im1 != 0 && size_diag_im1 != 0) {
		for (size_t j = 0; j < Nc_im1; j++) {
			// Get cols
			vectordb col_W_j = W_i.getcol(list_active_index_c_im1[j]);

			// Assign
			subdiag_i.setcol(size_diag_im1 - Nc_im1 + j, col_W_j);
		}
	}

	// Assign
	list_diag.push_back(diag_i);
	list_subdiag.push_back(subdiag_i);
	
	return pair<
		vector<matrixdb>,
		vector<matrixdb>>(list_diag, list_subdiag);
}
