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
	X_U_(0), EQ_INEQ_(0), der_EQ_INEQ_(0) {}

// Constructors
PNSolver::PNSolver(AULSolver const& AULsolver) : AULsolver_(AULsolver),
	list_x_(AULsolver.list_x()), list_u_(AULsolver.list_u()),
	cost_(AULsolver.cost()),
	list_der_cost_(vector<vector<matrixdb>>(AULsolver.list_eq().size() + 1)),
	X_U_(), EQ_INEQ_(), der_EQ_INEQ_() {
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
	X_U_ = vectordb(N*(Nx + Nu));
	EQ_INEQ_ = vectordb(N*(Nx + Neq + Nineq) + Nteq + Ntineq);
	der_EQ_INEQ_ = vector<matrixdb>(N*6 + 4);

	// First control
	for (size_t j=0; j<Nu; j++)  {
		X_U_[j] = list_u_[0][j];
	}

	// Loop
	for (size_t i=1; i<N; i++) {
		for (size_t j=0; j<Nx; j++)  {
			X_U_[(i - 1)*(Nx + Nu) + Nu + j] = list_x_[i][j];
		}
		for (size_t j=0; j<Nu; j++)  {
			X_U_[i*(Nx + Nu) + j] = list_u_[i][j];
		}
	}

	// Last state
	for (size_t j=0; j<Nx; j++)  {
		X_U_[(N - 1)*(Nx + Nu) + Nu + j] = list_x_[N][j];
	}
}

// Copy constructor
PNSolver::PNSolver(
	PNSolver const& solver) : AULsolver_(solver.AULsolver_),
	list_x_(solver.list_x_), list_u_(solver.list_u_),
	cost_(solver.cost_),
	list_der_cost_(vector<vector<matrixdb>>(solver.AULsolver_.list_eq().size() + 1)),
	X_U_(solver.X_U_), EQ_INEQ_(solver.EQ_INEQ_),
	der_EQ_INEQ_(solver.der_EQ_INEQ_) {}

// Destructors
PNSolver::~PNSolver() {}

// Getters
const AULSolver PNSolver::AULsolver() const { return AULsolver_; }
const DDPSolver PNSolver::DDPsolver() const { return AULsolver_.DDPsolver(); }
const vector<vectordb> PNSolver::list_x() const { return list_x_; }
const vector<vectordb> PNSolver::list_u() const { return list_u_; }
const double PNSolver::cost() const { return cost_; }

// Setters
void PNSolver::set_list_x_u() {
	// Unpack
	DDPSolver ddp_solver = DDPsolver();
	SolverParameters solver_parameters = ddp_solver.solver_parameters();
	unsigned int N = solver_parameters.N();
	unsigned int Nx = solver_parameters.Nx();
	unsigned int Nu = solver_parameters.Nu();


	// Loop for x
	for (size_t i=1; i< N + 1; i++) {
		for (size_t j=0; j<Nx; j++)  {
			list_x_[i][j] = X_U_[(i - 1)*(Nx + Nu) + Nu + j];
		}
	}

	// Loop for u
	for (size_t i=0; i<N; i++) {
		for (size_t j=0; j<Nu; j++)  {
			list_u_[i][j] = X_U_[i*(Nx + Nu) + j];
		}
	}
}

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
	double prev_violation = get_max_constraint_(EQ_INEQ_);

	// Init loop
	double violation(prev_violation);
	for (size_t i = 0; i < max_iter; i++) {
		// Output
		if (verbosity < 1) {
			// if (i % 5 == 0)
				cout << i << " - " << prev_violation << " - " << X_U_[X_U_.size() - 2] * constants.massu() << endl;
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
	set_list_x_u();
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

	// Init loop
	double violation(violation_0); vectordb d(d_0);
	for (size_t i = 0; i < 10; i++) {
		// Solve Sigma * z = d <=> L * L^t * z = d using
		// block tridiagonal Cholesky factorisation.
		vectordb z = solve_cholesky_(tridiag_L, d);
		mat_z.setcol(0, z);	
		vectordb correction(N * (Nu + Nx));
 
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
				correction[counter_corr + k] += alpha * corr_j[k];
			}
			counter_z += N_j;
			counter_corr += corr_j.size() - Nx;
		}

		// Comoute X_U
		vectordb X_U = X_U_ - correction;

		// Evaluate constraints
		vectordb EQ_INEQ = update_constraints_double_( // TO DO : DA method
			x_goal, X_U, correction);

		// TO DO use DA expansion

		// Get the max constraint
		violation = get_max_constraint_(EQ_INEQ);

		// Check exit
		if (violation <= violation_0) {
			// Update state and constraints
			X_U_ = X_U;
			EQ_INEQ_ = EQ_INEQ;
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
	DACE::vectordb const& EQ_INEQ) {
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
	double maximum = -1e15;
	for (size_t i = 0; i < N; i++) {

		// Find maximum
		for (size_t j = 0; j < Neq + Nx; j++) {
			maximum = max(maximum, abs(EQ_INEQ[i*(Neq + Nineq + Nx) + j]));
		}
		for (size_t j = 0; j < Nineq; j++) {
			maximum = max(maximum, EQ_INEQ[i*(Neq + Nineq + Nx) + Neq + Nx + j]);
		}
	}

	// Terminal constraints
	for (size_t j = 0; j < Nteq; j++) {
		maximum = max(maximum, abs(EQ_INEQ[N*(Neq + Nineq + Nx) + j]));
	}
	for (size_t j = 0; j < Ntineq; j++) {
		maximum = max(maximum, EQ_INEQ[N*(Neq + Nineq + Nx) + Nteq + j]);
	}

	return maximum;
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
	vectorDA x_DA, u_DA;
	vectordb x_kp1;
	for (size_t i = 0; i < N; i++) {

		// Get DA x, u
		if (i == 0)
			x_DA = id_vector(list_x_[0], 0, 0, Nx - 1);
		else 
			x_DA = id_vector(
				X_U_.extract(
					Nu + (i - 1)*(Nx + Nu),
					Nu + Nx - 1 + (i - 1)*(Nx + Nu)),
				0, 0, Nx - 1);
		u_DA = id_vector(
			X_U_.extract(
				i*(Nx + Nu),
				Nu - 1 + i*(Nx + Nu)),
			0, Nx, Nu + Nx);
		x_kp1 = X_U_.extract(
					Nu + i*(Nx + Nu),
					Nu + Nx - 1 + i*(Nx + Nu));

		// Constraints evaluations
		vectorDA eq_eval = dynamics.equality_constraints()(
			x_DA, u_DA, spacecraft_parameters, dynamics.constants(), solver_parameters);
		vectorDA ineq_eval = dynamics.inequality_constraints()(
			x_DA, u_DA, spacecraft_parameters, dynamics.constants(), solver_parameters);

		// Continuity constraints TO DO DA
		vectorDA x_kp1_eval = dynamics.dynamic()(
			x_DA, u_DA,
			spacecraft_parameters, dynamics.constants(), solver_parameters); // TO DO : test radius + evaluate 
		list_dynamics_.push_back(x_kp1_eval);
		x_kp1_eval -= x_kp1;

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
		der_EQ_INEQ_[6*i] = der_eq[0];
		der_EQ_INEQ_[6*i + 1] = der_eq[1];
		der_EQ_INEQ_[6*i + 2] = id_Nx;
		der_EQ_INEQ_[6*i + 3] = der_ineq[0];
		der_EQ_INEQ_[6*i + 4] = der_ineq[1];
		der_EQ_INEQ_[6*i + 5] = null_Nx;

		// Assign
		for (size_t k = 0; k < Neq; k++) {
			EQ_INEQ_[i*(Nx + Neq + Nineq) + k] = eq_eval[k].cons();}
		for (size_t k = 0; k < Nx; k++) {
			EQ_INEQ_[i*(Nx + Neq + Nineq) + Neq + k] = x_kp1_eval[k].cons();}
		for (size_t k = 0; k < Nineq; k++) {
			EQ_INEQ_[i*(Nx + Neq + Nineq) + Neq + Nx + k] = ineq_eval[k].cons();}
	}

	// Update terminal constraints

	// Get DA x, u
	x_DA = id_vector(
		X_U_.extract(
			Nu + (N - 1)*(Nx + Nu),
			Nu + Nx - 1 + (N - 1)*(Nx + Nu)),
		0, 0, Nx - 1);

	// Constraints evaluations TO DO DA
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
	for (size_t k = 0; k < Nteq; k++) {
		EQ_INEQ_[N*(Nx + Neq + Nineq) + k] = teq_eval.cons()[k];}
	for (size_t k = 0; k < Ntineq; k++) {
		EQ_INEQ_[N*(Nx + Neq + Nineq) + Nteq + Nx + k] = tineq_eval.cons()[k];}
	der_EQ_INEQ_[6*N] = der_teq[0];
	der_EQ_INEQ_[6*N + 1] = der_teq[1];
	der_EQ_INEQ_[6*N + 2] = der_tineq[0];
	der_EQ_INEQ_[6*N + 3] = der_tineq[1];
}

// Computes the new constraints given states and controls without DA
// Return the list of equalities, and inequalities
vectordb PNSolver::update_constraints_double_(
	vectordb const& x_goal,
	vectordb const& X_U,
	vectordb const& correction) {
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
	vectordb EQ_INEQ(N*(Neq + Nx + Nineq) + Nteq + Ntineq);

	// Update path constraints

	// Loop on all steps
	vectordb x(Nx), xp1(Nx), u(Nu);
	for (size_t i = 0; i < N; i++) {

		// Get x, u, xp1
		for (size_t j=0; j<Nu; j++) {
			u[j] = X_U[i*(Nu + Nx) + j];
		}
		if (i == 0)
			x = list_x_[0]; // Never changes
		else {
			for (size_t j=0; j<Nx; j++) {
				x[j] = X_U[(i - 1)*(Nu + Nx) + Nu + j];
			}
		}
		for (size_t j=0; j<Nx; j++) {
			xp1[j] = X_U[i*(Nu + Nx) + Nu + j];
		}
	
		/**/
		// Constraints evaluations
		vectordb eq_eval = dynamics.equality_constraints_db()(
			x, u, spacecraft_parameters, dynamics.constants(), solver_parameters);
		vectordb ineq_eval = dynamics.inequality_constraints_db()(
			x, u, spacecraft_parameters, dynamics.constants(), solver_parameters);
		eq_eval.reserve(Nx);

		// Continuity constraints TO DO DA
		vectordb x_kp1_eval = (dynamics.dynamic_db()(
			x, u, spacecraft_parameters, dynamics.constants(), solver_parameters) - xp1); // TO DO : use DA, test

		// Assign
		for (size_t k = 0; k < Neq; k++) {
			EQ_INEQ[i*(Nx + Neq + Nineq) + k] = eq_eval[k];}
		for (size_t k = 0; k < Nx; k++) {
			EQ_INEQ[i*(Nx + Neq + Nineq) + Neq + k] = x_kp1_eval[k];}
		for (size_t k = 0; k < Nineq; k++) {
			EQ_INEQ[i*(Nx + Neq + Nineq) + Neq + Nx + k] = ineq_eval[k];}
	}

	// Update terminal constraints

	// Get x_f
	for (size_t j=0; j<Nx; j++) {
		x[j] = X_U[(N - 1)*(Nu + Nx) + Nu + j];
	}

	// Constraints evaluations TO DO DA
	vectordb teq_eval = dynamics.terminal_equality_constraints_db()(
		x, x_goal, spacecraft_parameters, dynamics.constants(), solver_parameters);
	vectordb tineq_eval = dynamics.terminal_inequality_constraints_db()(
		x, x_goal, spacecraft_parameters, dynamics.constants(), solver_parameters);

	// Assign
	for (size_t k = 0; k < Nteq; k++) {
		EQ_INEQ[N*(Nx + Neq + Nineq) + k] = teq_eval[k];}
	for (size_t k = 0; k < Ntineq; k++) {
		EQ_INEQ[N*(Nx + Neq + Nineq) + Nteq + Nx + k] = tineq_eval[k];}


	return EQ_INEQ;
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

		// Equality constraints
		for (size_t j = 0; j < Neq; j++) {
			double abs_eq_i_j = abs(EQ_INEQ_[i*(Neq + Nx + Nineq) + j]);
			if (abs_eq_i_j > active_constraint_tol) {
				// Assign to d
				d.push_back(abs_eq_i_j);

				// Save index
				list_active_eq_index.push_back(j);
			}
		}

		// Inequality constraints
		for (size_t j = 0; j < Nineq; j++) {
			double ineq_i_j = EQ_INEQ_[i*(Neq + Nx + Nineq) + Neq + Nx + j];
			if (ineq_i_j > active_constraint_tol) {
				// Assign to d
				d.push_back(ineq_i_j);

				// Save index
				list_active_ineq_index.push_back(j);
			}
		}

		// Continuity constraints
		for (size_t j = 0; j < Nx; j++) {
			double eq_i_j = EQ_INEQ_[i*(Neq + Nx + Nineq) + Neq + j];
			if (abs(eq_i_j) > active_constraint_tol) {
				// Assign to d
				d.push_back(eq_i_j);

				// Save index
				list_active_c_index.push_back(Neq + j);
			}
		}

		// Build gradient
		size_t Neq_active(list_active_eq_index.size());
		size_t Nineq_active(list_active_ineq_index.size());
		size_t Nc_active(list_active_c_index.size());
		matrixdb D_i(
			Neq_active + Nineq_active + Nc_active,
			(2 * Nx + Nu));
		matrixdb der_eq_x = der_EQ_INEQ_[6*i];
		matrixdb der_eq_u = der_EQ_INEQ_[6*i + 1];
		matrixdb der_eq_x_kp1 = der_EQ_INEQ_[6*i + 2];
		matrixdb der_ineq_x = der_EQ_INEQ_[6*i + 3];
		matrixdb der_ineq_u = der_EQ_INEQ_[6*i + 4];
		matrixdb der_ineq_x_kp1 = der_EQ_INEQ_[6*i + 5];

		// Equality constraints
		for (size_t j = 0; j < list_active_eq_index.size(); j++) {
			// Get derivatives
			size_t active_eq_index_j = list_active_eq_index[j];
			vectordb dx(der_eq_x.getrow(active_eq_index_j));
			vectordb du(der_eq_u.getrow(active_eq_index_j));
			vectordb dxkp1(der_eq_x_kp1.getrow(active_eq_index_j));

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
			vectordb dx(der_ineq_x.getrow(active_ineq_index_j));
			vectordb du(der_ineq_u.getrow(active_ineq_index_j));
			vectordb dxkp1(der_ineq_x_kp1.getrow(active_ineq_index_j));

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
			vectordb dx(der_eq_x.getrow(active_c_index_j));
			vectordb du(der_eq_u.getrow(active_c_index_j));
			vectordb dxkp1(der_eq_x_kp1.getrow(active_c_index_j));

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
		double teq_j = EQ_INEQ_[N*(Neq + Nx + Nineq) + j];
		if (abs(teq_j) > active_constraint_tol) {
			// Assign to d
			d.push_back(teq_j);

			// Save index
			list_active_teq_index.push_back(j);
		}
	}

	// Inequality constraints
	for (size_t j = 0; j < Ntineq; j++) {
		double tineq_j = EQ_INEQ_[N*(Neq + Nx + Nineq) + Nteq + j];
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
	matrixdb der_teq_x = der_EQ_INEQ_[6*N];
	matrixdb der_tineq_x = der_EQ_INEQ_[6*N + 2];
	for (size_t j = 0; j < list_active_teq_index.size(); j++) {
		vectordb dx(der_teq_x.getrow(list_active_teq_index[j]));
		D_i.setrow(j, dx);
	}

	// Inequality constraints
	for (size_t j = 0; j < list_active_tineq_index.size(); j++) {
		vectordb dx(der_tineq_x.getrow(list_active_tineq_index[j]));
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
