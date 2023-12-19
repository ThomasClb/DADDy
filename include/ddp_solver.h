/**
	ddp_solver.h

	Purpose: Implementation of the DDPSolver class.

	@author Thomas Caleb

	@version 1.0 17/11/2023
*/

#ifndef DEF_DDP
#define DEF_DDP

#pragma once

#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include <dace/dace_s.h>

#include "parameters.h"
#include "derivation.h"
#include "dynamics.h"
#include "linalg.h"
#include "settings.h"

class DDPSolver {

// Attributes
protected:
	SolverParameters solver_parameters_; // Solver parameters
	SpacecraftParameters spacecraft_parameters_; // Spacecraft parameters
	Dynamics dynamics_; // Dynamics
	std::vector<DACE::vectordb> list_x_; // List of states
	std::vector<DACE::vectordb> list_u_; // List of controls
	std::vector<DACE::vectordb> list_eq_; // List of equality constraints
	std::vector<DACE::vectordb> list_ineq_; // List of inequality constraints
	DACE::vectordb teq_; // List of terminal equality constraints
	DACE::vectordb tineq_; // List of terminal inequality constraints
	double cost_; // Output cost [-]

	// Looping variables
	unsigned int n_iter_;
	double rho_;
	double alpha_;
	double d_rho_;
	std::vector<DACE::vectorDA> list_dynamic_eval_;
	std::vector<DACE::DA> list_ctg_eval_;
	DACE::DA tc_eval_;
	std::vector<DACE::matrixdb> list_Qu_;
	std::vector<DACE::matrixdb> list_k_;
	std::vector<DACE::matrixdb> list_K_;

// Methods
public:
	// Empty constructor
	DDPSolver();

	// Constructor
	DDPSolver(
		SolverParameters const& solver_parameters,
		SpacecraftParameters const& spacecraft_parameters,
		Dynamics const& dynamics);
	
	// Copy constructor
	DDPSolver(DDPSolver const& solver);

	// Destructors
	~DDPSolver();

	// Getters
	SolverParameters solver_parameters();
	const SpacecraftParameters spacecraft_parameters() const;
	const Dynamics dynamics() const;
	const std::vector<DACE::vectordb> list_x() const;
	const std::vector<DACE::vectordb> list_u() const;
	const std::vector<DACE::vectordb> list_eq() const;
	const std::vector<DACE::vectordb> list_ineq() const;
	const DACE::vectordb teq() const;
	const DACE::vectordb tineq() const;
	const double cost() const;
	const DACE::DA tc_eval() const;

	// Setters
	void set_list_lambda(std::vector<DACE::vectordb> const& list_lambda);
	void set_list_mu(std::vector<DACE::vectordb> const& list_mu);
	void set_ToF(double const& ToF);
	void set_homotopy_coefficient(double const& homotopy_coefficient);
	void set_huber_loss_coefficient(double const& huber_loss_coefficient);

	// Returns the Augmented lagrangian cost-to-go:
	// AUL_ctg = ctg + Lambda^T * c(x) + 0.5 * c(x)^T * I_mu * c(x)
	// DA version
	DACE::vectorDA get_AUL_cost_to_go(
		DACE::vectorDA const& x_star_DA, DACE::vectorDA const& u_star_DA, std::size_t const& index);

	// Returns the Augmented lagrangian terminal cost:
	// AUL_tc = tc + Lambda^T * c(x) + 0.5 * c(x)^T * I_mu * c(x)
	// Double version
	DACE::vectordb get_AUL_terminal_cost(
		DACE::vectordb const& x_star_DA, DACE::vectordb const& x_goal);

	// Returns the Augmented lagrangian terminal cost:
	// AUL_tc = tc + Lambda^T * c(x) + 0.5 * c(x)^T * I_mu * c(x)
	// DA version
	DACE::vectorDA get_AUL_terminal_cost(
		DACE::vectorDA const& x_star_DA, DACE::vectordb const& x_goal);

	// Returns the Augmented lagrangian cost-to-go:
	// AUL_ctg = ctg + Lambda^T * c(x) + 0.5 * c(x)^T * I_mu * c(x)
	// Double version
	DACE::vectordb get_AUL_cost_to_go(
		DACE::vectordb const& x_star_DA, DACE::vectordb const& u_star_DA, std::size_t const& index);


	// Increases the regulation to ensure Quu + rho*I
	// is symeteric positive definite.
	// Inspired from ALTRO (Julia).
	// See: https://github.com/RoboticExplorationLab/Altro.jl
	void increase_regularisation_();

	// Decreases the regulation to ensure Quu + rho*I
	// is symeteric positive definite.
	// Inspired from ALTRO (Julia).
	// See: https://github.com/RoboticExplorationLab/Altro.jl
	void decrease_regularisation_();

	// Computes the expected cost after backward sweep for linesearch.
	// Inspired from ALTRO (Julia).
	// See: https://github.com/RoboticExplorationLab/Altro.jl
	double expected_cost_(
		double const& alpha);

	// Evaluates the convergence of DDP optimisation
	bool evaluate_convergence_(
		double const& d_cost);

	// Compute the value of the max constraint.
	double get_max_constraint_();

	// Performs the DDP backward sweep, that consists in the computation
	// of the gains corrections.
	// Inspired from ALTRO (Julia).
	// See: https://github.com/RoboticExplorationLab/Altro.jl
	void backward_sweep_();

	// Performs the DDP backward sweep, that consists in the computation
	// of the gains corrections, with hessians.
	// Inspired from ALTRO (Julia).
	// See: https://github.com/RoboticExplorationLab/Altro.jl
	void backward_sweep_hessian_();

	// Performs the DDP backward sweep, that consists in the computation
	// of the gains corrections, Q is evaluated using DA, thus, with hessian.
	// Inspired from ALTRO (Julia).
	// See: https://github.com/RoboticExplorationLab/Altro.jl
	void backward_sweep_DA_Q_();

	// Performs the DDP forward pass, that consists in the computation
	// of the new states and control after correction.
	// Inspired from ALTRO (Julia).
	// See: https://github.com/RoboticExplorationLab/Altro.jl
	void forward_pass_(
		std::vector<DACE::vectordb> const& list_x,
		std::vector<DACE::vectordb> const& list_u,
		DACE::vectordb const& x_goal);

	// Performs the DDP forward pass, that consists in the computation
	// of the new states and control after correction using the DA mapping
	// Inspired from ALTRO (Julia).
	// See: https://github.com/RoboticExplorationLab/Altro.jl
	void forward_pass_convRadius_(
		std::vector<DACE::vectordb> const& list_x,
		std::vector<DACE::vectordb> const& list_u,
		DACE::vectordb const& x_goal);

	// Performs the DDP forward pass, that consists in the computation
	// of the new states and control after correction using the DA mapping
	// The linesearch is tweaked to implement a memory from one iteration to the other.
	// Inspired from ALTRO (Julia).
	// See: https://github.com/RoboticExplorationLab/Altro.jl
	void forward_pass_convRadius_ls_(
		std::vector<DACE::vectordb> const& list_x,
		std::vector<DACE::vectordb> const& list_u,
		DACE::vectordb const& x_goal);

	// Performs the DDP forward pass, that consists in the computation
	// of the new states and control after correction using the DA mapping
	// The linesearch is tweaked to implement a memory from one iteration to the other.
	// The linesearch computation are done with floats (in DA vectors), the DA mappings are computed
	// only when the linesearch is ended
	// Inspired from ALTRO (Julia).
	// See: https://github.com/RoboticExplorationLab/Altro.jl
	void forward_pass_convRadius_ls_Spencer_(
		std::vector<DACE::vectordb> const& list_x,
		std::vector<DACE::vectordb> const& list_u,
		DACE::vectordb const& x_goal);

	// Performs DDP solving given a starting point,
	// initial controls and a final state.
	// Inspired from ALTRO (Julia).
	// See: https://github.com/RoboticExplorationLab/Altro.jl
	void solve(
		DACE::vectordb const& x0,
		std::vector<DACE::vectordb> const& list_u_init,
		DACE::vectordb const& x_goal);
};

// Evaluates the convergence radius of a DA vector.
double convRadius(DACE::vectorDA const& x, double const& tol);

#endif