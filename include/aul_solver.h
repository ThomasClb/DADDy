/**
	aul_solver.h

	Purpose: Implementation of the AULSolver class.

	@author Thomas Caleb

	@version 1.0 17/11/2023
*/

#ifndef DEF_AUL
#define DEF_AUL

#pragma once

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <chrono>

#include <dace/dace_s.h>

#include "ddp_solver.h"

class AULSolver {

// Attributes
protected:
	DDPSolver DDPsolver_; // Solver parameters
	std::vector<DACE::vectordb> list_x_; // List of states
	std::vector<DACE::vectordb> list_u_; // List of controls
	double cost_; // Output cost [-]
	std::vector<DACE::vectordb> list_eq_; // List of equality constraints
	std::vector<DACE::vectordb> list_ineq_; // List of inequality constraints
	DACE::vectordb teq_; // List of terminal equality constraints
	DACE::vectordb tineq_; // List of terminal inequality constraints
	std::vector<DACE::vectordb> list_lambda_; // List of Lagrange multiplicator
	std::vector<DACE::vectordb> list_mu_; // List of penalty factors


	// Iterations
	unsigned int AUL_n_iter_; // Number of AUL iterations
	unsigned int DDP_n_iter_; // Total number of DDP iterations

// Methods
public:
	// Empty constructor
	AULSolver();

	// Constructor
	AULSolver(
		SolverParameters const& solver_parameters,
		SpacecraftParameters const& spacecraft_parameters,
		Dynamics const& dynamics);

	// Copy constructor
	AULSolver(AULSolver const& solver);

	// Destructors
	~AULSolver();

	// Getters
	const DDPSolver DDPsolver() const;
	const std::vector<DACE::vectordb> list_x() const;
	const std::vector<DACE::vectordb> list_u() const;
	const double cost() const;
	const std::vector<DACE::vectordb> list_eq() const;
	const std::vector<DACE::vectordb> list_ineq() const;
	const DACE::vectordb teq() const;
	const DACE::vectordb tineq() const;
	const std::vector<DACE::vectordb> list_lambda() const;
	const std::vector<DACE::vectordb> list_mu() const;
	const unsigned int AUL_n_iter() const;
	const unsigned int DDP_n_iter() const;

	// Setters
	void set_ToF(double const& ToF);
	void set_homotopy_coefficient(double const& homotopy_coefficient);
	void set_huber_loss_coefficient(double const& huber_loss_coefficient);

	// Compute the value of the max constraint.
	double get_max_constraint_();

	// Update dual state in Augmented Lagrangian formulation.
	void update_lambda_();

	// Update penalty in Augmented Lagrangian formulation.
	void update_mu_();

	// Performs AUL solving given a starting point,
	// initial controls and a final state.
	// Inspired from ALTRO (Julia).
	// See: https://github.com/RoboticExplorationLab/Altro.jl
	void solve(
		DACE::vectordb const& x0,
		std::vector<DACE::vectordb> const& list_u_init,
		DACE::vectordb const& x_goal);
};

#endif