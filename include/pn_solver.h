/**
	pn_solver.h

	Purpose: Implementation of the PNSolver class.

	@author Thomas Caleb

	@version 1.0 17/11/2023
*/

#ifndef DEF_PN
#define DEF_PN

#pragma once

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <chrono>

#include <dace/dace_s.h>

#include "aul_solver.h"

// Stores the list of active index
// N * [active_index_eq, active_index_ineq, active_index_cont]
// [active_index_teq, active_index_tineq, [] ]
using tri_vector_size_t = std::vector
	<std::vector<std::vector<std::size_t>>>;


class PNSolver {

	// Attributes
protected:
	AULSolver AULsolver_;
	std::vector<DACE::vectordb> list_x_; // List of states
	std::vector<DACE::vectordb> list_u_; // List of controls
	double cost_; // Output cost [-]
	std::vector<std::vector<DACE::matrixdb>> list_der_cost_; // Output cost derivatives [-]
	std::vector<DACE::vectordb> list_eq_; // List of equality constraints
	std::vector<DACE::vectordb> list_ineq_; // List of inequality constraints
	DACE::vectordb teq_; // List of terminal equality constraints
	DACE::vectordb tineq_; // List of terminal inequality constraints
	std::vector<std::vector<DACE::matrixdb>> list_der_eq_; // List of equality constraints derivatives
	std::vector<std::vector<DACE::matrixdb>> list_der_ineq_; // List of inequality constraints derivatives
	std::vector<DACE::matrixdb> der_teq_; // List of terminal equality constraints derivatives
	std::vector<DACE::matrixdb> der_tineq_; // List of terminal inequality constraints derivatives
	std::vector<DACE::vectorDA> list_dynamics_; // List of dynamics evaluations

	// Looping attributes
	DACE::vectordb  X_U_;
	// MAKE vector X_U
 	// Make vector EQ_INEQ 

// Methods
public:
	// Empty constructors
	PNSolver();

	// Constructors
	PNSolver(AULSolver const& AULSolver);

	// Copy constructor
	PNSolver(PNSolver const& solver); 

	// Destructors
	~PNSolver();

	// Getters
	const AULSolver AULsolver() const;
	const DDPSolver DDPsolver() const;
	const std::vector<DACE::vectordb> list_x() const;
	const std::vector<DACE::vectordb> list_u() const;
	const double cost() const;
	const std::vector<DACE::vectordb> list_eq() const;
	const std::vector<DACE::vectordb> list_ineq() const;
	const DACE::vectordb teq() const;
	const DACE::vectordb tineq() const;

	// Setters
	void set_list_x_u();


	// Solves the optimisation problem with a projected Newton method
	// Inspired from ALTRO (Julia).
	// See: https://github.com/RoboticExplorationLab/Altro.jl
	void solve(DACE::vectordb const& x_goal);
	
	// Iterative line search for PN.
	// Inspired from ALTRO (Julia).
	// See: https://github.com/RoboticExplorationLab/Altro.jl
	double line_search_(
		DACE::vectordb const& x_goal,
		sym_tridiag_matrixdb const& tridiag_L,
		std::vector<DACE::matrixdb> const& block_D,
		DACE::vectordb const& d_0,
		double const& violation_0);

	// Computes the maximum constraints given eq in ineq constraints
	double get_max_constraint_(
		std::vector<DACE::vectordb> const& list_eq,
		std::vector<DACE::vectordb> const& list_ineq);

	// Computes the maximum constraints using the attributes
	double get_max_constraint_();

	// Updates the values of controls and states given corrections
	// Return the list of equalities, and inequalities
	std::pair<
		std::vector<DACE::vectordb>,
		std::vector<DACE::vectordb>> update_list_x_u_(
			DACE::vectordb const& correction);

	// Computes the new constraints given states and controls
	void update_constraints_(DACE::vectordb const& x_goal);

	// Computes the new constraints given states and controls without DA
	// Return the list of equalities, and inequalities
	std::pair<
		std::vector<DACE::vectordb>,
		std::vector<DACE::vectordb>> update_constraints_double_(
			DACE::vectordb const& x_goal,
			DACE::vectordb const& X_U,
			DACE::vectordb const& correction);

	// Assign the constraints without DA
	// Given the list of equalities, and inequalities
	void assign_constraints_double_(
		std::pair<
		std::vector<DACE::vectordb>,
		std::vector<DACE::vectordb>> const& list_eq_ineq);

	// Returns the vector of active constraints and their gradients
	// first it the active constraints vector
	// second is a pair with the list of gradients of constraints first
	// second.second is the list of active constraints.
	std::pair<
		DACE::vectordb,
		std::pair<
		std::vector<DACE::matrixdb>,
		tri_vector_size_t> > get_d_block_D_();

	// Return the matrix Sigma = D_a * D_a^t 
	// Where is D_a without the active constraints.
	// Using tridiagonal symetric block computation.
	// TO DO: add reference.
	sym_tridiag_matrixdb get_block_sigma_sq_(
		std::vector<DACE::matrixdb> const& block_Delta,
		tri_vector_size_t const& list_active_index);
};


#endif