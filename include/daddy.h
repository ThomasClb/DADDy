/**
	daddy.h

	Purpose: Implementation of the DADDy solver class.

	@author Thomas Caleb

	@version 1.0 08/03/2024
*/

#ifndef DEF_DADDY
#define DEF_DADDY

#pragma once

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <chrono>

#include <dace/dace_s.h>

#include "settings.h"
#include "constants.h"
#include "dynamics.h"
#include "IO.h"
#include "aul_solver.h"
#include "pn_solver.h"

class DADDy {

// Attributes
protected:
	AULSolver AULsolver_; // Augmented Lagrangian solver
	SolverParameters solver_parameters_;
	SpacecraftParameters spacecraft_parameters_;
	Dynamics dynamics_;
	PNSolver PNsolver_; // Projected Newton solver

	// Iterations
	double PN_runtime_; // Runtime of PN
	double AUL_runtime_; // Runtime of AUL
	double runtime_; // Runtime
	std::size_t PN_n_iter_; // Number of PN iterations
	std::size_t AUL_n_iter_; // Number of AUL iterations
	std::size_t DDP_n_iter_; // Total number of DDP iterations

// Methods
public:
	// Empty constructor
	DADDy();

	// Constructor
	DADDy(
		SolverParameters const& solver_parameters,
		SpacecraftParameters const& spacecraft_parameters,
		Dynamics const& dynamics);

	// Copy constructor
	DADDy(DADDy const& solver);

	// Destructors
	~DADDy();

	// Getters
	const AULSolver AULsolver() const;
	const PNSolver PNsolver() const;
	const std::vector<DACE::vectordb> list_x() const;
	const std::vector<DACE::vectordb> list_u() const;
	const double cost() const;
	const double PN_runtime() const;
	const double AUL_runtime() const;
	const double runtime() const;
	const std::size_t PN_n_iter() const;
	const std::size_t DDP_n_iter() const;
	const std::size_t AUL_n_iter() const;

	// Computes the constraints with propagated dynamics.
	double real_constraints(
		DACE::vectordb const& x_goal);

	// Performs solving given a starting point,
	// initial controls and a final state.
	// Inspired from ALTRO (Julia).
	// See: https://github.com/RoboticExplorationLab/Altro.jl
	void solve(
		DACE::vectordb const& x0,
		std::vector<DACE::vectordb> const& list_u_init,
		DACE::vectordb const& x_goal,
		bool const& fuel_optimal,
		bool const& pn_solving);
};

#endif