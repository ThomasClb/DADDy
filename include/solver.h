/**
	solver.h

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