/**
	test_cases.h

	Purpose: Implementation of the test cases for DADDy.

	@author Thomas Caleb

	@version 1.0 07/12/2023
*/

#ifndef DEF_TEST_CASES
#define DEF_TEST_CASES

#pragma once

#include <dace/dace_s.h>
#include <chrono>

#include "settings.h"
#include "constants.h"
#include "dynamics.h"
#include "ddp_solver.h"
#include "aul_solver.h"
#include "pn_solver.h"
#include "visualisation.h"

// Test cases
void tbp_SUN_low_thrust_earth_to_mars(bool const& plot_graphs);
void cr3bp_EARTH_MOON_low_thrust_haloL2_to_haloL1(bool const& plot_graphs);

// Verification function
double real_constraints(
	DACE::vectordb const& x_goal,
	PNSolver const& PNsolver);

void run_test_cases(int argc, char** argv);

#endif