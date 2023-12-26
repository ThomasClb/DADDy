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

void low_thrust_2bp_SUN();
void low_thrust_cr3bp();


double real_constraints(
	DACE::vectordb const& x_goal,
	PNSolver const& PNsolver);

#endif