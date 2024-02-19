/**
	runtime_analysis.h

	Purpose: Preforms a runtime analysis on test cases.

	@author Thomas Caleb

	@version 1.0 19/02/2024
*/

#ifndef DEF_RUNTIME_ANALYSIS
#define DEF_RUNTIME_ANALYSIS

#pragma once

#include <dace/dace_s.h>
#include <chrono>

#include "test_cases.h"

void run_runtime_analysis(int argc, char** argv);

#endif