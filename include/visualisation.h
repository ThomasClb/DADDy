/**
	visualisation.h

	Purpose: Implementation of the visulasation tools.

	@author Thomas Caleb

	@version 1.0 17/11/2023
*/

#ifndef DEF_VISUALISATION
#define DEF_VISUALISATION

#pragma once

#include <vector>
#include <string>

#include <dace/dace_s.h>

// Axis label constants
const std::vector<std::string> axis_label{"x [LU]", "y [LU]", "z [LU]" };

// The following functions are wrappers of the library matplotlibcpp
// See: https://github.com/gmrukwa/matplotlib-cpp
void figure();
void show();
void legend();
void save(std::string name);
void xlabel(std::string name);
void ylabel(std::string name);
void plot(std::vector<double> x, std::vector<double> y);
void plot(std::vector<double> x, std::vector<double> y, std::string s);
void named_plot(std::string name, std::vector<double> x, std::vector<double> y);

#endif