/**
	visualisation.cpp

	Purpose: Implementation of the visulasation tools.

	@author Thomas Caleb

	@version 1.0 17/11/2023
*/

#include "visualisation.h"

#include <matplotlibcpp.h> // Matplotlib must be in .cpp (IDK why)

using namespace DACE;
using namespace std;
namespace plt = matplotlibcpp;

// The following functions are wrappers of the library matplotlibcpp
// See: https://github.com/gmrukwa/matplotlib-cpp
void figure() { plt::figure(); }
void show() { plt::show(); }
void save(string name) { plt::save(name); }
void legend() { plt::legend();}
void xlabel(string name) { plt::xlabel(name);}
void ylabel(string name) { plt::ylabel(name);}
void plot(std::vector<double> x, std::vector<double> y) { plt::plot(x, y); }
void plot(std::vector<double> x, std::vector<double> y, std::string s) { plt::plot(x, y, s); }
void named_plot(std::string name, std::vector<double> x, std::vector<double> y) { plt::named_plot(name, x, y); }
