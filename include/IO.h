/**
	IO.h

	Purpose: Implementation of the input and outputs of data

	@author Thomas Caleb

	@version 1.0 10/01/2024
*/

#ifndef DEF_IO
#define DEF_IO

#pragma once

#include <vector>
#include <string>

#include <dace/dace_s.h>
#include "parameters.h"
#include "dynamics.h"

// Function to print a dataset at a given name in order to
// produce python visuals
void print_dataset(
	std::string const& file_name,
	std::string const& system_name,
	SpacecraftParameters const& spacecraft_parameters,
	std::vector<std::vector<std::string>> const& list_title,
	std::vector<std::vector<DACE::vectordb>> const& list_data);

// Function to propagate a vector without control
std::vector<DACE::vectordb> get_reference_trajectory(
	DACE::vectordb const& x_0,
	Dynamics const& dynamics,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters,
	int const& nb_point);

// Prints a transfer dataset on a standardised format
// with reference orbits
void print_transfer_dataset(
	std::string const& file_name,
	std::string const& system_name,
	std::vector<DACE::vectordb> const& list_x,
	std::vector<DACE::vectordb> const& list_u,
	DACE::vectordb const& x_0, DACE::vectordb const& x_f,
	double const& ToF,
	Dynamics const& dynamics,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters);

#endif