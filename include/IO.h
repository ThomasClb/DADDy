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

// Function to print a dataset at a given name in order to
// produce python visuals
void print_dataset(
	std::string const& file_name,
	std::string const& system_name,
	SpacecraftParameters const& spacecraft_parameters,
	std::vector<std::vector<std::string>> const& list_title,
	std::vector<std::vector<DACE::vectordb>> const& list_data);

#endif