/**
	derivation.h

	Purpose: Implementation of the DA-based
	automatic differentiation methods.

	@author Thomas Caleb

	@version 1.0 14/11/2023
*/

#ifndef DEF_DERIVATION
#define DEF_DERIVATION

#pragma once

#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include <dace/dace_s.h>

// DACE vector conversion

// Turns a vectordb into a vectorDA without perturbation
DACE::vectorDA db_2_DA(
	DACE::vectordb const& vector_db);

// Turns a vectordb into a vectorDA with identity perturbation
DACE::vectorDA id_vector(
	DACE::vectordb const& vector_db,
	unsigned int const& begin_index,
	unsigned int const& begin_da_index,
	unsigned int const& end_da_index);

// Differentiation functions

// Differentiates with respect to x, and u
std::vector<DACE::matrixdb> deriv_xu(
	DACE::vectorDA const& f_eval,
	unsigned int const& Nx, unsigned int const& Nu,
	bool const& hessian_computation);

// Differentiates with respect to x
std::vector<DACE::matrixdb> deriv_x(
	DACE::vectorDA const& f_eval,
	unsigned int const& Nx,
	bool const& hessian_computation);

#endif