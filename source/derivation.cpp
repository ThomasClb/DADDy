/**
	derivation.cpp

	Purpose: Implementation of the DA-based
	automatic differentiation methods.

	@author Thomas Caleb

	@version 1.0 14/11/2023
*/

#include "derivation.h"

using namespace DACE;
using namespace std;

// DACE vector conversion

// Turns a vectordb into a vectorDA without perturbation
vectorDA db_2_DA(
	vectordb const& vector_db) {
	size_t n = vector_db.size();
	vectorDA vector_DA; vector_DA.reserve(n);
	for (size_t i = 0; i <n; i++) {
		vector_DA.push_back(vector_db[i]);
	}
	return vector_DA;
}

// Turns a vectordb into a vectorDA with identity perturbation
vectorDA id_vector(
	vectordb const& vector_db,
	unsigned int const& begin_index,
	unsigned int const& begin_da_index,
	unsigned int const& end_da_index) {
	
	
	size_t n = vector_db.size();
	vectorDA vector_DA; vector_DA.reserve(vector_db.size());
	for (unsigned int i = 0; i < n; i++) {

		// Copy constant part
		DA coord_i = vector_db[i];

		// Add linear terms starting from begin_index
		// From variable begin_da_index to end_da_index
		if (i >= begin_index && i < begin_index + end_da_index - begin_da_index)
			coord_i += DA(begin_da_index + i - begin_index + 1);

		// Assign
		vector_DA.push_back(coord_i);
	}

	return vector_DA;
}

// Differentiation functions

// Differentiates with respect to x, and u
vector<matrixdb> deriv_xu(
	vectorDA const& f_eval,
	unsigned int const& Nx, unsigned int const& Nu,
	bool const& hessian_computation) {
	// Unpack
	size_t n = f_eval.size();

	if (!hessian_computation) {
		if (n == 0) {
			return vector<matrixdb>(2);
		}

		// Get global jacobian
		matrixdb jac_f = f_eval.linear();

		// First order derivatives
		matrixdb fx = jac_f.submat(0, 0, n - 1, Nx - 1);
		matrixdb fu = jac_f.submat(0, Nx, n - 1, Nx + Nu - 1);

		return vector<matrixdb>{ fx, fu };
	}
	else {
		if (f_eval.size() == 0) {
			return vector<matrixdb>(2);
		}

		// Make output
		vector<matrixdb> output(2 + 3 * n);

		// Build DA jacobian
		matrixDA df_eval(n, Nx + Nu);
		for (size_t i = 0; i < Nx + Nu; i++) {
			df_eval.setcol(i, f_eval.deriv(i + 1));
		}
		matrixdb jac_f = df_eval.cons();

		// First order derivatives
		matrixdb fx = jac_f.submat(0, 0, n - 1, Nx - 1);
		matrixdb fu = jac_f.submat(0, Nx, n - 1, Nx + Nu - 1);

		// Assign
		output[0] = fx; output[1] = fu;

		// Second order derivatives
		for (size_t i = 0; i < n; i++) {
			matrixdb hessian = vectorDA(df_eval.getrow(i)).linear();

			// Slice hesssians
			output[2 + i] = hessian.submat(0, 0, Nx - 1, Nx - 1); // fxx
			output[2 + i + n] = hessian.submat(Nx, Nx, Nx + Nu - 1, Nx + Nu - 1); // fuu
			output[2 + i + 2*n] = hessian.submat(Nx, 0, Nx + Nu - 1, Nx - 1); // fux
		}

		// Output
		return output;
	}
}

// Differentiates with respect to x
vector<matrixdb> deriv_x(
	vectorDA const& f_eval,
	unsigned int const& Nx,
	bool const& hessian_computation) {
	// Unpack
	size_t n = f_eval.size();

	if (!hessian_computation) {
		// Get global jacobian
		matrixdb jac_f = f_eval.linear();

		// First order derivatives
		matrixdb fx = jac_f.submat(0, 0, n - 1, Nx - 1);
		return vector<matrixdb>{fx};
	}
	else {
		// Init
		vector<matrixdb> output; 
		output.reserve(1 + n);

		// Build DA jacobian
		matrixDA df_eval(n, Nx);
		for (size_t i = 0; i < Nx; i++) {
			df_eval.setcol(i, f_eval.deriv(i + 1));
		}
		matrixdb jac_f = df_eval.cons();

		// First order derivatives
		output.emplace_back(jac_f);

		// Add hessians
		for (size_t i = 0; i < n; i++) {
			output.emplace_back(vectorDA(
				df_eval.getrow(i)
			).linear().submat(Nx - 1, Nx - 1));
		}
		return output;
	}
}
