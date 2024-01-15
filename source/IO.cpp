/**
	IO.cpp

	Purpose: Implementation of the input and outputs of data

	@author Thomas Caleb

	@version 1.0 10/01/2024
*/

#include "IO.h"

using namespace DACE;
using namespace std;

// Function to print a dataset at a given name in order to
// produce python visuals
void print_dataset(
	string const& file_name,
	string const& system_name,
	SpacecraftParameters const& spacecraft_parameters,
	vector<vector<string>> const& list_title,
	vector<vector<vectordb>> const& list_data) {

	// Open file
	ofstream ofs(file_name);

	// Print header
	ofs << file_name << endl;
	ofs << system_name << endl;
	ofs << spacecraft_parameters;

	// Print lists
	
	size_t nb_lists = list_data.size();
	ofs << nb_lists << endl; // Number of lists
	
	for (size_t i = 0; i < nb_lists; i++) {
		// Unpack
		vector<string> list_title_i = list_title[i];
		vector<vectordb> list_data_i = list_data[i];
		size_t nb_rows = list_data_i.size();
		size_t nb_colomns = list_data_i[0].size();

		// Print

		// Headers
		ofs << list_title_i[0] << endl; // Series title
		ofs << nb_rows << ", " << nb_colomns << endl; // Dataset size
		
		for (size_t l = 0; l < nb_colomns; l++) {
			ofs << list_title_i[1 + l] << ", ";
		}
		ofs << endl;
		
		// Data
		for (size_t k = 0; k < nb_rows; k++) {
			for (size_t l = 0; l < nb_colomns; l++) {
				ofs << list_data_i[k][l] << ", ";
			}
			ofs << endl;
		}
	}
	ofs.close();
	return;
}

// Function to propagate a vector without control
vector<vectordb> get_reference_trajectory(
	vectordb const& x_0,
	Dynamics const& dynamics,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters,
	int const& nb_point) {
	// Unpack
	int Nu = solver_parameters.Nu();
	int Nx = solver_parameters.Nx();
	double period = x_0[Nx - 1];

	// Init
	vector<vectordb> output;
	vectordb x_i = x_0;
	double dt = period / (nb_point - 1.0);
	x_i[Nx - 1] = dt;
	vectordb null_control(Nu, 0.0);
	
	// Loop
	for (size_t i = 0; i < nb_point; i++) {
		output.push_back(x_i);
		x_i = dynamics.dynamic_db()(
			x_i, null_control,
			spacecraft_parameters, constants, solver_parameters);
	}

	return output;
}

// Prints a transfer dataset on a standardised format
// with reference orbits
void print_transfer_dataset(
	string const& file_name,
	string const& system_name,
	vector<vectordb> const& list_x, vector<vectordb> const& list_u,
	vectordb const& x_0, vectordb const& x_f,
	Dynamics const& dynamics,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {
	// Unpack
	int Nu = solver_parameters.Nu();
	int Nx = solver_parameters.Nx();
	int N = solver_parameters.N();

	// Get reference trajectories
	vector<vectordb> list_departure = get_reference_trajectory(
		x_0, dynamics,
		spacecraft_parameters, constants, solver_parameters,
		2*N);
	vector<vectordb> list_arrival = get_reference_trajectory(
		x_f, dynamics,
		spacecraft_parameters, constants, solver_parameters,
		2*N);

	// Make lists
	vector<string> title_state{
		"State",
		"x [LU]", "y [LU]", "z [LU]",
		"vx [VU]", "vy [VU]", "vz [VU]",
		"mass [MASSU]", "dt [TU]" };
	vector<string> title_control{
		"Control",
		"ux [THRUSTU]", "uy [THRUSTU]", "uz [THRUSTU]"};
	vector<string> title_departure(title_state), title_arrival(title_state);
	title_departure[0] = "Departure orbit";
	title_arrival[0] = "Arrival orbit";
	vector<vector<string>> list_title{
		title_state , title_control, title_departure, title_arrival };
	vector<vector<vectordb>> list_data{
		list_x, list_u,  list_departure, list_arrival };

	print_dataset(
		file_name, system_name,
		spacecraft_parameters,
		list_title, list_data);
}