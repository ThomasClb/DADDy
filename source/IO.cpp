/**
	IO.cpp

	Purpose: Implementation of the input and outputs of data.

	@author Thomas Caleb

	@version 1.0 10/01/2024
*/

#include "IO.h"

using namespace DACE;
using namespace std;

// Splits a string into substring.
vector<string> split(string s, string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    string token;
    vector<string> res;
    while ((pos_end = s.find(delimiter, pos_start)) != string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }
    res.push_back (s.substr (pos_start));
    return res;
}


// Function to print a dataset at a given name in order to
// produce python visuals.
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

// Reads a dataset
matrixdb read_dataset(ifstream& ifs) {
	// Init
	string buffer_str;
	string delimiter(", ");

	// Skip header
	getline(ifs, buffer_str);

	// Get n_rows, n_cols
	getline(ifs, buffer_str);
	vector<string> list_str = split(buffer_str, delimiter);
	unsigned int n_rows(stoi(list_str[0])), n_cols(stoi(list_str[1]));

	// Skip header
	getline(ifs, buffer_str);

	// Init
	matrixdb output(n_rows, n_cols);

	// Get nominal controls
	vectordb row(n_cols);
	for (size_t i=0; i<n_rows; i++) {
		// Get string
		getline(ifs, buffer_str);
		list_str = split(buffer_str, delimiter);

		// Assign
		for (size_t j=0; j<n_cols; j++) {
			row[j] = stod(list_str[j]);
			
		}
		output.setrow(i, row);
	}
	return output;
}

// Function to propagate a vector without control.
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
// with reference orbits.
void print_transfer_dataset(
	string const& file_name,
	string const& system_name,
	vector<vectordb> const& list_x, vector<vectordb> const& list_u,
	vectordb const& x_0, vectordb const& x_f,
	double const& ToF,
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

	// Convert to cartesian if needed
	if (system_name.find("EQUINOCTIAL") != std::string::npos) {
		
		// Get mu
		double mu(1);
		if (system_name.find("EARTH") != std::string::npos) {
			mu = MU_EARTH;
		}
		else if (system_name.find("SUN") != std::string::npos) {
			mu = MU_SUN;
		}
		mu /= constants.mu(); // Normalise
			
		for (size_t i = 0; i < list_departure.size(); i++) {
			// Unpack
			vectordb x_d_i = list_departure[i];
			vectordb x_a_i = list_arrival[i];

			// Convert to cartesian
			x_d_i = equi_2_kep(x_d_i); // Convert to Keplerian
			x_d_i = kep_2_cart(x_d_i, mu); // Convert to Cartesian
			x_a_i = equi_2_kep(x_a_i); // Convert to Keplerian
			x_a_i = kep_2_cart(x_a_i, mu); // Convert to Cartesian

			// Assign
			list_departure[i] = x_d_i;
			list_arrival[i] = x_a_i;
		}
	}

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


	// Make file name
	int power = 9;
	double thrust_mass = ((spacecraft_parameters.thrust() * spacecraft_parameters.constants().thrustu()
		/ (spacecraft_parameters.initial_mass() * spacecraft_parameters.constants().massu())
		) * pow(10.0, power));
	int exposant = log10(thrust_mass);
	int mantisse = static_cast<int>(thrust_mass) / static_cast<int>(pow(10, exposant));
	string str_T2m = to_string(mantisse) + "e" + to_string(exposant - power);
	string file_name_ = file_name + "_"
		+ str_T2m + "_"
		+ to_string((int)(ToF*spacecraft_parameters.constants().tu()*SEC2DAYS)) + "_"
		+ to_string(solver_parameters.DDP_type())
		+ ".dat";

	print_dataset(
		file_name_, system_name,
		spacecraft_parameters,
		list_title, list_data);
}
