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
