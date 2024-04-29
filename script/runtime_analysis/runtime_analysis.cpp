/**
	runtime_analysis.cpp

	Purpose: Preforms a runtime analysis on test cases.

	@author Thomas Caleb

	@version 1.0 19/02/2024
*/

#include "runtime_analysis.h"

using namespace DACE;
using namespace std::chrono;
using namespace std;

void run_runtime_analysis(int argc, char** argv) {

	// Input check
	if (argc < 11) {
		cout << "Wrong number of arguments." << endl;
		cout << "Requested number : 9 + 1" << endl;
		cout << "0-8 - Test case parameters." << endl;
		cout << "9 - Number of iterations." << endl;
		return;
	}

	// Unpack inputs
	string spacecraft_parameters_file = argv[2];
	unsigned int DDP_type = atoi(argv[3]);
	unsigned int N = atoi(argv[4]);
	double ToF = atof(argv[5]);
	bool fuel_optimal = false;
	bool pn_solving = false;
	bool save_results = false;
	int verbosity = atoi(argv[9]);
	if (atoi(argv[6]) == 1) { fuel_optimal = true; }
	if (atoi(argv[7]) == 1) { pn_solving = true; }
	if (atoi(argv[8]) == 1) { save_results = true; }
	unsigned int nb_iteration = atoi(argv[10]);

	// Spacecraft parameters
	SpacecraftParameters spacecraft_parameters(spacecraft_parameters_file);

	// Normalisation cosntants
	Constants constants(spacecraft_parameters.constants());
	double lu = constants.lu();
	double massu = constants.massu();
	double tu = constants.tu();
	double thrustu = constants.thrustu();
	double vu = constants.vu();
	
	// Output

	// Execute correct test case
	vectordb list_rt(nb_iteration, 0.0);
	auto start = high_resolution_clock::now();
	for (size_t i = 0; i < nb_iteration; i++) {
		auto start_i = high_resolution_clock::now();
			run_test_cases(argc, argv);
		auto stop_i = high_resolution_clock::now();
		auto duration_i = duration_cast<microseconds>(stop_i - start_i);
		double duration_db_i = static_cast<double>(duration_i.count()) / (1e6);
		list_rt[i] = duration_db_i;
	}
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	double duration_db = static_cast<double>(duration.count()) / (1e6);
	double mean_rt = duration_db / nb_iteration;

	// Get std and min, max
	double std_db = 0.0;
	double min_rt(1e15), max_rt(0.0);
	for (size_t i = 0; i < nb_iteration; i++) {
		std_db += sqr(list_rt[i] - mean_rt);

		if (list_rt[i] > max_rt) {
			max_rt = list_rt[i];
		}
		if (list_rt[i] < min_rt) {
			min_rt = list_rt[i];
		}
	}
	std_db = sqrt(std_db/nb_iteration);

	// Set double precision
	typedef std::numeric_limits<double> dbl;
	cout.precision(7);


	// Output

	// ID
	cout << atoi(argv[1]) << ", ";
	cout << spacecraft_parameters.thrust()*thrustu/(spacecraft_parameters.initial_mass()*massu) << ", ";
	cout << ToF << ", ";
	cout << DDP_type << ", ";

	// Data
	cout << nb_iteration << ", "; // Nb of interation
	cout << duration_db << ", "; // Total run time
	cout << mean_rt << ", "; // mean run time
	cout << std_db << ", "; // Std
	cout << min_rt << ", "; // Min run time
	cout << max_rt << endl; // Max run time
}