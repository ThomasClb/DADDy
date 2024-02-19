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
	unsigned int test_case = atoi(argv[1]);
	unsigned int nb_iteration = atoi(argv[10]);

	// Output
	cout << "Test case : " + to_string(test_case) << endl;
	cout << "Number of computations : " + to_string(nb_iteration) << endl;

	// Execute correct test case
	auto start = high_resolution_clock::now();
	for (size_t i = 0; i < nb_iteration; i++) {
		run_test_cases(argc, argv);
	}
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	double duration_db = static_cast<double>(duration.count()) / (1e6);

	// Output
	cout << "Total runtime : " + to_string(duration_db) + "s" << endl;
	cout << "Average runtime : " + to_string(duration_db / nb_iteration ) + "s" << endl;
}