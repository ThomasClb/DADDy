/**
	test_cases.cpp

	Purpose: Executes a test case.

	@author Thomas Caleb

	@version 1.0 07/12/2023
*/

#include "test_cases.h"

using namespace DACE;
using namespace std::chrono;
using namespace std;

// Runs the selected test case.
void run_test_cases(int argc, char** argv) {

	// Input check
	if (argc < 2) {
		cout << "Wrong number of arguments." << endl;
		cout << "Requested number : 1" << endl;
		cout << "0 - Test case number from 0 to 8." << endl;
		return;
	}

	// Unpack inputs
	unsigned int test_case = atoi(argv[1]);

	// Excute correct test case

	// Double integrator
	if (test_case == 0) {}

	// TBP
	
	// SUN centered

	// Earth-mars transfer
	else if (test_case == 1)
		tbp_SUN_lt_earth_to_mars(argc, argv);

	// Earth centered

	// GTO to GEO
	else if (test_case == 2) {
		tbp_EARTH_lt_leo_to_geo(argc, argv);
	}

	// LEO to LEO
	else if (test_case == 3) {
		tbp_EARTH_lt_leo_to_leo(argc, argv);
	}

	// MEO to MEO
	else if (test_case == 4) {
		tbp_EARTH_lt_meo_to_meo(argc, argv);
	}

	// CR3BP

	// Earth-Moon

	// Halo L2 to Halo L1
	else if (test_case == 5)
		cr3bp_EARTH_MOON_lt_haloL2_to_haloL1(argc, argv);

	// Halo L2 (Gateway) to DRO
	else if (test_case == 6) {
		cr3bp_EARTH_MOON_lt_nrho_to_dro(argc, argv);
	}

	// Lyapunov L2 to Lyapunov L1
	else if (test_case == 7) {
		cr3bp_EARTH_MOON_lt_lyapunovL1_to_lyapunovL2(argc, argv);
	}

	// Manyrev DRO to DRO
	else if (test_case == 8) {
		cr3bp_EARTH_MOON_lt_dro_to_dro(argc, argv);
	}

	// Manyrev GSO to NRHO
	else if (test_case == 9) {}
}