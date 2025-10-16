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

// Misc
void save_control(
	string const& file_name,
	vector<vectordb> const& list_u) {
	// Open file
	string file_name_(file_name);
	ofstream ofs(file_name_);

	// Set double precision
	typedef std::numeric_limits<double> dbl;
	ofs.precision(dbl::max_digits10);


	// Store the object to file
	size_t N = list_u.size();
	ofs << N << endl;

	for (size_t i=0; i<N; i++) {
		ofs << list_u[i];
	}

	ofs.close();
	return;
}
vector<vectordb> load_control(std::string const& file_name) {

	// Open file
	string file_name_(file_name);
	ifstream ifs(file_name_);

	// Load data
	string N_str; size_t N = 1;
	getline(ifs, N_str); 
	istringstream(N_str) >> N;

	vector<vectordb> list_u(N);
	for (size_t i=0; i<N; i++) {
		ifs >> list_u[i];
	}

	ifs.close();
	return list_u;
}

pair<vector<vectordb>, vector<vectordb>> load_dataset(
	string const& file_name,
	double const& ToF,
	Dynamics const& dynamics,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {

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

	cout << file_name_ << endl;

	// Open file
	ifstream ifs(file_name_);

	// Skip preamble
	string buffer_str;
	for (size_t i=0; i<13; i++) {getline(ifs, buffer_str);}

	// Skip states
	matrixdb buff_mat;
	vector<vectordb> list_x, list_u;

	// Read state
	buff_mat = read_dataset(ifs);
	for (size_t i=0; i<buff_mat.nrows(); i++) {
		list_x.emplace_back(buff_mat.getrow(i));
	}

	// Read control
	buff_mat = read_dataset(ifs);
	for (size_t i=0; i<buff_mat.nrows(); i++) {
		list_u.emplace_back(buff_mat.getrow(i));
	}

	ifs.close();
	return pair<vector<vectordb>, vector<vectordb>>(list_x, list_u);
}

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
	if (test_case == 0) {
		double_integrator(argc, argv);
	}

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
