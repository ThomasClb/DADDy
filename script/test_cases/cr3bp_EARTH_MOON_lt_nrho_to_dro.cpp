/**
	cr3bp_EARTH_MOON_lt_nrho_to_dro.cpp.cpp

	Purpose: Low-thrust L2 NRHO to DRO transfer execution script.
	In the Earth-Moon system.

	@author Thomas Caleb

	@version 1.0 07/12/2023
*/

#include "test_cases.h"

using namespace DACE;
using namespace std::chrono;
using namespace std;

SolverParameters get_SolverParameters_cr3bp_EARTH_MOON_lt_nrho_to_dro(
	unsigned int const& N, unsigned int const& DDP_type,
	unsigned int verbosity) {
	// Solver parameters
	unsigned int Nx = (SIZE_VECTOR + 1) + 1;
	unsigned int Nu = SIZE_VECTOR / 2;
	unsigned int Neq = 0;
	unsigned int Nineq =  4;
	unsigned int Nteq = 6;
	unsigned int Ntineq = 0;
	bool with_J2 = false;
	double cost_to_go_gain = 1e-4;
	double terminal_cost_gain = 1e6;
	double mass_leak = 1e-8;
	double homotopy_coefficient = 0.0;
	double huber_loss_coefficient = 5e-4;
	vectordb homotopy_sequence{0, 0.75, 0.95, 0.999};
	vectordb huber_loss_coefficient_sequence{ 1e-2, 5e-3, 5e-3, 1e-4};
	double DDP_tol = 1e-4;
	double AUL_tol = 1e-6;
	double PN_tol = 1e-10;
	double PN_active_constraint_tol = 1e-13;
	unsigned int max_iter = 10000;
	unsigned int DDP_max_iter = 200;
	unsigned int AUL_max_iter = max_iter / DDP_max_iter;
	unsigned int PN_max_iter = 400;
	vectordb lambda_parameters{0.0, 1e8};
	vectordb mu_parameters{1, 1e8, 10};
	vectordb line_search_parameters{1e-8, 10.0, 0.5, 20};
	bool backward_sweep_regulation = true;
	vectordb backward_sweep_regulation_parameters{0, 1e-8, 1e8, 1.6};
	double PN_regularisation(1e-8);
	double PN_cv_rate_threshold(1.1);
	double PN_alpha(1.0); double PN_gamma(0.5);
	unsigned int saving_iterations = 0;

	return SolverParameters(
		N, Nx, Nu,
		Neq, Nineq,
		Nteq, Ntineq, with_J2,
		cost_to_go_gain, terminal_cost_gain, mass_leak,
		homotopy_coefficient, huber_loss_coefficient,
		homotopy_sequence,
		huber_loss_coefficient_sequence,
		DDP_type,
		DDP_tol, AUL_tol, PN_tol,
		DDP_max_iter, AUL_max_iter, PN_max_iter,
		line_search_parameters,
		backward_sweep_regulation,
		backward_sweep_regulation_parameters,
		lambda_parameters, mu_parameters,
		PN_regularisation, PN_active_constraint_tol,
		PN_cv_rate_threshold, PN_alpha, PN_gamma,
		verbosity, saving_iterations);
}

void cr3bp_EARTH_MOON_lt_nrho_to_dro(int argc, char** argv) {
	// Input check
	if (argc < 10) {
		cout << "Wrong number of arguments." << endl;
		cout << "Requested number : 9" << endl;
		cout << "0 - Test case number." << endl;
		cout << "1 - SpacecraftParameter adress." << endl;
		cout << "2 - DDP type [0-7]." << endl;
		cout << "3 - Number of nodes [-]." << endl;
		cout << "4 - Time of flight [days]." << endl;
		cout << "5 - Perform fuel optimal optimisation [0/1]." << endl;
		cout << "6 - Perform projected Newton solving [0/1]." << endl;
		cout << "7 - Save results [0/1]." << endl;
		cout << "8 - Verbosity [0-2]." << endl;
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

	// Set double precision
	typedef std::numeric_limits<double> dbl;
	cout.precision(7);

	// Get dynamics
	Dynamics dynamics = get_cr3bp_EARTH_MOON_lt_dynamics();

	// Normalisation cosntants
	Constants constants(dynamics.constants());
	double lu = constants.lu();
	double massu = constants.massu();
	double tu = constants.tu();
	double thrustu = constants.thrustu();
	double vu = constants.vu();

	// Spacecraft parameters
	SpacecraftParameters spacecraft_parameters(spacecraft_parameters_file);

	// Init solver parameters
	SolverParameters solver_parameters = get_SolverParameters_cr3bp_EARTH_MOON_lt_nrho_to_dro(
		N, DDP_type, verbosity);

	// Solver parameters
	unsigned int Nx = solver_parameters.Nx();
	unsigned int Nu = solver_parameters.Nu();

	// Initial conditions [3*LU, 3*VU, MASSU, TU]
	// From [Boone MacMahon 2024]
	ToF = ToF / SEC2DAYS / tu; // [TU]
	double dt = ToF / N; // [TU]
	vectordb x_departure{ 
		1.021968177, 0,
		-0.18206, 0,
		-0.103140143, 0,
		spacecraft_parameters.initial_mass(),
		2.5748200748171399 };
	vectordb x_arrival{
		0.983368093, 0.259208967, 0,
		0.351341295, -0.008333464, 0,
		spacecraft_parameters.dry_mass(),
		3.2746644337639852 };
	vectordb x0 = x_departure; x0[Nx - 1] = dt; // Time step
	vectordb x_goal = x_arrival; x_goal[Nx - 1] = ToF; // ToF

	// First guess command
	vectordb u_init(Nu, 1e-4); // [VU]
	vector<vectordb> list_u_init(N, u_init);

	// Solver
	DADDy solver(solver_parameters, spacecraft_parameters, dynamics);
	solver.solve(x0, list_u_init, x_goal, fuel_optimal, pn_solving);
	
	// Main outputs
	if (verbosity == 3) {
		// ID
		cout << atoi(argv[1]) << ", ";
		cout << spacecraft_parameters.thrust()*thrustu/(spacecraft_parameters.initial_mass()*massu) << ", ";
		cout << ToF*tu*SEC2DAYS << ", ";
		cout << DDP_type << ", ";

		// Data

		// Results
		cout << solver.list_x()[N][SIZE_VECTOR]*massu << ", ";
		cout << solver.list_x()[N][SIZE_VECTOR]/spacecraft_parameters.initial_mass() << ", ";
		cout << solver.real_constraints(x_goal) << ", ";

		// Convergence metrics
		cout << solver.AUL_runtime() << ", ";
		cout << solver.PN_runtime() << ", ";
		cout << solver.runtime() << ", ";
		cout << solver.DDP_n_iter() << ", ";
		cout << solver.AUL_n_iter() << ", ";
		cout << solver.PN_n_iter() << endl;
	}
	
	// Print datasets
	if (save_results) {
		string file_name = "./data/datasets/cr3bp_EARTH_MOON_lt_nrho_to_dro";
		string system_name = "CR3BP EARTH-MOON CARTESIAN LT";
		print_transfer_dataset(
			file_name, system_name,
			solver.list_x(), solver.list_u(),
			x_departure, x_arrival, ToF,
			dynamics, spacecraft_parameters, constants, solver_parameters);
	}
}
