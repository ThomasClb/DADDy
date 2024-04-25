/**
	tbp_EARTH_lt_meo_to_meo.cpp

	Purpose: Low-thrust MEO to MEO transfer execution script.

	@author Thomas Caleb

	@version 1.0 29/01/2023
*/

#include "test_cases.h"

using namespace DACE;
using namespace std::chrono;
using namespace std;

SolverParameters get_SolverParameters_tbp_EARTH_lt_meo_to_meo(
	unsigned int const& N, unsigned int const& DDP_type,
	unsigned int verbosity) {
	// Solver parameters
	unsigned int Nx = (SIZE_VECTOR + 1) + 1;
	unsigned int Nu = SIZE_VECTOR / 2;
	unsigned int Neq = 0;
	unsigned int Nineq = 2;
	unsigned int Nteq = 5;
	unsigned int Ntineq = 0;
	bool with_J2 = false;
	double cost_to_go_gain = 1e-5;
	double terminal_cost_gain = 1e5;
	double mass_leak = 1e-8;
	double homotopy_coefficient = 0.0;
	double huber_loss_coefficient = 5e-3;
	vectordb homotopy_sequence{0, 0.5, 0.99};
	vectordb huber_loss_coefficient_sequence{1e-2, 1e-2, 5e-3};
	double DDP_tol = 1e-4;
	double AUL_tol = 1e-6; 
	double PN_tol = 1e-12;
	double PN_active_constraint_tol = 1e-13;
	unsigned int DDP_max_iter = 200;
	unsigned int AUL_max_iter = 200;
	unsigned int PN_max_iter = 200;
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

void tbp_EARTH_lt_meo_to_meo(int argc, char** argv) {
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
	cout.precision(5);

	// Set dynamics
	Dynamics dynamics = get_tbp_EARTH_lt_dynamics();
	
	// Normalisation constants
	Constants constants(dynamics.constants());
	double lu = constants.lu();
	double mu = constants.mu();
	double massu = constants.massu();
	double tu = constants.tu();
	double thrustu = constants.thrustu();
	double vu = constants.vu();

	// Spacecraft parameters
	SpacecraftParameters spacecraft_parameters(spacecraft_parameters_file);

	// Init solver parameters
	SolverParameters solver_parameters = get_SolverParameters_tbp_EARTH_lt_meo_to_meo(
		N, DDP_type, verbosity);

	// Solver parameters
	unsigned int Nx = solver_parameters.Nx();
	unsigned int Nu = solver_parameters.Nu();

	// Initial conditions [Equinoctial elements, MASSU, TU]
	ToF = ToF / SEC2DAYS / tu; // [TU]
	double dt = ToF / N; // [TU]
	double altitude = 28000;
	double r_d = R_EARTH + altitude;
	double r_a = r_d + 0;
	vectordb x_departure{ // Kep coordinates
		r_d / lu, 0,
		60 * DEG_2_RAD, 180 * DEG_2_RAD,
		0 * DEG_2_RAD, 0 * DEG_2_RAD,
		spacecraft_parameters.initial_mass(), 2 * PI * sqrt(pow(r_d, 3) / MU_EARTH) / tu };
	vectordb x_arrival{ // Kep coordinates
		r_a / lu, 0,
		60 * DEG_2_RAD, 155 * DEG_2_RAD,
		0 * DEG_2_RAD, 0 * DEG_2_RAD,
		spacecraft_parameters.dry_mass(), 2 * PI * sqrt(pow(r_a, 3) / MU_EARTH) / tu };
	x_departure = kep_2_equi(x_departure); // Equinoctial coordinates
	x_arrival = kep_2_equi(x_arrival);
	vectordb x0 = x_departure; x0[Nx - 1] = dt; // Time step
	vectordb x_goal = x_arrival; x_goal[Nx - 1] = ToF; // ToF

	// First guess command
	vectordb u_init(Nu, 1e-6 / thrustu); // [VU]
	// u_init[1] = 0.1 * T;
	vector<vectordb> list_u_init(N, u_init);

	// Solver
	DADDy solver(solver_parameters, spacecraft_parameters, dynamics);
	auto start = high_resolution_clock::now();
	solver.solve(x0, list_u_init, x_goal, fuel_optimal, pn_solving);
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	
	// Main outputs
	if (verbosity == 3) {
		cout << atoi(argv[1]) << " - ";
		cout << DDP_type << " - ";
		cout << spacecraft_parameters.thrust()*thrustu/(spacecraft_parameters.initial_mass()*massu) << " - ";
		cout << ToF*tu*SEC2DAYS << " - ";
		cout << solver.list_x()[N][SIZE_VECTOR]*massu << " - ";
		cout << solver.list_x()[N][SIZE_VECTOR]/spacecraft_parameters.initial_mass() << " - ";
		cout << solver.real_constraints(x_goal) << " - ";
		cout << to_string(static_cast<double>(duration.count()) / 1e6) << endl;
	}

	// Unpack
	vector<vectordb> list_x = solver.list_x();
	vector<vectordb> list_u = solver.list_u();

	// Print datasets
	if (save_results) {

		// Convert to cart
		mu = mu / constants.mu();
		for (size_t i = 0; i < list_u.size(); i++) {
			// Unpack
			vectordb x_i = list_x[i];
			vectordb u_i = list_u[i];

			// Convert to cartesian
			x_i = equi_2_kep(x_i); // Convert to Keplerian
			x_i = kep_2_cart(x_i, mu); // Convert to Cartesian
			u_i = RTN_2_cart(u_i, x_i);

			// Assign
			list_u[i] = u_i;
			list_x[i] = x_i;
		}
		list_x[list_u.size()] = kep_2_cart(equi_2_kep(list_x[list_u.size()]), mu);

		string file_name = "./data/datasets/tbp_EARTH_lt_meo_to_meo";
		string system_name = "TBP EARTH EQUINOCTIAL LT";
		print_transfer_dataset(
			file_name, system_name,
			list_x, list_u,
			x_departure, x_arrival,
			dynamics, spacecraft_parameters, constants, solver_parameters);
	}
}
