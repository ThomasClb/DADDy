/**
	tbp_EARTH_lt_leo_to_geo.cpp

	Purpose: Low-thrust LEO to GEO transfer execution script.

	@author Thomas Caleb

	@version 1.0 23/01/2023
*/

#include "test_cases.h"

using namespace DACE;
using namespace std::chrono;
using namespace std;

SolverParameters get_SolverParameters_tbp_EARTH_lt_leo_to_geo(
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
	double cost_to_go_gain = 1e-4; // 1e-5
	double terminal_cost_gain = 1e5; // 1e-5
	double mass_leak = 1e-8;
	double homotopy_coefficient = 0.0;
	double huber_loss_coefficient = 5e-3;
	vectordb homotopy_sequence{0, 0.5, 0.9, 0.99};
	vectordb huber_loss_coefficient_sequence{1e-2, 1e-2, 5e-3, 1e-3};

	homotopy_sequence = vectordb{0, 0.5};
	huber_loss_coefficient_sequence = vectordb{1e-2, 1e-2};

	double DDP_tol = 1e-4;
	double AUL_tol = 1e-6; 
	double PN_tol = 1e-10;
	double PN_active_constraint_tol = 1e-11;
	unsigned int DDP_max_iter = 400;
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

void tbp_EARTH_lt_leo_to_geo(int argc, char** argv) {
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
	SolverParameters solver_parameters = get_SolverParameters_tbp_EARTH_lt_leo_to_geo(
		N, DDP_type, verbosity);

	// Solver parameters
	unsigned int Nx = solver_parameters.Nx();
	unsigned int Nu = solver_parameters.Nu();

	// Initial conditions [Equinoctial elements, MASSU, TU]
	ToF = ToF / SEC2DAYS / tu; // [TU]
	double dt = ToF / N; // [TU]
	double altitude = 400;


	double r_p = R_EARTH + altitude;

	r_p = 6846.8;

	vectordb x_departure{ // Kep coordinates
		(lu + r_p) / 2.0 / lu, (lu - r_p) / (lu + r_p),
		0 * DEG_2_RAD, 135 * DEG_2_RAD,
		0 * DEG_2_RAD, 0 * DEG_2_RAD,
		spacecraft_parameters.initial_mass(), 2 * PI};
	vectordb x_arrival{ // Kep coordinates
		lu / lu, 0,
		0 * DEG_2_RAD, 0 * DEG_2_RAD,
		0 * DEG_2_RAD, 0 * DEG_2_RAD,
		spacecraft_parameters.dry_mass(), 2*PI};

	x_departure[2] = 7 * DEG_2_RAD;
	x_departure[3] = 0 * DEG_2_RAD;

	x_departure = kep_2_equi(x_departure); // Equinoctial coordinates
	x_arrival = kep_2_equi(x_arrival);
	vectordb x0 = x_departure; x0[Nx - 1] = dt; // Time step
	vectordb x_goal = x_arrival; x_goal[Nx - 1] = ToF; // ToF

	// First guess command
	vectordb u_init(Nu, 1e-6 / thrustu); // [VU]
	vector<vectordb> list_u_init(N, u_init);
	
	string file_name = "./data/datasets/tbp_EARTH_lt_leo_to_geo";
	/*
	pair<vector<vectordb>, vector<vectordb>> traj = load_dataset(
		file_name, ToF,
		dynamics, spacecraft_parameters, constants, solver_parameters);

	vector<vectordb> list_x_ = traj.first;
	vector<vectordb> list_u_ = traj.second;
	for (size_t i = 0; i < list_u_init.size(); i++) {
		// Unpack
		vectordb x_i = list_x_[i];
		vectordb u_i = list_u_[i];

		// Convert to cartesian
		x_i = cart_2_kep(x_i, mu); // Convert to Keplerian
		x_i = kep_2_equi(x_i); // Convert to Equinoctial
		u_i = cart_2_RTN(u_i, list_x_[i]);

		// Assign
		list_u_[i] = u_i;
		list_x_[i] = x_i;
	}
	list_u_init = list_u_;*/

	// Set double precision
	typedef std::numeric_limits<double> dbl;
	cout.precision(5);

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

	// Unpack
	vector<vectordb> list_x = solver.list_x();
	vector<vectordb> list_u = solver.list_u();
	
	// Print datasets
	if (save_results) {
		// Convert to cart
		mu = MU_EARTH / constants.mu();
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

		
		string system_name = "TBP EARTH EQUINOCTIAL LT";
		print_transfer_dataset(
			file_name, system_name,
			list_x, list_u,
			x_departure, x_arrival, ToF,
			dynamics, spacecraft_parameters, constants, solver_parameters);
	}

}
