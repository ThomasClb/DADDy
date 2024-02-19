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
	unsigned int Nineq = 3;
	unsigned int Nteq = 5;
	unsigned int Ntineq = 0;
	double homotopy_coefficient = 0.0;
	double cost_to_go_gain = 1e-3;
	double terminal_cost_gain = 1e9;
	double huber_loss_coefficient = 5e-3;
	double DDP_tol = 1e-4;
	double AUL_tol = 1e-6; 
	double PN_tol = 1e-12;
	double PN_active_constraint_tol = 1e-13;
	unsigned int max_iter = 10000;
	unsigned int DDP_max_iter = 100;
	unsigned int AUL_max_iter = max_iter / DDP_max_iter;
	unsigned int PN_max_iter = 50;
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
		Nteq, Ntineq,
		cost_to_go_gain, terminal_cost_gain,
		homotopy_coefficient, huber_loss_coefficient,
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
	double m_0 = 1000.0 / massu; // [MASSU]
	double dry_mass = m_0/2; // [MASSU]
	double T = 1 / thrustu; // [THRUSTU]
	double Isp = 2000.0 / tu; // [TU]
	SpacecraftParameters spacecraft_parameters(
		dynamics.constants(),
		m_0, dry_mass, T, Isp);
	// spacecraft_parameters.save(spacecraft_parameters_file);

	// Init solver parameters
	SolverParameters solver_parameters = get_SolverParameters_tbp_EARTH_lt_leo_to_geo(
		N, DDP_type, verbosity);

	// Solver parameters
	unsigned int Nx = solver_parameters.Nx();
	unsigned int Nu = solver_parameters.Nu();

	// Init DACE
	DA::init(2, Nx + Nu);
	DA::setEps(1e-90);

	// Initial conditions [Equinoctial elements, MASSU, TU]
	ToF = ToF / SEC2DAYS / tu; // [TU]
	double dt = ToF / N; // [TU]
	double altitude = 5000;
	double r_p = 6371 + altitude;
	vectordb x_departure{ // Kep coordinates
		lu / lu, 0,
		0 * DEG_2_RAD, 150 * DEG_2_RAD,
		0 * DEG_2_RAD, 0 * DEG_2_RAD,
		m_0, 2 * PI * sqrt(pow(lu, 3) / mu) / tu };
	vectordb x_arrival{ // Kep coordinates
		(lu + r_p) / 2.0 / lu, (lu - r_p) / (lu + r_p),
		0 * DEG_2_RAD, 150 * DEG_2_RAD,
		0 * DEG_2_RAD, 0 * DEG_2_RAD,
		dry_mass, 2 * PI * sqrt(pow((lu + r_p) / 2.0, 3) / mu) / tu   };
	x_departure = kep_2_equi(x_departure); // Equinoctial coordinates
	x_arrival = kep_2_equi(x_arrival);
	vectordb x0 = x_departure; x0[Nx - 1] = dt; // Time step
	vectordb x_goal = x_arrival; x_goal[Nx - 1] = ToF; // ToF

	// First guess command
	vectordb u_init(Nu, 1e-6 / thrustu); // [VU]
	vector<vectordb> list_u_init(N, u_init);

	// AULSolver
	AULSolver solver(solver_parameters, spacecraft_parameters, dynamics);

	// Run DDP
	auto start = high_resolution_clock::now();

	solver.set_homotopy_coefficient(0.0);
	solver.solve(x0, list_u_init, x_goal);
	vectordb homotopy_sequence{1, 1.0};
	vectordb huber_loss_coefficient_sequence{1e-2, 1e-3};
	if (fuel_optimal) {
		for (size_t i = 0; i < homotopy_sequence.size(); i++) {
			solver.set_huber_loss_coefficient(huber_loss_coefficient_sequence[i]);
			solver.set_homotopy_coefficient(homotopy_sequence[i]);
			solver.solve(x0, solver.list_u(), x_goal);
		}
	}
	// Set DACE at order 1 (No Hessian needed)
	DA::setTO(1);

	// PN test
	auto start_inter = high_resolution_clock::now();
	PNSolver pn_solver(solver);
	pn_solver.solve(x_goal);
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	auto duration_AUL = duration_cast<microseconds>(start_inter - start);
	auto duration_PN = duration_cast<microseconds>(stop - start_inter);

	// Get data
	vector<vectordb> list_x(pn_solver.list_x()), list_u(pn_solver.list_u());
	double cost(pn_solver.cost());

	// Compute errors
	vectordb loss = (list_x[N] - x_goal).extract(0, 5);

	// Get final mass
	double final_mass = list_x[N][6]; // [kg]

	// Output
	cout << endl;
	cout << "Optimised" << endl;
	cout << "	Total runtime : " + to_string(static_cast<double>(duration.count()) / 1e6) + "s" << endl;
	cout << "	AUL solver runtime : " + to_string(static_cast<double>(duration_AUL.count()) / 1e6) + "s" << endl;
	cout << "	PN solver runtime : " + to_string(static_cast<double>(duration_PN.count()) / 1e6) + "s" << endl;
	cout << "	FINAL MASS [kg] : " << massu * final_mass << endl;
	cout << "	FINAL ERROR [-] : " << real_constraints(x_goal, pn_solver) << endl;

	// Print datasets

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

	string file_name = "./data/datasets/tbp_EARTH_lt_leo_to_geo.dat";
	string system_name = "TBP EARTH EQUINOCTIAL LT";
	print_transfer_dataset(
		file_name, system_name,
		list_x, list_u,
		x_departure, x_arrival,
		dynamics, spacecraft_parameters, constants, solver_parameters);
}
