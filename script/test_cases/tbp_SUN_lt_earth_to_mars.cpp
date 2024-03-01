/**
	tbp_SUN_lt_earth_to_mars.cpp

	Purpose: Low-thrust Earth-Mars transfer execution script.

	@author Thomas Caleb

	@version 1.0 07/12/2023
*/

#include "test_cases.h"

using namespace DACE;
using namespace std::chrono;
using namespace std;

SolverParameters get_SolverParameters_tbp_SUN_lt_earth_to_mars(
	unsigned int const& N, unsigned int const& DDP_type,
	unsigned int verbosity) {
	// Solver parameters
	unsigned int Nx = (SIZE_VECTOR + 1) + 1;
	unsigned int Nu = SIZE_VECTOR / 2;
	unsigned int Neq = 0;
	unsigned int Nineq = 3;
	unsigned int Nteq = 6;
	unsigned int Ntineq = 0;
	double homotopy_coefficient = 0.0;
	double cost_to_go_gain = 1e-2;
	double terminal_cost_gain = 1e4;
	double huber_loss_coefficient = 5e-3;
	double DDP_tol = 1e-4;
	double AUL_tol = 1e-4; 
	double PN_tol = 1e-10;
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

void tbp_SUN_lt_earth_to_mars(int argc, char** argv) {
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
	Dynamics dynamics = get_tbp_SUN_lt_dynamics();
	
	// Normalisation constants
	Constants constants(dynamics.constants());
	double lu = constants.lu();
	double massu = constants.massu();
	double tu = constants.tu();
	double thrustu = constants.thrustu();
	double vu = constants.vu();

	// Spacecraft parameters
	SpacecraftParameters spacecraft_parameters(spacecraft_parameters_file);

	// Init solver parameters
	SolverParameters solver_parameters = get_SolverParameters_tbp_SUN_lt_earth_to_mars(
		N, DDP_type, verbosity);

	// Solver parameters
	unsigned int Nx = solver_parameters.Nx();
	unsigned int Nu = solver_parameters.Nu();
	N = solver_parameters.N();

	// Init DACE
	DA::init(2, Nx + Nu);
	DA::setEps(1e-90);

	// Initial conditions [3*LU, 3*VU, MASSU, TU]
	ToF /=  (SEC2DAYS * tu); // [TU]
	double dt = ToF / N; // [TU]
	vectordb x_departure{
		-140699693 / lu, -51614428 / lu, 980 / lu,
		9.774596 / vu, -28.07828 / vu, 4.337725e-4 / vu,
		spacecraft_parameters.initial_mass(), 365.25/SEC2DAYS/tu };
	vectordb x_arrival{
		-172682023 / lu, 176959469 / lu, 7948912 / lu,
		-16.427384 / vu, -14.860506 / vu, 9.21486e-2 / vu,
		spacecraft_parameters.dry_mass(), 700 / SEC2DAYS / tu };
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
	vectordb homotopy_sequence{0.5, 0.9, 0.999};
	vectordb huber_loss_coefficient_sequence{1e-2, 2e-3, 1e-3};
	if (fuel_optimal) {
		for (size_t i = 0; i < homotopy_sequence.size(); i++) {
			solver.set_huber_loss_coefficient(huber_loss_coefficient_sequence[i]);
			solver.set_homotopy_coefficient(homotopy_sequence[i]);
			solver.solve(x0, solver.list_u(), x_goal);
		}
	}

	// PN test
	auto start_inter = high_resolution_clock::now();
	PNSolver pn_solver(solver);
	if (pn_solving)
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
	if (verbosity <= 1) {
		cout << endl;
		cout << "Optimised" << endl;
		cout << "	Total runtime : " + to_string(static_cast<double>(duration.count()) / 1e6) + "s" << endl;
		cout << "	AUL solver runtime : " + to_string(static_cast<double>(duration_AUL.count()) / 1e6) + "s" << endl;
		cout << "	PN solver runtime : " + to_string(static_cast<double>(duration_PN.count()) / 1e6) + "s" << endl;
		cout << "	FINAL MASS [kg] : " << massu * final_mass << endl;
		cout << "	FINAL ERROR [-] : " << real_constraints(x_goal, pn_solver) << endl;
	}

	// Print datasets
	if (save_results) {
		string file_name = "./data/datasets/tbp_SUN_lt_earth_to_mars";
		string system_name = "TBP SUN CARTESIAN LT";
		print_transfer_dataset(
			file_name, system_name,
			list_x, list_u,
			x_departure, x_arrival,
			dynamics, spacecraft_parameters, constants, solver_parameters);
	}
}
