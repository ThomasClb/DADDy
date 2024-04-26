/**
	double_integrator.cpp

	Purpose: Double_integrator execution script.

	@author Thomas Caleb

	@version 1.0 26/04/2023
*/

#include "test_cases.h"

using namespace DACE;
using namespace std::chrono;
using namespace std;

SolverParameters get_SolverParameters_double_integrator(
	unsigned int const& N, unsigned int const& DDP_type,
	unsigned int verbosity) {
	// Solver parameters
	unsigned int Nx = 7;
	unsigned int Nu = 3;
	unsigned int Neq = 0;
	unsigned int Nineq = 0;
	unsigned int Nteq = 0;
	unsigned int Ntineq = 0;
	bool with_J2 = false;
	double cost_to_go_gain = 1e-2;
	double terminal_cost_gain = 1e4;
	double mass_leak = 1e-8;
	double homotopy_coefficient = 0.0;
	double huber_loss_coefficient = 5e-3;
	vectordb homotopy_sequence{0};
	vectordb huber_loss_coefficient_sequence{1e-2};
	double DDP_tol = 1e-6;
	double AUL_tol = 1e-8; 
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

void double_integrator(int argc, char** argv) {
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
	Dynamics dynamics = get_double_integrator_dynamics();
	
	// Normalisation constants
	Constants constants(dynamics.constants());

	// Spacecraft parameters
	SpacecraftParameters spacecraft_parameters(spacecraft_parameters_file);

	// Init solver parameters
	SolverParameters solver_parameters = get_SolverParameters_double_integrator(
		N, DDP_type, verbosity);

	// Solver parameters
	unsigned int Nx = solver_parameters.Nx();
	unsigned int Nu = solver_parameters.Nu();
	N = solver_parameters.N();

	// Initial conditions [3*LU, 3*VU, MASSU, TU]
	vectordb x_departure{ 1.0, 1.0, 1.0, 1, 1, 1, 0};
	vectordb x_arrival{1.0, -1.0, 0.0, 0, 0, 0, 0};
	vectordb x0 = x_departure;
	vectordb x_goal = x_arrival;

	// First guess command
	vectordb u_init(Nu, 2e-5); // [VU]
	vector<vectordb> list_u_init(N, u_init);

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
		cout << DDP_type << ", ";
		cout << "1" << ", ";
		cout << N << ", ";

		// Data

		// Results

		cout << solver.cost() << ", ";
		cout << solver.cost() << ", ";
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
		string file_name = "./data/datasets/double_integrator";
		string system_name = "DOUBLE INTEGRATOR";
		print_transfer_dataset(
			file_name, system_name,
			solver.list_x(), solver.list_u(),
			x_departure, x_arrival, ToF,
			dynamics, spacecraft_parameters, constants, solver_parameters);
	}
}
