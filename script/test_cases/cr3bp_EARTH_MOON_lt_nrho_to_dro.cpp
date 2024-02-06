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

SolverParameters get_SolverParameters_cr3bp_EARTH_MOON_lt_nrho_to_dro(unsigned int const& N) {
	// Solver parameters
	unsigned int Nx = (SIZE_VECTOR + 1) + 1;
	unsigned int Nu = SIZE_VECTOR / 2;
	unsigned int Neq = 0;
	unsigned int Nineq = 3;
	unsigned int Nteq = 6;
	unsigned int Ntineq = 0;
	double cost_to_go_gain = 5e-4;
	double terminal_cost_gain = 1e9;
	double homotopy_coefficient = 0.0;
	double huber_loss_coefficient = 5e-4;
	unsigned int DDP_type = 3 + 0*4;
	double DDP_tol = 1e-4;
	double AUL_tol = 5e-6;
	double PN_tol = 1e-10;
	double PN_active_constraint_tol = 1e-13;
	unsigned int max_iter = 10000;
	unsigned int DDP_max_iter = 100;
	unsigned int AUL_max_iter = max_iter / DDP_max_iter;
	unsigned int PN_max_iter = 100;
	vectordb lambda_parameters{0.0, 1e8};
	vectordb mu_parameters{1, 1e8, 10};
	vectordb line_search_parameters{1e-8, 10.0, 0.5, 20};
	bool backward_sweep_regulation = true;
	vectordb backward_sweep_regulation_parameters{0, 1e-8, 1e8, 1.6};
	double PN_regularisation(1e-8);
	double PN_cv_rate_threshold(1.1);
	double PN_alpha(1.0); double PN_gamma(0.5);
	unsigned int verbosity = 0;
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

void cr3bp_EARTH_MOON_lt_nrho_to_dro(int argc, char** argv) {
	// Input check
	if (argc < 5) {
		cout << "Wrong number of arguments." << endl;
		cout << "Requested number : 4" << endl;
		cout << "0 - Test case number." << endl;
		cout << "1 - SpacecraftParameter adress." << endl;
		cout << "2 - Number of nodes [-]." << endl;
		cout << "3 - Time of flight [days]." << endl;
		return;
	}

	// Unpack inputs
	string spacecraft_parameters_file = argv[2];
	int N = atoi(argv[3]);
	double ToF = atof(argv[4]);
	bool fuel_optimal = true;

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
	SolverParameters solver_parameters = get_SolverParameters_cr3bp_EARTH_MOON_lt_nrho_to_dro(N);

	// Solver parameters
	unsigned int Nx = solver_parameters.Nx();
	unsigned int Nu = solver_parameters.Nu();

	// Init DACE
	DA::init(2, Nx + Nu);
	DA::setEps(1e-90);

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
	vectordb u_init(Nu, 1e-6 / thrustu); // [VU]
	vector<vectordb> list_u_init(N, u_init);

	// Output
	cout << "DEPARTURE : " << endl << x0.extract(0, Nx - 1 - 1) << endl;
	cout << "ARRIVAL : " << endl << x_goal.extract(0, Nx - 1 - 1) << endl;

	// AULSolver
	AULSolver solver(solver_parameters, spacecraft_parameters, dynamics);

	// Run DDP
	auto start = high_resolution_clock::now();
	solver.set_homotopy_coefficient(0.0);
	solver.solve(x0, list_u_init, x_goal);
	/*	*/
	vectordb huber_loss_coefficient_sequence{1e-2, 1e-2, 5e-3, 5e-3, 1e-4};
	vectordb homotopy_sequence{0.3, 0.66, 0.9, 0.99, 0.99};
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
	cout << "	Total runtime : " + to_string(static_cast<int>(duration.count()) / 1e6) + "s" << endl;
	cout << "	AUL solver runtime : " + to_string(static_cast<int>(duration_AUL.count()) / 1e6) + "s" << endl;
	cout << "	PN solver runtime : " + to_string(static_cast<int>(duration_PN.count()) / 1e6) + "s" << endl;
	cout << "	FINAL MASS [kg] : " << massu * final_mass << endl;
	cout << "	FINAL ERROR [-] : " << real_constraints(x_goal, pn_solver) << endl;

	// Print datasets
	string file_name = "./data/datasets/cr3bp_EARTH_MOON_lt_nrho_to_dro.dat";
	string system_name = "CR3BP EARTH-MOON CARTESIAN LT";
	print_transfer_dataset(
		file_name, system_name,
		list_x, list_u,
		x_departure, x_arrival,
		dynamics, spacecraft_parameters, constants, solver_parameters);
}
