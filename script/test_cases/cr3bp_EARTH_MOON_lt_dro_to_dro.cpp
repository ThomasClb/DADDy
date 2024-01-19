/**
	cr3bp_EARTH_MOON_lt_dro_to_dro.cpp

	Purpose: Low-thrust DRO to DRO transfer execution script.
	In the Earth-Moon system.

	@author Thomas Caleb

	@version 1.0 07/12/2023
*/

#include "test_cases.h"

using namespace DACE;
using namespace std::chrono;
using namespace std;

SolverParameters get_SolverParameters_cr3bp_EARTH_MOON_lt_dro_to_dro() {
	// Solver parameters
	unsigned int Nx = (SIZE_VECTOR + 1) + 1;
	unsigned int Nu = SIZE_VECTOR / 2;
	unsigned int Neq = 0;
	unsigned int Nineq = 3;
	unsigned int Nteq = 6;
	unsigned int Ntineq = 0;
	unsigned int N = 300;
	double cost_to_go_gain = 1e-3;
	double terminal_cost_gain = 1e7;
	double homotopy_coefficient = 0.0;
	double huber_loss_coefficient = 1e-4;
	unsigned int DDP_type = 3 + 0*4;
	double DDP_tol = 1e-4;
	double AUL_tol = 5e-6;
	double PN_tol = 1e-12;
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
		verbosity);
}

vector<vectordb> make_first_guess(
	double const& thrust,
	vectordb const& x_0,
	Dynamics const& dynamics,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {

	// Unpack
	int Nu = solver_parameters.Nu();
	int Nx = solver_parameters.Nx();
	int N = solver_parameters.N();
	double period = x_0[Nx - 1];

	// Init
	vector<vectordb> output;
	vectordb x_i = x_0;
	double dt = period / (N - 1.0);
	x_i[Nx - 1] = dt;
	vectordb null_control(Nu, 0.0);

	// Loop
	for (size_t i = 0; i < N; i++) {
		output.push_back(x_i.extract(3, 5).normalize()* thrust);
		x_i = dynamics.dynamic_db()(
			x_i, null_control,
			spacecraft_parameters, constants, solver_parameters);
	}
	return output;
}

void cr3bp_EARTH_MOON_lt_dro_to_dro(bool const& plot_graphs) {

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

	// Spacecraft parameters (GTOC 12)
	double m_0 = 1000 / massu; // [MASSU]
	double dry_mass = 500 / massu; // [MASSU]
	double T = 0.5 / thrustu; // [N]
	double Isp = 2000.0 / tu; // [s]
	SpacecraftParameters spacecraft_parameters(
		dynamics.constants(),
		m_0, dry_mass, T, Isp);

	// Init solver parameters
	SolverParameters solver_parameters = get_SolverParameters_cr3bp_EARTH_MOON_lt_dro_to_dro();

	// Solver parameters
	unsigned int Nx = solver_parameters.Nx();
	unsigned int Nu = solver_parameters.Nu();
	unsigned int N = solver_parameters.N();

	// Init DACE
	DA::init(2, Nx + Nu);
	DA::setEps(1e-90);

	// Initial conditions [3*LU, 3*VU, MASSU, TU]
	// From [Boone MacMahon 2024]
	int nb_revs = 2;
	double ToF = nb_revs *17.5 / SEC2DAYS / tu; // [TU]
	double dt = ToF / N; // [TU]
	vectordb x_departure{ 
		1.171359, 0, 0.0,
		0, -0.489458, 0.0,
		m_0, 13.4 /SEC2DAYS / tu };
	vectordb x_arrival{
		1.301844, 0, 0.0,
		0, -0.642177, 0.0,
		dry_mass, 21.6 / SEC2DAYS / tu };
	vectordb x0 = x_departure; x0[Nx - 1] = dt; // Time step
	vectordb x_goal = x_arrival; x_goal[Nx - 1] = ToF; // ToF

	// First guess command
	vectordb u_init(Nu, 1e-1 / thrustu); // [VU]
	vector<vectordb> list_u_init(N, u_init);
	list_u_init = make_first_guess(
		0.1 * T, x_departure,
		dynamics, spacecraft_parameters,
		constants, solver_parameters);

	// Output
	cout << "DEPARTURE : " << endl << x0.extract(0, Nx - 1 - 1) << endl;
	cout << "ARRIVAL : " << endl << x_goal.extract(0, Nx - 1 - 1) << endl;

	// AULSolver
	AULSolver solver(solver_parameters, spacecraft_parameters, dynamics);

	// Run DDP
	auto start = high_resolution_clock::now();
	solver.set_homotopy_coefficient(0.0);
	solver.solve(x0, list_u_init, x_goal);


	vectordb homotopy_sequence, huber_loss_coefficient_sequence;
	if (nb_revs == 1) {
		homotopy_sequence = vectordb{0.95, 0.95, 1};
		huber_loss_coefficient_sequence = vectordb{1e-2, 1e-3, 1e-3};
	}
	if (nb_revs == 2) {
		homotopy_sequence = vectordb{ 0.95, 0.99, 1-1e-3 };
		huber_loss_coefficient_sequence = vectordb{ 1e-2, 5e-3, 1e-3 };
	}
	if (nb_revs == 3) {
		homotopy_sequence = vectordb{ 0.95, 0.99, 1 - 1e-3 };
		huber_loss_coefficient_sequence = vectordb{ 1e-2, 5e-3, 1e-3 };
	}
	/**/
	// Set DACE at order 1 (No Hessian needed)
	DA::setTO(1);

	// PN test
	auto start_inter = high_resolution_clock::now();
	PNSolver pn_solver(solver);
	//pn_solver.solve(x_goal);
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
	string file_name = "./data/datasets/cr3bp_EARTH_MOON_lt_dro_to_dro.dat";
	string system_name = "CR3BP LT";
	print_transfer_dataset(
		file_name, system_name,
		list_x, list_u,
		x_departure, x_arrival,
		dynamics, spacecraft_parameters, constants, solver_parameters);
}
