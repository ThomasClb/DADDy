/**
	tbp_SUN_low_thrust_earth_to_mars.cpp

	Purpose: Low-thrust Earth-Mars transfer execution script.

	@author Thomas Caleb

	@version 1.0 07/12/2023
*/

#include "test_cases.h"

using namespace DACE;
using namespace std::chrono;
using namespace std;

SolverParameters get_SolverParameters_tbp_SUN() {
	// Solver parameters
	unsigned int Nx = (SIZE_VECTOR + 1) + 1;
	unsigned int Nu = SIZE_VECTOR / 2;
	unsigned int Neq = 0;
	unsigned int Nineq = 3;
	unsigned int Nteq = 6;
	unsigned int Ntineq = 0;
	unsigned int N = 40;
	double homotopy_coefficient = 0.0;
	double cost_to_go_gain = 1e-2;
	double terminal_cost_gain = 1e4;
	double huber_loss_coefficient = 1e-3;
	unsigned int DDP_type = 3 + 0*4;
	double DDP_tol = 1e-4;
	double AUL_tol = 1e-4;
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

void tbp_SUN_low_thrust_earth_to_mars(bool const& plot_graphs) {

	// Set double precision
	typedef std::numeric_limits<double> dbl;
	cout.precision(5);

	// Set dynamics
	Dynamics dynamics = get_tbp_SUN_low_thrust_dynamics();
	
	// Normalisation constants
	Constants constants(dynamics.constants());
	double lu = constants.lu();
	double massu = constants.massu();
	double tu = constants.tu();
	double thrustu = constants.thrustu();
	double vu = constants.vu();

	// Spacecraft parameters
	double m_0 = 1000.0 / massu; // [MASSU]
	double dry_mass = 500.0 / massu; // [MASSU]
	double T = 0.5 / thrustu; // [THRUSTU]
	double Isp = 2000.0 / tu; // [TU]
	SpacecraftParameters spacecraft_parameters(
		constants,
		m_0, dry_mass, T, Isp);

	// Init solver parameters
	SolverParameters solver_parameters = get_SolverParameters_tbp_SUN();

	// Solver parameters
	unsigned int Nx = solver_parameters.Nx();
	unsigned int Nu = solver_parameters.Nu();
	unsigned int N = solver_parameters.N();

	// Init DACE
	DA::init(2, Nx + Nu);
	DA::setEps(1e-90);

	// Initial conditions [3*LU, 3*VU, MASSU, TU]
	double ToF = 348.79 / SEC2DAYS / tu; // [TU]
	double dt = ToF / N; // [TU]
	vectordb x0{
		-140699693 / lu, -51614428 / lu, 980 / lu,
		9.774596 / vu, -28.07828 / vu, 4.337725e-4 / vu,
		m_0, dt };
	vectordb x_goal{
		-172682023 / lu, 176959469 / lu, 7948912 / lu,
		-16.427384 / vu, -14.860506 / vu, 9.21486e-2 / vu,
		dry_mass, ToF };

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
	vectordb homotopy_sequence{0.95, 1.0 - 1e-2};
	for (size_t i = 0; i < homotopy_sequence.size(); i++) {
		solver.set_homotopy_coefficient(homotopy_sequence[i]);
		solver.solve(x0, solver.list_u(), x_goal);
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

	if (plot_graphs) {
		/**/

		// Plot
		vectordb list_0(N + 1), list_1(N + 1), list_2(N + 1), list_m(N + 1), list_N(N + 1);
		for (size_t i = 0; i < N + 1; i++) {
			list_N[i] = i;
			list_0[i] = list_x[i][0];
			list_1[i] = list_x[i][1];
			list_2[i] = list_x[i][2];
			list_m[i] = list_x[i][6];
		}

		// Transfer x,y
		figure();
		plot(list_0, list_1);
		xlabel("X [AU]"); ylabel("Y [AU]");

		figure();
		plot(list_N * dt * SEC2DAYS * tu, massu * list_m);
		xlabel("Time [days]"); ylabel("Remaining mass [kg]");

		// Thrust
		list_0 = vectordb(N); list_1 = vectordb(N); list_2 = vectordb(N); list_N = vectordb(N); vectordb list_T(N);
		for (size_t i = 0; i < N; i++) {
			list_N[i] = i;
			list_0[i] = list_u[i][0];
			list_1[i] = list_u[i][1];
			list_2[i] = list_u[i][2];
			list_T[i] = list_u[i].vnorm();
		}
		figure();
		/*
		plot(list_N* dt* SEC2DAYS* tu, list_0);
		plot(list_N * dt * SEC2DAYS * tu, list_1);
		plot(list_N* dt* SEC2DAYS* tu, list_2);
		*/
		plot(list_N * dt * SEC2DAYS * tu, list_T* thrustu);
		plot(list_N * dt * SEC2DAYS * tu, list_T * 0 + T * thrustu);
		ylabel("T [N]"); xlabel("Time [days]");

		show();

	}

}
