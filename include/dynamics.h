/**
	dynamics.cpp

	Purpose: Implementation of the Dynamics class
	methods.

	@author Thomas Caleb

	@version 1.0 16/11/2023
*/

#ifndef DEF_DYNAMIC
#define DEF_DYNAMIC

#pragma once

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <functional>

#include <dace/dace_s.h>

#include "settings.h"
#include "constants.h"
#include "parameters.h"
#include "integration.h"


// Define the 7 types of functions

// vDA-vDA-SpacecraftParameters-SolverParameters->vDA

// Returns the next state given
// the current state, the control and the parameters.
using dynFunction = std::function<DACE::vectorDA(
	DACE::vectorDA const&,
	DACE::vectorDA const&, SpacecraftParameters const&,
	Constants const&, SolverParameters const&)>;

// Returns the path equality contraints given
// the current state, the control and the parameters.
using eqFunction = dynFunction;

// Returns the path inequality contraints given
// the current state, the control and the parameters.
using ineqFunction = dynFunction;

// vDA-vDA-SpacecraftParameters-SolverParameters->DA

// Returns the cost-to-go given
// the current state, the control and the parameters.
using ctgFunction = std::function<DACE::DA(
	DACE::vectorDA const&, DACE::vectorDA const&, 
	SpacecraftParameters const&,
	Constants const&,
	SolverParameters const&)>;

// vDA-vdb-SpacecraftParameters-SolverParameters->DA

// Returns the terminal cost given
// the current state, the target state and the parameters.
using tcFunction = std::function<DACE::DA(
	DACE::vectorDA const&, DACE::vectordb const&,
	SpacecraftParameters const&,
	Constants const&,
	SolverParameters const&)>;

// vDA-vdb-SpacecraftParameters-SolverParameters->vDA

// Returns the terminal equality constraints given
// the current state, the target state and the parameters.
using teqFunction = std::function<DACE::vectorDA(
	DACE::vectorDA const&, DACE::vectordb const&,
	SpacecraftParameters const&,
	Constants const&,
	SolverParameters const&)>;

// Returns the terminal inequality constraints given
// the current state, the target state and the parameters.
using tineqFunction = teqFunction;


// Double versions

// vdb-vdb-SpacecraftParameters-SolverParameters->vdb

// Returns the next state given
// the current state, the control and the parameters.
using dynFunction_db = std::function<DACE::vectordb(
	DACE::vectordb const&,
	DACE::vectordb const&, SpacecraftParameters const&,
	Constants const&,
	SolverParameters const&)>;

// Returns the path equality contraints given
// the current state, the control and the parameters.
using eqFunction_db = dynFunction_db;

// Returns the path inequality contraints given
// the current state, the control and the parameters.
using ineqFunction_db = dynFunction_db;

// vdb-vdb-SpacecraftParameters-SolverParameters->db

// Returns the cost-to-go given
// the current state, the control and the parameters.
using ctgFunction_db = std::function<double(
	DACE::vectordb const&, DACE::vectordb const&,
	SpacecraftParameters const&,
	Constants const&,
	SolverParameters const&)>;

// vdb-vdb-SpacecraftParameters-SolverParameters->db

// Returns the terminal cost given
// the current state, the target state and the parameters.
using tcFunction_db = std::function<double(
	DACE::vectordb const&, DACE::vectordb const&,
	SpacecraftParameters const&,
	Constants const&,
	SolverParameters const&)>;

// vdb-vdb-SpacecraftParameters-SolverParameters->vdb

// Returns the terminal equality constraints given
// the current state, the target state and the parameters.
using teqFunction_db = std::function<DACE::vectordb(
	DACE::vectordb const&, DACE::vectordb const&,
	SpacecraftParameters const&,
	Constants const&,
	SolverParameters const&)>;

// Returns the terminal inequality constraints given
// the current state, the target state and the parameters.
using tineqFunction_db = teqFunction_db;


/*

	DYNAMICS

*/

class Dynamics {

	// Attributes
protected:
	Constants constants_;

	// DA functions
	dynFunction dynamic_;
	ctgFunction cost_to_go_;
	eqFunction equality_constraints_;
	ineqFunction inequality_constraints_;
	tcFunction terminal_cost_;
	teqFunction terminal_equality_constraints_;
	tineqFunction terminal_inequality_constraints_;

	// Double versions
	dynFunction_db dynamic_db_;
	ctgFunction_db cost_to_go_db_;
	eqFunction_db equality_constraints_db_;
	ineqFunction_db inequality_constraints_db_;
	tcFunction_db terminal_cost_db_;
	teqFunction_db terminal_equality_constraints_db_;
	tineqFunction_db terminal_inequality_constraints_db_;

// Methods
public:
	// Constructors

	// Default constructors
	Dynamics();

	// Constructor
	Dynamics(
		Constants const& constants,
		dynFunction const& dynamic,
		ctgFunction const& cost_to_go,
		eqFunction const& equality_constraints,
		ineqFunction const& inequality_constraints,
		tcFunction const& terminal_cost,
		teqFunction const& terminal_equality_constraints,
		tineqFunction const& terminal_inequality_constraints,
		dynFunction_db const& dynamic_db,
		ctgFunction_db const& cost_to_go_db,
		eqFunction_db const& equality_constraints_db,
		ineqFunction_db const& inequality_constraints_db,
		tcFunction_db const& terminal_cost_db,
		teqFunction_db const& terminal_equality_constraints_db,
		tineqFunction_db const& terminal_inequality_constraints_db);

	// Copy constructor
	Dynamics(Dynamics const& dynamics);

	// Destructors
	~Dynamics();

	// Getters
	const Constants constants() const;

	const dynFunction dynamic() const;
	const ctgFunction cost_to_go() const;
	const eqFunction equality_constraints() const;
	const ineqFunction inequality_constraints() const;
	const tcFunction terminal_cost() const;
	const teqFunction terminal_equality_constraints() const;
	const tineqFunction terminal_inequality_constraints() const;

	const dynFunction_db dynamic_db() const;
	const ctgFunction_db cost_to_go_db() const;
	const eqFunction_db equality_constraints_db() const;
	const ineqFunction_db inequality_constraints_db() const;
	const tcFunction_db terminal_cost_db() const;
	const teqFunction_db terminal_equality_constraints_db() const;
	const tineqFunction_db terminal_inequality_constraints_db() const;
};


/*

	FUNCTIONS

*/

// Wraps a value between 0 and mod > 0
double wrap_mod(double const& value, double const& mod);
DACE::DA wrap_mod(DACE::DA const& value, double const& mod); // DA version

// Acceleration functions

// Computes the derivaties in an Sun-centered 2-body problem.
// With 3D continuous thrust.
// It takes at input [3*LU, 3*VU, MASSU, TU], [3*N], TU, SpacecraftParameters.
// It returns [3*VU, 3*VU/TU, MASSU/TU, 1].
template<typename T>
DACE::AlgebraicVector<T> acceleration_tbp_SUN_low_thrust(
	DACE::AlgebraicVector<T> const& x,
	DACE::AlgebraicVector<T> const& u, double const& t,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants) {
	// Unpack
	double v_e = spacecraft_parameters.ejection_velocity(); // [VU]

	// Init output
	DACE::AlgebraicVector<T> output(SIZE_VECTOR + 2);
	output[0] = x[3];
	output[1] = x[4];
	output[2] = x[5];
	output[7] = 0.0; // ToF

	// Acceleration kepler
	double mu = MU_SUN/ constants.mu();
	DACE::AlgebraicVector<T> r = x.extract(0, 2);
	T r_2 = r.dot(r);
	DACE::AlgebraicVector<T> acc_kep = -mu * r * pow(r_2, -1.5);

	// Thrust
	T inv_mass = 1 / x[6];
	DACE::AlgebraicVector<T> acc_thrust = inv_mass * u;

	// Assign
	DACE::AlgebraicVector<T> acc = acc_kep + acc_thrust;

	// dm [MASSU/LU]
	T thrust_norm = u.vnorm();
	T m_p = thrust_norm * pow(-1.0*v_e, -1);

	// Assign to output
	output[3] = acc[0];
	output[4] = acc[1];
	output[5] = acc[2];
	output[6] = m_p;
	
	return x[7] * output;
}

// Computes the derivaties in an Earth-centered 2-body problem.
// With 3D continuous thrust using Gauss planetary equations.
// It takes at input [3*LU, 3*VU, MASSU, TU], [3*N], TU, SpacecraftParameters.
// It returns [3*VU, 3*VU/TU, MASSU/TU, 1].
template<typename T>
DACE::AlgebraicVector<T> acceleration_tbp_EARTH_low_thrust(
	DACE::AlgebraicVector<T> const& x,
	DACE::AlgebraicVector<T> const& u, double const& t,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants) {
	// Unpack
	double v_e = spacecraft_parameters.ejection_velocity(); // [VU]
	double mu = MU_EARTH / constants.mu();
	T a = x[0]; // Equinoctial elements
	T P_1 = x[1];
	T P_2 = x[2];
	T Q_1 = x[3];
	T Q_2 = x[4];
	T L = x[5];
	T u_R = u[0]; 
	T u_T = u[1];
	T u_N = u[2];

	// Init output
	DACE::AlgebraicVector<T> output(SIZE_VECTOR + 2);
	output[7] = 0.0; // ToF

	// Compute useful data
	T cos_L = cos(L);
	T sin_L = sin(L);
	T p = a * (1 - P_1*P_1 + P_2 * P_2); // Ellipse parameter
	T r = p/(1+ P_1*cos_L + P_2*sin_L); // Distance
	T h = sqrt(mu * a * p); // Angular momentum
	T n = sqrt(mu * pow(a, -3.0)); // Mean motion
	T inv_cos_sqr_inc = (1 + Q_1 * Q_1 + Q_2 * Q_2);
	T r_over_h = (r / h);
	T p_over_r = (p / r);
	T d_Q = 0.5 * r_over_h * inv_cos_sqr_inc * u_N;
	T Q_1_cos_L = Q_1 * cos_L;
	T Q_2_sin_L = Q_2 * sin_L;

	// Acceleration Gauss
	T d_a = 2 * (a * a / h) * (
		(P_2 * sin_L - P_1 * cos_L) * u_R
		+ p_over_r * u_T);
	T d_P_1 = r_over_h * (
		-p_over_r * cos_L * u_R
		+ (P_1 + (1 + p_over_r)*sin_L) * u_T
		- P_2 * (Q_1_cos_L - Q_2_sin_L) * u_N);
	T d_P_2 = r_over_h * (
		-p_over_r * cos_L * u_R
		+ (P_2 + (1 + p_over_r) * sin_L) * u_T
		- P_1 * (Q_1_cos_L - Q_2_sin_L) * u_N);
	T d_Q_1 = d_Q * sin_L;
	T d_Q_2 = d_Q * cos_L;
	T d_L = n - r_over_h * (Q_1_cos_L - Q_2_sin_L) * u_N;

	// Thrust
	T inv_mass = 1 / x[6];
	DACE::AlgebraicVector<T> acc_thrust = inv_mass * u;

	// dm [MASSU/LU]
	T thrust_norm = u.vnorm();
	T m_p = thrust_norm * pow(-1.0 * v_e, -1);

	// Assign to output
	output[0] = d_a;
	output[1] = d_P_1;
	output[2] = d_P_2;
	output[3] = d_Q_1;
	output[4] = d_Q_2;
	output[5] = d_L;
	output[6] = m_p; // Mass

	return x[7] * output;
}

// Computes the derivaties in an Earth-Moon circular restricted 3-body problem.
// With 3D continuous thrust.
// It takes at input [3*LU, 3*VU, MASSU, TU], [3*N], TU, SpacecraftParameters.
// It returns [3*VU, 3*VU/TU, MASSU/TU, 1].
template<typename T>
DACE::AlgebraicVector<T>  acceleration_cr3bp_low_thrust(
	DACE::AlgebraicVector<T>  const& state_vector,
	DACE::AlgebraicVector<T> const& u, double const& t,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants) {
	// Unpack
	T x = state_vector[0];
	T y = state_vector[1];
	T z = state_vector[2];
	T x_p = state_vector[3];
	T y_p = state_vector[4];
	T z_p = state_vector[5];
	T mass = state_vector[6];
	T s = state_vector[7]; // Period
	double v_e = spacecraft_parameters.ejection_velocity(); // [VU]

	// velocity
	DACE::AlgebraicVector<T>  output(SIZE_VECTOR + 2);
	output[0] = x_p;
	output[1] = y_p;
	output[2] = z_p;
	output[7] = 0.0; // Period is constant

	// Inversion
	double mu = constants.mu();
	T d_y_z_2 = y* y + z* z;
	T inv_r_1_3 = (1 - mu) * pow((x + mu)* (x + mu) + d_y_z_2, -1.5);
	T inv_r_2_3 = mu * pow(pow(x - (1 - mu), 2) + d_y_z_2, -1.5);

	// Thrust
	T inv_mass = 1 / (mass);
	DACE::AlgebraicVector<T> acc_thrust = inv_mass * u;

	// Acceleration
	output[3] = 2 * y_p + x - (x + mu) * inv_r_1_3 - (x - (1 - mu)) * inv_r_2_3 + acc_thrust[0];
	output[4] = -2 * x_p + y * (1 - inv_r_1_3 - inv_r_2_3) + acc_thrust[1];
	output[5] = -z * (inv_r_1_3 + inv_r_2_3) + acc_thrust[2];

	// dm [MASSU/TU]
	T thrust_norm = u.vnorm(); // [-]
	output[6] = thrust_norm *pow(-1.0*v_e, -1);
	return s * output;
}

// Returns the next state given
// the current state, the control and the parameters.
// With acceleration acceleration_2b_SUN.
template<typename T>
DACE::AlgebraicVector<T> dynamic_tbp_SUN_low_thrust(
	DACE::AlgebraicVector<T> const& x, DACE::AlgebraicVector<T>const& u,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {
	return RK78(acceleration_tbp_SUN_low_thrust, x, u, 0, 1.0,
		spacecraft_parameters, constants);
}

// Returns the next state given
// the current state, the control and the parameters.
// With acceleration acceleration_2b_EARTH.
template<typename T>
DACE::AlgebraicVector<T> dynamic_tbp_EARTH_low_thrust(
	DACE::AlgebraicVector<T> const& x, DACE::AlgebraicVector<T>const& u,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {
	DACE::AlgebraicVector<T> output = RK78(
		acceleration_tbp_EARTH_low_thrust, x, u, 0, 1.0,
		spacecraft_parameters, constants);

	// std::cout << output << std::endl;
	return output;
}

// Returns the next state given
// the current state, the control and the parameters.
// With acceleration acceleration_cr3bp.
template<typename T>
DACE::AlgebraicVector<T> dynamic_cr3bp_low_thrust(
	DACE::AlgebraicVector<T> const& x, DACE::AlgebraicVector<T>const& u,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {
	return RK78(acceleration_cr3bp_low_thrust, x, u, 0, 1.0,
		spacecraft_parameters, constants);
}

// Returns the cost-to-go given
// the current state, the control and the parameters.
// Homotopy between energy-optimal and fuel-optimal low thrust.
template<typename T>
T cost_to_go(
	DACE::AlgebraicVector<T> const& x, DACE::AlgebraicVector<T> const& u,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {
	// Unpack parameters
	double tol = solver_parameters.DDP_tol();
	unsigned int N = solver_parameters.N();
	double gain = solver_parameters.cost_to_go_gain(); // [-]
	double homotopy_coefficient = solver_parameters.homotopy_coefficient(); // [-]
	double delta = solver_parameters.huber_loss_coefficient(); // [-]
	double ToF = solver_parameters.ToF(); // [-]
	double T_max = spacecraft_parameters.thrust(); // [THRUSTU]

	// Get time ratio
	double dt_avg = ToF / N; // [TU]
	T dt = x[SIZE_VECTOR + 1];

	// NRJ cost to go
	DACE::AlgebraicVector<T>  u_norm = u * (dt * (1 / (dt_avg * T_max))); // [THRUSTU]
	T NRJ_ctg = u_norm.dot(u_norm); // [THRUSTU^2]

	// Return cost
	if (homotopy_coefficient == 0.0)
		return (0.5 * gain) * NRJ_ctg;

	// Pseudo-Huber loss
	NRJ_ctg = NRJ_ctg / pow(delta, 2.0); // Normalise NRJ optimal term
	double mass_leak = tol * tol / pow(delta, 2.0);
	T fuel_ctg = sqrt(1.0 + NRJ_ctg + mass_leak) - 1.0;

	// Homotopy
	double c_1 = 0.5 * gain * delta * delta * (1.0 - homotopy_coefficient);
	double c_2 = delta *  gain * homotopy_coefficient;
	return(c_1 * NRJ_ctg + c_2 * fuel_ctg);
}

// Returns the path equality contraints given
// the current state, the control and the parameters.
// Null.
template<typename T>
DACE::AlgebraicVector<T> equality_constraints(
	DACE::AlgebraicVector<T> const& x, DACE::AlgebraicVector<T> const& u,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {
	return DACE::AlgebraicVector<T>(0);
}

// Returns the path inequality contraints given
// the current state, the control and the parameters.
// Boarder constraints, and thrust constraints.
// For low-thrust
// Nineq = 3
template<typename T>
DACE::AlgebraicVector<T> inequality_constraints_low_thrust(
	DACE::AlgebraicVector<T> const& x, DACE::AlgebraicVector<T> const& u,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {
	// Unpack parameters
	double T_max = spacecraft_parameters.thrust(); // [THRUSTU]
	double dry_mass = spacecraft_parameters.dry_mass(); // [MASSU]
	double initial_mass = spacecraft_parameters.initial_mass(); // [MASSU]

	// Init
	DACE::AlgebraicVector<T> output; output.reserve(2 + 1);

	// Thrust (1)
	T T_const = u.dot(u) / (T_max * T_max) - 1.0; // [-]
	output.push_back(T_const);

	// Mass (2)
	output.push_back(dry_mass - x[SIZE_VECTOR]); output.push_back(x[SIZE_VECTOR] - initial_mass); // Mass
	
	return output;
}

// Returns the terminal cost given
// the current state, the target state and the parameters.
// Final difference.
template<typename T>
T terminal_cost(
	DACE::AlgebraicVector<T> const& x, DACE::vectordb const& x_goal,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {
	// Unpack
	double gain = solver_parameters.terminal_cost_gain(); // [-]

	// Get vectors
	DACE::AlgebraicVector<T> x_(x.extract(0, SIZE_VECTOR - 1));
	DACE::vectordb x_goal_(x_goal.extract(0, SIZE_VECTOR - 1));

	// Compute error
	DACE::AlgebraicVector<T> loss = x_ - x_goal_;
	T output = (0.5 * gain) * loss.dot(loss);

	return output;
}

// Returns the terminal cost given
// the current state, the target state and the parameters.
// Final difference, without equinoctial anomaly.
template<typename T>
T terminal_cost_equinoctial(
	DACE::AlgebraicVector<T> const& x, DACE::vectordb const& x_goal,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {
	// Unpack
	double gain = solver_parameters.terminal_cost_gain(); // [-]

	// Get vectors
	DACE::AlgebraicVector<T> x_(x.extract(0, SIZE_VECTOR - 2));
	DACE::vectordb x_goal_(x_goal.extract(0, SIZE_VECTOR - 2));

	// Compute error
	DACE::AlgebraicVector<T> loss = x_ - x_goal_;
	T output = (0.5 * gain) * loss.dot(loss);

	return output;
}

// Returns the terminal equality constraints given
// the current state, the target state and the parameters.
// Final difference.
template<typename T>
DACE::AlgebraicVector<T> terminal_equality_constraints(
	DACE::AlgebraicVector<T> const& x, DACE::vectordb const& x_goal,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {
	// Get vectors
	DACE::AlgebraicVector<T> x_(x.extract(0, SIZE_VECTOR - 1));
	DACE::vectordb x_goal_(x_goal.extract(0, SIZE_VECTOR - 1));

	// Compute error
	DACE::AlgebraicVector<T> loss = (x_ - x_goal_);

	return loss;
}

// Returns the terminal equality constraints given
// the current state, the target state and the parameters.
// Final difference, without equinoctial anomaly.
template<typename T>
DACE::AlgebraicVector<T> terminal_equality_constraints_equinoctial(
	DACE::AlgebraicVector<T> const& x, DACE::vectordb const& x_goal,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {
	// Get vectors
	DACE::AlgebraicVector<T> x_(x.extract(0, SIZE_VECTOR - 2));
	DACE::vectordb x_goal_(x_goal.extract(0, SIZE_VECTOR - 2));

	// Compute error
	DACE::AlgebraicVector<T> loss = (x_ - x_goal_);

	return loss;
}

// Returns the terminal inequality constraints given
// the current state, the target state and the parameters.
// Null.
template<typename T>
DACE::AlgebraicVector<T> terminal_inequality_constraints(
	DACE::AlgebraicVector<T> const& x, DACE::vectordb  const& x_goal,
	SpacecraftParameters const& spacecraft_parameters,
	Constants const& constants,
	SolverParameters const& solver_parameters) {
	return DACE::AlgebraicVector<T>(0);
}

// Transformations

// Transforms a keplerian state vector into a cartesian one 
DACE::vectordb kep_2_cart(
	DACE::vectordb const& kep_state_vector,
	double const& mu);

// Transforms a keplerian state vector into an equinoctial one 
DACE::vectordb equi_2_kep(DACE::vectordb const& equi_state_vector);

// Transforms a keplerian state vector into an equinoctial one 
DACE::vectordb kep_2_equi(DACE::vectordb const& kep_state_vector);

// Transforms RTN reference frame into cartesian coordinates
DACE::vectordb RTN_2_cart(
	DACE::vectordb const& RTN_vector,
	DACE::vectordb const& cart_state_vector);

// Returns dynamics with acceleration_2b_SUN as accelerations.
// Terminal constraints and thrust constraints.
Dynamics get_tbp_SUN_lt_dynamics();
Dynamics get_tbp_EARTH_lt_dynamics();
Dynamics get_cr3bp_EARTH_MOON_lt_dynamics();


#endif