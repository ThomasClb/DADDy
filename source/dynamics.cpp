/**
	dynamics.cpp

	Purpose: Implementation of the Dynamics class
	methods.

	@author Thomas Caleb

	@version 1.0 16/11/2023
*/

#include "dynamics.h"

using namespace DACE;
using namespace std;


/*

	DYNAMICS

*/


// Constructors

// Default constructors
Dynamics::Dynamics() {
	Constants constants; constants_ = constants;
	dynamic_ = dynFunction([](
		vectorDA const&, vectorDA const&,
		SpacecraftParameters const&, Constants const&, SolverParameters const&) {
			return vectorDA(); });
	cost_to_go_ = ctgFunction([](
		vectorDA const&, vectorDA const&,
		SpacecraftParameters const&, Constants const&, SolverParameters const&) {
			return 0.0; });
	equality_constraints_ = eqFunction([](
		vectorDA const&, vectorDA const&,
		SpacecraftParameters const&, Constants const&, SolverParameters const&) {
			return vectorDA(); });
	inequality_constraints_ = ineqFunction([](
		vectorDA const&, vectorDA const&,
		SpacecraftParameters const&, Constants const&, SolverParameters const&) {
			return vectorDA(); });
	terminal_cost_ = tcFunction([](
		vectorDA const&, vectordb const&,
		SpacecraftParameters const&, Constants const&, SolverParameters const&) {
			return 0.0; });
	terminal_equality_constraints_ = teqFunction([](
		vectorDA const&, vectordb const&,
		SpacecraftParameters const&, Constants const&, SolverParameters const&) {
			return vectorDA(); });
	terminal_inequality_constraints_ = tineqFunction([](
		vectorDA const&, vectordb const&,
		SpacecraftParameters const&, Constants const&, SolverParameters const&) {
			return vectorDA(); });
	dynamic_db_ = dynFunction_db([](
		vectordb const&, vectordb const&,
		SpacecraftParameters const&, Constants const&, SolverParameters const&) {
			return vectordb(); });
	cost_to_go_db_ = ctgFunction_db([](
		vectordb const&, vectordb const&,
		SpacecraftParameters const&, Constants const&, SolverParameters const&) {
			return 0.0; });
	equality_constraints_db_ = eqFunction_db([](
		vectordb const&, vectordb const&,
		SpacecraftParameters const&, Constants const&, SolverParameters const&) {
			return vectordb(); });
	inequality_constraints_db_ = ineqFunction_db([](
		vectordb const&, vectordb const&,
		SpacecraftParameters const&, Constants const&, SolverParameters const&) {
			return vectordb(); });
	terminal_cost_db_ = tcFunction_db([](
		vectordb const&, vectordb const&,
		SpacecraftParameters const&, Constants const&, SolverParameters const&) {
			return 0.0; });
	terminal_equality_constraints_db_ = teqFunction_db([](
		vectordb const&, vectordb const&,
		SpacecraftParameters const&, Constants const&, SolverParameters const&) {
			return vectordb(); });
	terminal_inequality_constraints_db_ = tineqFunction_db([](
		vectordb const&, vectordb const&,
		SpacecraftParameters const&, Constants const&, SolverParameters const&) {
			return vectordb(); });
}

// Constructor
Dynamics::Dynamics(
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
	tineqFunction_db const& terminal_inequality_constraints_db) :
	constants_(constants),
	dynamic_(dynamic),
	cost_to_go_(cost_to_go), terminal_cost_(terminal_cost),
	equality_constraints_(equality_constraints), inequality_constraints_(inequality_constraints),
	terminal_equality_constraints_(terminal_equality_constraints),
	terminal_inequality_constraints_(terminal_inequality_constraints),
	dynamic_db_(dynamic_db),
	cost_to_go_db_(cost_to_go_db), terminal_cost_db_(terminal_cost_db),
	equality_constraints_db_(equality_constraints_db), inequality_constraints_db_(inequality_constraints_db),
	terminal_equality_constraints_db_(terminal_equality_constraints_db),
	terminal_inequality_constraints_db_(terminal_inequality_constraints_db) {}

// Copy constructor
Dynamics::Dynamics(
	Dynamics const& dynamics) : 
	constants_(dynamics.constants_),
	dynamic_(dynamics.dynamic_),
	cost_to_go_(dynamics.cost_to_go_), terminal_cost_(dynamics.terminal_cost_),
	equality_constraints_(dynamics.equality_constraints_),
	inequality_constraints_(dynamics.inequality_constraints_),
	terminal_inequality_constraints_(dynamics.terminal_inequality_constraints_),
	terminal_equality_constraints_(dynamics.terminal_equality_constraints_),
	dynamic_db_(dynamics.dynamic_db_),
	cost_to_go_db_(dynamics.cost_to_go_db_), terminal_cost_db_(dynamics.terminal_cost_db_),
	equality_constraints_db_(dynamics.equality_constraints_db_),
	inequality_constraints_db_(dynamics.inequality_constraints_db_),
	terminal_inequality_constraints_db_(dynamics.terminal_inequality_constraints_db_),
	terminal_equality_constraints_db_(dynamics.terminal_equality_constraints_db_) {}

// Destructor
Dynamics::~Dynamics() {}

// Getters
const Constants Dynamics::constants() const { return constants_; }
const dynFunction Dynamics::dynamic() const { return dynamic_; }
const ctgFunction Dynamics::cost_to_go() const { return cost_to_go_; }
const eqFunction Dynamics::equality_constraints() const { return equality_constraints_; }
const ineqFunction Dynamics::inequality_constraints() const { return inequality_constraints_; }
const tcFunction Dynamics::terminal_cost() const { return terminal_cost_; }
const teqFunction Dynamics::terminal_equality_constraints() const {
	return terminal_equality_constraints_;
}
const tineqFunction Dynamics::terminal_inequality_constraints() const {
	return terminal_inequality_constraints_;
}

const dynFunction_db Dynamics::dynamic_db() const { return dynamic_db_; }
const ctgFunction_db Dynamics::cost_to_go_db() const { return cost_to_go_db_; }
const eqFunction_db Dynamics::equality_constraints_db() const { return equality_constraints_db_; }
const ineqFunction_db Dynamics::inequality_constraints_db() const { return inequality_constraints_db_; }
const tcFunction_db Dynamics::terminal_cost_db() const { return terminal_cost_db_; }
const teqFunction_db Dynamics::terminal_equality_constraints_db() const {
	return terminal_equality_constraints_db_;
}
const tineqFunction_db Dynamics::terminal_inequality_constraints_db() const {
	return terminal_inequality_constraints_db_;
}

/*

	FUNCTIONS

*/

// Transformations

// Transforms a keplerian state vector into an equinoctial one.
vectordb kep_2_equi(vectordb const& kep_state_vector) {
	// Unpack
	double sma(kep_state_vector[0]);
	double ecc(kep_state_vector[1]);
	double inc(kep_state_vector[2]);
	double RAAN(kep_state_vector[3]);
	double omega(kep_state_vector[4]);
	double M(kep_state_vector[5]);

	// Computing eccentric anomaly (Newton's method)
	double ecc_anomaly_init(0.0);
	double ecc_anomaly(M);

	// Pi is init value is e is too high
	if (ecc > 0.8) {
		ecc_anomaly = PI;
	}

	// Newton's method
	do {
		ecc_anomaly_init = ecc_anomaly;
		ecc_anomaly = M + ecc * sin(ecc_anomaly_init);
		ecc_anomaly -= ecc * ecc_anomaly_init * cos(ecc_anomaly_init);
		ecc_anomaly /= 1 - ecc * cos(ecc_anomaly_init);
	} while (abs(ecc_anomaly - ecc_anomaly_init) > EPS);

	// Compute elements
	double P_1 = ecc * sin(RAAN + omega);
	double P_2 = ecc * cos(RAAN + omega);
	double Q_1 = tan(inc / 2) * sin(RAAN);
	double Q_2 = tan(inc / 2) * cos(RAAN);
	double true_anomaly = (2 * atan(sqrt((1 + ecc) / (1 - ecc)) * tan(ecc_anomaly / 2)));
	double L = RAAN + omega + true_anomaly;
	
	// Store values in a vector
	vectordb output(kep_state_vector);
	output[0] = sma;
	output[1] = P_1;
	output[2] = P_2;
	output[3] = Q_1;
	output[4] = Q_2;
	output[5] = L;

	return output;
}

// Transforms a equinoctial state vector into an keplerian one.
vectordb equi_2_kep(vectordb const& equi_state_vector) {
	// Unpack
	double sma(equi_state_vector[0]);
	double P_1(equi_state_vector[1]);
	double P_2(equi_state_vector[2]);
	double Q_1(equi_state_vector[3]);
	double Q_2(equi_state_vector[4]);
	double L(equi_state_vector[5]);

	// Compute elements
	double ecc = sqrt(P_1*P_1 + P_2*P_2);
	double inc = 2 * atan(sqrt(Q_1 * Q_1 + Q_2 * Q_2));
	double RAAN = atan2(Q_1, Q_2);
	double omega = atan2(P_1, P_2) - RAAN;
	double true_anomaly = L - omega - RAAN;
	double ecc_anomaly = 2 * atan(sqrt((1 - ecc) / (1 + ecc)) * tan(true_anomaly / 2));
	double M = fmod(ecc_anomaly - ecc *sin(ecc_anomaly), 2*PI);

	// Store values in a vector
	vectordb output(equi_state_vector);
	output[0] = sma;
	output[1] = ecc;
	output[2] = inc;
	output[3] = RAAN;
	output[4] = omega;
	output[5] = M;

	return output;
}

// Transforms a keplerian state vector into a cartesian one.
vectordb kep_2_cart(vectordb const& kep_state_vector, double const& mu) {
	// Unpack
	double sma(kep_state_vector[0]);
	double ecc(kep_state_vector[1]);
	double inc(kep_state_vector[2]);
	double RAAN(kep_state_vector[3]);
	double omega(kep_state_vector[4]);
	double M(kep_state_vector[5]);

	// Computing eccentric anomaly (Newton's method)
	double ecc_anomaly_init(0.0);
	double ecc_anomaly(M);

	// Pi is init value is e is too high
	if (ecc > 0.8) {
		ecc_anomaly = PI;
	}

	// Newton's method
	do {
		ecc_anomaly_init = ecc_anomaly;
		ecc_anomaly = M + ecc * sin(ecc_anomaly_init);
		ecc_anomaly -= ecc * ecc_anomaly_init * cos(ecc_anomaly_init);
		ecc_anomaly /= 1 - ecc * cos(ecc_anomaly_init);
	} while (abs(ecc_anomaly - ecc_anomaly_init) > EPS);

	// Computing true anomaly
	ecc_anomaly_init = sqrt((1 + ecc) / (1 - ecc)) * tan(ecc_anomaly / 2);
	double v(2 * atan(ecc_anomaly_init));

	// Computing radius, ellipse parameter, and angular momentom
	double r(sma * (1.0 - ecc * cos(ecc_anomaly)));
	double p(sma * (1.0 - sqr(ecc)));

	// Get mu
	double h(sqrt(mu * p));

	// Computing positions
	double x(r * (cos(RAAN) * cos(omega + v) - sin(RAAN) * sin(omega + v) * cos(inc)));
	double y(r * (sin(RAAN) * cos(omega + v) + cos(RAAN) * sin(omega + v) * cos(inc)));
	double z(r * sin(omega + v) * sin(inc));

	// Computing velocities
	double d_x(h / r * (ecc * sin(v) * x / p -
		(cos(RAAN) * sin(omega + v) + sin(RAAN) * cos(omega + v) * cos(inc))));
	double d_y(h / r * (ecc * sin(v) * y / p -
		(sin(RAAN) * sin(omega + v) - cos(RAAN) * cos(omega + v) * cos(inc))));
	double d_z(h / r * (ecc * sin(v) * z / p + cos(omega + v) * sin(inc)));

	// Store values in a vector
	vectordb output(kep_state_vector);
	output[0] = x;
	output[1] = y;
	output[2] = z;
	output[3] = d_x;
	output[4] = d_y;
	output[5] = d_z;

	return output;
}

// Transforms a cartesian state vector into a keplerian one.
vectordb cart_2_kep(
	vectordb const& cart_state_vector,
	double const& mu) {
	// Unpack
	vectordb r(cart_state_vector.extract(0, SIZE_VECTOR/2 - 1));
	vectordb v(cart_state_vector.extract(SIZE_VECTOR/2, SIZE_VECTOR - 1));
	
	// Get moment and n
	vectordb h(r.cross(v));
	vectordb k{0, 0, 1};
	vectordb n(k.cross(h));

	// Compute norms
	double r_norm = r.vnorm();
	double v_norm = v.vnorm();
	double h_norm = h.vnorm();
	double n_norm = n.vnorm();

	// Specific NRJ spécifique ε = v²/2 – μ/r
    double epsilon = 0.5*v_norm*v_norm - mu/r_norm;

	// SMA
	double a = -mu/(2.0*epsilon);

	// Eccentricity
	double dot_r_v = r.dot(v);
	vectordb e_vec = (v_norm*v_norm - mu/r_norm)*r/mu - dot_r_v*v;
    double e = e_vec.vnorm();

	// Inclinaison
	double i = acos(h[2]/h_norm);

	// RAAN
	double Omega = 0.0;
    if (n_norm > 1e-12) {                    
        Omega = atan2(n[1], n[0]);
        Omega = wrap_mod(Omega, 2*PI);
    }

	// omega
	double omega = 0.0;
    if (n_norm > 1e-12 && e > 1e-12) {
        // ω = atan2( (h·(n×e))/h , n·e )
        vectordb n_cross_e = n.cross(e_vec);
        double sin_omega = n_cross_e.dot(h) / (h_norm * e);
        double cos_omega = n.dot(e_vec) / (n_norm * e);
        omega = atan2(sin_omega, cos_omega);
        omega = wrap_mod(omega, 2*PI);
    }

	// True anomaly
	double nu = 0.0;
    if (e > 1e-12) {
        double cos_nu = e_vec.dot(r) / (e * r_norm);
        double sin_nu = (e_vec.cross(r)).dot(h) / (h_norm * e * r_norm);
        nu = atan2(sin_nu, cos_nu);
        nu = wrap_mod(nu, 2*PI);
    } else {
        // circular case
        vectordb r_hat = r / r_norm;
        vectordb v_hat = v / v_norm;
        double cos_nu = r_hat.dot(v_hat);
        double sin_nu = (r_hat.cross(v_hat)).dot(h) / h_norm;
        nu = atan2(sin_nu, cos_nu);
        nu = wrap_mod(nu, 2*PI);
    }

	// Eccetric anomaly
	double E = 2*atan(sqrt((1-e)/(1+e))*tan(nu/2));

	// Mean anomaly
	double M = E - e*sin(E);

	return vectordb{a, e, i, Omega, omega, M};
}

// Transforms RTN reference frame into cartesian coordinates.
vectordb RTN_2_cart(
	vectordb const& RTN_vector,
	vectordb const& cart_state_vector) {
	// Unpack
	vectordb r = cart_state_vector.extract(0, 2);
	vectordb v = cart_state_vector.extract(3, 5);

	// Make RTN triade
	vectordb N = (r.cross(v)).normalize(); // Angular momentum
	vectordb T = v.normalize();
	vectordb R = (T.cross(N)).normalize();

	// Make matrix
	matrixdb matrix(R.size());
	matrix.setcol(0, R);
	matrix.setcol(1, T);
	matrix.setcol(2, N);

	return matrix * RTN_vector;
}

// Transforms cartesian coordinates into RTN reference frame.
vectordb cart_2_RTN(
	vectordb const& cart_vector,
	vectordb const& cart_state_vector) {

	// Unpack
	vectordb r = cart_state_vector.extract(0, 2);
	vectordb v = cart_state_vector.extract(3, 5);

	// Make RTN triade
	vectordb N = (r.cross(v)).normalize(); // Angular momentum
	vectordb T = v.normalize();
	vectordb R = (T.cross(N)).normalize();

	// Make matrix
	matrixdb matrix(R.size());
	matrix.setcol(0, R);
	matrix.setcol(1, T);
	matrix.setcol(2, N);

	return matrix.transpose() * cart_vector;
}

// Returns dynamics of double integrator problem.
// Terminal constraints and thrust constraints.
Dynamics get_double_integrator_dynamics() {
	Constants constants(1, 1, 1, 1);
	dynFunction dyn([](
		vectorDA const& a, vectorDA const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return dynamic_double_integrator(a, b, c, e, d); });
	ctgFunction ctg([](
		vectorDA const& a, vectorDA const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return cost_to_go_double_integrator(a, b, c, e, d); });
	eqFunction eq([](
		vectorDA const& a, vectorDA const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return equality_constraints_null(a, b, c, e, d); });
	ineqFunction ineq([](
		vectorDA const& a, vectorDA const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return inequality_constraints_null(a, b, c, e, d); });
	tcFunction tc([](
		vectorDA const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return terminal_cost_double_integrator(a, b, c, e, d); });
	teqFunction teq([](
		vectorDA const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return terminal_equality_constraints_null(a, b, c, e, d); });
	tineqFunction tineq([](
		vectorDA const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return terminal_inequality_constraints_null(a, b, c, e, d); });

	dynFunction_db dyn_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return dynamic_double_integrator(a, b, c, e, d); });
	ctgFunction_db ctg_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return cost_to_go_double_integrator(a, b, c, e, d); });
	eqFunction_db eq_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return equality_constraints_null(a, b, c, e, d); });
	ineqFunction_db ineq_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return inequality_constraints_null(a, b, c, e, d); });
	tcFunction_db tc_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return terminal_cost_double_integrator(a, b, c, e, d); });
	teqFunction_db teq_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return terminal_equality_constraints_null(a, b, c, e, d); });
	tineqFunction_db tineq_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return terminal_inequality_constraints_null(a, b, c, e, d); });
	return Dynamics(
		constants,
		dyn, ctg, eq, ineq, tc, teq, tineq,
		dyn_db, ctg_db, eq_db, ineq_db, tc_db, teq_db, tineq_db);
}

// Returns dynamics with acceleration_2bp_SUN as accelerations.
// Terminal constraints and thrust constraints.
// Anomaly is not free.
Dynamics get_tbp_SUN_lt_dynamics() {
	Constants constants(MU_SUN, SUN_EARTH_DISTANCE,
		sqrt((MU_SUN) / pow(SUN_EARTH_DISTANCE, 3)), 1000);
	dynFunction dyn([](
		vectorDA const& a, vectorDA const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return dynamic_tbp_SUN_lt(a, b, c, e, d); });
	ctgFunction ctg([](
		vectorDA const& a, vectorDA const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return cost_to_go(a, b, c, e, d); });
	eqFunction eq([](
		vectorDA const& a, vectorDA const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return equality_constraints_null(a, b, c, e, d); });
	ineqFunction ineq([](
		vectorDA const& a, vectorDA const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return inequality_constraints_tbp_lt(a, b, c, e, d); });
	tcFunction tc([](
		vectorDA const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return terminal_cost(a, b, c, e, d); });
	teqFunction teq([](
		vectorDA const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return terminal_equality_constraints(a, b, c, e, d); });
	tineqFunction tineq([](
		vectorDA const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return terminal_inequality_constraints_null(a, b, c, e, d); });

	dynFunction_db dyn_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return dynamic_tbp_SUN_lt(a, b, c, e, d); });
	ctgFunction_db ctg_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return cost_to_go(a, b, c, e, d); });
	eqFunction_db eq_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return equality_constraints_null(a, b, c, e, d); });
	ineqFunction_db ineq_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return inequality_constraints_tbp_lt(a, b, c, e, d); });
	tcFunction_db tc_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return terminal_cost(a, b, c, e, d); });
	teqFunction_db teq_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return terminal_equality_constraints(a, b, c, e, d); });
	tineqFunction_db tineq_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return terminal_inequality_constraints_null(a, b, c, e, d); });
	return Dynamics(
		constants,
		dyn, ctg, eq, ineq, tc, teq, tineq,
		dyn_db, ctg_db, eq_db, ineq_db, tc_db, teq_db, tineq_db);
}

// Returns dynamics with acceleration_2bp_EARTH as accelerations.
// Terminal constraints and thrust constraints.
// Anomaly is free.
Dynamics get_tbp_EARTH_lt_dynamics() {
	double a_GEO = pow(MU_EARTH * pow(SEC2DAYS, -2.0) / pow(2 * PI, 2), 1.0 / 3.0);
	Constants constants(MU_EARTH, a_GEO,
		sqrt((MU_EARTH) / pow(a_GEO, 3)), 1000);
	dynFunction dyn([](
		vectorDA const& a, vectorDA const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return dynamic_tbp_EARTH_lt(a, b, c, e, d); });
	ctgFunction ctg([](
		vectorDA const& a, vectorDA const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return cost_to_go(a, b, c, e, d); });
	eqFunction eq([](
		vectorDA const& a, vectorDA const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return equality_constraints_null(a, b, c, e, d); });
	ineqFunction ineq([](
		vectorDA const& a, vectorDA const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return inequality_constraints_tbp_lt(a, b, c, e, d); });
	tcFunction tc([](
		vectorDA const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return terminal_cost_equinoctial(a, b, c, e, d); });
	teqFunction teq([](
		vectorDA const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return terminal_equality_constraints_equinoctial(a, b, c, e, d); });
	tineqFunction tineq([](
		vectorDA const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return terminal_inequality_constraints_null(a, b, c, e, d); });

	dynFunction_db dyn_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return dynamic_tbp_EARTH_lt(a, b, c, e, d); });
	ctgFunction_db ctg_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return cost_to_go(a, b, c, e, d); });
	eqFunction_db eq_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return equality_constraints_null(a, b, c, e, d); });
	ineqFunction_db ineq_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return inequality_constraints_tbp_lt(a, b, c, e, d); });
	tcFunction_db tc_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return terminal_cost_equinoctial(a, b, c, e, d); });
	teqFunction_db teq_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return terminal_equality_constraints_equinoctial(a, b, c, e, d); });
	tineqFunction_db tineq_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return terminal_inequality_constraints_null(a, b, c, e, d); });
	return Dynamics(
		constants,
		dyn, ctg, eq, ineq, tc, teq, tineq,
		dyn_db, ctg_db, eq_db, ineq_db, tc_db, teq_db, tineq_db);
}

// Returns dynamics with acceleration_2bp_EARTH as accelerations.
// Terminal constraints and thrust constraints.
// Anomaly is free.
Dynamics get_tbp_EARTH_geo_lt_dynamics() {
	double a_GEO = pow(MU_EARTH * pow(SEC2DAYS, -2.0) / pow(2 * PI, 2), 1.0 / 3.0);
	Constants constants(MU_EARTH, a_GEO,
		sqrt((MU_EARTH) / pow(a_GEO, 3)), 1000);
	dynFunction dyn([](
		vectorDA const& a, vectorDA const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return dynamic_tbp_EARTH_lt(a, b, c, e, d); });
	ctgFunction ctg([](
		vectorDA const& a, vectorDA const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return cost_to_go(a, b, c, e, d); });
	eqFunction eq([](
		vectorDA const& a, vectorDA const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return equality_constraints_null(a, b, c, e, d); });
	ineqFunction ineq([](
		vectorDA const& a, vectorDA const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return inequality_constraints_tbp_lt(a, b, c, e, d); });
	tcFunction tc([](
		vectorDA const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return terminal_cost_geo(a, b, c, e, d); });
	teqFunction teq([](
		vectorDA const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return terminal_equality_constraints_geo(a, b, c, e, d); });
	tineqFunction tineq([](
		vectorDA const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return terminal_inequality_constraints_null(a, b, c, e, d); });

	dynFunction_db dyn_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return dynamic_tbp_EARTH_lt(a, b, c, e, d); });
	ctgFunction_db ctg_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return cost_to_go(a, b, c, e, d); });
	eqFunction_db eq_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return equality_constraints_null(a, b, c, e, d); });
	ineqFunction_db ineq_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return inequality_constraints_tbp_lt(a, b, c, e, d); });
	tcFunction_db tc_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return terminal_cost_geo(a, b, c, e, d); });
	teqFunction_db teq_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return terminal_equality_constraints_geo(a, b, c, e, d); });
	tineqFunction_db tineq_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return terminal_inequality_constraints_null(a, b, c, e, d); });
	return Dynamics(
		constants,
		dyn, ctg, eq, ineq, tc, teq, tineq,
		dyn_db, ctg_db, eq_db, ineq_db, tc_db, teq_db, tineq_db);
}


// Returns dynamics with acceleration_cr3bp as accelerations.
// Terminal constraints and thrust constraints.
Dynamics get_cr3bp_EARTH_MOON_lt_dynamics() {
	Constants constants(MU_MOON / (MU_MOON + MU_EARTH), EARTH_MOON_DISTANCE,
		sqrt((MU_MOON + MU_EARTH) / pow(EARTH_MOON_DISTANCE, 3)), 1000);
	dynFunction dyn([](
		vectorDA const& a, vectorDA const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return dynamic_cr3bp_lt(a, b, c, e, d); });
	ctgFunction ctg([](
		vectorDA const& a, vectorDA const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return cost_to_go(a, b, c, e, d); });
	eqFunction eq([](
		vectorDA const& a, vectorDA const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return equality_constraints_null(a, b, c, e, d); });
	ineqFunction ineq([](
		vectorDA const& a, vectorDA const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return inequality_constraints_cr3bp_EARTH_MOON_lt(a, b, c, e, d); });
	tcFunction tc([](
		vectorDA const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return terminal_cost(a, b, c, e, d); });
	teqFunction teq([](
		vectorDA const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return terminal_equality_constraints(a, b, c, e, d); });
	tineqFunction tineq([](
		vectorDA const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return terminal_inequality_constraints_null(a, b, c, e, d); });

	dynFunction_db dyn_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return dynamic_cr3bp_lt(a, b, c, e, d); });
	ctgFunction_db ctg_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return cost_to_go(a, b, c, e, d); });
	eqFunction_db eq_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return equality_constraints_null(a, b, c, e, d); });
	ineqFunction_db ineq_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return inequality_constraints_cr3bp_EARTH_MOON_lt(a, b, c, e, d); });
	tcFunction_db tc_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return terminal_cost(a, b, c, e, d); });
	teqFunction_db teq_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return terminal_equality_constraints(a, b, c, e, d); });
	tineqFunction_db tineq_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return terminal_inequality_constraints_null(a, b, c, e, d); });
	return Dynamics(
		constants,
		dyn, ctg, eq, ineq, tc, teq, tineq,
		dyn_db, ctg_db, eq_db, ineq_db, tc_db, teq_db, tineq_db);
}
