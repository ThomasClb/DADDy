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

// Returns dynamics with acceleration_2bp_SUN as accelerations.
// Terminal constraints and thrust constraints.
Dynamics get_low_trust_2bp_SUN_dynamics() {
	Constants constants(MU_SUN, SUN_EARTH_DISTANCE,
		sqrt((MU_SUN) / pow(SUN_EARTH_DISTANCE, 3)), 1000);
	dynFunction dyn([](
		vectorDA const& a, vectorDA const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return dynamic_2bp_SUN(a, b, c, e, d); });
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
			return equality_constraints(a, b, c, e, d); });
	ineqFunction ineq([](
		vectorDA const& a, vectorDA const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return inequality_constraints_2bp_SUN(a, b, c, e, d); });
	tcFunction tc([](
		vectorDA const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return terminal_cost_2bp_SUN(a, b, c, e, d); });
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
			return terminal_inequality_constraints(a, b, c, e, d); });

	dynFunction_db dyn_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return dynamic_2bp_SUN(a, b, c, e, d); });
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
			return equality_constraints(a, b, c, e, d); });
	ineqFunction_db ineq_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return inequality_constraints_2bp_SUN(a, b, c, e, d); });
	tcFunction_db tc_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return terminal_cost_2bp_SUN(a, b, c, e, d); });
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
			return terminal_inequality_constraints(a, b, c, e, d); });
	return Dynamics(
		constants,
		dyn, ctg, eq, ineq, tc, teq, tineq,
		dyn_db, ctg_db, eq_db, ineq_db, tc_db, teq_db, tineq_db);
}

// Returns dynamics with acceleration_cr3bp as accelerations.
// Terminal constraints and thrust constraints.
Dynamics get_low_trust_cr3bp_dynamics() {
	Constants constants(MU_MOON / (MU_MOON + MU_EARTH), EARTH_MOON_DISTANCE,
		sqrt((MU_MOON + MU_EARTH) / pow(EARTH_MOON_DISTANCE, 3)), 1000);
	dynFunction dyn([](
		vectorDA const& a, vectorDA const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return dynamic_cr3bp(a, b, c, e, d); });
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
			return equality_constraints(a, b, c, e, d); });
	ineqFunction ineq([](
		vectorDA const& a, vectorDA const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return inequality_constraints_cr3bp(a, b, c, e, d); });
	tcFunction tc([](
		vectorDA const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return terminal_cost_cr3bp(a, b, c, e, d); });
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
			return terminal_inequality_constraints(a, b, c, e, d); });

	dynFunction_db dyn_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return dynamic_cr3bp(a, b, c, e, d); });
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
			return equality_constraints(a, b, c, e, d); });
	ineqFunction_db ineq_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return inequality_constraints_cr3bp(a, b, c, e, d); });
	tcFunction_db tc_db([](
		vectordb const& a, vectordb const& b,
		SpacecraftParameters const& c,
		Constants const& e,
		SolverParameters const& d) {
			return terminal_cost_cr3bp(a, b, c, e, d); });
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
			return terminal_inequality_constraints(a, b, c, e, d); });
	return Dynamics(
		constants,
		dyn, ctg, eq, ineq, tc, teq, tineq,
		dyn_db, ctg_db, eq_db, ineq_db, tc_db, teq_db, tineq_db);
}
