/**
	constants.cpp

	Purpose: Implementation of the Constants classes.

	@author Thomas Caleb

	@version 1.0 24/12/2023
*/

#include "constants.h"

using namespace DACE;
using namespace std;


/*

	CONSTANTS

*/

// Constructors

// Default constructors
// Sun-centered 2-body problem
// Similar as in [Caleb et al. 2023]
// DOI: https://doi.org/10.1007/s11071-023-08375-0
Constants::Constants() :
	mu_(MU_SUN), lu_(SUN_EARTH_DISTANCE),
	wu_(sqrt((MU_SUN) / pow(SUN_EARTH_DISTANCE, 3))),
	massu_(1000) {
	tu_ = 1 / wu_;
	vu_ = lu_ * wu_;
	thrustu_ = 1000 * vu_ * massu_ * wu_;
}

// Constructor
// Similar as in [Caleb et al. 2023]
// DOI: https://doi.org/10.1007/s11071-023-08375-0
Constants::Constants(
	double const& mu,
	double const& lu,
	double const& wu,
	double const& massu) :
	mu_(mu), lu_(lu), wu_(wu), massu_(massu) {
	tu_ = 1 / wu_;
	vu_ = lu_ * wu_;
	thrustu_ = 1000 * vu_ * massu_ * wu_;
}

// Copy constructor
Constants::Constants(Constants const& constants) :
	mu_(constants.mu_), lu_(constants.lu_),
	wu_(constants.wu_), massu_(constants.massu_),
	tu_(constants.tu_),
	vu_(constants.vu_), thrustu_(constants.thrustu_) {}

// Destructor
Constants::~Constants() {}

// Getters
const double Constants::mu() const { return mu_; }
const double Constants::lu() const { return lu_; }
const double Constants::wu() const { return wu_; }
const double Constants::massu() const { return massu_; }
const double Constants::tu() const { return tu_; }
const double Constants::vu() const { return vu_; }
const double Constants::thrustu() const { return thrustu_; }

// IO operator
ostream& operator<<(ostream& os, const Constants& constants) {
	// Set double precision
	typedef std::numeric_limits<double> dbl;
	os.precision(dbl::max_digits10);

	// Write attributs
	os << constants.mu() << endl;
	os << constants.lu() << endl;
	os << constants.wu() << endl;
	os << constants.massu() << endl;

	return os;
}
istream& operator>>(istream& is, Constants& constants) {

	// Reading simple property from a line
	string mu_str, lu_str, wu_str, massu_str;
	getline(is, mu_str);
	getline(is, lu_str);
	getline(is, wu_str);
	getline(is, massu_str);
	istringstream(mu_str) >> constants.mu_;
	istringstream(lu_str) >> constants.lu_;
	istringstream(wu_str) >> constants.wu_;
	istringstream(massu_str) >> constants.massu_;

	// Compute others 
	constants.tu_ = 1 / constants.wu_;
	constants.vu_ = constants.lu_ * constants.wu_;
	constants.thrustu_ = 1000 * constants.vu_ * constants.massu_ * constants.wu_;

	return is;
}
