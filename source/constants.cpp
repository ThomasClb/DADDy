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
