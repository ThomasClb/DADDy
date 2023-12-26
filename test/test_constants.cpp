/**
	test_constants.cpp

	Purpose: Test of the implementation of the Constants class.

	@author Thomas Caleb

	@version 1.0 24/12/2023
*/

#include "pch.h"
#include "constants.cpp"

using namespace DACE;
using namespace std;

/*

	CONSTANTS

*/

// Constructors
TEST(TestConstants, EmptyConstructor) {
	// Init
	Constants constants;
	double mu = MU_SUN; // [km^3.s^-2]
	double lu = SUN_EARTH_DISTANCE; // [km]
	double wu = sqrt((MU_SUN) / pow(SUN_EARTH_DISTANCE, 3)); // [s^-1]
	double tu = 1/wu; // [s]
	double vu = lu*wu; // [m.s^-1]
	double massu = 1000; // [kg]
	double thrustu = 1000 * vu * massu * wu; // [N]

	// Tests
	EXPECT_EQ(constants.mu(), mu);
	EXPECT_EQ(constants.lu(), lu);
	EXPECT_EQ(constants.wu(), wu);
	EXPECT_EQ(constants.tu(), tu);
	EXPECT_EQ(constants.vu(), vu);
	EXPECT_EQ(constants.massu(), massu);
	EXPECT_EQ(constants.thrustu(), thrustu);
}
TEST(TestConstants, FilledConstructor) {
	// Init
	double mu = MU_MOON/(MU_EARTH + MU_MOON); // [-]
	double lu = EARTH_MOON_DISTANCE; // [km]
	double wu = sqrt((MU_EARTH + MU_MOON) / pow(EARTH_MOON_DISTANCE, 3)); // [s^-1]
	double tu = 1 / wu; // [s]
	double vu = lu * wu; // [m.s^-1]
	double massu = 1000; // [kg]
	double thrustu = 1000 * vu * massu * wu; // [N]
	Constants constants(mu, lu, wu, massu);

	// Tests
	EXPECT_EQ(constants.mu(), mu);
	EXPECT_EQ(constants.lu(), lu);
	EXPECT_EQ(constants.wu(), wu);
	EXPECT_EQ(constants.tu(), tu);
	EXPECT_EQ(constants.vu(), vu);
	EXPECT_EQ(constants.massu(), massu);
	EXPECT_EQ(constants.thrustu(), thrustu);
}
