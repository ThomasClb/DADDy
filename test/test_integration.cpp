/**
	test_integration.cpp

	Purpose: Test of the implementation of the integration
	methods.

	@author Thomas Caleb

	@version 1.0 16/11/2023
*/

#include "pch.h"
#include "integration.cpp"
#include "dynamics.h"

using namespace DACE;
using namespace std;

TEST(TestIntegration, MinMaxNorm) {
	// Init
	double a = 10;
	double b = -15;
	vector<double> vect{ a, b };

	// Tests
	EXPECT_EQ(min(a, b), b);
	EXPECT_EQ(max(a, b), a);
	EXPECT_EQ(normtmp(vect.size(), vect), abs(b));
}

TEST(TestIntegration, RK4) {
	// Init
	vectordb x0(SIZE_VECTOR + 2, 1.0);
	vectordb u(3, 1.0);
	double t = 0; double dt = 1.0;
	SpacecraftParameters sp;
	Constants constants;
	vectordb xf = RK4(acceleration_2bp_SUN, x0, u, t, dt, sp, constants);

	// Tests
	EXPECT_EQ(xf.size(), x0.size());
}

TEST(TestIntegration, RK78) {
	// Init DACE
	size_t nb_variables = 3;
	DA::init(2, nb_variables);

	// Init
	vectorDA x0(SIZE_VECTOR + 2, 1.0 + 1 * DA(2));
	vectorDA u(3, 1.0 + DA(1));
	double t = 0; double dt = 1.0;
	SpacecraftParameters sp;
	Constants constants;
	vectorDA xf = RK78(acceleration_2bp_SUN, x0, u, t, dt, sp, constants);

	// Tests
	EXPECT_EQ(xf.size(), x0.size());

	// For doubles
	/**/
	vectordb xf_db = RK78(acceleration_2bp_SUN, x0.cons(), u.cons(), t, dt, sp, constants);

	// Tests
	EXPECT_EQ(xf_db.size(), xf_db.size());
	EXPECT_EQ(xf_db, xf.cons());
}