/**
	parameters.h

	Purpose: Implementation of the SpacecraftParameters
	and SolverParameters classes.

	@author Thomas Caleb

	@version 1.0 14/11/2023
*/

#ifndef DEF_PARAMETERS
#define DEF_PARAMETERS

#pragma once

#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include <dace/dace_s.h>

#include "settings.h"
#include "constants.h"


/*

	SPACECRAFT PARAMETERS

*/
class SpacecraftParameters {

// Attributes
protected:
	Constants constants_; // Dynamics normalisation constants 
	double initial_mass_; // Spacecraft mass [MASSU]
	double dry_mass_; // Spacecraft dry mass [MASSU]
	double thrust_; // Spacecraft thrust [THRUSTU]
	double Isp_; // Spacecraft thrust Isp [TU]

// Methods
public:
	// Constructors

	// Default constructors
	SpacecraftParameters();

	// Constructors
	SpacecraftParameters(Constants const& constants);
	SpacecraftParameters(
		Constants const& constants,
		double const& initial_mass,
		double const& dry_mass,
		double const& thrust,
		double const& Isp);

	// Copy constructor
	SpacecraftParameters(
		SpacecraftParameters const& param);

	// Destructors
	~SpacecraftParameters();

	// Getters
	const Constants constants() const;
	const double initial_mass() const;
	const double dry_mass() const;
	const double initial_wet_mass() const;
	const double thrust() const;
	const double Isp() const;
	const double ejection_velocity() const;
	const double mass_flow() const;
};


/*

	SOLVER PARAMETERS

*/
class SolverParameters {

// Attributes
protected:
	unsigned int N_; // Number of steps [-]
	unsigned int Nx_; // Size of state vector including dt [-]
	unsigned int Nu_; // Size of control vector [-]
	unsigned int Neq_; // Size of equality constraints vector [-]
	unsigned int Nineq_; // Size of equality constraints vector [-]
	unsigned int Nteq_; // Size of terminal equality constraints vector [-]
	unsigned int Ntineq_; // Size of terminal equality constraints vector [-]
	double ToF_; // Time-of-flight, initialized at solving [TU]
	double cost_to_go_gain_; // ctg magnitude [-]
	double terminal_cost_gain_; // tc magnitude [-]
	double homotopy_coefficient_; // Homotpy coefficient between the NRJ optimal and the fuel optimal [0, 1]
	double huber_loss_coefficient_; // coefficient of the Huber loss regularisation
	unsigned int DDP_type_; // Choice of the DDP method 0 = classic, 1 = DA-based
	double DDP_tol_; // DDPSolver tolerance [-]
	double AUL_tol_; // AULSolver tolerance [-]
	double PN_tol_; // PNSolver tolerance [-]
	unsigned int DDP_max_iter_; // Maximum number of iterations for DDP [-]
	unsigned int AUL_max_iter_; // Maximum number of iterations for AUL [-]
	unsigned int PN_max_iter_; // Maximum number of iterations for PN [-]
	std::vector<DACE::vectordb> list_lambda_; // List of Lagrange multiplicator
	std::vector<DACE::vectordb> list_mu_; // List of penalty factors
	DACE::vectordb line_search_parameters_; // [alpha_lb, alpha_ub, alpha_factor, max_iter]
	bool backward_sweep_regulation_; // Regulation on back 
	DACE::vectordb backward_sweep_regulation_parameters_; // [rho_initial_value, rho_lb, rho_ub, rho_factor]
	DACE::vectordb lambda_parameters_; // [lambda_initial_value, lambda_ub]
	DACE::vectordb mu_parameters_; // [mu_initial_value, mu_ub, mu_factor]
	double PN_regularisation_; // PN regulation
	double PN_active_constraint_tol_; // PN active constraint tol
	double PN_cv_rate_threshold_; // PN convergence rate threshold
	double PN_alpha_; // PN intial line search parameter
	double PN_gamma_; // PN line search reduction factor
	unsigned int verbosity_; // Quantity of display data, 0=full, 1=no DDP.

// Methods
public:
	// Constructors

	// Default constructor
	SolverParameters();

	// Constructor
	SolverParameters(
		unsigned int const& N,
		unsigned int const& Nx, unsigned int const& Nu,
		unsigned int const& Neq, unsigned int const& Nineq,
		unsigned int const& Nteq, unsigned int const& Ntineq,
		double const& cost_to_go_gain, double const& terminal_cost_gain,
		double const& homotopy_coefficient, double const& huber_loss_coefficient,
		unsigned int const& DDP_type,
		double const& DDP_tol, double const& AUL_tol, double const& PN_tol,
		unsigned int const& DDP_max_iter, unsigned int const& AUL_max_iter,
		unsigned int const& PN_max_iter,
		DACE::vectordb const& line_search_parameters,
		bool const& backward_sweep_regulation,
		DACE::vectordb const& backward_sweep_regulation_parameters,
		DACE::vectordb const& lambda_parameters, DACE::vectordb const& mu_parameters,
		double const& PN_regularisation, double const& PN_active_constraint_tol,
		double const& PN_cv_rate_threshold, double const& PN_alpha,
		double const& PN_gamma, unsigned int const& verbosity);

	// Copy constructor
	SolverParameters(SolverParameters const& param);

	// Destructors
	~SolverParameters();

	// Getters
	const unsigned int N() const;
	const unsigned int Nx() const;
	const unsigned int Nu() const;
	const unsigned int Neq() const;
	const unsigned int Nineq() const;
	const unsigned int Nteq() const;
	const unsigned int Ntineq() const;
	const double ToF() const;
	const double cost_to_go_gain() const;
	const double terminal_cost_gain() const;
	const double homotopy_coefficient() const;
	const double huber_loss_coefficient() const;
	const unsigned int DDP_type() const;
	const double DDP_tol() const;
	const double AUL_tol() const;
	const double PN_tol() const;
	const unsigned int DDP_max_iter() const;
	const unsigned int AUL_max_iter() const;
	const unsigned int PN_max_iter() const;
	const std::vector<DACE::vectordb> list_lambda() const;
	const std::vector<DACE::vectordb> list_mu() const;
	const DACE::vectordb line_search_parameters() const;
	const bool backward_sweep_regulation() const;
	const DACE::vectordb backward_sweep_regulation_parameters() const;
	const DACE::vectordb lambda_parameters() const;
	const DACE::vectordb mu_parameters() const;
	const double PN_regularisation() const;
	const double PN_active_constraint_tol() const;
	const double PN_cv_rate_threshold() const;
	const double PN_alpha() const;
	const double PN_gamma() const;
	const unsigned int verbosity() const;

	// Setters
	void set_homotopy_coefficient(double const& homotopy_coefficient);
	void set_huber_loss_coefficient(double const& huber_loss_coefficient);
	void set_ToF(double const& ToF);
	void set_list_lambda(std::vector<DACE::vectordb> const& list_lambda);
	void set_list_mu(std::vector<DACE::vectordb> const& list_mu);

};
#endif