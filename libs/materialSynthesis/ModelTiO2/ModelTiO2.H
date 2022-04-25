/*-----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
|   License                                                               |
|                                                                         |
|   Copyright(C) 2022 Alberto Cuoci                                       |
|   OpenSMOKE++ is free software: you can redistribute it and/or modify   |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKE++ is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKE++. If not, see <http://www.gnu.org/licenses/>.   |
|                                                                         |
\*-----------------------------------------------------------------------*/

#pragma once

#include <string>
#include <vector>

class ModelTiO2 
{
public:

	/**
	*@brief Default constructor
	*/
	ModelTiO2();


	/**
	*@brief Summary on screen
	*/
	void SummaryOnScreen();


	/**
	*@brief Calculates the nucleation source terms for the 3-eq model: Ti(OH)4 + Ti(OH)4 => 2TiO2 + 4H2O
	*@param T temperature (K)
	*@param rho density (kg/m3)
	*@param cTiOH4 concentration of Ti(OH)4 (kmol/m3)
	*@param r reaction rate (kmol/m3/s)
	*@param OmegaYs source term: mass fraction of particles (1/s)
	*@param OmegaN source term: particle number density (#/m3/s)
	*@param OmegaS source term: total surface area concentration (m2/m3/s)
	*/
	void NucleationSourceTerms(	const double T, const double rho, const double cTiOH4,
					double& r, 
					double& OmegaYs, double& OmegaN, double& OmegaS	);

	/**
	*@brief Calculates the coagulation source terms for the 3-eq model
	*@param T temperature (K)
	*@param rho density (kg/m3)
	*@param mu dynamic viscosity (kg/m/s)
	*@param Ys mass fraction of particles (-)
	*@param N particle number density (#/m3)
	*@param S total surface area concentration (m2/m3)
	*@param OmegaYs source term: mass fraction of particles (1/s)
	*@param OmegaN source term: particle number density (#/m3/s)
	*@param OmegaS source term: total surface area concentration (m2/m3/s)
	*/
	void CoagulationSourceTerms( 	const double T, const double rho, const double mu, 
					const double Ys, const double N, const double S, 
					double& OmegaYs, double& OmegaN, double& OmegaS );

	/**
	*@brief Calculates the sintering source terms for the 3-eq model
	*@param T temperature (K)
	*@param rho density (kg/m3)
	*@param Ys mass fraction of particles (-)
	*@param N particle number density (#/m3)
	*@param S total surface area concentration (m2/m3)
	*@param OmegaYs source term: mass fraction of particles (1/s)
	*@param OmegaN source term: particle number density (#/m3/s)
	*@param OmegaS source term: total surface area concentration (m2/m3/s)
	*/
	void SinteringSourceTerms( 	const double T, const double rho,
					const double Ys, const double N, const double S, 
					double& OmegaYs, double& OmegaN, double& OmegaS );

	/**
	*@brief Calculates the coagulation source terms for the 3-eq model
	*@param T temperature (K)
	*@param rho density (kg/m3)
	*@param Ys mass fraction of particles (-)
	*@param N particle number density (#/m3)
	*@param S total surface area concentration (m2/m3)
	*/
	void Properties(	const double T, const double rho,
				const double Ys, const double N, const double S,
				double& fv, double& dp, double& dc, double& da, double& ssph, double& np, double& Taus );


	/**
	*@brief Calculates the mass diffusion coefficient
	*@param T temperature (K)
	*@param rho density (kg/m3)
	*@param mu dynamic viscosity (kg/m/s)
	*@param Sc Schmidt number (-)
	*@param Ys mass fraction of particles (-)
	*@param N particle number density (#/m3)
	*@param S total surface area concentration (m2/m3)
	*@param Gamma mass diffusion coefficient (kg/m/s)
	*/
	void DiffusionCoefficient(	const double T, const double rho, const double mu, const double Sc,
					const double Ys, const double N, const double S,
					double& Gamma );

	/**
	*@brief Set collision enhancement factor
	*@param epsilon collision enhancement factor (default: 2.2)
	*/
	void SetNucleationCollisionEnhancementFactor(const double epsilon);

	/**
	*@brief Set nucleation activation temperature
	*@param Ta activation temperature (default: 0 K)
	*/
	void SetNucleationActivationTemperature(const double Ta);

	/**
	*@brief Set nucleation temperature exponent
	*@param na temperature exponent (default: 0.50)
	*/
	void SetNucleationTemperatureExponent(const double na);

	/**
	*@brief Set sintering frequency factor
	*@param As frequency factor (default: 7.44e16 1/s/K)
	*/
	void SetSinteringFrequencyFactor(const double As);

	/**
	*@brief Set sintering activation temperature
	*@param Ts activation temperature (default: -31000 K)
	*/
	void SetSinteringActivationTemperature(const double Ts);

	/**
	*@brief Set sintering temperature exponent
	*@param na temperature exponent (default: 1.00)
	*/
	void SetSinteringTemperatureExponent(const double ns);

	/**
	*@brief Set number of particles
	*@param n number of particles (default: 5)
	*/
	void SetNumberOfParticles(const unsigned int n);

	/**
	*@brief Turn on/off the nucleation process
	*@param flag if true, nucleation is turned on (default: true)
	*/
	void SetNucleation(const bool flag);

	/**
	*@brief Turn on/off the coagulation process
	*@param flag if true, coagulation is turned on (default: true)
	*/
	void SetCoagulation(const bool flag);

	/**
	*@brief Turn on/off the sintering process
	*@param flag if true, sintering is turned on (default: true)
	*/
	void SetSintering(const bool flag);

	/**
	*@brief Set the threshold for N for the calculation of coagulation/sintering contributions
	*@param N_threshold threshold (#/m3)
	*/	
	void SetNThreshold(const double N_threshold);

	/**
	*@brief Set the threshold for volume fraction for calculation of properties
	*@param fv_threshold threshold (-)
	*/	
	void SetFvThreshold(const double fv_threshold);


protected:

	/**
	*@brief Constant variable depending on user-defined parameters only
	*/
	void Precalculations();

	/**
	*@brief Fatal error message
	*/
	void FatalErrorMessage(const std::string message);


private:

	// Nucleation parameters
	unsigned int n0_;		// number of TiO2 molecules in a monomer (default: 5)
	double epsilon_;		// nucleation collision enhancement factor (-)
	double na_;			// nucleation temperature exponent (-)
	double Ta_;			// nucleation activation temperature (K)

	// Sintering parameters	
	double As_;			// sintering frequency factor (s,K)
	double ns_;			// sintering temperature exponent (-)
	double Ts_;			// sintering activation temperature (K)

	// Activation of processes
	bool is_nucleation_;		// if true, nucleation is turned on
	bool is_coagulation_;		// if true, coagulation is turned on
	bool is_sintering_;		// if true, sintering is turned on

	// Numerical parameters
	double N_threshold_;		// N threshold for calculation of coagulation/sintering source terms (#/m3)
	double fv_threshold_;		// fv threshold for calculation of properties (-)

	// Constant variable depending on user-defined parameters only
	double v0_;			// monomer volume (m3)
	double s0_;			// monomer surface (m2)
	double d0_;			// monomer diameter (m)
	double alpha_;			// kinetic coefficient (m3/s/K^(1/2))


private:

	static const double pi_;	// PI
	static const double kB_;	// Boltzmann's constant (m2kg/s2/K)
	static const double Nav_mol_;	// Avogadro's constant (#/mol)
	static const double Nav_kmol_;	// Avogadro's constant (#/kmol)

	static const double W_TiO2_;	// TiO2 molecular weight (kg/kmol or g/mol)
	static const double rhos_;	// solid density (kg/m3) 
	static const double m_TiO2_;	// mass of TiO2 molecule (kg)
	static const double d_TiO2_;	// diameter of TiO2 molecule (m)

	static const double m_TiOH4_;	// mass of TiOH4 molecule (kg)
	static const double d_TiOH4_;	// diameter of TiOH4 molecule (m)
};
