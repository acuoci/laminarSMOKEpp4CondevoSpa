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

class ModelNOxPostProcessor 
{

public:

	enum O_Radical_Model { O_RADICAL_EQUILIBRIUM, O_RADICAL_PARTIAL_EQUILIBRIUM, O_RADICAL_KINETICS } ;
	enum OH_Radical_Model { OH_RADICAL_EXCLUSION, OH_RADICAL_PARTIAL_EQUILIBRIUM, OH_RADICAL_KINETICS } ;

public:

	/**
	*@brief Default constructor
	*/
	ModelNOxPostProcessor();


	/**
	*@brief Summary on screen
	*/
	void SummaryOnScreen();


	/**
	*@brief Formation rate of NO from the thermal NOx path
	*@param T temperature (K)
	*@param concN2 concentration of N2 (kmol/m3)
	*@param concO2 concentration of O2 (kmol/m3)
	*@param concH2O concentration of H2O (kmol/m3)
	*@param concO concentration of O (kmol/m3)
	*@param concOH concentration of OH (kmol/m3)
	*@param concNO concentration of NO (kmol/m3)
	*@param OmegaNO formation rate (kg/m3/s)
	*/
	void ThermalNOx(	const double T,
				const double concN2, const double concO2, const double concH2O, const double concO, const double concOH, const double concNO, 
				double& OmegaNO	);

	/**
	*@brief Formation rate of NO from the prompt NOx path
	*@param T temperature (K)
	*@param P pressure (Pa)
	*@param concN2 concentration of N2 (kmol/m3)
	*@param concO2 concentration of O2 (kmol/m3)
	*@param concFuel concentration of FUEL (kmol/m3)
	*@param OmegaNO formation rate (kg/m3/s)
	*/
	void PromptNOx(	const double T, const double P, 
			const double concN2, const double concO2, const double concFuel, 
			double& OmegaNO );

	/**
	*@brief Formation rate of NO from the intermediate N@O path
	*@param T temperature (K)
	*@param P pressure (Pa)
	*@param concN2 concentration of N2 (kmol/m3)
	*@param concO2 concentration of O2 (kmol/m3)
	*@param concFuel concentration of O (kmol/m3)
	*@param concFuel concentration of NO (kmol/m3)
	*@param OmegaNO formation rate (kg/m3/s)
	*/
	void IntermediateN2O(	const double T, const double P,
				const double concN2, const double concO2, const double concO, const double concNO, 
				double& OmegaNO	);


	/**
	*@brief Turns on the thermal NOx path
	*@param flag if true, the thermal NOx path is turned on
	*/
	void SetThermalNOx(const bool flag);

	/**
	*@brief Turns on the prompt NOx path
	*@param flag if true, the prompt NOx path is turned on
	*/
	void SetPromptNOx(const bool flag);

	/**
	*@brief Turns on the intermediate N2O path
	*@param flag if true, the intermediate N2O path is turned on
	*/
	void SetIntermediateN2O(const bool flag);

	/**
	*@brief Selects the model to estiamte the O radical concentration
	*@param model the O radical model
	*/
	void SetORadicalMode(const ModelNOxPostProcessor::O_Radical_Model model);

	/**
	*@brief Selects the model to estiamte the OH radical concentration
	*@param model the OH radical model
	*/
	void SetOHRadicalMode(const ModelNOxPostProcessor::OH_Radical_Model model);


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

	// Sub-models
	bool is_thermal_;			//!< thermal NOx on/off
	bool is_prompt_;			//!< prompt NOx on/off
	bool is_fuel_nox_;			//!< fuel NOx on/off
	bool is_intermediate_n2o_;		//!< formation of NOx from intermediate N2O
	bool is_reduction_by_reburning_;	//!< reduction of NOx by reburning
	bool is_reduction_by_sncr_;		//!< reduction of NOx by SNCR	


	// Radical models
	O_Radical_Model		o_radical_model_;
	OH_Radical_Model	oh_radical_model_;


	// Thermal NOx

	// Reaction: O+N2 = N+NO (extended Zeldovich mechanism, NOx from molecular nitrogen)
	double Af1_;		// (m, mol, s)
	double nf1_;		// (-)
	double Tf1_;		// (K)
	double Ar1_;		// (m, mol, s)
	double nr1_;		// (-)
	double Tr1_;		// (K)

	// Reaction: N+O2 = O+NO (extended Zeldovich mechanism, NOx from molecular nitrogen)
	double Af2_ ;		// (m, mol, s)
	double nf2_;		// (-)
	double Tf2_;		// (K)
	double Ar2_;		// (m, mol, s)
	double nr2_;		// (-)
	double Tr2_;		// (K)

	// Reaction: N+OH = H+NO (important near stoichiometric conditions and in fuel-rich mixtures)
	double Af3_;		// (m, mol, s)
	double nf3_;		// (-)
	double Tf3_;		// (K)
	double Ar3_;		// (m, mol, s)
	double nr3_;		// (-)
	double Tr3_;		// (K)

	// O radical (equilibrium approach)
	double AOe_ ;		// (m,mol)
	double nOe_;		// (-)
	double TOe_;		// (K)

	// O radical (partial equilibrium approach, O2+M = O+O+M)
	double AOpe_;		// (m,mol)
	double nOpe_;		// (-)
	double TOpe_;		// (K)

	// OH radical (partial equilibrium approach)
	double AOHpe_;		// (m,mol)
	double nOHpe_;		// (-)
	double TOHpe_;		// (K)


	// Prompt mechanism
	double APrompt_ ; 	// (1/s)		
	double APrimePrompt_;	// (1/s)
	double EPrompt_;	// (cal/mol)
	double EPrimePrompt_;	// (cal/mol)
	bool prompt_fuel_correction_factor_;
	double nFuel_;
	double phi_;

	// Intermediate N2O  mechanism

	// Reaction: N+O+M = N2O+M
	double An2of1_;		// (m, mol, s)
	double nn2of1_;		// (-)
	double Tn2of1_;		// (K)
	double An2or1_;		// (m, mol, s)
	double nn2or1_;		// (-)
	double Tn2or1_;		// (K)

	// Reaction: N2O+O=2NO
	double An2of2_ ;	// (m, mol, s)
	double nn2of2_;		// (-)
	double Tn2of2_;		// (K)
	double An2or2_;		// (m, mol, s)
	double nn2or2_;		// (-)
	double Tn2or2_;		// (K)


private:

	static const double R_J_kmol_;	// R (J/kmol)
	static const double R_cal_mol_;	// R (cal/mol)
	static const double MW_NO_;	// NO molecular weight (kg/kmol)
	static const double MW_HCN_;	// HCN molecular weight (kg/kmol)
	static const double MW_NH3_;	// NH3 molecular weight (kg/kmol)
	static const double MW_N2O_;	// N2O molecular weight (kg/kmol)
};
