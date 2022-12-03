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

#include "ModelNOxPostProcessor.H"

#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <numeric>
#include <iomanip>

const double ModelNOxPostProcessor::R_J_kmol_ = 8314.4621;	// R (J/kmol)
const double ModelNOxPostProcessor::R_cal_mol_ = 1.987000;	// R (cal/mol)
const double ModelNOxPostProcessor::MW_NO_ = 30.006000;		// molecular weight of NO (kg/kmol)
const double ModelNOxPostProcessor::MW_HCN_ = 27.026000;	// molecular weight of HCN (kg/kmol)
const double ModelNOxPostProcessor::MW_NH3_ = 17.031000;	// molecular weight of NH3 (kg/kmol)
const double ModelNOxPostProcessor::MW_N2O_ = 44.013000;	// molecular weight of N2O (kg/kmol)

ModelNOxPostProcessor::ModelNOxPostProcessor()
{
	// User-defined parameters 
	is_thermal_ = true;			//!< thermal NOx on/off
	is_prompt_ = true;			//!< prompt NOx on/off
	is_intermediate_n2o_ = true;		//!< formation of NOx from intermediate N2O

	is_fuel_nox_ = false;			//!< fuel NOx on/off
	is_reduction_by_reburning_ = false;	//!< reduction of NOx by reburning
	is_reduction_by_sncr_ = false;		//!< reduction of NOx by SNCR


	// Radical models
	o_radical_model_ = ModelNOxPostProcessor::O_RADICAL_PARTIAL_EQUILIBRIUM;
	oh_radical_model_ = ModelNOxPostProcessor::OH_RADICAL_PARTIAL_EQUILIBRIUM;	


	// Thermal NOx parameters
	// Hanson and Salimian, "Survey of rate constants in H/N/O systems", in "Combustion Chemistry", 361 (1984)
	
	// Reaction: O+N2 = N+NO (extended Zeldovich mechanism, NOx from molecular nitrogen)
	Af1_ = 1.8e8;		// (m, mol, s)
	nf1_ = 0.;		// (-)
	Tf1_ = 38370.;		// (K)
	Ar1_ = 3.80e7;		// (m, mol, s)
	nr1_ = 0.;		// (-)
	Tr1_ = 425.;		// (K)

	// Reaction: N+O2 = O+NO (extended Zeldovich mechanism, NOx from molecular nitrogen)
	Af2_ = 1.8e4;		// (m, mol, s)
	nf2_ = 1.;		// (-)
	Tf2_ = 4680.;		// (K)
	Ar2_ = 3.81e3;		// (m, mol, s)
	nr2_ = 1.;		// (-)
	Tr2_ = 20820.;		// (K)

	// Reaction: N+OH = H+NO (important near stoichiometric conditions and in fuel-rich mixtures)
	Af3_ = 7.1e7;		// (m, mol, s)
	nf3_ = 0.;		// (-)
	Tf3_ = 450.;		// (K)
	Ar3_ = 1.70e8;		// (m, mol, s)
	nr3_ = 0.;		// (-)
	Tr3_ = 24560.;		// (K)

	// O radical (equilibrium approach)
	// Westenberg, Cpmbustion Science and Technology, 4, 59 (1971)
	AOe_ = 3.97e5;		// (m, mol,K)
	nOe_ = 0.50;		// (-)
	TOe_ = 31090;		// (K)

	// O radical (partial equilibrium approach, O2+M = O+O+M)
	// Warnatz, "NOx formation in high temperature processes", University of Stuttgart, Germany
	AOpe_ = 36.64;		// (m, mol, K)
	nOpe_ = 0.50;		// (-)
	TOpe_ = 27123.;		// (K)

	// OH radical (partial equilibrium approach)
	// Baulch et al., "Evaluated kinetic data for combustion modeling", J. Physical and Chemical Reference Data, 21(3) (1992)
	// Westbrook and Dryer, "Chemical kinetic modelling of hydrocarbon combustion", Progress in Energy Combustion Science 1, (1984)
	AOHpe_ = 2.129e2;	// (m, mol, K)
	nOHpe_ = -0.57;		// (-)
	TOHpe_ = 4595.;		// (K)


	// Prompt NOx
	APrompt_ = 1.2e7;			// 1/s
	APrimePrompt_ = 6.4e6;			// 1/s
	EPrompt_ = 60000.;			// cal/mol
	EPrimePrompt_ = 72500.;			// cal/mol
	prompt_fuel_correction_factor_ = false;
	nFuel_ = 2.0;
	phi_ = 1.0;


	// Intermediate N2O
	// Melte and Pratt, "Measurement of atomic oxygen and nitrogen oxides in jet stirred combustion", 15th Symposium on Combustion, 1061-1070 (1974)
	// Reaction: N+O+M = N2O+M
	An2of1_ = 4.44e32;		// (m, mol, s)
	nn2of1_ = -8.358;		// (-)
	Tn2of1_ = 28234.;		// (K)
	An2or1_ = 4.00e8;		// (m, mol, s)
	nn2or1_ = 0.;			// (-)
	Tn2or1_ = 28234.;		// (K)

	// Reaction: N2O+O=2NO
	An2of2_ = 2.90e7;		// (m, mol, s)
	nn2of2_ = 0.;			// (-)
	Tn2of2_ = 11651.;		// (K)
	An2or2_ = 1.453e-29;		// (m, mol, s)
	nn2or2_ = 9.259;		// (-)
	Tn2or2_ = 11651.;		// (K)

	// Precalculation of constant variables dependent on user defined parameters
	Precalculations();
}


void ModelNOxPostProcessor::Precalculations()
{

}

void ModelNOxPostProcessor::SetThermalNOx(const bool flag)
{
	is_thermal_ = flag;
	Precalculations();
}

void ModelNOxPostProcessor::SetPromptNOx(const bool flag)
{
	is_prompt_ = flag;
	Precalculations();
}

void ModelNOxPostProcessor::SetIntermediateN2O(const bool flag)
{
	is_intermediate_n2o_ = flag;
	Precalculations();
}

void ModelNOxPostProcessor::SetORadicalMode(const ModelNOxPostProcessor::O_Radical_Model model)
{
	o_radical_model_ = model;
	Precalculations();
}

void ModelNOxPostProcessor::SetOHRadicalMode(const ModelNOxPostProcessor::OH_Radical_Model model)
{
	oh_radical_model_ = model;
	Precalculations();
}

void ModelNOxPostProcessor::ThermalNOx(	const double T, const double concN2, const double concO2, const double concH2O, const double concO, const double concOH,
					const double concNO, 
					double& SNO)
{
	const double cN2 = std::max(concN2, 1.e-16) * 1000.;	// (mol/m3)
	const double cO2 = std::max(concO2, 1.e-16) * 1000.;	// (mol/m3)
	const double cH2O = std::max(concH2O, 0.) * 1000.;	// (mol/m3)
	const double cNO = std::max(concNO, 0.) * 1000.;	// (mol/m3)

	double cO = std::max(concO, 0.) * 1000.;	// (mol/m3)
	double cOH = std::max(concOH, 0.) * 1000.;	// (mol/m3)


	SNO = 0.; // (kg/m3/s)

	if (is_thermal_ == true)
	{		
		const double kf1 = Af1_*std::pow(T, nf1_)*std::exp(-Tf1_/T);	// (m3/mol/s)
		const double kr1 = Ar1_*std::pow(T, nr1_)*std::exp(-Tr1_/T);	// (m3/mol/s)

		const double kf2 = Af2_*std::pow(T, nf2_)*std::exp(-Tf2_/T);	// (m3/mol/s)
		const double kr2 = Ar2_*std::pow(T, nr2_)*std::exp(-Tr2_/T);	// (m3/mol/s)

		const double kf3 = Af3_*std::pow(T, nf3_)*std::exp(-Tf3_/T);	// (m3/mol/s)
		//const double kr3 = Ar3_*std::pow(T, nr3_)*std::exp(-Tr3_/T);	// (m3/mol/s)

		// O radical concentration
		if (o_radical_model_ == ModelNOxPostProcessor::O_RADICAL_EQUILIBRIUM)
		{
			const double kOe = AOe_*std::pow(T,nOe_)*std::exp(-TOe_/T);	// (mol/m3)^(1/2)
			cO = kOe*std::sqrt(cO2);					// (mol/m3)
		}
		else if (o_radical_model_ == ModelNOxPostProcessor::O_RADICAL_PARTIAL_EQUILIBRIUM)
		{
			const double kOpe = AOpe_*std::pow(T,nOpe_)*std::exp(-TOpe_/T);	// (mol/m3)^(1/2)
			cO = kOpe*std::sqrt(cO2);					// (mol/m3)
		}
		else if (o_radical_model_ == ModelNOxPostProcessor::O_RADICAL_KINETICS)
		{
			cO = cO;	// (mol/m3)
		}

		// OH radical concentration
		if (oh_radical_model_ == ModelNOxPostProcessor::OH_RADICAL_EXCLUSION)
		{
			cOH = 0.;	// (mol/m3)
		}
		else if (oh_radical_model_ == ModelNOxPostProcessor::OH_RADICAL_PARTIAL_EQUILIBRIUM)
		{
			const double kOHpe = AOHpe_*std::pow(T,nOHpe_)*std::exp(-TOHpe_/T);	// (-)
			cOH = kOHpe*std::sqrt(cO)*std::sqrt(cH2O);				// (mol/m3)
		}
		else if (oh_radical_model_ == ModelNOxPostProcessor::OH_RADICAL_KINETICS)
		{
			cOH = cOH;	// (mol/m3)
		}

		// Source term
		SNO = 2.*kf1*cO*cN2*(1.-(kr1*kr2*cNO*cNO)/(kf1*cN2*kf2*cO2))/(1.+(kr1*cNO)/(kf2*cO2+kf3*cOH));	// (mol/m3/s)
		SNO /= 1000.;	// (kmol/m3/s)
		SNO *= MW_NO_;	// (kg/m3/s)
	}
}

void ModelNOxPostProcessor::PromptNOx(	const double T, const double P, 
					const double concN2, const double concO2, const double concFuel, 
					double& SNO	)
{
	const double cN2 = std::max(concN2, 0.);		// (kmol/m3)
	const double cO2 = std::max(concO2, 1.e-16);		// (kmol/m3)
	const double cFuel = std::max(concFuel, 0.);		// (kmol/m3)
	const double cTot = P/R_J_kmol_/T;			// (kmol/m3)
	const double xO2 = cO2/cTot;				// (-)

	
	SNO = 0.; // (kg/m3/s)
	if (is_prompt_ == true)
	{		
		// Oxygen reaction order
		// De Soete, "Overall reaction rates of NO and N2 formation from fuel nitrogen", 15th Symposium International on Combustion, 1093-1102 (1975)
		double aPrompt = 1.;
		     if (xO2<=4.09e-3)			aPrompt = 1.;
		else if (xO2>=4.09e-3 && xO2<=1.11e-2)	aPrompt = -3.95-0.9*std::log(xO2);
		else if (xO2>=1.11e-2 && xO2<=0.03)	aPrompt = -0.25-0.1*std::log(xO2);
		else if (xO2>=0.03)			aPrompt = 0.;
	
		// Correction factor 
		if (prompt_fuel_correction_factor_ == true)
		{
			// Backmier, Eberius and Jus, Combustion Science and Technology 7, 77 (1973)
			double f = 4.75+0.0819*nFuel_-23.2*phi_+32.*(phi_*phi_)-12.2*(phi_*phi_*phi_);
			if (phi_<=0.60)		f = 4.75+0.0819*nFuel_-23.2*0.60+32.*(0.60*0.60)-12.2*(0.60*0.60*0.60);
			else if (phi_>=1.60)	f = 4.75+0.0819*nFuel_-23.2*1.60+32.*(1.60*1.60)-12.2*(1.60*1.60*1.60);
			
			// Kinetic constant
			const double kPrimePrompt = APrimePrompt_ * std::pow(1./cTot,aPrompt+1.);	// (m3/kmol)^(a+1)
		
			// Source term
			SNO = f*kPrimePrompt*std::exp(-EPrimePrompt_/R_cal_mol_/T) * std::pow(cO2,aPrompt)*cN2*cFuel;	// (kmol/m3/s)
			SNO *= MW_NO_;	// (kg/m3/s)
		}
		else
		{
			// Kinetic constant
			const double kPrompt = APrompt_ * std::pow(1./cTot,aPrompt+1.);	// (m3/kmol)^(a+1)
		
			// Source term
			SNO = kPrompt*std::exp(-EPrompt_/R_cal_mol_/T) * std::pow(cO2,aPrompt)*cN2*cFuel;	// (kmol/m3/s)
			SNO *= MW_NO_;	// (kg/m3/s)
		}
	}
}

void ModelNOxPostProcessor::IntermediateN2O(	const double T, const double P,
						const double concN2, const double concO2, const double concO, const double concNO, 
						double& SNO	)
{
	const double cN2 = std::max(concN2, 0.) * 1000.;	// (mol/m3)
	const double cO2 = std::max(concO2, 0.) * 1000.;	// (mol/m3)
	const double cNO = std::max(concNO, 0.) * 1000.;	// (mol/m3)
	const double cM = P/R_J_kmol_/T * 1000.;		// (mol/m3)
	double cO = std::max(concO, 0.) * 1000.;		// (mol/m3)

	// Source term
	SNO = 0.; // (kg/m3/s)

	if (is_intermediate_n2o_ == true)
	{		
		const double kf1 = An2of1_*std::pow(T, nn2of1_)*std::exp(-Tn2of1_/T);	// (m3/mol/s)
		const double kr1 = An2or1_*std::pow(T, nn2or1_)*std::exp(-Tn2or1_/T);	// (m3/mol/s)

		const double kf2 = An2of2_*std::pow(T, nn2of2_)*std::exp(-Tn2of2_/T);	// (m3/mol/s)
		const double kr2 = An2or2_*std::pow(T, nn2or2_)*std::exp(-Tn2or2_/T);	// (m3/mol/s)

		// O radical concentration
		if (o_radical_model_ == ModelNOxPostProcessor::O_RADICAL_EQUILIBRIUM)
		{
			const double kOe = AOe_*std::pow(T,nOe_)*std::exp(-TOe_/T);	// (mol/m3)^(1/2)
			cO = kOe*std::sqrt(cO2);					// (mol/m3)
		}
		else if (o_radical_model_ == ModelNOxPostProcessor::O_RADICAL_PARTIAL_EQUILIBRIUM)
		{
			const double kOpe = AOpe_*std::pow(T,nOpe_)*std::exp(-TOpe_/T);	// (mol/m3)^(1/2)
			cO = kOpe*std::sqrt(cO2);					// (mol/m3)
		}
		else if (o_radical_model_ == ModelNOxPostProcessor::O_RADICAL_KINETICS)
		{
			cO = cO;	// (mol/m3)
		}

		// N2O concentration (quasi-steady-state assumption)
		const double cN2O = (kf1*cN2*cO*cM+kr2*cNO*cNO)/(kr1*cM+kf2*cO);	// (mol/m3/s)

		// Source term
		SNO = 2.*(kf2*cN2O*cO-kr2*cNO*cNO);	// (mol/m3/s)
		SNO /= 1000.;	// (kmol/m3/s)
		SNO *= MW_NO_;	// (kg/m3/s)
	}
}

void ModelNOxPostProcessor::SummaryOnScreen()
{	
	std::cout << std::endl;
	std::cout << "------------------------------------------------------------------------------------------" << std::endl; 
	std::cout << "                         NOx post-processor Model Summary                                 " << std::endl;
	std::cout << "------------------------------------------------------------------------------------------" << std::endl; 
	std::cout << " * Processes" << std::endl;
	std::cout << "    + Thermal NOx:            " << is_thermal_ << std::endl;
	std::cout << "    + Prompt NOx:             " << is_prompt_ << std::endl;
	std::cout << "    + Intermediate N2O:       " << is_intermediate_n2o_ << std::endl;
	std::cout << "    + Fuel NOx:               " << is_fuel_nox_ << std::endl;
	std::cout << "    + Reduction by reburning: " << is_reduction_by_reburning_ << std::endl;
	std::cout << "    + Reduction by SNCR:      " << is_reduction_by_sncr_ << std::endl;
	std::cout << std::endl;
	std::cout << " * Radical models" << std::endl;
	std::cout << "    + O:                      " << o_radical_model_ << std::endl;
	std::cout << "    + OH:                     " << oh_radical_model_ << std::endl;
	std::cout << std::endl;
	std::cout << " * Thermal NOx" << std::endl;
	std::cout << "    + Reaction 1: O+N2 = N+NO" << std::endl;
	std::cout << "       - Forward:             " << Af1_ << ", " << nf1_ << ", " << Tf1_ << std::endl;
	std::cout << "       - Backward:            " << Ar1_ << ", " << nr1_ << ", " << Tr1_ << std::endl;
	std::cout << "    + Reaction 2: N+O2 = O+NO" << std::endl;
	std::cout << "       - Forward:             " << Af2_ << ", " << nf2_ << ", " << Tf2_ << std::endl;
	std::cout << "       - Backward:            " << Ar2_ << ", " << nr2_ << ", " << Tr2_ << std::endl;
	std::cout << "    + Reaction 3: N+OH = H+NO" << std::endl;
	std::cout << "       - Forward:             " << Af3_ << ", " << nf3_ << ", " << Tf3_ << std::endl;
	std::cout << "       - Backward:            " << Ar3_ << ", " << nr3_ << ", " << Tr3_ << std::endl;
	std::cout << std::endl;
	std::cout << "------------------------------------------------------------------------------------------" << std::endl; 
	std::cout << std::endl;
}

void ModelNOxPostProcessor::FatalErrorMessage(const std::string message)
{
	std::cout << "ModelNOxPostProcessor: Fatal error" << std::endl;
	std::cout << "Error message: " << message << std::endl;
	std::cout << "Press enter to exit..." << std::endl;
	getchar();
	exit(-1);
}
