/*-----------------------------------------------------------------------*\
|                                                                         |
|   ╭╮╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╭━━━┳━╮╭━┳━━━┳╮╭━┳━━━╮                               |
|   ┃┃╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱┃╭━╮┃┃╰╯┃┃╭━╮┃┃┃╭┫╭━━╯                               |
|   ┃┃╭━━┳╮╭┳┳━╮╭━━┳━┫╰━━┫╭╮╭╮┃┃╱┃┃╰╯╯┃╰━━┳╮╱╭╮                           |
|   ┃┃┃╭╮┃╰╯┣┫╭╮┫╭╮┃╭┻━━╮┃┃┃┃┃┃┃╱┃┃╭╮┃┃╭━┳╯╰┳╯╰╮                          |
|   ┃╰┫╭╮┃┃┃┃┃┃┃┃╭╮┃┃┃╰━╯┃┃┃┃┃┃╰━╯┃┃┃╰┫╰━┻╮╭┻╮╭╯                          |
|   ╰━┻╯╰┻┻┻┻┻╯╰┻╯╰┻╯╰━━━┻╯╰╯╰┻━━━┻╯╰━┻━━━┻╯╱╰╯                           |
|                                                                         |
|   Authors: Alberto Cuoci                                                |
|                                                                         |
|   Contacts: Alberto Cuoci                                               |
|   email: alberto.cuoci@polimi.it                                        |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano (Italy)                      |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of laminarSMOKE++ solver.                           |
|                                                                         |
|   License                                                               |
|                                                                         |
|   Copyright(C) 2022 Alberto Cuoci                                       |
|   laminarSMOKE++ is free software: you can redistribute it and/or       |
|   modify it under the terms of the GNU General Public License           |
|   as published by the Free Software Foundation, either version 3 of     |
|   the License, or (at your option) any later version.                   |
|                                                                         |
|   laminarSMOKE++ is distributed in the hope that it will be useful,     |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with laminarSMOKE++.                                            |
|   If not, see <http://www.gnu.org/licenses/>.                           |
|                                                                         |
\*-----------------------------------------------------------------------*/

#ifndef BatchReactorHomogeneousConstantPressure_NOxPP_H
#define	BatchReactorHomogeneousConstantPressure_NOxPP_H

#include "ModelNOxPostProcessor.H"

class BatchReactorHomogeneousConstantPressure_NOxPP
{
public:

	//- Default constructor
	BatchReactorHomogeneousConstantPressure_NOxPP(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, 
							OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
							ModelNOxPostProcessor& NOxPP);

	//- Set reactor
	void SetReactor( const double T0, const double P0, const double rho0);

	//- Set reactor mass fractions
	void SetMassFractions( const Eigen::VectorXd& omega );

	//- Return the total number of equations
	unsigned int NumberOfEquations() const;

	//- Returns the differential equations
	int Equations(const double t, const Eigen::VectorXd& y, Eigen::VectorXd& dy);

	//- Print function
	int Print(const double t, const Eigen::VectorXd& y);

	//- Integrate NOxPP equations (only in case the soot volume fraction is sufficiently large)
	void IntegrateNOxPP(const bool flag) { integrate_NOxPP_ = flag; }

	//- Set the NOxPP options
	void SetNOxPPOptions(const int index_OH, const int index_O, const int index_CH, const int index_Fuel);

	//- Enable the debug mode
	void SetDebug(const int flag) { __debug__ = flag; }

private:

	OpenSMOKE::ThermodynamicsMap_CHEMKIN&	thermodynamicsMap_;		//!< thermodynamic map
	OpenSMOKE::KineticsMap_CHEMKIN&		kineticsMap_;			//!< kinetic map

	double T0_;			//!< temperature (in K)
	double P0_;			//!< pressure (in Pa)
	double rho0_;			//!< density (in kg/m3)
	double MW0_;			//!< molecular weight (in kg/kmol)
	double cTot0_;			//!< total concentration (in kmol/m3)

	unsigned int NC_;		//!< number of species
	unsigned int NR_;		//!< number of reactions

	double YN2_;
	double YO2_;
	double YH2O_;
	double YFuel_;
	double YO_;
	double YOH_;

	bool debug_;			//!< if true, additional debug (time-consuming) operations are carried out

	int index_OH_;			//!< index of OH species (0-based)
	int index_O_;			//!< index of O species (0-based)
	int index_CH_;			//!< index of CH species (0-based)
	int index_Fuel_;		//!< index of Fuel species (0-based)

	bool integrate_NOxPP_;			//!< if true, the chemical step is solved

	ModelNOxPostProcessor& NOxPP_;			//!< NOxPP model

	int __debug__;				//!< debug mode (if __debug__>= 0)	

private:

	static const double MW_NO_;	// NO molecular weight (kg/kmol)
	static const double MW_HCN_;	// HCN molecular weight (kg/kmol)
	static const double MW_NH3_;	// NH3 molecular weight (kg/kmol)
	static const double MW_N2O_;	// N2O molecular weight (kg/kmol)

	static const double MW_N2_;	// N2 molecular weight (kg/kmol)
	static const double MW_O2_;	// O2 molecular weight (kg/kmol)
	static const double MW_H2O_;	// H2O molecular weight (kg/kmol)
	static const double MW_O_;	// O molecular weight (kg/kmol)
	static const double MW_OH_;	// OH molecular weight (kg/kmol)
		     double MW_Fuel_;	// FUEL molecular weight (kg/kmol)
		
};

const double BatchReactorHomogeneousConstantPressure_NOxPP::MW_NO_ = 30.006000;		// molecular weight of NO (kg/kmol)
const double BatchReactorHomogeneousConstantPressure_NOxPP::MW_HCN_ = 27.026000;	// molecular weight of HCN (kg/kmol)
const double BatchReactorHomogeneousConstantPressure_NOxPP::MW_NH3_ = 17.031000;	// molecular weight of NH3 (kg/kmol)
const double BatchReactorHomogeneousConstantPressure_NOxPP::MW_N2O_ = 44.013000;	// molecular weight of N2O (kg/kmol)

const double BatchReactorHomogeneousConstantPressure_NOxPP::MW_N2_ = 28.014000;		// molecular weight of N2O (kg/kmol)
const double BatchReactorHomogeneousConstantPressure_NOxPP::MW_O2_ = 31.997999;		// molecular weight of N2O (kg/kmol)
const double BatchReactorHomogeneousConstantPressure_NOxPP::MW_H2O_ = 18.015000;	// molecular weight of N2O (kg/kmol)
const double BatchReactorHomogeneousConstantPressure_NOxPP::MW_O_ = 15.999000;		// molecular weight of N2O (kg/kmol)
const double BatchReactorHomogeneousConstantPressure_NOxPP::MW_OH_ = 17.007000;		// molecular weight of N2O (kg/kmol)


unsigned int BatchReactorHomogeneousConstantPressure_NOxPP::NumberOfEquations() const
{
	return 1;
}

void BatchReactorHomogeneousConstantPressure_NOxPP::SetNOxPPOptions(const int index_OH, const int index_O, const int index_CH, const int index_Fuel)
{
	index_OH_    	= index_OH;
	index_O_	= index_O;
	index_CH_    	= index_CH;
	index_Fuel_    	= index_Fuel;	

	MW_Fuel_ = thermodynamicsMap_.MWs()[index_Fuel_];
}

BatchReactorHomogeneousConstantPressure_NOxPP::BatchReactorHomogeneousConstantPressure_NOxPP
(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap, ModelNOxPostProcessor& NOxPP ):
	thermodynamicsMap_(thermodynamicsMap), 
	kineticsMap_(kineticsMap),
	NOxPP_(NOxPP)
{
	NC_ = thermodynamicsMap_.NumberOfSpecies();
	NR_ = kineticsMap_.NumberOfReactions();

	debug_ = false;
	integrate_NOxPP_ = true;

	__debug__ = -1;
}

void BatchReactorHomogeneousConstantPressure_NOxPP::SetReactor( const double T0, const double P0, const double rho0 )
{
	T0_    = T0;					// (K)
	P0_    = P0;					// (Pa)
	rho0_  = rho0;					// (kg/m3)
	cTot0_ = P0/PhysicalConstants::R_J_kmol/T0_;	// (kmol/m3)
	MW0_ = rho0_/cTot0_;				// (kg/kmol)
}

void BatchReactorHomogeneousConstantPressure_NOxPP::SetMassFractions( const Eigen::VectorXd& omega )
{
	YN2_ = omega(0);
	YO2_ = omega(1);
	YH2O_ = omega(2);
	YFuel_ = omega(3);

	YO_ = omega(4);
	YOH_ = omega(5);
}

int BatchReactorHomogeneousConstantPressure_NOxPP::Equations(const double t, const Eigen::VectorXd& y, Eigen::VectorXd& dy)
{
	// Recover NOxPP variables	
	const double YNO = y(0); 	// mass fraction
//	const double YHCN = y(1); 	// mass fraction
//	const double YNH3 = y(2);	// mass fraction
//	const double YN2O = y(3);	// mass fraction

	// Concentration of relevant variables (kmol/m3)
	const double cNO = cTot0_*(YNO*MW0_/MW_NO_);
//	const double cHCN = cTot0_*(YHCN*MW0_/MW_NHCN_);
//	const double cNH3 = cTot0_*(YNH3*MW0_/MW_NH3_);
//	const double cN2O = cTot0_*(YN2O*MW0_/MW_N2O_);

	const double cFuel = cTot0_*(YFuel_*MW0_/MW_Fuel_);
	const double cN2 = cTot0_*(YN2_*MW0_/MW_N2_);
	const double cO2 = cTot0_*(YO2_*MW0_/MW_O2_);
	const double cH2O = cTot0_*(YH2O_*MW0_/MW_H2O_);
	const double cO = cTot0_*(YO_*MW0_/MW_O_);
	const double cOH = cTot0_*(YOH_*MW0_/MW_OH_);


	// Debug info on the screen
	if (__debug__ != -1)
	{
		std::cout << "Cell: " << __debug__ << " T: " << T0_ << " P: " << P0_ << " MW: " << MW0_ << " rho: " << rho0_ << " cTot: " << cTot0_ << std::endl;
	//	std::cout << " + Y: " << YNO << " " << YHCN << " " << YNH3 << " " << YN2O << std::endl;
		std::cout << " + Y: " << YNO << std::endl;
		std::cout << " + Y: " << YN2_ << " " << YO2_ << " " << YH2O_ << " " << YO_ << " " << YOH_ << " " << YFuel_ << std::endl;
		std::cout << " + c: " << cN2 << " " << cO2 << " " << cH2O << " " << cO << " " << cOH << " " << cNO << " " << cFuel << " " << MW_Fuel_ << std::endl;
	}

	// NOxPP Equations
	dy.setZero();
	if (integrate_NOxPP_ == true)
	{
		// Thermal NOx
		double OmegaNOThermal = 0.;
		NOxPP_.ThermalNOx( T0_, cN2, cO2, cH2O, cO, cOH, cNO, OmegaNOThermal);

		// Prompt NOx
		double OmegaNOPrompt = 0.;
		NOxPP_.PromptNOx( T0_, P0_, cN2, cO2, cFuel, OmegaNOPrompt);

		// Intermediate N2O
		double OmegaNOIntermediateN2O = 0.;
		NOxPP_.IntermediateN2O( T0_, P0_, cN2, cO2, cO, cNO, OmegaNOIntermediateN2O);
		
		dy(0) = (OmegaNOThermal+OmegaNOPrompt+OmegaNOIntermediateN2O)/rho0_;
	//	dy(1) = 0.;	
	//	dy(2) = 0.;	
	//	dy(3) = 0.;	

		// Debug info on the screen
		if (__debug__ != -1)
		{
			std::cout << " + Omega: " << OmegaNOThermal << " " << OmegaNOPrompt << " " << OmegaNOIntermediateN2O << std::endl;	
		}
	}

	return 0;
}

int BatchReactorHomogeneousConstantPressure_NOxPP::Print(const double t, const Eigen::VectorXd& y)
{
	return 0;
}

#endif // BatchReactorHomogeneousConstantPressure_NOxPP_H
