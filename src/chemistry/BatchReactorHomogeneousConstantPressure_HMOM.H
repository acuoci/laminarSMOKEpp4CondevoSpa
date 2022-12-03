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
|   Copyright(C) 2020, 2021 Alberto Cuoci                                 |
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

#ifndef BatchReactorHomogeneousConstantPressure_HMOM_H
#define	BatchReactorHomogeneousConstantPressure_HMOM_H

#include "utilities/soot/hmom/HMOM.h"

class BatchReactorHomogeneousConstantPressure_HMOM
{
public:

	//- Default constructor
	BatchReactorHomogeneousConstantPressure_HMOM(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, 
							OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
							OpenSMOKE::HMOM& hmom);

	//- Set reactor pressure (in Pa)
	void SetReactor( const double P0 );

	//- Set reactor mass fractions
	void SetMassFractions( const Eigen::VectorXd& omega );
	
	//- Enable the solution of energy equation
	void SetEnergyEquation(const bool flag) { energyEquation_ = flag; }
	
	//- Enable mass fraction checkings (default: false)
	void SetCheckMassFractions( const bool flag ) { checkMassFractions_ = flag; };

	//- Returns the reactor temperature
	double T() const;

	//- Returns the vector of formation rates of species (in kmol/m3/s)
	const Eigen::VectorXd& R() const { return R_; }

	//- Returns the reaction heat (in W/m3)
	double QR() const { return QR_; }

	//- Return the total number of equations
	unsigned int NumberOfEquations() const;

	//- Returns the differential equations
	int Equations(const double t, const Eigen::VectorXd& y, Eigen::VectorXd& dy);

	//- Print function
	int Print(const double t, const Eigen::VectorXd& y);

	//- Integrate HMOM equations (only in case the soot volume fraction is sufficiently large)
	void IntegrateHMOM(const bool flag) { integrate_hmom_ = flag; }

	//- Set the HMOM options
	void SetHMOMOptions(const int index_H, const int index_OH, const int index_H2, const int index_H2O, const int index_C2H2, const int index_O2, const std::vector<int> index_PAH);

	//- Enable the PAH consumption
	void SetPAHGasConsumption(const bool pah_gas_consumption);


private:

	OpenSMOKE::ThermodynamicsMap_CHEMKIN&	thermodynamicsMap_;		//!< thermodynamic map
	OpenSMOKE::KineticsMap_CHEMKIN&		kineticsMap_;			//!< kinetic map

	double T_;			//!< temperature (in K)
	double P0_;			//!< pressure (in Pa)
	
	double rho_;			//!< density (in kg/m3)
	double cTot_;			//!< total concentration (in kmol/m3)
	double MW_;			//!< molecular weight (in kg/kmol)
	double CpMixMass_;		//!< constant pressure specific heat (in J/kg/K)
	double QR_;			//!< heat release (in W/m3)

	unsigned int NC_;		//!< number of species
	unsigned int NR_;		//!< number of reactions
	unsigned int NHMOM_;		//!< number of moments 

	Eigen::VectorXd omega_;		//!< mass fractions
	Eigen::VectorXd x_;		//!< mole fractions
	Eigen::VectorXd c_;		//!< concentrations (in kmol/m3)
	Eigen::VectorXd R_;		//!< formation rates (in kmol/m3/s)

	bool checkMassFractions_;	//!< if true, mass fractions are normalized
	bool energyEquation_;		//!< if true, the energy equation is solved
	bool debug_;			//!< if true, additional debug (time-consuming) operations are carried out

	int hmom_index_H_;
	int hmom_index_OH_;
	int hmom_index_H2_;
	int hmom_index_H2O_;
	int hmom_index_C2H2_;
	int hmom_index_O2_;
	std::vector<int> hmom_index_PAH_;

	bool pah_gas_consumption_;
	bool integrate_hmom_;

	OpenSMOKE::HMOM& hmom_;


};

unsigned int BatchReactorHomogeneousConstantPressure_HMOM::NumberOfEquations() const
{
	return NC_ + 1 + NHMOM_;
}

void BatchReactorHomogeneousConstantPressure_HMOM::SetHMOMOptions(const int index_H, const int index_OH, const int index_H2, const int index_H2O, const int index_C2H2, const int index_O2, const std::vector<int> index_PAH)
{
	hmom_index_H_    = index_H;
	hmom_index_OH_   = index_OH;
	hmom_index_H2_   = index_H2;
	hmom_index_H2O_  = index_H2O;
	hmom_index_C2H2_ = index_C2H2;
	hmom_index_O2_   = index_O2;
	hmom_index_PAH_  = index_PAH;
}

void BatchReactorHomogeneousConstantPressure_HMOM::SetPAHGasConsumption(const bool pah_gas_consumption)
{
	pah_gas_consumption_ = pah_gas_consumption;
}

BatchReactorHomogeneousConstantPressure_HMOM::BatchReactorHomogeneousConstantPressure_HMOM
(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap, OpenSMOKE::HMOM& hmom ):
	thermodynamicsMap_(thermodynamicsMap), 
	kineticsMap_(kineticsMap),
	hmom_(hmom)
{
	NC_ = thermodynamicsMap_.NumberOfSpecies();
	NR_ = kineticsMap_.NumberOfReactions();
	NHMOM_ = 4;

	QR_ = 0.;

	omega_.resize(NC_); omega_.setZero();
	x_.resize(NC_); x_.setZero();
	c_.resize(NC_); c_.setZero();
	R_.resize(NC_); R_.setZero();

	checkMassFractions_ = false;
	energyEquation_ = true;
	debug_ = false;
	pah_gas_consumption_ = true;
	integrate_hmom_ = true;
}

void BatchReactorHomogeneousConstantPressure_HMOM::SetReactor( const double P0 )
{
	P0_    = P0;
}

void BatchReactorHomogeneousConstantPressure_HMOM::SetMassFractions( const Eigen::VectorXd& omega )
{
	omega_ = omega;
}

double BatchReactorHomogeneousConstantPressure_HMOM::T() const
{
	return T_;
}

int BatchReactorHomogeneousConstantPressure_HMOM::Equations(const double t, const Eigen::VectorXd& y, Eigen::VectorXd& dy)
{
	// Recover mass fractions
	if (checkMassFractions_ == true)
	{	
		for(unsigned int i=0;i<NC_;i++)
			omega_(i) = std::max(y(i), 0.);
	}
	else
	{
		for(unsigned int i=0;i<NC_;i++)
			omega_(i) = y(i);
	}

	// Recover temperature
	T_ = y(NC_); 

	// Recover normalized moments	
	const double HMOM_Norm_M00_ = y(NC_+1); 
	const double HMOM_Norm_M10_ = y(NC_+2); 
	const double HMOM_Norm_M01_ = y(NC_+3); 
	const double HMOM_Norm_N0_  = y(NC_+4);


	// Calculates the pressure and the concentrations of species
	thermodynamicsMap_.MoleFractions_From_MassFractions(x_.data(), MW_, omega_.data());
	cTot_ = P0_/PhysicalConstants::R_J_kmol/T_;
	c_ = cTot_*x_;
    	rho_ = cTot_*MW_;

	// Calculates thermodynamic properties
	thermodynamicsMap_.SetTemperature(T_);
	thermodynamicsMap_.SetPressure(P0_);

	// Calculates kinetics
	kineticsMap_.SetTemperature(T_);
	kineticsMap_.SetPressure(P0_);
	kineticsMap_.KineticConstants();
	kineticsMap_.ReactionRates(c_.data());
	kineticsMap_.FormationRates(R_.data());

	// Species equations
	for (unsigned int i=0;i<NC_;++i)	
		dy(i) = thermodynamicsMap_.MW(i)*R_(i)/rho_;

    	// Energy equation
    	dy(NC_) = 0.;     
    	if (energyEquation_ == true)
    	{
		CpMixMass_ = thermodynamicsMap_.cpMolar_Mixture_From_MoleFractions(x_.data()) / MW_;
		QR_ = kineticsMap_.HeatRelease(R_.data());
	
		dy(NC_)  = QR_ / (rho_*CpMixMass_);
	}

	// HMOM Equations
	if (integrate_hmom_ == true)
	{
		const double viscosity = 1.716e-5*Foam::pow(T_/273.15,1.5)*(273.15+110.4)/(T_+110.4);		// [kg/m/s]

		// Concentrations of PAHs
		Eigen::VectorXd conc_PAH(hmom_.pah_species().size());
		for (unsigned int j = 0; j<hmom_.pah_species().size(); j++)
			conc_PAH(j) = cTot_*x_(hmom_index_PAH_[j]);

		// Setup
		hmom_.SetNormalizedMoments(HMOM_Norm_M00_, HMOM_Norm_M10_, HMOM_Norm_M01_, HMOM_Norm_N0_);
		hmom_.SetTemperatureAndPressure(T_, P0_);
		hmom_.SetMassFractions( omega_(hmom_index_OH_), omega_(hmom_index_H_) );
		hmom_.SetConcentrations("kmol/m3", c_(hmom_index_OH_), c_(hmom_index_H_), c_(hmom_index_H2O_), c_(hmom_index_H2_), c_(hmom_index_C2H2_), c_(hmom_index_O2_), conc_PAH);
		hmom_.SetViscosity(viscosity);
		
		// Calculates
		hmom_.CalculateSourceMoments();

		// Moments equations [mol/m3/s]
		dy(NC_+1) = hmom_.SourceM00();
		dy(NC_+2) = hmom_.SourceM10();
		dy(NC_+3) = hmom_.SourceM01();
		dy(NC_+4) = hmom_.SourceN0();

		// PAH consumption
		if (pah_gas_consumption_ == true)
		{
			if (conc_PAH.sum() > 1.e-64)
			{

				double Omega_PAH = 0.;	// [kg/m3]
				for (unsigned int j = 0; j < hmom_.pah_species().size(); j++)
				{
					const double R_PAH = hmom_.PAHConsumptionRate(j) / 1000.;			// [kmol/m3/s]
					const double omega = thermodynamicsMap_.MW( hmom_index_PAH_[j] ) * R_PAH;	// [kg/m3/s]

					dy(hmom_index_PAH_[j]) -= omega/rho_;
					Omega_PAH += omega;
				}

				const int index_N2  = thermodynamicsMap_.IndexOfSpecies("N2") - 1;
				dy(index_N2)  += Omega_PAH/rho_;
			}
		}
	}
	else
	{
		dy(NC_+1) = 0.;
		dy(NC_+2) = 0.;
		dy(NC_+3) = 0.;
		dy(NC_+4) = 0.;
	}

	return 0;
}

int BatchReactorHomogeneousConstantPressure_HMOM::Print(const double t, const Eigen::VectorXd& y)
{
	return 0;
}

#endif // BatchReactorHomogeneousConstantPressure_HMOM_H
