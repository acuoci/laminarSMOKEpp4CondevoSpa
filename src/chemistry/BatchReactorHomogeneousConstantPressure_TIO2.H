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

#ifndef BatchReactorHomogeneousConstantPressure_TIO2_H
#define	BatchReactorHomogeneousConstantPressure_TIO2_H

#include "ModelTiO2.H"

class BatchReactorHomogeneousConstantPressure_TIO2
{
public:

	//- Default constructor
	BatchReactorHomogeneousConstantPressure_TIO2(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, 
							OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
							ModelTiO2& TIO2);

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

	//- Integrate TIO2 equations (only in case the soot volume fraction is sufficiently large)
	void IntegrateTIO2(const bool flag) { integrate_TIO2_ = flag; }

	//- Set the TIO2 options
	void SetTIO2Options(const int index_TiO2, const int index_TiOH4, const int index_H2O);

	//- Enable the gaseous species consumption because of nucleation
	void SetGaseousFeedback(const bool gaseous_feedback);

	//- Enable the debug mode
	void SetDebug(const int flag) { __debug__ = flag; }

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

	Eigen::VectorXd omega_;		//!< mass fractions
	Eigen::VectorXd x_;		//!< mole fractions
	Eigen::VectorXd c_;		//!< concentrations (in kmol/m3)
	Eigen::VectorXd R_;		//!< formation rates (in kmol/m3/s)

	bool checkMassFractions_;	//!< if true, mass fractions are normalized
	bool energyEquation_;		//!< if true, the energy equation is solved
	bool debug_;			//!< if true, additional debug (time-consuming) operations are carried out

	int index_TiO2_;			//!< index of TiO2 species (0-based)
	int index_TiOH4_;			//!< index of TiOH4 species (0-based)
	int index_H2O_;				//!< index of H2O species (0-based)

	bool gaseous_feedback_;			//!< if true, gaseous species are consumed because of nucleation
	bool integrate_TIO2_;			//!< if true, the chemical step is solved

	ModelTiO2& TIO2_;			//!< TiO2 model

	double scaling_factor_;			//!< scaling factor for N (default: 6.02214086e26 #/kmol)

	int __debug__;				//!< debug mode (if __debug__>= 0)			
};

unsigned int BatchReactorHomogeneousConstantPressure_TIO2::NumberOfEquations() const
{
	return (NC_+1) + 3;
}

void BatchReactorHomogeneousConstantPressure_TIO2::SetTIO2Options(const int index_TiO2, const int index_TiOH4, const int index_H2O)
{
	index_TiO2_    	= index_TiO2;
	index_TiOH4_	= index_TiOH4;
	index_H2O_    	= index_H2O;	
}

void BatchReactorHomogeneousConstantPressure_TIO2::SetGaseousFeedback(const bool gaseous_feedback)
{
	gaseous_feedback_ = gaseous_feedback;
}

BatchReactorHomogeneousConstantPressure_TIO2::BatchReactorHomogeneousConstantPressure_TIO2
(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap, ModelTiO2& TIO2 ):
	thermodynamicsMap_(thermodynamicsMap), 
	kineticsMap_(kineticsMap),
	TIO2_(TIO2)
{
	NC_ = thermodynamicsMap_.NumberOfSpecies();
	NR_ = kineticsMap_.NumberOfReactions();

	QR_ = 0.;

	omega_.resize(NC_); omega_.setZero();
	x_.resize(NC_); x_.setZero();
	c_.resize(NC_); c_.setZero();
	R_.resize(NC_); R_.setZero();

	checkMassFractions_ = false;
	energyEquation_ = true;
	debug_ = false;
	gaseous_feedback_ = true;
	integrate_TIO2_ = true;
	scaling_factor_ = 6.02214086e26;

	__debug__ = -1;
}

void BatchReactorHomogeneousConstantPressure_TIO2::SetReactor( const double P0 )
{
	P0_    = P0;
}

void BatchReactorHomogeneousConstantPressure_TIO2::SetMassFractions( const Eigen::VectorXd& omega )
{
	omega_ = omega;
}

double BatchReactorHomogeneousConstantPressure_TIO2::T() const
{
	return T_;
}

int BatchReactorHomogeneousConstantPressure_TIO2::Equations(const double t, const Eigen::VectorXd& y, Eigen::VectorXd& dy)
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

	// Recover TiO2 variables	
	const double TIO2_Ys = y(NC_+1); 			// mass fraction
	const double TIO2_N = y(NC_+2) * scaling_factor_; 	// particle number density (#/m3)
	const double TIO2_S = y(NC_+3);				// surface concentration (m2/m3)

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

	// Debug info on the screen
	if (__debug__ != -1)
	{
		std::cout << "Cell: " << __debug__ << " T: " << T_ << " P: " << P0_ << " MW: " << MW_ << " rho: " << rho_ << " cTot: " << cTot_ << std::endl;
		std::cout << " + Ys, N, S: " << TIO2_Ys << " " << TIO2_N << " " << TIO2_S << std::endl;
		std::cout << " + RH2O RTiOH4 RTiO2: " 	<< thermodynamicsMap_.MW(index_H2O_)*R_(index_H2O_) << " " 
							<< thermodynamicsMap_.MW(index_TiOH4_)*R_(index_TiOH4_) << " "
							<< thermodynamicsMap_.MW(index_TiO2_)*R_(index_TiO2_) << std::endl;
	}

	// TIO2 Equations
	if (integrate_TIO2_ == true)
	{
		const double mu = 1.716e-5*Foam::pow(T_/273.15,1.5)*(273.15+110.4)/(T_+110.4);	// [kg/m/s]
		const double cTiOH4 = cTot_*x_(index_TiOH4_);					// [kmol/m3]

		// Source terms: nucleation
		double OmegaYs_Nucleation = 0.;		// (1/s)
		double OmegaN_Nucleation = 0.;		// (1/m3/s)
		double OmegaS_Nucleation = 0.;		// (1/m/s)
		double r_Nucleation = 0.;		// (kmol/m3/s)
		TIO2_.NucleationSourceTerms(T_, rho_, cTiOH4, r_Nucleation, OmegaYs_Nucleation, OmegaN_Nucleation, OmegaS_Nucleation );

		// Source terms: coagulation
		double OmegaYs_Coagulation = 0.;	// (1/s)
		double OmegaN_Coagulation = 0.;		// (1/m3/s)
		double OmegaS_Coagulation = 0.;		// (1/m/s)
		TIO2_.CoagulationSourceTerms(T_, rho_, mu, TIO2_Ys, TIO2_N, TIO2_S, OmegaYs_Coagulation, OmegaN_Coagulation, OmegaS_Coagulation);

		// Source terms: sintering
		double OmegaYs_Sintering = 0.;	// (1/s)
		double OmegaN_Sintering = 0.;		// (1/m3/s)
		double OmegaS_Sintering = 0.;		// (1/m/s)
		TIO2_.SinteringSourceTerms(T_, rho_, TIO2_Ys, TIO2_N, TIO2_S, OmegaYs_Sintering, OmegaN_Sintering, OmegaS_Sintering);

		// 3eq-model equations
		dy(NC_+1) = OmegaYs_Nucleation + OmegaYs_Coagulation + OmegaYs_Sintering;			// (1/s)
		dy(NC_+2) = (OmegaN_Nucleation + OmegaN_Coagulation + OmegaN_Sintering)/scaling_factor_;	// (1/m3/s)
		dy(NC_+3) = OmegaS_Nucleation + OmegaS_Coagulation+ OmegaS_Sintering;				// (1/m/s)

		// Debug info on the screen
		if (__debug__ != -1)
		{
			std::cout << " + mu, x_(index_TiOH4_), cTiOH4: " << OmegaYs_Nucleation << " " << x_(index_TiOH4_) << " " << cTot_*x_(index_TiOH4_) << std::endl;
			std::cout << " + Omega nucleation: " << OmegaYs_Nucleation << " " << OmegaN_Nucleation << " " << OmegaS_Nucleation << std::endl;
			std::cout << " + Omega coagulation: " << OmegaYs_Coagulation << " " << OmegaN_Coagulation << " " << OmegaS_Coagulation << std::endl;
			std::cout << " + Omega sintering: " << OmegaYs_Sintering << " " << OmegaN_Sintering << " " << OmegaS_Sintering << std::endl;	
		}

		// Consumption
		if (gaseous_feedback_ == true)
		{
			const double OmegaTiOH4 = -2.*r_Nucleation * thermodynamicsMap_.MW(index_TiOH4_);	// [kg/m3/s]
			const double OmegaTiO2 = 2.*r_Nucleation * thermodynamicsMap_.MW(index_TiO2_);		// [kg/m3/s]
			const double OmegaH2O = 4.*r_Nucleation * thermodynamicsMap_.MW(index_H2O_) ;		// [kg/m3/s]

			dy(index_TiOH4_) += OmegaTiOH4/rho_;
			dy(index_TiO2_) += OmegaTiO2/rho_;
			dy(index_H2O_) += OmegaH2O/rho_;
		}
	}
	else
	{
		dy(NC_+1) = 0.;
		dy(NC_+2) = 0.;
		dy(NC_+3) = 0.;
	}

	return 0;
}

int BatchReactorHomogeneousConstantPressure_TIO2::Print(const double t, const Eigen::VectorXd& y)
{
	return 0;
}

#endif // BatchReactorHomogeneousConstantPressure_TIO2_H
