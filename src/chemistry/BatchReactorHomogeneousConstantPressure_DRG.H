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

#ifndef BatchReactorHomogeneousConstantPressure_DRG_H
#define	BatchReactorHomogeneousConstantPressure_DRG_H

class BatchReactorHomogeneousConstantPressure_DRG
{
public:

	//- Default constructor
	BatchReactorHomogeneousConstantPressure_DRG(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, 
							OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap);

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

	//- Set the DRG options
	void SetDRGOptions(const Eigen::VectorXi& indices_important_species, const Eigen::VectorXi& indices_unimportant_reactions, const Eigen::VectorXi& important_reaction_steps);


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

	unsigned int NC_;		//!< number of species (total)
	unsigned int NR_;		//!< number of reactions (total)

	unsigned int nis_;						//!< number of important species
	unsigned int nur_;						//!< number of unimportant reactions
	Eigen::VectorXi indices_unimportant_reactions_;	//!< indices of unimportant reactions (zero-based)
	Eigen::VectorXi indices_important_species_;		//!< indices of important species (zero-based)
	Eigen::VectorXi important_reaction_steps_;

	Eigen::VectorXd omega_;		//!< mass fractions
	Eigen::VectorXd x_;		//!< mole fractions
	Eigen::VectorXd c_;		//!< concentrations (in kmol/m3)
	Eigen::VectorXd R_;		//!< formation rates (in kmol/m3/s)
	Eigen::VectorXd r_;		//!< reaction rates (in kmol/m3/s)

	bool checkMassFractions_;	//!< if true, mass fractions are normalized
	bool energyEquation_;		//!< if true, the energy equation is solved
	bool debug_;			//!< if true, additional debug (time-consuming) operations are carried out

};

unsigned int BatchReactorHomogeneousConstantPressure_DRG::NumberOfEquations() const
{
	return nis_+1;
}

BatchReactorHomogeneousConstantPressure_DRG::BatchReactorHomogeneousConstantPressure_DRG
(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap):
	thermodynamicsMap_(thermodynamicsMap), 
	kineticsMap_(kineticsMap)
{
	NC_ = thermodynamicsMap_.NumberOfSpecies();
	NR_ = kineticsMap_.NumberOfReactions();
	nis_ = NC_;
	nur_ = 0;

	QR_ = 0.;

	omega_.resize(NC_); omega_.setZero();
	x_.resize(NC_); x_.setZero();
	c_.resize(NC_); c_.setZero();
	R_.resize(NC_); R_.setZero();
	r_.resize(NR_); r_.setZero();

	checkMassFractions_ = false;
	energyEquation_ = true;
	debug_ = false;
}

void BatchReactorHomogeneousConstantPressure_DRG::SetReactor( const double P0 )
{
	P0_    = P0;
}

void BatchReactorHomogeneousConstantPressure_DRG::SetMassFractions( const Eigen::VectorXd& omega )
{
	omega_ = omega;
}

double BatchReactorHomogeneousConstantPressure_DRG::T() const
{
	return T_;
}

void BatchReactorHomogeneousConstantPressure_DRG::SetDRGOptions(const Eigen::VectorXi& indices_important_species, const Eigen::VectorXi& indices_unimportant_reactions, const Eigen::VectorXi& important_reaction_steps)
{
	nis_ = indices_important_species.size();
	nur_ = indices_unimportant_reactions.size();

	indices_important_species_ = indices_important_species;
	indices_unimportant_reactions_ = indices_unimportant_reactions;
	important_reaction_steps_ = important_reaction_steps;
}

int BatchReactorHomogeneousConstantPressure_DRG::Equations(const double t, const Eigen::VectorXd& y, Eigen::VectorXd& dy)
{
	// Recover mass fractions
	if (checkMassFractions_ == true)
	{	
		for (unsigned int i=0;i<nis_;++i)	
		{
			const unsigned int j = indices_important_species_(i);
			omega_(j) = max(y(i), 0.);
		}
	}
	else
	{
		for (unsigned int i=0;i<nis_;++i)	
		{
			const unsigned int j = indices_important_species_(i);
			omega_(j) = y(i);
		}
	}

	// In case there are important reactions
	if (nur_ < NR_)
	{
		// Recover temperature
		const unsigned int index_T = nis_;
		T_ = y(index_T);

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

		// Calculates reaction rates
		kineticsMap_.KineticConstants(important_reaction_steps_);
		kineticsMap_.ReactionRates(c_.data());

		// Remove useless reactions
		kineticsMap_.GiveMeReactionRates(r_.data());
		for (unsigned int i=0;i<nur_;++i)
			r_(indices_unimportant_reactions_(i)) = 0.;

		// Formation rates
		kineticsMap_.stoichiometry().FormationRatesFromReactionRates(R_.data(), r_.data());

		// Species equations
		for (unsigned int i=0;i<nis_;++i)	
		{
			const unsigned int j = indices_important_species_(i);
			dy(i) = thermodynamicsMap_.MW(j)*R_(j)/rho_;
		}

	    	// Energy equation
	    	dy(index_T) = 0.;     
	    	if (energyEquation_ == true)
	    	{
			CpMixMass_ = thermodynamicsMap_.cpMolar_Mixture_From_MoleFractions(x_.data()) / MW_;
			QR_ = kineticsMap_.HeatRelease(R_.data());
	
			dy(index_T)  = QR_ / (rho_*CpMixMass_);
		}
	}
	// If there are no important reactions
	else
	{
		dy.setZero();
		QR_ = 0.;
	}

	return 0;
}

int BatchReactorHomogeneousConstantPressure_DRG::Print(const double t, const Eigen::VectorXd& y)
{
	return 0;
}

#endif // BatchReactorHomogeneousConstantPressure_DRG_H
