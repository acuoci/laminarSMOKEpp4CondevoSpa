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

#ifndef polimiSootModel_H
#define polimiSootModel_H

// Include standard OpenFOAM files
#include "volFields.H"
#include "typeInfo.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "wordIOList.H"
#include "hashedWordList.H"

// OpenSMOKE++ Polimi Soot Analyzer
#include "utilities/soot/polimi/OpenSMOKE_PolimiSoot_Analyzer.h"

namespace Foam
{

class PolimiSootModel : public IOdictionary
{

public: 

	//- Construct from mesh
	PolimiSootModel(const fvMesh& mesh, OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap, OpenSMOKE::TransportPropertiesMap_CHEMKIN& transportMap);

	//- Destructor
	virtual ~PolimiSootModel();

	//- Preparation of Polimi Soot Model
	void Setup(const Foam::fvMesh& mesh);

	//- Update the relevant properties having impact on the simulation
	void UpdateProperties(const volScalarField& p, const volScalarField& T, const PtrList<volScalarField>& Y);

	//- On-the-fly analysis
	void OnTheFlyAnalysis(const Foam::fvMesh& mesh, const volScalarField& rho, const PtrList<volScalarField>& Y);

	//- Reconstruction of soot fields (fv, N, partial density, mass and mole fractions)
	void FieldsReconstruction(const Foam::fvMesh& mesh, const volScalarField& p, const volScalarField& T, const PtrList<volScalarField>& Y);

	//- Name of main dictionary
	static const word dictName;

	//- Returns true if physical diffusivity is turned on
	Switch physicalDiffusion() const { return physicalDiffusion_; }

	//- Returns true if thermophoretic effect is turned on
	Switch thermophoreticEffect() const { return thermophoreticEffect_; }

	//- Data for physical diffusion of soot particles
	label physicalDiffusionReferenceIndex() const { return sootAnalyzer_->bin_physical_diffusion_reference_species(); }
	
	//- Data for physical diffusion of soot particles
	const std::vector<double>& physicalDiffusionCorrection() const { return sootAnalyzer_->bin_physical_diffusion_correction_factors(); }

	//- Data for physical diffusion of soot particles
	const std::vector<unsigned int> physicalDiffusionCorrectionIndex() const { return sootAnalyzer_->bin_indices(); }

	//- Returns the indices of species to which apply the thermophoretic diffusion
	const std::vector<unsigned int> thermophoretic_effect_list() const { return sootAnalyzer_->bin_indices_thermophoresis(); }
	
private:

	//- Thermodynamics map
	OpenSMOKE::ThermodynamicsMap_CHEMKIN&		thermodynamicsMap_; 

	//- Kinetics map
	OpenSMOKE::KineticsMap_CHEMKIN& 		kineticsMap_;

	//- Transport properties map
	OpenSMOKE::TransportPropertiesMap_CHEMKIN& 	transportMap_;

	//- Soot fields: 0=volume fraction, 1=density, 2=number of particles, 3=mole fraction, 4=mass fraction, 5=formation rate
	PtrList<volScalarField> sootFields_;

	//- PAH fields: 0=mass fraction 1-2 rings, 1=mass fraction 3-4 rings, 2=mass fraction large precursors
	PtrList<volScalarField> pahFields_;

	//- Thermophoretic effect tag
	Switch thermophoreticEffect_;

	//- Data for physical diffusion of soot particles
	Switch physicalDiffusion_;

	//- Data for writing report (on-the-fly analysis)
	Switch ontheflyAnalysis_;
	label frequency_;
	label frequencyCounter_;

	//- Polimi Soot Analyzer (from OpenSMOKE++ libraries)
	OpenSMOKE::PolimiSoot_Analyzer* sootAnalyzer_;

	//- Log file
	autoPtr<std::ofstream> fSootIntegrals_;


	//- Returns the integral mass of soot particles (large BINS) [kg] 
	double largeBinMass(const volScalarField& rho, const PtrList<volScalarField>& Y);

	//- Returns the integral volume of soot particles (large BINs) [m3] 
	double largeBinVolume(const volScalarField& rho, const PtrList<volScalarField>& Y);

	//- Returns the integral number of soot particles (large BINs) [#] 
	double largeBinNumberParticles(const volScalarField& rho, const PtrList<volScalarField>& Y);

	//- Returns the integral mass of soot particles (small BINS) [kg] 
	double smallBinMass(const volScalarField& rho, const PtrList<volScalarField>& Y);

	//- Returns the integral volume of soot particles (small BINs) [m3] 
	double smallBinVolume(const volScalarField& rho, const PtrList<volScalarField>& Y);

	//- Returns the integral number of soot particles (small BINs) [#] 
	double smallBinNumberParticles(const volScalarField& rho, const PtrList<volScalarField>& Y);

	//- Returns the integral mass of dimer particles [kg] 
	double sootDimerMass(const volScalarField& rho, const PtrList<volScalarField>& Y);

	//- Returns the integral mass of PAHs with 1 or 2 aromatic rings (340 nm) [kg] 
	double pah_1_2_rings_Mass(const volScalarField& rho, const PtrList<volScalarField>& Y);

	//- Returns the integral mass of PAHs with 3 or 4 aromatic rings (400 nm) [kg] 
	double pah_3_4_rings_Mass(const volScalarField& rho, const PtrList<volScalarField>& Y);

	//- Returns the integral mass of PAHs with more than 4 rings (500 nm) [kg] 
	double pah_more_than_4_rings_Mass(const volScalarField& rho, const PtrList<volScalarField>& Y);

	//- Returns the total formation rate of soot particles (large BINs) [kg/m3/s]
	double largeBinMassFormationRates(const Eigen::VectorXd& R_times_W);

	//- Returns the total formation rate of soot particles (small BINs) [kg/m3/s]
	double smallBinMassFormationRates(const Eigen::VectorXd& R_times_W);

	//- Returns the total formation rate of PAHs with 1 or 2 aromatic rings (340 nm) [kg/m3/s] 
	double pah12RingsFormationRates(const Eigen::VectorXd& R_times_W);

	//- Returns the total formation rate of PAHs with 3 or 4 aromatic rings (400 nm) [kg/m3/s] 
	double pah34RingsFormationRates(const Eigen::VectorXd& R_times_W);

	//- Returns the total formation rate of PAHs with more that aromatic rings (500 nm) [kg/m3/s] 
	double pahMoreThan4RingsFormationRates(const Eigen::VectorXd& R_times_W);

};

} // end namespace Foam

#include "PolimiSootModel.HPP"

#endif	

