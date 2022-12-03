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

#ifndef HMOMModel_H
#define HMOMModel_H

// Include standard OpenFOAM files
#include "volFields.H"
#include "typeInfo.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "wordIOList.H"
#include "hashedWordList.H"

// Soot HMOM in OpenSMOKE++
#include "utilities/soot/hmom/HMOM.h"

namespace Foam
{

class HMOMModel : public IOdictionary
{

public: 

	//- Construct from mesh
	HMOMModel(const fvMesh& mesh, OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap, OpenSMOKE::TransportPropertiesMap_CHEMKIN& transportMap);

	//- Destructor
	virtual ~HMOMModel();

	//- Preparation of Polimi Soot Model
	void Setup(const Foam::fvMesh& mesh);

	//- Update the relevant properties having impact on the simulation
	void UpdateProperties(const volScalarField& p, const volScalarField& T);

	//- On-the-fly analysis
	void OnTheFlyAnalysis(const Foam::fvMesh& mesh, const volScalarField& rho, const PtrList<volScalarField>& Y);

	//- Reconstruction of soot fields (fv, N, partial density, mass and mole fractions)
	void FieldsReconstruction(const Foam::fvMesh& mesh, const volScalarField& p, const volScalarField& T);

	//- Solve transport equations
	void SolveEquations(const Foam::fvMesh& mesh, const volScalarField& p, const volScalarField& T, const surfaceScalarField& phi, const volScalarField& rho, const volScalarField& mu, const Switch thermophoretic_effect);

	//- Name of main dictionary
	static const word dictName;

	//- Return the HMOM analyzer
	OpenSMOKE::HMOM& hmomAnalyzer() { return *hmomAnalyzer_; }

	//- Moment N0 [mol/m3]
	const volScalarField& Norm_N0() const;

	//- Moment N0 [mol/m3]
	//  Non-const access allowed for transport equations
	volScalarField& Norm_N0();

	//- Moment M00 [mol/m3]
	const volScalarField& Norm_M00() const;

	//- Moment M00 [mol/m3]
	//  Non-const access allowed for transport equations
	volScalarField& Norm_M00();

	//- Moment M10 [mol/m3]
	const volScalarField& Norm_M10() const;

	//- Moment M10 [mol/m3]
	//  Non-const access allowed for transport equations
	volScalarField& Norm_M10();

	//- Moment M01 [mol/m3]
	const volScalarField& Norm_M01() const;

	//- Moment M01 [mol/m3]
	//  Non-const access allowed for transport equations
	volScalarField& Norm_M01();

	//- Soot volume fraction
	const volScalarField& fv() const;

	//- Soot volume fraction
	//  Non-const access allowed for transport equations
	volScalarField& fv();

	//- Should the PAH consumption be enabled?
	Switch pah_gas_consumption() const { return hmom_pah_gas_consumption_; }

	//- This is the minimu soot volume fraction for integrating the moment equations
	double minimum_fv_for_integration() const { return hmom_minimum_fv_for_integration_; }

	//- Indices of PAHs (0-index based)
	const std::vector<int> pah_species_indices() const { return hmom_pah_species_indices_; }

	//- Indices of relevant species (0-index based)
	int index_H() const { return hmom_index_H_; }
	int index_OH() const { return hmom_index_OH_; }
	int index_H2O() const { return hmom_index_H2O_; }
	int index_H2() const { return hmom_index_H2_; }
	int index_O2() const { return hmom_index_O2_; }
	int index_C2H2() const { return hmom_index_C2H2_; }


private:

	//- Thermodynamics map
	OpenSMOKE::ThermodynamicsMap_CHEMKIN&		thermodynamicsMap_; 

	//- Kinetics map
	OpenSMOKE::KineticsMap_CHEMKIN& 		kineticsMap_;

	//- Transport properties map
	OpenSMOKE::TransportPropertiesMap_CHEMKIN& 	transportMap_;

	//- Moment N0 [mol/m3]
	volScalarField Norm_N0_;
		
	//- Moment M00 [mol/m3]
	volScalarField Norm_M00_;
		
	//- Moment M10 [mol/m3]
	volScalarField Norm_M10_;
		
	//- Moment M01 [mol/m3]
	volScalarField Norm_M01_;
		
	//- Soot volume fraction fv [-]
	volScalarField fv_;
		
	//- Soot number of particles N [#/m3]
	volScalarField N_;
		
	//- Primary particle diameter dp [m]
	volScalarField dp_;
		
	//- Xollisional diameter dc [m]
	volScalarField dc_;
		
	//- Number o primary particles np [-]
	volScalarField np_;

	//- Soot volume fraction fv [-] (needed by the radiation model)
	volScalarField soot_fv_;

	//- Diffusion coefficient [kg/m/s]
	volScalarField gamma_;

	//- Data for writing report (on-the-fly analysis)
	Switch ontheflyAnalysis_;
	label frequency_;
	label frequencyCounter_;

	//- Polimi Soot Analyzer (from OpenSMOKE++ libraries)
	OpenSMOKE::HMOM* hmomAnalyzer_;

	//- Log file
	autoPtr<std::ofstream> fSootIntegrals_;

	//- This is the minimum allowed value for normalized moments
	scalar hmom_minimum_normalized_value_;

	//- This is the minimu soot volume fraction for integrating the moment equations
	scalar hmom_minimum_fv_for_integration_;

	//- Should the PAH consumption be enabled?
	Switch hmom_pah_gas_consumption_;

	//- Indices of PAHs
	std::vector<int> hmom_pah_species_indices_;

	//- Indices of relevant species (0-index based)
	int hmom_index_H_ ;
	int hmom_index_OH_;
	int hmom_index_H2O_;
	int hmom_index_H2_;
	int hmom_index_O2_;
	int hmom_index_C2H2_;
};

} // end namespace Foam

#include "HMOMModel.HPP"

#endif	
