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

#ifndef TiO2_H
#define TiO2_H

// Include standard OpenFOAM files
#include "volFields.H"
#include "typeInfo.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "wordIOList.H"
#include "hashedWordList.H"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

// Clustering algorithms
#include "ModelTiO2.H"

namespace Foam
{

//!  A class to apply the TiO2 model based on 3 equations
/*!
This class provides the toolsto apply the TiO2 model based on 3 equations
*/

class TiO2 : public IOdictionary
{

private:

	enum DiffusionType { FIXED_SCHMIDT_NUMBER, COLLISION_DIAMETER };

public:

	/**
	*@brief Default constructor
	*@param mesh mesh
	*@param thermodynamicsMapXML thermodynamic map
	*@param kineticsMapXML kinetics map
	*/
	TiO2(const fvMesh& mesh, OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap);

	/**
	*@brief Default destructor
	*/
	virtual ~TiO2();

	/**
	*@brief Name of main dictionary
	*/
	static const word dictName;

	/**
	*@brief Additional operations for preparing the TiO2 model
	*@param mesh mesh
	*/
	void Setup(const Foam::fvMesh& mesh);

	//- Update the relevant properties having impact on the simulation
	void UpdateProperties(const volScalarField& p, const volScalarField& T, const volScalarField& rho, const volScalarField& mu);

	//- On-the-fly analysis
	void OnTheFlyAnalysis(const Foam::fvMesh& mesh, const volScalarField& rho, const PtrList<volScalarField>& Y);

	//- Reconstruction of TiO2 fields (fv, N, partial density, mass and mole fractions)
	void FieldsReconstruction(const Foam::fvMesh& mesh, const volScalarField& rho, const volScalarField& T);

	//- Solve transport equations
	void SolveEquations(const Foam::fvMesh& mesh, const volScalarField& p, const volScalarField& T, const surfaceScalarField& phi, const volScalarField& rho, const volScalarField& mu);

	//- Mass fraction Ys (-)
	const volScalarField& Ys() const;

	//- Mass fraction Ys (-)
	//  Non-const access allowed for transport equations
	volScalarField& Ys();

	//- Surface per unit of volume (m2/m3)
	const volScalarField& S() const;

	//- Surface per unit of volume (m2/m3)
	//  Non-const access allowed for transport equations
	volScalarField& S();

	//- Normalized particle number density N (kmol/m3)
	const volScalarField& Norm_N() const;

	//- Normalized particle number density N (kmol/m3)
	//  Non-const access allowed for transport equations
	volScalarField& Norm_N();

	//- Return the HMOM analyzer
	ModelTiO2& tio2Model() { return *model_tio2_; }

	//- Indices of relevant species (0-index based)
	int index_TiO2() const { return index_tio2_; }
	int index_TiOH4() const { return index_tioh4_; }
	int index_H2O() const { return index_h2o_; }
	
	//- Gaseous feedback
	Switch is_gaseous_feedback() const { return is_gaseous_feedback_; }


private:

	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap_;		//!< reference to the thermodynamic map
	OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap_;				//!< reference to the kinetics map 

	ModelTiO2* model_tio2_;

	//- Mass fraction of TiO2 [-]
	volScalarField Ys_;

	//- Normalized particle number density [kmol/m3]
	volScalarField Norm_N_;

	//- Particle number density [#/m3]
	volScalarField N_;

	//- Surface concentration [m2/m3]
	volScalarField S_;

	//- Diffusion coefficient [kg/m/s]
	volScalarField gamma_;

	//- Soot volume fraction fv [-]
	volScalarField fv_; 

	//- Primary particle diameter dp [m]
	volScalarField dp_;

	//- Collision diameter dc [m]
	volScalarField dc_; 

	//- Aggregate diameter da [m]
	volScalarField da_;

	//- Surface area of the completely fused particles [m2]
	volScalarField ss_;

 	//- Number o primary particles np [-]
	volScalarField np_;

	//- Characteristic sintering time for TiO2 [s]
	volScalarField taus_;			

	unsigned int NR_;			//!< total number of reactions
	unsigned int NS_;			//!< total number of species

	double schmidt_number_;			//!< Schmidt number for TiO2 particles
	DiffusionType diffusion_type_; 		//!< Diffsuion type: FIXED_SCHMIDT_NUMBER | COLLISION_DIAMETER

	int index_tio2_;			//!< index of TiO2 species in the kineti mechanism (0-based)
	int index_tioh4_;			//!< index of TiOH4 species in the kineti mechanism (0-based)
	int index_h2o_;				//!< index of H2O species in the kineti mechanism (0-based)
	
	Switch is_nucleation_;			//!< nucleation on/off
	Switch is_coagulation_;			//!< coagulation on/off
	Switch is_sintering_;			//!< sintering on/off
	
	Switch is_thermophoretic_effect_;	//!< thermophoretic effect on/off
	Switch is_gaseous_feedback_; 		//!< consumption of gaseous species because of nucleation on/off

	//- Data for writing report (on-the-fly analysis)
	Switch ontheflyAnalysis_;
	label frequency_;
	label frequencyCounter_;
};

} // end namespace Foam

#include "TiO2.HPP"

#endif	
