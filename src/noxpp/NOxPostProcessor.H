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

#ifndef NOxPostProcessor_H
#define NOxPostProcessor_H

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
#include "ModelNOxPostProcessor.H"

namespace Foam
{

//!  A class to estimate NOx emissions (derived from the Ansys Fluent NOx post-processor)
/*!
This class provides the tools to estimate NOx emissions (derived from the Ansys Fluent NOx post-processor)
*/

class NOxPostProcessor : public IOdictionary
{

private:

//	enum DiffusionType { FIXED_SCHMIDT_NUMBER, COLLISION_DIAMETER };

public:

	/**
	*@brief Default constructor
	*@param mesh mesh
	*@param thermodynamicsMapXML thermodynamic map
	*@param kineticsMapXML kinetics map
	*/
	NOxPostProcessor(const fvMesh& mesh, OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap);

	/**
	*@brief Default destructor
	*/
	virtual ~NOxPostProcessor();

	/**
	*@brief Name of main dictionary
	*/
	static const word dictName;

	/**
	*@brief Additional operations for preparing the NOxPostProcessor model
	*@param mesh mesh
	*/
	void Setup(const Foam::fvMesh& mesh);

	//- Update the relevant properties having impact on the simulation
	void UpdateProperties(const volScalarField& p, const volScalarField& T, const volScalarField& rho, const PtrList<volScalarField>& Dmix);

	//- On-the-fly analysis
	void OnTheFlyAnalysis(const Foam::fvMesh& mesh, const volScalarField& rho, const PtrList<volScalarField>& Y);

	//- Reconstruction of NOxPostProcessor fields (fv, N, partial density, mass and mole fractions)
	void FieldsReconstruction(const Foam::fvMesh& mesh, const volScalarField& rho, const volScalarField& p, const volScalarField& T, const PtrList<volScalarField>& Y);

	//- Solve transport equations
	void SolveEquations(const Foam::fvMesh& mesh, const volScalarField& p, const volScalarField& T, const surfaceScalarField& phi, const volScalarField& rho, const PtrList<volScalarField>& Y);

	//- Mass fraction of NO (-)
	const volScalarField& YNO() const;

	//- Mass fraction of NO (-)
	//  Non-const access allowed for transport equations
	volScalarField& YNO();
/*
	//- Mass fraction of HCN (-)
	const volScalarField& YHCN() const;

	//- Mass fraction of HCN (-)
	//  Non-const access allowed for transport equations
	volScalarField& YHCN();

	//- Mass fraction of NH3 (-)
	const volScalarField& YNH3() const;

	//- Mass fraction of NH3 (-)
	//  Non-const access allowed for transport equations
	volScalarField& YNH3();

	//- Mass fraction of N2O (-)
	const volScalarField& YN2O() const;

	//- Mass fraction of N2O (-)
	//  Non-const access allowed for transport equations
	volScalarField& YN2O();
*/
	//- Return the NOx analyzer
	ModelNOxPostProcessor& NOxPostProcessorModel() { return *model_NOxPostProcessor_; }

	//- Indices of relevant species (0-index based)
	int index_Fuel() const { return index_fuel_; }
	int index_N2() const { return index_n2_; }
	int index_O2() const { return index_o2_; }
	int index_H2O() const { return index_h2o_; }
	int index_OH() const { return index_oh_; }
	int index_O() const { return index_o_; }
	int index_CH() const { return index_ch_; }

private:

	void reactionSourceTerms( const Eigen::VectorXd& Y, const double T0, const double P0, const double rho0, const Eigen::VectorXd& Y0, Eigen::VectorXd& S );
	void reactionJacobian( const Eigen::VectorXd& y, const double T0, const double P0, const double rho0, const Eigen::VectorXd& y0, Eigen::VectorXd &J);

private:

	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap_;		//!< reference to the thermodynamic map
	OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap_;				//!< reference to the kinetics map 

	ModelNOxPostProcessor* model_NOxPostProcessor_;

	unsigned int ne_;
	Switch is_implicit_mode_;

	//- Mass fraction of NO [-]
	volScalarField YNO_;
/*
	//- Mass fraction of HCN [-]
	volScalarField YHCN_;

	//- Mass fraction of NH3 [-]
	volScalarField YNH3_;

	//- Mass fraction of N2O [-]
	volScalarField YN2O_;
*/
	//- Formation rate of NO from thermal mechanism (kg/m3/s)
	volScalarField SNO_thermal_;

	//- Formation rate of NO from prompt mechanism (kg/m3/s)
	volScalarField SNO_prompt_;

	//- Formation rate of NO from intermediate N2O mechanism (kg/m3/s)
	volScalarField SNO_intermediate_n2o_;

	//- Source terms for steady-state solver (kg/m3/s)
	volScalarField sourceImplicitNO_;
	volScalarField sourceExplicitNO_;

	//- Diffusion coefficient [kg/m/s]
	volScalarField gamma_;			

	unsigned int NR_;			//!< total number of reactions
	unsigned int NS_;			//!< total number of species

	int index_fuel_;			//!< index of Fuel species in the kinetic mechanism (0-based)
	int index_n2_;				//!< index of N2 species in the kinetic mechanism (0-based)
	int index_o2_;				//!< index of O2 species in the kinetic mechanism (0-based)
	int index_h2o_;				//!< index of H2O species in the kinetic mechanism (0-based)

	int index_oh_;				//!< index of OH species in the kinetic mechanism (0-based)
	int index_o_;				//!< index of O species in the kinetic mechanism (0-based)
	int index_ch_;				//!< index of CH species in the kinetic mechanism (0-based)
	
	Switch is_thermal_;			//!< thermal NOx on/off
	Switch is_prompt_;			//!< prompt NOx on/off
	Switch is_intermediate_n2o_;		//!< formation of NOx from intermediate N2O

	Switch is_fuel_nox_;			//!< fuel NOx on/off
	Switch is_reduction_by_reburning_;	//!< reduction of NOx by reburning
	Switch is_reduction_by_sncr_;		//!< reduction of NOx by SNCR

	//- Data for writing report (on-the-fly analysis)
	Switch ontheflyAnalysis_;
	label frequency_;
	label frequencyCounter_;

private:

	static const double MW_NO_;	// NO molecular weight (kg/kmol)
	static const double MW_N2_;	// N2 molecular weight (kg/kmol)
	static const double MW_O2_;	// O2 molecular weight (kg/kmol)
	static const double MW_H2O_;	// H2O molecular weight (kg/kmol)
	static const double MW_O_;	// O molecular weight (kg/kmol)
	static const double MW_OH_;	// OH molecular weight (kg/kmol)
		     double MW_Fuel_;	// Fuel molecular weight (kg/kmol)

};

} // end namespace Foam

#include "NOxPostProcessor.HPP"

#endif	
