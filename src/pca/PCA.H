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
|   Copyright(C) 2021 Alberto Cuoci                                       |
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

#ifndef PCA_H
#define PCA_H

// Include standard OpenFOAM files
#include "volFields.H"
#include "typeInfo.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "wordIOList.H"
#include "hashedWordList.H"

// PCA algorithms
#include "pcaEigen.H"

namespace Foam
{

//!  A class to perform PCA of cells on the basis of composition and temperature
/*!
This class provides the tools to perform PCA of cells on the basis of composition and temperature
*/

class PCA : public IOdictionary
{
public:

	enum SamplingType { SAMPLING_SIMPLE_STEP };

public:

	/**
	*@brief Default constructor
	*@param mesh mesh
	*@param thermodynamicsMapXML thermodynamic map
	*/
	PCA(const fvMesh& mesh);

	/**
	*@brief Default destructor
	*/
	virtual ~PCA();

	/**
	*@brief Name of main dictionary
	*/
	static const word dictName;

	/**
	*@brief Additional operations for preparing the PCA model
	*@param mesh mesh
	*/
	void Setup(const Foam::fvMesh& mesh, const OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap);

	/**
	*@brief Perform the PCA
	*@param mesh mesh
	*/
	void Analysis(const volScalarField& T, const PtrList<volScalarField>& Y);

	/**
	*@brief Perform the PCA without calculating the scores
	*@param mesh mesh
	*/
	void Analysis(const volScalarField& T, const PtrList<volScalarField>& Y, Eigen::VectorXd& explained, Eigen::MatrixXd& cosines);

	/**
	*@brief Returns the list of enabled species indices
	*/
	const std::vector<int>& enabled_species_indices() const { return enabled_species_indices_; }


private:

	//- Principal Component Fields
	PtrList<volScalarField> pc_;

	//- Number of retained PCs
	unsigned int n_;

	//- Sampling type
	SamplingType samplingType_;

	//- Sampling step: 1 means every cell, 2 means every 2 cells, etc.
	unsigned int samplingStep_;

	//- Verbosity level
	label verbosity_;

	//- List of species on which to apply the PCA
	List<word> enabled_species_names_;

	//- Indices (0-based) of species on which to apply the PCA
	std::vector<int> enabled_species_indices_;

	//- Enable calculation of PCs on boundary faces (for graphical reasons only)
	Switch boundaryFaces_;
};

} // end namespace Foam

#include "PCA.HPP"

#endif	
