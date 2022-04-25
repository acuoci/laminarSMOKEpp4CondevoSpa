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
|   Copyright(C) 2020 Alberto Cuoci                                       |
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

#ifndef Clustering_H
#define Clustering_H

// Include standard OpenFOAM files
#include "volFields.H"
#include "typeInfo.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "wordIOList.H"
#include "hashedWordList.H"

// Clustering algorithms
#include "BinningAlgorithm.H"
#include "DMZAlgorithm.H"

namespace Foam
{

//!  A class to perform clustering of computational cells based on arbitrary fields provided by the user
/*!
This class provides the tools to perform clustering of computational cells based on arbitrary fields provided by the user
*/

class Clustering : public IOdictionary
{
public:

	enum SpeciesBackwardMapping { SPECIES_LINEAR_INTERPOLATION_WITH_CORRECTION, SPECIES_LINEAR_INTERPOLATION } ;

	enum TemperatureBackwardMapping { TEMPERATURE_LINEAR_INTERPOLATION, TEMPERATURE_RECONSTRUCTED } ;
public:

	/**
	*@brief Default constructor
	*@param mesh mesh
	*@param thermodynamicsMapXML thermodynamic map
	*/
	Clustering(const fvMesh& mesh);

	/**
	*@brief Default destructor
	*/
	virtual ~Clustering();

	/**
	*@brief Name of main dictionary
	*/
	static const word dictName;

	/**
	*@brief Additional operations for preparing the Clustering model
	*@param mesh mesh
	*/
	void Setup(const Foam::fvMesh& mesh);

	/**
	*@brief Perform the clustering operations
	*@param mesh mesh
	*/
	void Analysis(const Foam::fvMesh& mesh);

	/**
	*@brief Sets the key-fields 
	*@param names_key_fields vector containing the names of key-fields
	*/
	void SetKeyFields(List<word>& fields) { fields_ = fields; }

	
	unsigned int nc() const { return nc_; }
	const std::vector<int>& gj() const { return gj_; }
	const std::vector< std::vector<int> >& lj() const { return lj_; }

	/**
	*@brief Returns the backward mapping policy adopted for species
	*/
	SpeciesBackwardMapping species_backward_mapping() const { return species_backward_mapping_; }

	/**
	*@brief Returns the backward mapping policy adopted for temperature
	*/
	TemperatureBackwardMapping temperature_backward_mapping() const { return temperature_backward_mapping_; }

	/**
	*@brief Returns the index of cluster to which the cell belongs
	*/
	const volScalarField& cluster_index() const;

	/**
	*@brief Returns the index of cluster to which the cell belongs (non-const access)
	*/
	volScalarField& cluster_index();

	/**
	*@brief Returns the size of the cluster to which the cell belongs
	*/
	const volScalarField& cluster_size() const;

	/**
	*@brief Returns the size of the cluster to which the cell belongs (non-const access)
	*/
	volScalarField& cluster_size();


private:

	List<double>  epsilons_;
	List<word>    fields_;
	List<bool>    epsilons_rel_;

	label D_;
	enum {CLUSTERING_DMZ, CLUSTERING_BINNING} algorithm_;
	label verbosity_;
	SpeciesBackwardMapping species_backward_mapping_; 
	TemperatureBackwardMapping temperature_backward_mapping_;

	unsigned int nc_;
	std::vector<int> gj_;
	std::vector< std::vector<int> > lj_;

	volScalarField cluster_index_;
	volScalarField cluster_size_;
};

} // end namespace Foam

#include "Clustering.HPP"

#endif	
