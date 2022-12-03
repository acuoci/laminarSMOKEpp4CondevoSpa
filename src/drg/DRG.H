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

#ifndef DRG_H
#define DRG_H

// Include standard OpenFOAM files
#include "volFields.H"
#include "typeInfo.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "wordIOList.H"
#include "hashedWordList.H"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

namespace Foam
{

//!  A class to perform DRG (Direct Relation Graph)
/*!
This class provides the tools to perform DRG on detailed kinetic mechanisms with arbitrary number of species and reactions
*/

class DRG : public IOdictionary
{
public:

	/**
	*@brief Default constructor
	*@param mesh mesh
	*@param thermodynamicsMapXML thermodynamic map
	*@param kineticsMapXML kinetics map
	*/
	DRG(const fvMesh& mesh, OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap);

	/**
	*@brief Default destructor
	*/
	virtual ~DRG();

	/**
	*@brief Name of main dictionary
	*/
	static const word dictName;

	/**
	*@brief Additional operations for preparing the DRG model
	*@param mesh mesh
	*/
	void Setup(const Foam::fvMesh& mesh);

	/**
	*@brief Sets the key-species (i.e. the target species which are considered important)
	*@param names_key_species vector containing the names of target or key species (0-index based)
	*/
	void SetKeySpecies(const std::vector<std::string> names_key_species);

	/**
	*@brief Sets the key-species (i.e. the target species which are considered important)
	*@param names_key_species vector containing the indices of target or key species (0-index based)
	*/
	void SetKeySpecies(const Eigen::VectorXi& indices_key_species);

	/**
	*@brief Performs the DRG analysis for the given conditions
	*@param T temperature in K
	*@param P_Pa pressure in Pa
	*@param c vector of concentrations in kmol/m3
	*/
	void Analysis(const double T, const double P_Pa, const Eigen::VectorXd& c);

	/**
	*@brief Returns the number of key species
	*/
	unsigned int number_key_species() const { return index_key_species_.size(); }

	/**
	*@brief Returns the number of important species
	*/
	unsigned int number_important_species() const { return number_important_species_; }

	/**
	*@brief Returns the number of important reactions
	*/
	unsigned int number_important_reactions() const { return number_important_reactions_; }

	/**
	*@brief Returns the indices of important species (zero-based)
	*/
	const Eigen::VectorXi& indices_important_species() const { return indices_important_species_; }

	/**
	*@brief Returns the indices of important reactions (zero-based)
	*/
	const Eigen::VectorXi& indices_important_reactions() const { return indices_important_reactions_; }

	/**
	*@brief Returns the indices of unimportant reactions (zero-based)
	*/
	const Eigen::VectorXi& indices_unimportant_reactions() const { return indices_unimportant_reactions_; }

	/**
	*@brief Returns the indices of important reactions steps (zero-based)
	*/
	const Eigen::VectorXi& important_reactions_steps() const { return important_reactions_steps_; }

	/**
	*@brief Returns number of important species
	*/
	const volScalarField& important_species() const;

	/**
	*@brief Returns number of important species (non-const access)
	*/
	volScalarField& important_species();

	/**
	*@brief Returns number of important reactions
	*/
	const volScalarField& important_reactions() const;

	/**
	*@brief Returns number of important reactions (non-const access)
	*/
	volScalarField& important_reactions();

	/**
	*@brief Returns CPU time for solving chemistry
	*/
	const volScalarField& cpu_time() const;

	/**
	*@brief Returns CPU time for solving chemistry (non-const access)
	*/
	volScalarField& cpu_time();
        
private:

	/**
	*@brief Builds the pair wise error matrix
	*@param T temperature in K
	*@param P_Pa pressure in Pa
	*@param c vector of concentrations in kmol/m3 
	*@param epsilon threshold
	*/
	void PairWiseErrorMatrix(const double T, const double P_Pa, const Eigen::VectorXd& c);

	/**
	*@brief Analyzes the pair wise error matrix and calculates (internally) the important_species vector
	*/
	void ParsePairWiseErrorMatrix(const double epsilon, const double T);

private:

	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap_;		/**< reference to the thermodynamic map */
	OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap_;				/**< reference to the kinetics map */

	Eigen::VectorXi index_key_species_;					/**< indices of target or key species (zero-based) */

	Eigen::VectorXd rNet_;							/**< vector containing the net reaction rates */
	unsigned int NR_;							/**< total number of reactions */
	unsigned int NS_;							/**< total number of species */

	unsigned int number_important_species_;					/**< current number of important species */
	unsigned int number_important_reactions_;				/**< number of unimportant reactions */
	unsigned int number_unimportant_reactions_;				/**< number of unimportant reactions */

	Eigen::VectorXi indices_important_reactions_;				/**< indices of important reactions (zero-based) */
	Eigen::VectorXi indices_unimportant_reactions_;				/**< indices of unimportant reactions (zero-based) */
	Eigen::VectorXi indices_important_species_;				/**< indices of important species (zero-based) */
	Eigen::VectorXi important_reactions_steps_;				/**< indices of important reactions steps (zero-based) */

	Eigen::MatrixXd					r_;			/**< full matrix, pair wise error matrix, (NS x NS) */
	Eigen::SparseMatrix<double> 			nu_sparse_;		/**< sparse matrix containing the net stoichiometric coefficients, (NS x NR) */
	std::vector< Eigen::SparseMatrix<double> >	nu_times_delta_;	/**< vector of sparse matrices (NS x NS) */
	Eigen::SparseMatrix<double> 			delta_sparse_;		/**< delta matrix (NS x NR) */

	volScalarField important_species_;
	volScalarField important_reactions_;
	volScalarField cpu_time_;

	List<double>  epsilons_;
	List<double>  temperatures_;
	double threshold_temperature_;
	double formation_rate_threshold_;

	Eigen::VectorXd Rp_;
	Eigen::VectorXd Rd_;
	Eigen::VectorXd numerator_;
	Eigen::VectorXd denominator_;
	

};

} // end namespace Foam

#include "DRG.HPP"

#endif	
