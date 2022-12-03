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

#ifndef chemistryLinearModel_H
#define chemistryLinearModel_H

// Include standard OpenFOAM files
#include "volFields.H"
#include "typeInfo.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "wordIOList.H"
#include "hashedWordList.H"

#if CHEMKINBASED == 1
#include "ckwyp.H"
#endif

namespace Foam
{

enum species_order_policy_enum {SPECIES_ORDER_POLICY_CONSTANT, SPECIES_ORDER_POLICY_SWEEP, SPECIES_ORDER_POLICY_RANDOM_SHUFFLE, SPECIES_ORDER_POLICY_ROTATE };

class ChemistryLinearModel : public IOdictionary
{
public: 

	/**
	*@brief Default constructor
	*@param mesh mesh
	*@param thermodynamicsMapXML thermodynamic map
	*@param kineticsMapXML kinetics map
	*/
	ChemistryLinearModel(const fvMesh& mesh, OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap);

	/**
	*@brief Default destructor
	*/
	virtual ~ChemistryLinearModel();

	/**
	*@brief Name of main dictionary
	*/
	static const word dictName;

	/**
	*@brief Returns the number of iterations at which the Jacobian matrix has to be updated
	*/
	label jacobianUpdate() const { return jacobianUpdate_; }

	/**
	*@brief Returns the number of iterations at which the mixture properties have to be updated
	*/
	label propertiesUpdate() const { return propertiesUpdate_; }

	/**
	*@brief Returns the current number of iterations for tracking the update of Jacobian matrix
	*/
	label jacobianCounter() const { return jacobianCounter_; }

	/**
	*@brief Returns the current number of iterations for tracking the update of mixture properties
	*/
	label propertiesCounter() const { return propertiesCounter_; }

	/**
	*@brief Returns true if the chemical source terms are treated implicitly (i.e. if the Jacobian calculation is required)
	*/
	Switch implicitSourceTerm() const { return implicitSourceTerm_; }

	/**
	*@brief Returns true if the sparsity of kinetic mechanism is exploited for building the Jacobian matrix
	*/
	Switch sparseJacobian() const { return sparseJacobian_; }

	/**
	*@brief Updates the counter for updating the Jacobian matrix
	*/
	void updateJacobianCounter();

	/**
	*@brief Updates the counter for updating the mixture properties
	*/
	void updatePropertiesCounter();

	/**
	*@brief Updates the sequence of species to be solved
	*/
	void updateSpeciesOrderPolicy();

	/**
	*@brief Returns the index of species to be solved according to the current policy
	*@param j the position in the current policy
	*/
	label species_order(const label j) const { return species_order_[j]; } 

	/**
	*@brief Calculates the chemical source terms explicitly
	*@param y the current vector of species mass fractions and temperature
	*@param P0 pressure (in Pa)
	*@param S the vector of source terms (mass fractions and temperature)
	*/
	void reactionSourceTerms( const Eigen::VectorXd& y, const double P0, Eigen::VectorXd& S );

	/**
	*@brief Calculates the Jacobian elements without exploiting the sparsity of kinetic mechanism
	*@param y the current vector of species mass fractions and temperature
	*@param P0 pressure (in Pa)
	*@param J the Jacobian elements along the main diagonal
	*/	
	void reactionJacobian( const Eigen::VectorXd& y, const double P0, Eigen::VectorXd &J);

	/**
	*@brief Calculates the Jacobian elements by exploiting the sparsity of kinetic mechanism
	*@param y the current vector of species mass fractions and temperature
	*@param P0 pressure (in Pa)
	*@param J the Jacobian elements along the main diagonal
	*@param energyEquation true in case the energy equation is turned on
	*@param energyEquation true in case the species equations are turned on
	*/
	void reactionJacobianSparse( const Eigen::VectorXd& y, const double P0, Eigen::VectorXd &J, const bool energyEquation, const bool speciesEquations);

private:

	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap_;		/**< reference to the thermodynamic map */
	OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap_;				/**< reference to the kinetics map */

	unsigned int NC_;			/**< number of species */
	unsigned int NE_; 			/**< number of equations */

	Eigen::VectorXd omega_;			/**< mass fractions */
	Eigen::VectorXd x_;			/**< mole fractions */
	Eigen::VectorXd R_;			/**< formation rates */
	Eigen::VectorXd c_;			/**< concentrations */

	Eigen::VectorXd y_plus_;		/**< auxiliary vector for storing the variables during the Jacobian calculation */
     	Eigen::VectorXd dy_plus_;		/**< auxiliary vector for storing the variables during the Jacobian calculation */
     	Eigen::VectorXd dy_original_;		/**< auxiliary vector for storing the variables during the Jacobian calculation */

	Eigen::VectorXd Jdiagonal_;		/**< auxiliary vector for storing the variables during the Jacobian calculation */

	label jacobianUpdate_;			/**< number of iterations at which the Jacobian matrix has to be updated */
	label propertiesUpdate_;		/**< number of iterations at which the mixture properties to be updated */
	label jacobianCounter_;			/**< current number of iterations for tracking the update of Jacobian matrix */
	label propertiesCounter_;		/**< current number of iterations for tracking the update of mixture properties */

	Switch implicitSourceTerm_;		/**< true if the chemical source terms are treated implicitly (i.e. if the Jacobian calculation is required) */
	Switch sparseJacobian_;			/**< true if the sparsity of kinetic mechanism is exploited for building the Jacobian matrix */

	species_order_policy_enum species_order_policy_;	/**< policy for defining the sequence of species to be solved */
	std::vector<int> species_order_;			/**< sequence of species to be solved */

};

} // end namespace Foam

#include "ChemistryLinearModel.HPP"

#endif	
