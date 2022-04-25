/*-----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
|   License                                                               |
|                                                                         |
|   Copyright(C) 2022  Alberto Cuoci                                      |
|   OpenSMOKE++ is free software: you can redistribute it and/or modify   |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKE++ is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKE++. If not, see <http://www.gnu.org/licenses/>.   |
|                                                                         |
\*-----------------------------------------------------------------------*/

#ifndef OpenSMOKE_LookUpTable_FPI_H
#define OpenSMOKE_LookUpTable_FPI_H

// External parameters
#include <Eigen/Dense>

// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"


namespace OpenSMOKE
{
	//!  A class to manage lookup tables based on the mixture-fraction/progress-variable approach
	/*!
	A class to manage lookup tables based on the mixture-fraction/progress-variable approach
	*/

	class LookUpTable_FPI
	{

	private:

		enum LookUpTable_FPI_Type { UNIFORM_Y };
		enum Interpolated_Types { TEMPERATURE, DENSITY, SPECIFIC_HEAT, THERMAL_DIFFUSIVITY, DYNAMIC_VISCOSITY, PLANCK_ABS_COEFFICIENT, MASS_FRACTION };

	public:

		/**
		*@brief Default constructor
		*@param thermodynamicsMap	reference to the thermodynamic map
		*/
		LookUpTable_FPI(OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap);

		/**
		*@brief Setup from a dictionary
		*@param dictionary dictionary name
		*/
		void SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary);

		/**
		*@brief Setup from user provided path
		*@param path_to_xml path to the CSV lookup table
		*@param is_apriori if true, a priori analysis is performed
		*/
		void Setup(const boost::filesystem::path& path_to_csv, const bool is_apriori);

		/**
		*@brief Interpolate
		*@param c progress-variable (-)
		*/
		void Interpolate(const double c);

		/**
		*@brief Interpolate from normalized progress variable
		*@param ctilde normalized progress-variable (-)
		*/
		void InterpolateFromNormalizedProgressVariable(const double ctilde);

		/**
		*@brief Returns the j-th interpolated mass fraction (according to the order reported in the list_fields variable)
		*@param j index of required interpolated value
		*@return the j-th interpolated mass fraction
		*/
		double Y(const unsigned int j) const { return interpolated_Y_[j]; }

		/**
		*@brief Returns the interpolated variables
		*@return the interpolated variable
		*/
		double rho() const { return interpolated_rho_; }
		double Cp() const { return interpolated_Cp_; }
		double lambda() const { return interpolated_lambda_; }
		double mu() const { return interpolated_mu_; }
		double MW() const { return interpolated_MW_; }
		double Q() const { return interpolated_Q_; }
		double a() const { return interpolated_a_; }
		double T() const { return interpolated_T_; }
		double OmegaC() const { return interpolated_OmegaC_; }
		double OmegaCtilde() const { return interpolated_OmegaCtilde_; }

	
	private:

		/**
		*@brief Read lookup table from XML file
		*/
		void ReadTableFromCSVFile();

		/**
		*@brief Check input options
		*/
		void CheckInput();

		/**
		*@brief Summary on the screen
		*/
		void Summary();


	private:

		bool is_active_;		//!< true only if the table has been activated
		bool is_apriori_;		//!< a priori analysis

		OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap_;		//!< thermodynamics map

		LookUpTable_FPI_Type type_;		//!< type of dynamic boundary condition

		boost::filesystem::path path_output_;	//!< path to output folder

		unsigned int ns_;
		unsigned int nc_;				//<! number of points along progress variable field

		Eigen::VectorXd c_;				//<! progress-variable coordinates
		Eigen::VectorXd ctilde_;			//<! normalized progress-variable coordinates
		Eigen::VectorXd omega_c_;			//<! progress-variable coordinates
		Eigen::VectorXd omega_ctilde_;		//<! normalized progress-variable coordinates

		double c_min_;			//<! minimum progress variable at each mixture fraction
		double c_max_;			//<! maximum progress variable at each mixture fraction

		double dc_;					//<! progress-variable spacing (only in case of uniform grid)
		double dctilde_;				//<! normalized progress-variable spacing (only in case of uniform grid)

		Eigen::VectorXd rho_;
		Eigen::VectorXd Cp_;
		Eigen::VectorXd lambda_;
		Eigen::VectorXd mu_;
		Eigen::VectorXd MW_;
		Eigen::VectorXd Q_;
		Eigen::VectorXd a_;
		Eigen::VectorXd T_;

		Eigen::MatrixXd Y_;	//<! lookup tables

		double query_c_;		//!< current progress-variable for which the interpolation was carried out

		unsigned int jC_;		//!< index of C atom in the matrix containing atomic composition of species
		unsigned int jO_;		//!< index of O atom in the matrix containing atomic composition of species
		unsigned int jH_;		//!< index of H atom in the matrix containing atomic composition of species

		double WC_;		//!< molecular weight of C atom (in kg/kmol)
		double WO_;		//!< molecular weight of O atom (in kg/kmol)
		double WH_;		//!< molecular weight of H atom (in kg/kmol)

		double interpolated_rho_;
		double interpolated_Cp_;
		double interpolated_lambda_;
		double interpolated_mu_;
		double interpolated_MW_;
		double interpolated_Q_;
		double interpolated_a_;
		double interpolated_T_;
		double interpolated_OmegaC_;
		double interpolated_OmegaCtilde_;

		Eigen::VectorXd interpolated_Y_;
		
	};
}

#include "LookUpTable_FPI.hpp"

#endif /* OpenSMOKE_LookUpTable_FPI_H */

