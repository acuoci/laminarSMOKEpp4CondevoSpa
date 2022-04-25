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
|   Copyright(C) 2020 Alberto Cuoci                                       |
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

#pragma once

#include "ClusteringAlgorithms.H"

class BinningAlgorithm : public ClusteringAlgorithms
{
public:

	/**
	*@brief Default constructor
	*/
	BinningAlgorithm();

	/**
	*@brief Set the number of features
	*@param D number of features (min 1, max 6)
	*/
	virtual void SetNumberOfFeatures(const unsigned int D);

	/**
	*@brief Set the vector of data corresponding to one of the features on which to cluster
	*@param index index of the feature (from 0)
	*@param psi vector of data
	*@param epsilon accuracy level
	*@param relative flag to specify if accuracy is given in relative (true) or absolute (false) units
	*/
	virtual void SetData(const unsigned int index, const std::vector<double>& psi, const double epsilon, const bool relative);

	/**
	*@brief Cluster the data provided by the user via SetData function
	*/
	virtual void ClusterData();


protected:

	/**
	*@brief Check the data provided by the user
	*/
	virtual void CheckInputData();

	/**
	*@brief Fatal error message
	*/
	virtual void FatalErrorMessage(const std::string message);


private:

	std::vector<int> L_;								//!< number of intervals along each direction
	std::vector< std::vector<double> > psi_tilde_;		//!< normalized feature
	std::vector<double> psi_min_;						//!< minimum value of (un-normalized) feature
	std::vector<double> psi_max_;						//!< maximum value of (un-normalized) feature
	std::vector< std::vector<int> >	idx_;				//!< 1D interval binning index
};
