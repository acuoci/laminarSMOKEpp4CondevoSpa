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
|   Copyright(C) 2021 Alberto Cuoci                                       |
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

// This class is adapted from the EigenPCA library by 
// ihar.safonau@gmail.com available on github at the following 
// address: https://github.com/ihar/EigenPCA

#pragma once

#include <vector>
#include <Eigen/Dense>

class PCAEigen
{

public:

	enum PCAMethod { PCA_SVD };

public:

 	/**
	*@brief 	Default constructor
	*/
  	PCAEigen();


 	/**
	*@brief 	Default destructor
	*/
  	~PCAEigen();


	/**
	*@brief 	Calculate PCA
	*@param 	X data matrix (m x n), where m is the number of observations and n is the number of variables
	*@return	returns 0 is calculation was OK
	*/
	int Calculate(const Eigen::MatrixXd& X);


	/**
	*@brief 	Calculate scores (i.e. PCs)
	*@param 	X data matrix (M x N), where M is the number of observations and N is the number of variables
	*@param 	n the number of retained variables (n<=N)
	*@param 	the score matrix, with size equal to (M x n)
	*/
	void Scores(Eigen::VectorXd& X, const unsigned int n, Eigen::VectorXd& scores) const;


	/**
	*@brief 	Calculate the cosines
	*@param 	Calculates the cosines
	*/
	void CalculateCosines();


	/**
	*@brief 	If variables are centered around the mean
	*@param 	flag: if true, variables are centered around the mean
	*/
  	void is_center(const bool flag);


	/**
	*@brief 	If variables are scaled with respect to the std deviation
	*@param 	flag: if true, variables are scaled with respect to the std deviation
	*/
  	void is_scale(const bool flag);


	/**
	*@brief 	If variables are scaled
	*@return 	true: variables are scaled, false: scaled
	*/
  	bool is_scale() const;	


	/**
	*@brief 	If variables are centered
	*@return 	true: variables are centered, false: otherwise
	*/
  	bool is_center() const;


	/**
	*@brief 	Method for calculation of principal components
			There are different methods used. The most used is SVD. 
			But in some cases it may be correlation or covariance matrices.
	*@return 	"svd" - PCA with singular value decomposition
	*/
  	PCAMethod method() const;


	/**
	*@brief 	Returns numbers of eliminated columns
			If standard deviation of a column is equal to 0, the column shoud be rejected, or PCA will fail.
	*@return 	Numbers of eliminated columns, empty vector otherwise
	*/
  	const Eigen::VectorXi& eliminated_columns() const;

	/**
	*@brief 	Mean of each variable
	*@return 	Vector of mean for each variable
	*/
  	const Eigen::VectorXd& mu() const;  

	/**
	*@brief 	Standard deviation of each variable
	*@return 	Vector of standard deviation for each variable
	*/
  	const Eigen::VectorXd& sigma() const;    


	/**
	*@brief 	Proportion of variance
	*@return 	Vector of variances for each component
	*/
  	const Eigen::VectorXd& explained() const;


	/**
	*@brief 	Cosines between the columns of weights matrix and the variables
	*@return 	Returns the cosines between the columns of weights matrix and the variables
			The row refers to the principal component, the column to the variable
	*/
  	const Eigen::MatrixXd& cosines() const;


	/**
	*@brief 	Returns the PCA weight matrix
			Return the weight matrix 
	*@return 	weight matrix 
	*/
  	const Eigen::MatrixXd& weights() const;

	/**
	*@brief 	Returns the minimum number of variables to explain variance (95%)
	*@return 	minimum number of variables to explain variance (95%)
	*/
	unsigned int explained_095() const ;

	/**
	*@brief 	Returns the minimum number of variables to explain variance (98%)
	*@return 	minimum number of variables to explain variance (98%)
	*/
	unsigned int explained_098() const ;

	/**
	*@brief 	Returns the minimum number of variables to explain variance (99%)
	*@return 	minimum number of variables to explain variance (99%)
	*/
	unsigned int explained_099() const ;


private:
  
  	bool is_center_;         	//!< Whether the variables should be shifted to be zero centered
        bool is_scale_;          	//!< Whether the variables should be scaled to have unit variance
  
	PCAMethod method_;        			//!< PCA method: PCA_SVD
  	Eigen::VectorXi eliminated_columns_;  		//!< Numbers of eliminated columns
	Eigen::VectorXd mu_;  				//!< Mean of each component
  	Eigen::VectorXd sigma_;  			//!< Standard deviation of each component
        Eigen::VectorXd explained_;      		//!< Expalained variance
        Eigen::MatrixXd weights_;        		//!< PCA weight matrix
	Eigen::MatrixXd cosines_;      			//!< Cosines
		
	unsigned int explained_095_;		//!< minimum number of variables to explain variance (95%)
	unsigned int explained_098_;		//!< minimum number of variables to explain variance (98%)
	unsigned int explained_099_;		//!< minimum number of variables to explain variance (99%)
};

