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

#define NODEBUG_PCAEIGEN

#include <iostream>
#ifdef DEBUG_PCAEIGEN
#include <iterator>
#endif
#include <Eigen/SVD>
#include "pcaEigen.H"


PCAEigen::PCAEigen() 
{
	method_ = PCA_SVD;
	is_center_ = true;
	is_scale_ = true;
}

PCAEigen::~PCAEigen() 
{ 
}

void PCAEigen::is_center(const bool flag)
{
	is_center_ = flag;
}

void PCAEigen::is_scale(const bool flag)
{
	is_scale_ = flag;
}

bool PCAEigen::is_scale() const
{  
	return is_scale_; 
}

bool PCAEigen::is_center() const
{ 
	return is_center_; 
}

const Eigen::VectorXd& PCAEigen::mu() const
{
	return mu_; 
}

const Eigen::VectorXd& PCAEigen::sigma() const
{ 
	return sigma_; 
}
   
const Eigen::VectorXd& PCAEigen::explained() const 
{ 
	return explained_; 
}

const Eigen::MatrixXd& PCAEigen::cosines() const 
{ 
	return cosines_; 
}

const Eigen::MatrixXd& PCAEigen::weights() const 
{ 
	return weights_; 
}

const Eigen::VectorXi& PCAEigen::eliminated_columns() const
{ 
	return eliminated_columns_; 
}

PCAEigen::PCAMethod PCAEigen::method() const
{ 
	return method_; 
}

unsigned int PCAEigen::explained_095() const
{
	return explained_095_;
}

unsigned int PCAEigen::explained_098() const
{
	return explained_098_;
}

unsigned int PCAEigen::explained_099() const
{
	return explained_099_;
}

int PCAEigen::Calculate(const Eigen::MatrixXd& X)
{
	// Number of variables
	unsigned int ncols = X.cols();

	// Number of observations
	unsigned int nrows = X.rows();

	// Mean for each column
	mu_.resize(ncols);
	mu_ = X.colwise().mean();

	// Standard deviation for each column
	sigma_.resize(ncols);
	unsigned int number_zero_sigma = 0;
	const double denom = static_cast<double>((nrows > 1)? nrows-1: 1);
	std::vector<unsigned int> zero_sigma_variables;
	for (unsigned int i = 0; i < ncols; ++i) 
	{
		Eigen::VectorXd col = Eigen::VectorXd::Constant(nrows, mu_(i)); 	// mu for column x
		col = X.col(i) - col; 							// x - mu
		col = col.array().square(); 						// (x-mu)^2  
		sigma_(i) = std::sqrt( (col.sum())/denom );
		if (std::fabs(sigma_(i)) <= 1.e-16) 
		{
			number_zero_sigma++;
			zero_sigma_variables.push_back(i+1);
		}
	}

	// Check the number of columns having zero variance
	if (number_zero_sigma > 0)
	{
		std::cout << "WARNING: number of variables having zero variance is: " << number_zero_sigma << std::endl;
		for (unsigned int i=0;i<zero_sigma_variables.size();i++)
			std::cout << " * Variable " << zero_sigma_variables[i] << " (starting from 1)" << std::endl;
	}

	// Allocate memory for normalized data matrix
	Eigen::MatrixXd Xstar(nrows, ncols);

	// Shift to zero
	if (is_center_ == true) 
	{
		for (unsigned int i = 0; i < ncols; ++i) 
			Xstar.col(i) = X.col(i) - Eigen::VectorXd::Constant(nrows, mu_(i));
	}

	// Scale to unit variance
	if ( is_scale_ == true ) 
	{
		for (unsigned int i = 0; i < ncols; ++i) 
			Xstar.col(i) /= std::sqrt(Xstar.col(i).array().square().sum()/denom);
	}

	//SVD decomposition
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(Xstar, Eigen::ComputeThinV);

	// Calculate explained variance
	{
		Eigen::VectorXd s = svd.singularValues();
		Eigen::VectorXd s2 = s.array().square();
		const double sum_s2 = s2.sum();
		explained_ = s2/sum_s2;

		explained_095_ = 0;
		explained_098_ = 0;
		explained_099_ = 0;

		double sum = 0.;
		for (unsigned int i=0;i<explained_.size();i++)
		{	
			sum += explained_(i);

			if (sum >= 0.95 && explained_095_ == 0)
				explained_095_ = i;
			if (sum >= 0.98 && explained_098_ == 0)
				explained_098_ = i;
			if (sum >= 0.99 && explained_099_ == 0)
			{
				explained_099_ = i;
				break;
			}
		}
	}

	// Get weights
	weights_ = svd.matrixV();

	return 0;
}

void PCAEigen::Scores(Eigen::VectorXd& X, const unsigned int n, Eigen::VectorXd& scores) const
{
	// Number of variables
	unsigned int ncols = X.size();

	// Shift and scaling (most common choice)
	if (is_center_ == true && is_scale_ == true)
	{
		for (unsigned int i = 0; i < ncols; ++i) 
			X(i) = (X(i)-mu_(i))/sigma_(i);

	}
	else
	{
		// Shift to zero
		if (is_center_ == true) 
		{
			for (unsigned int i = 0; i < ncols; ++i) 
				X(i) -= mu_(i);
		}

		// Scale to unit variance
		if (is_scale_ == true) 
		{
			for (unsigned int i = 0; i < ncols; ++i) 
				X(i) /= sigma_(i);
		}
	}

	// Calculates scores
	scores = X.transpose()*weights_.leftCols(n);
}

void PCAEigen::CalculateCosines()
{
	const unsigned int ncols = weights_.cols();
	cosines_.resize(ncols, ncols); 
	for (unsigned int i = 0; i < ncols; ++i)
	{
		const double norm = weights_.col(i).norm();
		for (unsigned int j = 0; j < ncols; ++j)
			cosines_(i,j) = weights_(j,i)/norm;
	}
}

/*
int PCAEigen::Calculate(const std::vector<double> &x, const unsigned int &nrows, const unsigned int &ncols, const bool is_corr, const bool is_center, const bool is_scale) 
{
	_ncols = ncols;
	_nrows = nrows;
	_is_corr = is_corr;
	_is_center = is_center;
	_is_scale = is_scale;

	if (x.size()!= _nrows*_ncols) 
	{
		return -1;
	}
	if ((1 == _ncols) || (1 == nrows)) 
	{
		return -1;
	}
	
	// Convert vector to Eigen 2-dimensional matrix
	//Map<Eigen::MatrixXd> _xXd(x.data(), _nrows, _ncols);
	_xXd.resize(_nrows, _ncols);
	for (unsigned int i = 0; i < _nrows; ++i) 
	{
		for (unsigned int j = 0; j < _ncols; ++j) 
		{
			_xXd(i, j) = x[j + i*_ncols];
		}
	}

	// Mean and standard deviation for each column
	Eigen::VectorXd mean_vector(_ncols);
	mean_vector = _xXd.colwise().mean();
	Eigen::VectorXd sd_vector(_ncols);
	unsigned int zero_sd_num = 0;
	const double denom = static_cast<double>((_nrows > 1)? _nrows - 1: 1);
	for (unsigned int i = 0; i < _ncols; ++i) 
	{
		Eigen::VectorXd curr_col  = Eigen::VectorXd::Constant(_nrows, mean_vector(i)); // mean(x) for column x
		curr_col = _xXd.col(i) - curr_col; // x - mean(x)
		curr_col = curr_col.array().square(); // (x-mean(x))^2  
		sd_vector(i) = sqrt((curr_col.sum())/denom);
		if (0 == sd_vector(i)) 
		{
			zero_sd_num++;
		}
	}

	// If colums with sd == 0 are too many, don't continue calculations
	if (1 > _ncols-zero_sd_num) 
	{
		return -1;
	}

	// Delete columns where sd == 0
	Eigen::MatrixXd tmp(_nrows, _ncols-zero_sd_num);
	Eigen::VectorXd tmp_mean_vector(_ncols-zero_sd_num);
	unsigned int curr_col_num = 0;
	for (unsigned int i = 0; i < _ncols; ++i) 
	{
		if (0 != sd_vector(i)) 
		{
			tmp.col(curr_col_num) = _xXd.col(i);
			tmp_mean_vector(curr_col_num) = mean_vector(i);
			curr_col_num++;
		} 
		else 
		{
			_eliminated_columns.push_back(i);
		}
	}

	_ncols -= zero_sd_num;
	_xXd = tmp;
	mean_vector = tmp_mean_vector;
	tmp.resize(0, 0); tmp_mean_vector.resize(0);
	
	// Shift to zero
	if (true == _is_center) 
	{
		for (unsigned int i = 0; i < _ncols; ++i) 
		{
			_xXd.col(i) -= Eigen::VectorXd::Constant(_nrows, mean_vector(i));
		}
	}

	// Scale to unit variance
	if ( (false == _is_corr) || (true == _is_scale)) 
	{
		for (unsigned int i = 0; i < _ncols; ++i) 
		{
			_xXd.col(i) /= sqrt(_xXd.col(i).array().square().sum()/denom);
		}
	}

	#ifdef DEBUG_PCAEIGEN
		std::cout << "\nScaled matrix:\n";
		std::cout << _xXd << std::endl;
		std::cout << "\nMean before scaling:\n" << mean_vector.transpose();
		std::cout << "\nStandard deviation before scaling:\n" << sd_vector.transpose();
	#endif

	// When _nrows < _ncols then svd will be used.
	// If corr is true and _nrows > _ncols then will be used correlation matrix
	// (TODO): What about covariance?
	if ( (_nrows < _ncols) || (false == _is_corr)) 
	{ 
		// Singular Value Decomposition is on
		_method = "svd";
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(_xXd, Eigen::ComputeThinV);
		Eigen::VectorXd eigen_singular_values = svd.singularValues();
		Eigen::VectorXd tmp_vec = eigen_singular_values.array().square();
		double tmp_sum = tmp_vec.sum();
		tmp_vec /= tmp_sum;
		
		// PC's standard deviation and
		// PC's proportion of variance
		_kaiser = 0;
		unsigned int lim = (_nrows < _ncols)? _nrows : _ncols;
		for (unsigned int i = 0; i < lim; ++i) 
		{
			_sd.push_back(eigen_singular_values(i)/sqrt(denom));
			if (_sd[i] >= 1) 
			{
				_kaiser = i + 1;
			}
			_prop_of_var.push_back(tmp_vec(i));
		}

		#ifdef DEBUG_PCAEIGEN
			std::cout << "\n\nStandard deviations for PCs:\n";
			copy(_sd.begin(), _sd.end(),std::ostream_iterator<double>(std::cout," "));  
			std::cout << "\n\nKaiser criterion: PC #" << _kaiser << std::endl;
		#endif

		tmp_vec.resize(0);
		
		// PC's cumulative proportion
		_thresh95 = 1;
		_cum_prop.push_back(_prop_of_var[0]); 
		for (unsigned int i = 1; i < _prop_of_var.size(); ++i) 
		{
			_cum_prop.push_back(_cum_prop[i-1]+_prop_of_var[i]);
			if (_cum_prop[i] < 0.95) 
			{
				_thresh95 = i+1;
			}
		}

		#ifdef DEBUG_PCAEIGEN
			std::cout << "\nCumulative proportion:\n";
			copy(_cum_prop.begin(), _cum_prop.end(),std::ostream_iterator<double>(std::cout," "));  
			std::cout << "\n\nThresh95 criterion: PC #" << _thresh95 << std::endl;
		#endif

		// Scores
		Eigen::MatrixXd eigen_scores = _xXd * svd.matrixV();

		#ifdef DEBUG_PCAEIGEN
			std::cout << "\n\nRotated values (scores):\n" << eigen_scores;
		#endif

		_scores.reserve(lim*lim);
		for (unsigned int i = 0; i < lim; ++i) 
		{
			for (unsigned int j = 0; j < lim; ++j) 
			{
				_scores.push_back(eigen_scores(i, j));
			}
		}

		eigen_scores.resize(0, 0);

		#ifdef DEBUG_PCAEIGEN
			std::cout << "\n\nScores in vector:\n";
			copy(_scores.begin(), _scores.end(),std::ostream_iterator<double>(std::cout," "));  
			std::cout << "\n";  
		#endif
	} 

	else 
	{ 
		// COR OR COV MATRICES ARE HERE
		_method = "cor";
		
		// Calculate covariance matrix
		Eigen::MatrixXd eigen_cov; // = Eigen::MatrixXd::Zero(_ncols, _ncols);
		Eigen::VectorXd sds;
		
		// (TODO) Should be weighted cov matrix, even if is_center == false
		eigen_cov = (1.0 /(_nrows-1 )) * _xXd.transpose() * _xXd;
		sds = eigen_cov.diagonal().array().sqrt();
		Eigen::MatrixXd outer_sds = sds * sds.transpose();
		eigen_cov = eigen_cov.array() / outer_sds.array();
		outer_sds.resize(0, 0);

		// ?If data matrix is scaled, covariance matrix is equal to correlation matrix
		
		#ifdef DEBUG_PCAEIGEN
			std::cout << eigen_cov << std::endl;
		#endif

		Eigen::EigenSolver<Eigen::MatrixXd> edc(eigen_cov);
		Eigen::VectorXd eigen_eigenvalues = edc.eigenvalues().real();
		
		#ifdef DEBUG_PCAEIGEN
			std::cout << endl << eigen_eigenvalues.transpose() << std::endl;
		#endif
		
		Eigen::MatrixXd eigen_eigenvectors = edc.eigenvectors().real();	

		#ifdef DEBUG_PCAEIGEN
			std::cout << endl << eigen_eigenvectors << std::endl;
		#endif

		// The eigenvalues and eigenvectors are not sorted in any particular order.
		// So, we should sort them
		typedef std::pair<double, int> eigen_pair;
		std::vector<eigen_pair> ep;	
		for (unsigned int i = 0 ; i < _ncols; ++i) 
		{
			ep.push_back(std::make_pair(eigen_eigenvalues(i), i));
		}
		sort(ep.begin(), ep.end()); // Ascending order by default

		// Sort them all in descending order
		Eigen::MatrixXd eigen_eigenvectors_sorted = Eigen::MatrixXd::Zero(eigen_eigenvectors.rows(), eigen_eigenvectors.cols());
		Eigen::VectorXd eigen_eigenvalues_sorted = Eigen::VectorXd::Zero(_ncols);
		int colnum = 0;
		int i = ep.size()-1;
		for (; i > -1; i--) 
		{
			eigen_eigenvalues_sorted(colnum) = ep[i].first;
			eigen_eigenvectors_sorted.col(colnum++) += eigen_eigenvectors.col(ep[i].second);
		}

		#ifdef DEBUG_PCAEIGEN
			std::cout << endl << eigen_eigenvalues_sorted.transpose() << std::endl;
			std::cout << endl << eigen_eigenvectors_sorted << std::endl;
		#endif  

		// We don't need not sorted arrays anymore
		eigen_eigenvalues.resize(0);
		eigen_eigenvectors.resize(0, 0);
		
		_sd.clear(); 
		_prop_of_var.clear(); 
		_kaiser = 0;
		double tmp_sum = eigen_eigenvalues_sorted.sum();
		for (unsigned int i = 0; i < _ncols; ++i) 
		{
			_sd.push_back(sqrt(eigen_eigenvalues_sorted(i)));
			if (_sd[i] >= 1) 
			{
				_kaiser = i + 1;
			}
			_prop_of_var.push_back(eigen_eigenvalues_sorted(i)/tmp_sum);
		}

		#ifdef DEBUG_PCAEIGEN
			std::cout << "\nStandard deviations for PCs:\n";
			copy(_sd.begin(), _sd.end(), std::ostream_iterator<double>(std::cout," "));  
			std::cout << "\nProportion of variance:\n";
			copy(_prop_of_var.begin(), _prop_of_var.end(), std::ostream_iterator<double>(std::cout," ")); 
			std::cout << "\nKaiser criterion: PC #" << _kaiser << std::endl;
		#endif

		// PC's cumulative proportion
		_cum_prop.clear(); _thresh95 = 1;
		_cum_prop.push_back(_prop_of_var[0]);
		for (unsigned int i = 1; i < _prop_of_var.size(); ++i) 
		{
			_cum_prop.push_back(_cum_prop[i-1]+_prop_of_var[i]);
			if (_cum_prop[i] < 0.95) 
			{
				_thresh95 = i+1;
			}
		}  

		#ifdef DEBUG_PCAEIGEN
			std::cout << "\n\nCumulative proportions:\n";
			copy(_cum_prop.begin(), _cum_prop.end(), std::ostream_iterator<double>(std::cout," "));  
			std::cout << "\n\n95% threshold: PC #" << _thresh95 << std::endl;
		#endif

		// Scores for PCA with correlation matrix
		// Scale before calculating new values
		for (unsigned int i = 0; i < _ncols; ++i) 
		{
			_xXd.col(i) /= sds(i);
		}

		sds.resize(0);
		Eigen::MatrixXd eigen_scores = _xXd * eigen_eigenvectors_sorted;
		
		#ifdef DEBUG_PCAEIGEN
			std::cout << "\n\nRotated values (scores):\n" << eigen_scores;
		#endif
		
		_scores.clear();
		_scores.reserve(_ncols*_nrows);
		for (unsigned int i = 0; i < _nrows; ++i) 
		{
			for (unsigned int j = 0; j < _ncols; ++j) 
			{
				_scores.push_back(eigen_scores(i, j));
			}
		}

		eigen_scores.resize(0, 0);

		#ifdef DEBUG_PCAEIGEN
			std::cout << "\n\nScores in vector:\n";
			copy(_scores.begin(), _scores.end(), std::ostream_iterator<double>(std::cout," "));  
			std::cout << "\n";  
		#endif
	}

	return 0;
}
*/







