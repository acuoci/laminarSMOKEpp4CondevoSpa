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

#include "BinningAlgorithm.H"

#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <numeric>
#include <iomanip>

BinningAlgorithm::BinningAlgorithm() :
	ClusteringAlgorithms()
{
}


void BinningAlgorithm::SetNumberOfFeatures(const unsigned int D)
{
	if (D >= 6)
		FatalErrorMessage("The maximum number of features is 6");

	D_ = D;

	for (unsigned int i = 0; i < D_; i++)
	{
		L_.resize(D_);
		epsilon_.resize(D_);
		psi_.resize(D_);
		psi_tilde_.resize(D_);
		psi_min_.resize(D_);
		psi_max_.resize(D_);
		idx_.resize(D_);
		C_.resize(D_);
		sigma_.resize(D_);
		C_avg_.resize(D_);
		sigma_avg_.resize(D_);
	}
}

void BinningAlgorithm::SetData(const unsigned int index, const std::vector<double>& psi, const double epsilon, const bool relative)
{
	if (verbosity_>0)
		std::cout << " * Setting variable #" << index << std::endl;

	if (index >= D_)
		FatalErrorMessage("The provided number of features is not sufficient to accomodate the input vector of data");
		
	n_ = psi.size();
	idx_[index].resize(n_);
	psi_[index] = psi;

	const double threshold = 1e-6;
	auto result  = std::minmax_element(std::begin(psi), std::end(psi));
	psi_min_[index] = *result.first;
	psi_max_[index] = *result.second;
	const double psi_min = *result.first - threshold;
	const double psi_max = *result.second + threshold;

	// Check if the accuracy is provided in relative or absolute units
	// The binning algorithm works with relative accuracy
	if (relative == false)
		epsilon_[index] = epsilon / (psi_max_[index] - psi_min_[index]);
	else
		epsilon_[index] = epsilon;

	// Correct epsilon to be sure to have an integer number of intervals
	L_[index] = static_cast<int>(1. / epsilon_[index]);
	epsilon_[index] = static_cast<double>(1. / L_[index]);

	// Dimensionless maps
	psi_tilde_[index].resize(n_);
	for (unsigned int i = 0; i < n_; i++)
		psi_tilde_[index][i] = (psi[i] - psi_min) / (psi_max - psi_min);

	for (unsigned int i = 0; i < n_; i++)
		idx_[index][i] = static_cast<int>( std::ceil(psi_tilde_[index][i] / epsilon_[index]) );

	if (verbosity_ > 0)
	{
		std::cout << "   - Min: " << psi_min_[index] << " Max: " << psi_max_[index] << " delta: " << psi_max_[index] - psi_min_[index] << std::endl;
		std::cout << "   - Epsilon: " << epsilon_[index] << " Absolute epsilon: " << epsilon_[index]*(psi_max_[index] - psi_min_[index]) << std::endl;
	}
}

void BinningAlgorithm::ClusterData()
{
	if (verbosity_ > 0)
		std::cout << " * Clustering data via binning algorithm..." << std::endl;

	CheckInputData();

	auto cpu_start = std::chrono::steady_clock::now();

	std::vector<int> j_(n_);

	for (unsigned int k = 0; k < n_; k++)
		j_[k] = idx_[0][k];

	if (D_ >= 2)
		for (unsigned int k = 0; k < n_; k++)
			j_[k] += idx_[1][k] * L_[0];

	if (D_ >= 3)
		for (unsigned int k = 0; k < n_; k++)
			j_[k] += idx_[2][k] * L_[0] * L_[1];

	if (D_ >= 4)
		for (unsigned int k = 0; k < n_; k++)
			j_[k] += idx_[3][k] * L_[0] * L_[1] * L_[2];

	if (D_ >= 5)
		for (unsigned int k = 0; k < n_; k++)
			j_[k] += idx_[4][k] * L_[0] * L_[1] * L_[2] * L_[3];

	if (D_ >= 6)
		for (unsigned int k = 0; k < n_; k++)
			j_[k] += idx_[5][k] * L_[0] * L_[1] * L_[2] * L_[3] * L_[4];


	// Identify unique elements
	std::vector<int> h = j_;
	std::sort(h.begin(), h.end());
	h.erase(std::unique(h.begin(), h.end()), h.end());

	// Number of clusters
	nc_ = h.size();

	// Identify members of clusters
	gj_.resize(n_);
	lj_.resize(nc_);
	for (unsigned int i = 0; i < nc_; i++)
		lj_[i].reserve(static_cast<int>(5 * n_ / nc_));
	for (unsigned int k = 0; k < n_; k++)
		for (unsigned int i = 0; i < nc_; i++)
		{
			if (j_[k] == h[i])
			{
				lj_[i].push_back(k);
				gj_[k] = i;
				break;
			}
		}

	// Number of elements in each cluster
	ne_.resize(nc_);
	for (unsigned int i = 0; i < nc_; i++)
		ne_[i] = lj_[i].size();

	auto cpu_end = std::chrono::steady_clock::now();

	std::chrono::duration<double> cpu_clustering = cpu_end - cpu_start;
	cpu_clustering_ = cpu_clustering.count();
}

void BinningAlgorithm::CheckInputData()
{
	for (unsigned int i = 0; i < D_; i++)
		if (psi_tilde_[i].size() != n_)
			FatalErrorMessage("Wrong input data: the sizes of input vectors are not consistent");
}

void BinningAlgorithm::FatalErrorMessage(const std::string message)
{
	std::cout << "BinningAlgorithm: Fatal error" << std::endl;
	std::cout << "Error message: " << message << std::endl;
	std::cout << "Press enter to exit..." << std::endl;
	getchar();
	exit(-1);
}
