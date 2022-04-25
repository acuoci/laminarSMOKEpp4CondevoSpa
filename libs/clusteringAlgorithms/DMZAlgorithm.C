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

#include "DMZAlgorithm.H"

#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <numeric>
#include <array>
#include <iomanip>

#ifdef __clang__
#pragma clang diagnostic ignored "-Wmissing-braces"
#endif

#include "dkm.hpp"

void BisectionAlgorithm(	double& C, double& sigma, std::vector<double>& y, std::vector<int>& lj,
							double& C_new, double& sigma_new, std::vector<double>& y_new, std::vector<int>& lj_new);

unsigned int PartitionSingleFeature(std::vector<std::vector<double>>& y, std::vector<std::vector<int>>& lj, const double epsilon);


DMZAlgorithm::DMZAlgorithm() : 
	ClusteringAlgorithms()
{
	nc_init_ = 5;
}

void DMZAlgorithm::SetNumberOfFeatures(const unsigned int D)
{
	D_ = D;

	for (unsigned int i = 0; i < D_; i++)
	{
		epsilon_.resize(D_);
		psi_.resize(D_);
		C_.resize(D_);
		sigma_.resize(D_);
		C_avg_.resize(D_);
		sigma_avg_.resize(D_);
	}
}

void DMZAlgorithm::SetData(const unsigned int index, const std::vector<double>& psi, const double epsilon, const bool relative)
{
	if (verbosity_ > 0)
		std::cout << " * Setting variable #" << index << std::endl;

	if (index >= D_)
		FatalErrorMessage("The provided number of features is not sufficient to accomodate the input vector of data");

	n_ = psi.size();
	psi_[index] = psi;
	epsilon_[index] = static_cast<double>(epsilon);
	
	// Check if the accuracy is provided in relative or absolute units
	// The binning algorithm works with relative accuracy
	if (relative == true)
	{
		auto result = std::minmax_element(std::begin(psi), std::end(psi));
		const double psi_min = *result.first;
		const double psi_max = *result.second;
		epsilon_[index] *= (psi_max - psi_min);
	}

	if (verbosity_ > 0)
	{
		auto result = std::minmax_element(std::begin(psi), std::end(psi));
		const double psi_min = *result.first;
		const double psi_max = *result.second;
		const int ne = std::ceil((psi_max - psi_min) / epsilon_[index]);

		std::cout << "   - Min: " << psi_min << " Max: " << psi_max << " delta: " << psi_max - psi_min << std::endl;
		std::cout << "   - Epsilon: " << epsilon_[index] << " Equivalent intervals: " << ne << std::endl;
		if (ne < 10)	std::cout << "   - WARNING: the number of equivalent intervals is small. Better to consider a smaller epsilon coefficient" << std::endl;
	}

	if (index == 0)
	{

		const double threshold = 1.e-6;
		auto result = std::minmax_element(std::begin(psi), std::end(psi));
		const double psi_min = *result.first - threshold;
		const double psi_max = *result.second + threshold;

		std::vector<double> psi_tilde(n_);
		for (unsigned int i = 0; i < n_; i++)
			psi_tilde[i] = (psi[i] - psi_min) / (psi_max - psi_min);

		std::vector<unsigned int> idx(n_);
		for (unsigned int i = 0; i < n_; i++)
			idx[i] = static_cast<unsigned int>(std::ceil(psi_tilde[i] * nc_init_)) - 1;

		// Calculate the number of elements in each bin and allocate memory
		{
			std::vector<unsigned int> ne(nc_init_, 0);
			for (unsigned int i = 0; i < n_; i++)
				ne[idx[i]]++;
			const unsigned int ne_tot = std::accumulate(ne.begin(), ne.end(), 0);
			if (ne_tot != n_)
				FatalErrorMessage("Something wrong in the preliminary binning operation along the first feature");

			// Memory allocation
			y0_.resize(nc_init_);
			lj_.resize(nc_init_);
			for (unsigned int i = 0; i < nc_init_; i++)
			{
				y0_[i].resize(ne[i]);
				lj_[i].resize(ne[i]);
			}
		}

		// Initialize the bin containers
		{
			std::vector<unsigned int> ne(nc_init_, 0);
			for (unsigned int i = 0; i < n_; i++)
			{
				const unsigned int j = idx[i];
				y0_[j][ne[j]] = psi_[index][i];
				lj_[j][ne[j]] = i;
				ne[j]++;
			}
		}
	}
}

void DMZAlgorithm::ClusterData()
{
	std::cout << " * Clustering data DMZ algorithm..." << std::endl;

	CheckInputData();

	auto cpu_start = std::chrono::steady_clock::now();

	// Index 0: main variable
	nc_ = PartitionSingleFeature(y0_, lj_, epsilon_[0]);

	if (verbosity_ > 0)
		std::cout << " * Total number of clusters (stage 1): " << nc_ << std::endl;

	// Index 1
	for (unsigned int j=1;j<D_;j++)
	{
		std::vector< std::vector<int> >		ljtmp_;
		ljtmp_.reserve(nc_);

		for (unsigned int i = 0; i < nc_; i++)
		{
			// Number of elements in the sub-cluster
			const unsigned int n = lj_[i].size();
			
			// Reusing always the same objects for multiple features
			std::vector< std::vector<double> >  yloc_(1);
			ljloc_.resize(1);

			// Fill the local objects 
			ljloc_[0] = lj_[i];
			yloc_[0].resize(n);
			for (unsigned k = 0; k < n; k++)
				yloc_[0][k] = psi_[j][ljloc_[0][k]];

			// Partition on sub-cluster
			const unsigned int nc_split = PartitionSingleFeature(yloc_, ljloc_, epsilon_[j]);

			// Move the local clustering results to a temporary object
			for (unsigned int k = 0; k < nc_split; k++)
				ljtmp_.push_back(ljloc_[k]);
		}

		// Update the total number of clusters
		nc_ = ljtmp_.size();
		
		// Copy the temporary object to the global clustering object
		lj_.resize(0);
		lj_.resize(nc_);
		for (unsigned int k = 0; k < nc_; k++)
			lj_[k] = ljtmp_[k];

		if (verbosity_ > 0)
			std::cout << " * Total number of clusters (stage " << j+1 << "): " << nc_ << std::endl;
	}

	// Global indices
	gj_.resize(n_);
	for (unsigned int k = 0; k < nc_; k++)
	{
		unsigned int n = lj_[k].size();
		for (unsigned int i = 0; i < n; i++)
			gj_[lj_[k][i]] = k;
	}

	// Number of elements in each cluster
	ne_.resize(nc_);
	for (unsigned int i = 0; i < nc_; i++)
		ne_[i] = lj_[i].size();

	auto cpu_end = std::chrono::steady_clock::now();

	std::chrono::duration<double> cpu_clustering = cpu_end - cpu_start;
	cpu_clustering_ = cpu_clustering.count();
}

void DMZAlgorithm::CheckInputData()
{
	for (unsigned int i = 0; i < D_; i++)
		if (psi_[i].size() != n_)
			FatalErrorMessage("Wrong input data: the sizes of input vectors are not consistent");
}

void DMZAlgorithm::FatalErrorMessage(const std::string message)
{
	std::cout << "BinningAlgorithm: Fatal error" << std::endl;
	std::cout << "Error message: " << message << std::endl;
	std::cout << "Press enter to exit..." << std::endl;
	getchar();
	exit(-1);
}

unsigned int PartitionSingleFeature(std::vector<std::vector<double>>& y, std::vector<std::vector<int>>& lj, const double epsilon)
{
	unsigned int nc = y.size();

	std::vector<double> C(nc);
	std::vector<double> sigma(nc);
	for (unsigned int i = 0; i < nc; i++)
	{
		C[i] = std::accumulate(y[i].begin(), y[i].end(), 0.) / static_cast<double>(y[i].size());
		sigma[i] = dispersion_coefficient(y[i]);
	}

	bool convergence = false;
	unsigned int loop_number = 0;

	while (convergence == false)
	{
		loop_number++;
	
		// Count the number of clusters to be split
		unsigned int added_cluster = 0;
		for (unsigned int k = 0; k < nc; k++)
			if (sigma[k] > epsilon)
				added_cluster++;

		// Split the clusters according to the bisection algorithm
		if (added_cluster != 0)
		{
			// Reserve additional memory
			C.resize(nc + added_cluster);
			sigma.resize(nc + added_cluster);
			y.resize(nc + added_cluster);
			lj.resize(nc + added_cluster);
			
			std::vector< std::vector<unsigned int> > couples;

			added_cluster = 0;
			for (unsigned int k = 0; k < nc; k++)
			{
				if (sigma[k] > epsilon)
				{
					added_cluster++;

					const unsigned int index = nc + added_cluster - 1;

					BisectionAlgorithm(C[k], sigma[k], y[k], lj[k], C[index], sigma[index], y[index], lj[index]);
					
					couples.push_back(std::vector<unsigned int>({ k, index }));
				}
			}

			nc = C.size();
			
			// Apply the k-means algorithm on the pairs of bins coming from the splitting operation
			{
				for (unsigned int i = 0; i < added_cluster; i++)
				{
					const unsigned int idx0 = couples[i][0];
					const unsigned int idx1 = couples[i][1];

					const unsigned int n = y[idx0].size() + y[idx1].size();
					std::vector< std::array<double, 1> > data(n);
					for (unsigned int j = 0; j < y[idx0].size(); j++)
						data[j] = { y[idx0][j] };
					
					for (unsigned int j = 0; j < y[idx1].size(); j++)
						data[y[idx0].size()+j] = { y[idx1][j] };
					
					std::vector<int> indices = lj[idx0];
					indices.insert(indices.end(), lj[idx1].begin(), lj[idx1].end());

					auto cluster_data = dkm::kmeans_lloyd(data, 2);
					
					// Recognize the number of elements in each cluster and reallocate memory
					{
						unsigned int label1 = 0;
						for (const auto& label : std::get<1>(cluster_data))
							label1 += label;
						unsigned int label0 = n - label1;

						y[idx0].resize(label0);
						y[idx1].resize(label1);
						lj[idx0].resize(label0);
						lj[idx1].resize(label1);
					}

					// Reconstruct cluster labels
					{
						

						unsigned int label0 = 0;
						unsigned int label1 = 0;
						unsigned int j = 0;
						for (const auto& label : std::get<1>(cluster_data))
						{
							if (label == 0)
							{
								lj[idx0][label0] = indices[j];
								y[idx0][label0] = data[j][0];
								label0++;
								j++;
							}
							else
							{
								lj[idx1][label1] = indices[j];
								y[idx1][label1] = data[j][0];
								label1++;
								j++;
							}
						}
					}

					// Update centroids and dispersion
					{
						C[idx0] = std::accumulate(y[idx0].begin(), y[idx0].end(), 0.) / static_cast<double>(y[idx0].size());
						C[idx1] = std::accumulate(y[idx1].begin(), y[idx1].end(), 0.) / static_cast<double>(y[idx1].size());
						sigma[idx0] = dispersion_coefficient(y[idx0]);
						sigma[idx1] = dispersion_coefficient(y[idx1]);
					}
				}
			}
		}
		else
		{
			convergence = true;
		}
	}

	return nc;
}

void BisectionAlgorithm(	double& C, double& sigma, std::vector<double>& y, std::vector<int>& lj,
							double& C_new, double& sigma_new, std::vector<double>& y_new, std::vector<int>& lj_new )
{
	const unsigned int n = y.size();
	std::vector<double> d(n);
	for (unsigned int i = 0; i < n; i++)
		d[i] = std::fabs(C - y[i]);
	auto result = std::max_element(std::begin(d), std::end(d));
	auto imax = std::distance(d.begin(), result);


	// New cluster
	C_new = y[imax];
	sigma_new = 0.;
	y_new.resize(1, y[imax]);
	lj_new.resize(1, lj[imax]);

	// Update original cluster
	C = (C * n - C_new) / static_cast<double>(n - 1);
	y.erase(y.begin() + imax);
	lj.erase(lj.begin() + imax);
	sigma = dispersion_coefficient(y);
}

