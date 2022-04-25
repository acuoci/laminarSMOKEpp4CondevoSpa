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

#include "ClusteringAlgorithms.H"

#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <numeric>
#include <array>
#include <iomanip>

double centroid(const std::vector<double>& y)
{
	return accumulate(y.begin(), y.end(), 0.) / static_cast<double>(y.size());
}

double dispersion_coefficient(const std::vector<double>& y)
{
	auto result = std::minmax_element(std::begin(y), std::end(y));
	return (*result.second - *result.first);
}

ClusteringAlgorithms::ClusteringAlgorithms()
{
	D_ = 0;
	n_ = 0;
	nc_ = 0;
	cpu_clustering_ = 0.;
	cpu_centroids_ = 0.;
	cpu_dispersion_ = 0.;
	verbosity_ = 1;
}

void ClusteringAlgorithms::SetVerbosityLevel(const unsigned int verbosity)
{
	verbosity_ = verbosity;
}

void ClusteringAlgorithms::SetNumberOfFeatures(const unsigned int D)
{
	D_ = D;

	for (unsigned int i = 0; i < D_; i++)
	{
		epsilon_.resize(D_);
		C_.resize(D_);
		sigma_.resize(D_);
		C_avg_.resize(D_);
		sigma_avg_.resize(D_);
	}
}

void  ClusteringAlgorithms::GetData(unsigned int& nc, std::vector<int>& gj, std::vector< std::vector<int> >& lj)
{
	nc = nc_;
	gj = gj_;
	lj = lj_;
}

void ClusteringAlgorithms::CalculateCentroids()
{
	auto cpu_start = std::chrono::steady_clock::now();

	for (unsigned int k = 0; k < D_; k++)
	{
		C_[k].resize(nc_);
		for (unsigned int i = 0; i < n_; i++)
			C_[k][gj_[i]] += psi_[k][i];
		for (unsigned int i = 0; i < nc_; i++)
			C_[k][i] /= static_cast<double>(lj_[i].size());

		C_avg_[k] = std::accumulate(C_[k].begin(), C_[k].end(), 0.0) / static_cast<double>(C_[k].size());
	}

	auto cpu_end = std::chrono::steady_clock::now();

	std::chrono::duration<double> cpu_centroids = cpu_end - cpu_start;
	cpu_centroids_ = cpu_centroids.count();
}

void ClusteringAlgorithms::CalculateDispersion()
{
	auto cpu_start = std::chrono::steady_clock::now();

	// Resize
	for (unsigned int j = 0; j < D_; j++)
		for (unsigned int k = 0; k < nc_; k++)
			sigma_[j].resize(nc_);

	// Fill
	for (unsigned int k = 0; k < nc_; k++)
	{
		for (unsigned int j = 0; j < D_; j++)
		{
			unsigned int n = lj_[k].size();
			std::vector<double> y(n);
			for (unsigned int i = 0; i < n; i++)
				y[i] = psi_[j][lj_[k][i]];

			sigma_[j][k] = dispersion_coefficient(y);
		}
	}

	// Average values
	for (unsigned int k = 0; k < D_; k++)
		sigma_avg_[k] = std::accumulate(sigma_[k].begin(), sigma_[k].end(), 0.0) / static_cast<double>(sigma_[k].size());

	auto cpu_end = std::chrono::steady_clock::now();

	std::chrono::duration<double> cpu_dispersion = cpu_end - cpu_start;
	cpu_dispersion_ = cpu_dispersion.count();
}

void ClusteringAlgorithms::CheckInputData()
{
	for (unsigned int i = 0; i < D_; i++)
		if (psi_[i].size() != n_)
			FatalErrorMessage("Wrong input data: the sizes of input vectors are not consistent");
}

void ClusteringAlgorithms::Print()
{
	auto ne_minmax = std::minmax_element(std::begin(ne_), std::end(ne_));
	const double ne_avg = std::accumulate(ne_.begin(), ne_.end(), 0) / static_cast<double>(nc_);

	std::cout << std::endl;
	std::cout << "------------------------------------------------------------------------------------" << std::endl;
	std::cout << " Summary of clustering algorithm                                                    " << std::endl;
	std::cout << "------------------------------------------------------------------------------------" << std::endl;
	std::cout << " - total number of points:   " << n_ << std::endl;
	std::cout << " - total number of clusters: " << nc_ << " (" << (100. * nc_) / static_cast<double>(n_) << "%)" << std::endl;
	std::cout << " - cluster size min/max/avg: " << *ne_minmax.first << "/" << *ne_minmax.second << "/" << ne_avg << std::endl;
	std::cout << " - clustering cpu time (ms): " << cpu_clustering_ * 1000. << std::endl;

	if (verbosity_ > 0)
	{
		std::cout << std::endl;
		std::cout << "------------------------------------------------------------------------------------" << std::endl;
		std::cout << " Centroids and dispersion coefficients                                              " << std::endl;
		std::cout << "------------------------------------------------------------------------------------" << std::endl;
		std::cout << std::setw(6) << std::left << "#";
		std::cout << std::setw(8) << std::left << "ncells";
		for (unsigned int k = 0; k < D_; k++)
		{
			std::string centroid = "C(" + std::to_string(k) + ")";
			std::cout << std::setw(11) << std::setprecision(4) << std::left << centroid;
			std::string dispersion = "sigma(" + std::to_string(k) + ")";
			std::cout << std::setw(11) << std::setprecision(4) << std::left << dispersion;
		}
		std::cout << std::endl;

		if (verbosity_ > 2)
		{
			for (unsigned int i = 0; i < nc_; i++)
			{
				std::cout << std::setw(6) << std::left << i;
				std::cout << std::setw(8) << std::left << lj_[i].size();
				for (unsigned int k = 0; k < D_; k++)
				{
					std::cout << std::setw(11) << std::setprecision(4) << std::left << C_[k][i];
					std::cout << std::setw(11) << std::setprecision(4) << std::left << sigma_[k][i];
				}
				std::cout << std::endl;
			}
		}

		std::cout << std::setw(6) << std::left << "Max";
		std::cout << std::setw(8) << std::left << n_;
		for (unsigned int k = 0; k < D_; k++)
		{
			auto max_C = std::max_element(std::begin(C_[k]), std::end(C_[k]));
			std::cout << std::setw(11) << std::setprecision(4) << std::left << *max_C;
			auto max_sigma = std::max_element(std::begin(sigma_[k]), std::end(sigma_[k]));
			std::cout << std::setw(11) << std::setprecision(4) << std::left << *max_sigma;
		}
		std::cout << std::endl;

		std::cout << std::setw(6) << std::left << "Avg";
		std::cout << std::setw(8) << std::left << n_;
		for (unsigned int k = 0; k < D_; k++)
		{
			std::cout << std::setw(11) << std::setprecision(4) << std::left << C_avg_[k];
			std::cout << std::setw(11) << std::setprecision(4) << std::left << sigma_avg_[k];
		}
		std::cout << std::endl;
	}

	// Statistical distribution of cluster size
	if (verbosity_ > 1)
	{
		std::vector<int> v = ne_;
		std::sort(v.begin(), v.end());

		std::vector<int> threshold{ 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 100000 };
		std::vector<int> cdf(threshold.size());
		for (unsigned int i = 0; i < threshold.size(); i++)
			cdf[i] = std::upper_bound(v.begin(), v.end(), threshold[i]) - v.begin();

		std::cout << std::endl;
		std::cout << "------------------------------------------------------------------------------------" << std::endl;
		std::cout << " Cumulative Distribution Function                                                   " << std::endl;
		std::cout << "------------------------------------------------------------------------------------" << std::endl;
		std::cout << std::setw(12) << std::left << "ncells";
		std::cout << std::setw(12) << std::left << "occur.";
		std::cout << std::setw(12) << std::left << "cdf(%)";
		std::cout << std::endl;
		for (unsigned int i = 0; i < threshold.size(); i++)
		{
			std::cout << std::setw(12) << std::left << threshold[i];
			std::cout << std::setw(12) << std::left << cdf[i];
			std::cout << std::setw(12) << std::left << cdf[i] / static_cast<double>(nc_) * 100.;
			std::cout << std::endl;
		}

		std::cout << std::endl;
		std::cout << "------------------------------------------------------------------------------------" << std::endl;
		std::cout << " Probability Distribution Function                                                  " << std::endl;
		std::cout << "------------------------------------------------------------------------------------" << std::endl;
		std::cout << std::setw(16) << std::left << "ncells";
		std::cout << std::setw(12) << std::left << "center";
		std::cout << std::setw(12) << std::left << "occur.";
		std::cout << std::setw(12) << std::left << "pdf(%)";
		std::cout << std::endl;
		std::cout << std::setw(16) << std::left << threshold[0];
		std::cout << std::setw(12) << std::left << threshold[0];
		std::cout << std::setw(12) << std::left << cdf[0];
		std::cout << std::setw(12) << std::left << cdf[0] / static_cast<double>(nc_) * 100.;
		std::cout << std::endl;
		std::cout << std::setw(16) << std::left << threshold[1];
		std::cout << std::setw(12) << std::left << threshold[1];
		std::cout << std::setw(12) << std::left << cdf[1] - cdf[0];
		std::cout << std::setw(12) << std::left << (cdf[1] - cdf[0]) / static_cast<double>(nc_) * 100.;
		std::cout << std::endl;
		for (unsigned int i = 2; i < threshold.size(); i++)
		{
			std::string interval = std::to_string(threshold[i - 1] + 1) + "-" + std::to_string(threshold[i]);
			std::cout << std::setw(16) << std::left << interval;
			std::cout << std::setw(12) << std::left << (threshold[i - 1] + 1 + threshold[i]) / 2.;
			std::cout << std::setw(12) << std::left << (cdf[i] - cdf[i - 1]);
			std::cout << std::setw(12) << std::left << (cdf[i] - cdf[i - 1]) / static_cast<double>(nc_) * 100.;
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
}
