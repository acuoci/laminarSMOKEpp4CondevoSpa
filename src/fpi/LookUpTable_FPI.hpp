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

#include "Grammar_LookUpTable_FPI.h"

#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <string_view>

class CSVRow
{
    public:

        std::string operator[](std::size_t index) const
        {
            return std::string(&m_line[m_data[index] + 1], m_data[index + 1] -  (m_data[index] + 1));
        }
        std::size_t size() const
        {
            return m_data.size() - 1;
        }
        void readNextRow(std::istream& str)
        {
            std::getline(str, m_line);

            m_data.clear();
            m_data.emplace_back(-1);
            std::string::size_type pos = 0;
            while((pos = m_line.find(',', pos)) != std::string::npos)
            {
                m_data.emplace_back(pos);
                ++pos;
            }
            // This checks for a trailing comma with no data after it.
            pos   = m_line.size();
            m_data.emplace_back(pos);
        }
    private:
        std::string         m_line;
        std::vector<int>    m_data;
};

std::istream& operator>>(std::istream& str, CSVRow& data)
{
    data.readNextRow(str);
    return str;
}   


namespace OpenSMOKE
{
	LookUpTable_FPI::LookUpTable_FPI(OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap) :
		thermodynamicsMap_(thermodynamicsMap)
	{
		is_active_ = false;

		boost::filesystem::path path_output_;

		type_ = LookUpTable_FPI::UNIFORM_Y;
		is_apriori_ = true;

		nc_ = 0;
		dc_ = 0.;
		dctilde_ = 0.;
		c_min_ = 1.e16;
		c_max_ = -1.e16;
		query_c_ = -1.e16;

		ns_ = thermodynamicsMap_.NumberOfSpecies();

		jC_ = thermodynamicsMap_.IndexOfElementWithoutError("C") - 1;
		jO_ = thermodynamicsMap_.IndexOfElementWithoutError("O") - 1;
		jH_ = thermodynamicsMap_.IndexOfElementWithoutError("H") - 1;

		WC_ = OpenSMOKE::AtomicWeights["C"];
		WO_ = OpenSMOKE::AtomicWeights["O"];
		WH_ = OpenSMOKE::AtomicWeights["H"];
	}

	void LookUpTable_FPI::SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary)
	{
		Grammar_LookUpTable_FPI grammar;
		dictionary.SetGrammar(grammar);

		if (dictionary.CheckOption("@Path") == true)
			dictionary.ReadPath("@Path", path_output_);

		if (dictionary.CheckOption("@Apriori") == true)
			dictionary.ReadBool("@Apriori", is_apriori_);

	//	if (dictionary.CheckOption("@Fields") == true)
	//		dictionary.ReadOption("@Fields", list_fields_);

		ReadTableFromCSVFile();

		Summary();

		CheckInput();
	}

	void LookUpTable_FPI::Setup(const boost::filesystem::path& path_to_csv, const bool is_apriori)
	{
		path_output_ = path_to_csv;
		//list_fields_ = list_fields;
		is_apriori_ = is_apriori;

		ReadTableFromCSVFile();

		Summary();

		CheckInput();
	}

	void LookUpTable_FPI::ReadTableFromCSVFile()
	{
		std::vector<std::string> header_tags;
		std::vector<std::vector<double>> data_matrix;

		{
			std::cout << "Open XML file: " << path_output_.string() << std::endl;
    			std::ifstream fCSV(path_output_.string());

   			CSVRow	row;
			int line = 0;
    			while(fCSV >> row)
    			{
				std::cout << "Reading line " << line << std::endl;
				if (line == 0)
				{
					for (unsigned int i=0;i<row.size();i++)
						header_tags.push_back(row[i]);
				}
				else
				{
					std::vector<double> line_data(row.size());
					for (unsigned int i=0;i<row.size();i++)
						line_data[i] = std::stod(row[i]);
					data_matrix.push_back(line_data);
				}
				
				line++;
			}
		}

		// Check consistency with thermodynamics
		for (unsigned int j=0;j<ns_;j++)
		{
			Info << "Table: " << header_tags[12+j] << "  Thermodynamics: " << thermodynamicsMap_.NamesOfSpecies()[j] << endl;
			if (header_tags[12+j] != thermodynamicsMap_.NamesOfSpecies()[j])
				FatalError << "The thermodynamic map is not consistent with the FPI table" << ::Foam::exit(FatalError);
		}

		nc_ = data_matrix.size();

		c_.resize(nc_);
		ctilde_.resize(nc_);
		omega_c_.resize(nc_);
		omega_ctilde_.resize(nc_);

		rho_.resize(nc_);
		Cp_.resize(nc_);
		lambda_.resize(nc_);
		mu_.resize(nc_);
		MW_.resize(nc_);
		Q_.resize(nc_);
		a_.resize(nc_);
		T_.resize(nc_);
		Y_.resize(nc_, ns_);
		interpolated_Y_.resize(ns_);

		for (unsigned int i=0;i<nc_;i++)
		{
			c_(i) = data_matrix[i][0];
			ctilde_(i) = data_matrix[i][1];
			omega_c_(i) = data_matrix[i][2];
			omega_ctilde_(i) = data_matrix[i][3];

			rho_(i) = data_matrix[i][4];
			Cp_(i) = data_matrix[i][5];
			lambda_(i) = data_matrix[i][6];
			mu_(i) = data_matrix[i][7];
			MW_(i) = data_matrix[i][8];
			Q_(i) = data_matrix[i][9];
			a_(i) = data_matrix[i][10];
			T_(i) = data_matrix[i][11];

			for (unsigned int j=0;j<ns_;j++)
				Y_(i,j) = data_matrix[i][12+j];
		}

		dc_ = c_(1) - c_(0);
		dctilde_ = ctilde_(1) - ctilde_(0);
		c_min_ = c_(0);
		c_max_ = c_(nc_-1);

		is_active_ = true;
	}

	void LookUpTable_FPI::CheckInput()
	{
		std::cout << "Testing lookup table..." << std::endl;

		// Internal points
		{
			const unsigned int np = 5;
			for (unsigned int i = 0; i < np; i++)
			{
				
				const double c = c_min_ + (c_max_ - c_min_) / static_cast<double>(np - 1) * static_cast<double>(i);

				Interpolate(c);

				std::cout << std::left << std::setw(16) << c << std::endl;
				std::cout << std::left << std::setw(16) << interpolated_rho_ << std::endl;
				std::cout << std::left << std::setw(16) << interpolated_Cp_ << std::endl;
				std::cout << std::left << std::setw(16) << interpolated_lambda_ << std::endl;
				std::cout << std::left << std::setw(16) << interpolated_mu_ << std::endl;
				std::cout << std::left << std::setw(16) << interpolated_MW_ << std::endl;
				std::cout << std::left << std::setw(16) << interpolated_Q_ << std::endl;			
				std::cout << std::left << std::setw(16) << interpolated_a_ << std::endl;
				std::cout << std::left << std::setw(16) << interpolated_T_ << std::endl;
				std::cout << std::left << std::setw(16) << interpolated_OmegaC_ << std::endl;
				std::cout << std::left << std::setw(16) << interpolated_OmegaCtilde_ << std::endl;
				for (unsigned int k = 0; k < 8; k++)
					std::cout << std::left << std::setw(16) << interpolated_Y_(k);
				std::cout << std::endl;
			}
		}
	}

	void LookUpTable_FPI::Summary()
	{
		std::cout << "Progress variable definition..." << std::endl;

		std::cout << std::endl;
		std::cout << "Progress variable space: " << c_min_ << " " << c_max_ << " (" << nc_ << ")" << std::endl;

		std::cout << "Min/Max/Mean Values (additional fields)" << std::endl;
		std::cout << " * Density (kg/m3):                " << rho_.minCoeff() << " " << rho_.maxCoeff() << " " << rho_.mean() << std::endl;
		std::cout << " * Specific heat (J/kg/K):         " << Cp_.minCoeff() << " " << Cp_.maxCoeff() << " " << Cp_.mean() << std::endl;
		std::cout << " * Thermal conductivity (W/m/K):   " << lambda_.minCoeff() << " " << lambda_.maxCoeff() << " " << lambda_.mean() << std::endl;
		std::cout << " * Viscosity (kg/m/s):             " << mu_.minCoeff() << " " << mu_.maxCoeff() << " " << mu_.mean() << std::endl;
		std::cout << " * Molecular weight (kg/kmol):     " << MW_.minCoeff() << " " << MW_.maxCoeff() << " " << MW_.mean() << std::endl;
		std::cout << " * Heat release rate (W/m3):       " << Q_.minCoeff() << " " << Q_.maxCoeff() << " " << Q_.mean() << std::endl;
		std::cout << " * Planck absorption (1/m):        " << a_.minCoeff() << " " << a_.maxCoeff() << " " << a_.mean() << std::endl;
		std::cout << " * Temperature (kg/m3):            " << T_.minCoeff() << " " << T_.maxCoeff() << " " << T_.mean() << std::endl;
		std::cout << " * Source term (kg/m3/s):          " << omega_c_.minCoeff() << " " << omega_c_.maxCoeff() << " " << omega_c_.mean() << std::endl;
		std::cout << " * Source term (norm.) (kg/m3/s):  " << omega_ctilde_.minCoeff() << " " << omega_ctilde_.maxCoeff() << " " << omega_ctilde_.mean() << std::endl;

		std::cout << "Min/Max/Mean Values (mass fractions)" << std::endl;
		for (unsigned int k = 0; k < ns_; k++)
			std::cout << " * " << thermodynamicsMap_.NamesOfSpecies()[k] << ": " << Y_.col(k).minCoeff() << " " << Y_.col(k).maxCoeff() << " " << Y_.col(k).mean() << std::endl;
	}

	void LookUpTable_FPI::InterpolateFromNormalizedProgressVariable(const double ctilde)
	{
		// Minimum/Maximum local progress variable
		const double Cmin = c_min_;
		const double Cmax = c_max_;
		const double c = ctilde * (Cmax - Cmin) + Cmin;

		Interpolate(c);
	}

	void LookUpTable_FPI::Interpolate(const double c)
	{
		//std::cout << "Query " << c << " " << query_c_ << std::endl;

		if (c == query_c_)
			return;
		query_c_ = c;

		{
			// Normlized progress variable
			const double Cmin = c_min_;
			const double Cmax = c_max_;
			const double query_ctilde = std::max(std::min((query_c_ - Cmin) / (Cmax - Cmin), 1.), 0.);

			//std::cout << "Interpolate " << c << " " << query_c_ << " " << Cmin << " " << Cmax << " " << query_ctilde << std::endl;

			// Regular points
			if (query_ctilde <= 0.)
			{
				interpolated_rho_ = rho_(0) ;
				interpolated_Cp_ = Cp_(0) ;
				interpolated_lambda_ = lambda_(0);
				interpolated_mu_ = mu_(0);
				interpolated_MW_ = MW_(0);
				interpolated_Q_ = Q_(0);
				interpolated_a_ = a_(0);
				interpolated_T_ = T_(0);

				interpolated_OmegaC_ = omega_c_(0);
				interpolated_OmegaCtilde_ = omega_ctilde_(0);

				for (unsigned int k=0; k<ns_; k++)
					interpolated_Y_(k) = Y_(0,k);
			}
			else if (query_ctilde >= 1.)
			{
				const unsigned int j = nc_-1;

				interpolated_rho_ = rho_(j) ;
				interpolated_Cp_ = Cp_(j) ;
				interpolated_lambda_ = lambda_(j);
				interpolated_mu_ = mu_(j);
				interpolated_MW_ = MW_(j);
				interpolated_Q_ = Q_(j);
				interpolated_a_ = a_(j);
				interpolated_T_ = T_(j);

				interpolated_OmegaC_ = omega_c_(j);
				interpolated_OmegaCtilde_ = omega_ctilde_(j);

				for (unsigned int k=0; k<ns_; k++)
					interpolated_Y_(k) = Y_(j,k);
			}
			else
			{
				const unsigned int ic = std::floor((query_ctilde - 0.) / dctilde_);

				const double coeff = (query_ctilde - ctilde_(ic+1))/(ctilde_(ic+1)-ctilde_(ic));
				
				interpolated_rho_ = rho_(ic) + (rho_(ic+1)-rho_(ic)) * coeff;
				interpolated_Cp_ = Cp_(ic) + (Cp_(ic+1)-Cp_(ic)) * coeff;
				interpolated_lambda_ = lambda_(ic) + (lambda_(ic+1)-lambda_(ic)) * coeff;
				interpolated_mu_ = mu_(ic) + (mu_(ic+1)-mu_(ic)) * coeff;
				interpolated_MW_ = MW_(ic) + (MW_(ic+1)-MW_(ic)) * coeff;
				interpolated_Q_ = Q_(ic) + (Q_(ic+1)-Q_(ic)) * coeff;
				interpolated_a_ = a_(ic) + (a_(ic+1)-a_(ic)) * coeff;
				interpolated_T_ = T_(ic) + (T_(ic+1)-T_(ic)) * coeff;

				interpolated_OmegaC_ = omega_c_(ic) + (omega_c_(ic+1)-omega_c_(ic)) * coeff;
				interpolated_OmegaCtilde_ = omega_ctilde_(ic) + (omega_ctilde_(ic+1)-omega_ctilde_(ic)) * coeff;

				for (unsigned int k=0; k<ns_; k++)
					interpolated_Y_(k) = Y_(ic,k) + (Y_(ic + 1,k)-Y_(ic,k)) * coeff;
			}
		}
	}
}
