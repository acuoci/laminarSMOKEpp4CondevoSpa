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
|   Copyright(C) 2014, 2013, 2012  Alberto Cuoci                          |
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

namespace OpenSMOKE
{
	#include <string>
	#include <sstream>
	#include <algorithm>

	bool BIN_is_radical(const std::string name, const char label_radical)
	{
		if (name.back() == label_radical) return true;
		else return false;
	}

	std::string BIN_hydrogenation_level(const std::string name, const std::string label, const int section)
	{
		std::string provisional = name;

		// Check and modify if radical species
		if ( BIN_is_radical(provisional, 'J') == true)
			provisional.pop_back();

		// Remove the main label
		std::stringstream section_label;
		section_label << section;
		const std::string toErase = label + section_label.str();
		size_t pos = provisional.find(toErase);
		if (pos != std::string::npos)
			provisional.erase(pos, toErase.length());
		
		// Return the hydrogenation level

		return provisional;
	}

	unsigned int BIN_extract_section(const std::string name, const std::string label)
	{
		std::string provisional = name;
		size_t pos = provisional.find(label);
		if (pos != std::string::npos)
			provisional.erase(pos, label.length());

		unsigned int count = 0;
		for (unsigned int i = 0; i < provisional.size(); i++)
		{
			if (isdigit(provisional[i]) == true)
				count++;
			else
				break;
		}

		const std::string section_string = provisional.substr(0, count);
		return std::stoi(section_string);
	}

	template<class T>
	std::vector<T> extract_unique_values(const std::vector<T>& v)
	{
		std::vector<T> v_ = v;
		std::sort(v_.begin(), v_.end());
		typename std::vector<T>::iterator it = std::unique(v_.begin(), v_.end());
		v_.resize(std::distance(v_.begin(), it));
		return v_;
	}

	PolimiSoot_Analyzer::PolimiSoot_Analyzer(OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML) :
		thermo_(*thermodynamicsMapXML)
	{
		// Initialize default values
		Initialize();
	}

	PolimiSoot_Analyzer::PolimiSoot_Analyzer(OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML, OpenSMOKE::OpenSMOKE_Dictionary& dictionary) :
		thermo_(*thermodynamicsMapXML)
	{
		SetupFromDictionary(dictionary);
	}

	void PolimiSoot_Analyzer::SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary)
	{
		// Initialize default values
		Initialize();

		// Read from dictionary
		{
			Grammar_PolimiSoot_Analyzer grammar;
			dictionary.SetGrammar(grammar);

			if (dictionary.CheckOption("@FractalDimension") == true)
				dictionary.ReadDouble("@FractalDimension", Df_);

			if (dictionary.CheckOption("@SootLabel") == true)
				dictionary.ReadString("@SootLabel", bin_label_);

			if (dictionary.CheckOption("@SootMinimumSectionSphericalParticles") == true)
			{
				dictionary.ReadInt("@SootMinimumSectionSphericalParticles", bin_minimum_spheres_section_);
				
				std::stringstream number;
				number << bin_minimum_spheres_section_;
				bin_minimum_spheres_ = bin_label_ + number.str();
			}

			if (dictionary.CheckOption("@SootMinimumSectionAggregates") == true)
			{
				dictionary.ReadInt("@SootMinimumSectionAggregates", bin_minimum_aggregates_section_);

				std::stringstream number;
				number << bin_minimum_aggregates_section_;
				bin_minimum_aggregates_ = bin_label_ + number.str();
			}

			if (dictionary.CheckOption("@ThermophoreticEffect") == true)
				dictionary.ReadBool("@ThermophoreticEffect", thermophoretic_effect_);

			if (dictionary.CheckOption("@ThermophoreticEffectAmplificationFactor") == true)
				dictionary.ReadDouble("@ThermophoreticEffectAmplificationFactor", thermophoretic_effect_amplification_factor_);			

			if (dictionary.CheckOption("@ThermophoreticEffectEuckenApproximation") == true)
				dictionary.ReadBool("@ThermophoreticEffectEuckenApproximation", thermophoretic_effect_eucken_approximation_);

			if (dictionary.CheckOption("@ThermophoreticEffectInCorrectionVelocity") == true)
				dictionary.ReadBool("@ThermophoreticEffectInCorrectionVelocity", thermophoretic_effect_included_in_correction_);

			if (dictionary.CheckOption("@ThermophoreticEffectSmoothingTime") == true)
			{
				std::string units;
				dictionary.ReadMeasure("@ThermophoreticEffectSmoothingTime", thermophoretic_effect_smoothing_time_, units);

				if (units == "s")			thermophoretic_effect_smoothing_time_ = thermophoretic_effect_smoothing_time_;
				else if (units == "ms")		thermophoretic_effect_smoothing_time_ /= 1000.;
				else if (units == "min")	thermophoretic_effect_smoothing_time_ *= 60.;
				else OpenSMOKE::FatalErrorMessage("@ThermophoreticEffectSmoothingTime: check units of time");
			}

			if (dictionary.CheckOption("@ThermophoreticEffectInEnthalpyFluxes") == true)
			{
				std::string flag;
				dictionary.ReadString("@ThermophoreticEffectInEnthalpyFluxes", flag);

				if (flag == "DoNotExclude")
					thermophoretic_effect_in_enthalpy_fluxes_ = ThermophoreticEffectInEnthalpyFluxes::DO_NOT_EXCLUDE_SOOT_CONTRIBUTION;
				else if (flag == "Exclude")
					thermophoretic_effect_in_enthalpy_fluxes_ = ThermophoreticEffectInEnthalpyFluxes::EXCLUDE_TOTAL_SOOT_CONTRIBUTION;
				else if (flag == "ExcludeOnlyThermophoreticEffect")
					thermophoretic_effect_in_enthalpy_fluxes_ = ThermophoreticEffectInEnthalpyFluxes::EXCLUDE_THERMOPHORETIC_SOOT_CONTRIBUTION;
				else
					OpenSMOKE::FatalErrorMessage("@ThermophoreticEffectInEnthalpyFluxes: available options: DoNotExclude | Exclude | ExcludeOnlyThermophoreticEffect");
			}

			if (dictionary.CheckOption("@RadiativeHeatTransfer") == true)
				dictionary.ReadBool("@RadiativeHeatTransfer", radiative_heat_transfer_);

			if (dictionary.CheckOption("@PhysicalDiffusion") == true)
				dictionary.ReadBool("@PhysicalDiffusion", physical_diffusion_);

			if (dictionary.CheckOption("@PhysicalDiffusionReductionCoefficient") == true)
			{
				dictionary.ReadDouble("@PhysicalDiffusionReductionCoefficient", physical_diffusion_reduction_coefficient_);
			}

			if (dictionary.CheckOption("@PhysicalDiffusionBinToCut") == true)
			{
				dictionary.ReadInt("@PhysicalDiffusionBinToCut", physical_diffusion_bin_to_cut_);
			}

			if (dictionary.CheckOption("@PhysicalDiffusionBinToStart") == true)
			{
				dictionary.ReadInt("@PhysicalDiffusionBinToStart", physical_diffusion_bin_to_start_);
			}

			if (dictionary.CheckOption("@PhysicalDiffusionExp") == true)
			{
				dictionary.ReadDouble("@PhysicalDiffusionExp", physical_diffusion_exp_);
			}

			if (dictionary.CheckOption("@ThermophoreticEffectMinimumBin") == true)
			{
				dictionary.ReadInt("@ThermophoreticEffectMinimumBin", thermophoretic_effect_minimum_bin_);
			}
			
			if (dictionary.CheckOption("@PlanckCoefficient") == true)
			{
				std::string flag;
				dictionary.ReadString("@PlanckCoefficient", flag);
				SetPlanckAbsorptionCoefficient(flag);
			}

			if (dictionary.CheckOption("@WritePSDF") == true)
				dictionary.ReadBool("@WritePSDF", write_psdf_);

			if (dictionary.CheckOption("@ThresholdForPSDF") == true)
				dictionary.ReadDouble("@ThresholdForPSDF", threshold_for_psdf_);

			if (dictionary.CheckOption("@Density") == true)
			{
				std::vector<double> values;
				dictionary.ReadOption("@Density", values);

				if (values.size() != 4)
					OpenSMOKE::FatalErrorMessage("@Density wrong density coefficients. Expected format: @Density BINstart DENSITYstart BINend DENSITYend (density in kg/m3)");

				SetDensity(static_cast<int>(values[0]), static_cast<int>(values[2]), values[1], values[3]);
			}

			if (dictionary.CheckOption("@FormationRatesByClasses") == true)
			{
				dictionary.ReadBool("@FormationRatesByClasses", formation_rates_by_classes_);
			}

			if (dictionary.CheckOption("@ClassCorrectionCoefficients") == true)
			{
				std::vector<std::string> list_values;
				dictionary.ReadOption("@ClassCorrectionCoefficients", list_values);

				if (list_values.size() % 2 != 0)
					OpenSMOKE::FatalErrorMessage("@ClassCorrectionCoefficients option has a wrong number of elements.");

				const int n = static_cast<int>(list_values.size() / 2);
				for (int i = 0; i < n; i++)
				{
					list_class_corrections_names_.push_back(list_values[i * 2]);
					list_class_corrections_coefficients_.push_back(boost::lexical_cast<double>(list_values[i * 2 + 1]));
				}
			}
		}

		// Setup
		Setup();
	}

	void PolimiSoot_Analyzer::Initialize()
	{
		//	Default values
		bin_label_			= "BIN";

		bin_minimum_spheres_ = "BIN5";
		bin_minimum_aggregates_ = "BIN13";
		bin_minimum_spheres_section_ = 5;
		bin_minimum_aggregates_section_ = 13;

		bin_density_index_zero_		= 10;
		bin_density_index_final_	= 20;
		bin_density_value_zero_   = 1500.;
		bin_density_value_final_  = 1700.;
		Df_	= 1.8;

		physical_diffusion_ = true;
		physical_diffusion_exp_ = -0.640;
		physical_diffusion_reduction_coefficient_ = 1.;
		physical_diffusion_bin_to_start_ = 5;
		physical_diffusion_bin_to_cut_ = 10;
		bin_physical_diffusion_reference_species_ = -1;


		thermophoretic_effect_ = true;
		thermophoretic_effect_amplification_factor_ = 1.;
		thermophoretic_effect_eucken_approximation_ = true;
		thermophoretic_effect_included_in_correction_ = false;
		thermophoretic_effect_in_enthalpy_fluxes_ = ThermophoreticEffectInEnthalpyFluxes::DO_NOT_EXCLUDE_SOOT_CONTRIBUTION;
		thermophoretic_effect_smoothing_time_ = 0.;
		thermophoretic_effect_minimum_bin_ = 5;
		radiative_heat_transfer_ = true;
		soot_planck_coefficient_ = SootPlanckCoefficient::SOOT_PLANCK_COEFFICIENT_SMOOKE;

		write_psdf_ = true;
		threshold_for_psdf_ = 1e-11;

		soot_classes_ = new PolimiSootClasses();
		formation_rates_by_classes_ = false;
	}

	void PolimiSoot_Analyzer::ClassesFromXMLFile(const boost::filesystem::path& file_name)
	{
		soot_classes_->ReadXMLFile(file_name);
	}

	void PolimiSoot_Analyzer::SetLabel(const std::string label)
	{
		bin_label_ = label;
	}

	void PolimiSoot_Analyzer::SetFractalDiameter(const double Df)
	{
		Df_ = Df;
	}

	void PolimiSoot_Analyzer::SetPhysicalDiffusionBinToStart(const unsigned int physical_diffusion_bin_to_start)
	{
		physical_diffusion_bin_to_start_ = physical_diffusion_bin_to_start;
	}

	void PolimiSoot_Analyzer::SetPhysicalDiffusionBinToCut(const unsigned int physical_diffusion_bin_to_cut)
	{
		physical_diffusion_bin_to_cut_ = physical_diffusion_bin_to_cut;
	}

	void PolimiSoot_Analyzer::SetThermophoreticEffectMinimumBin(const unsigned int thermophoretic_effect_minimum_bin)
	{
		thermophoretic_effect_minimum_bin_ = thermophoretic_effect_minimum_bin;
	}

	void PolimiSoot_Analyzer::SetMinimumSectionSphericalParticles(const int minimum_section_spheres)
	{
		std::stringstream number;
		number << minimum_section_spheres;
		bin_minimum_spheres_section_ = minimum_section_spheres;
		bin_minimum_spheres_ = bin_label_ + number.str();
	}

	void PolimiSoot_Analyzer::SetMinimumSectionSphericalParticles(const std::string bin_minimum_spheres)
	{
		bin_minimum_spheres_section_ = BIN_extract_section(bin_minimum_spheres, "BIN");
		bin_minimum_spheres_ = bin_minimum_spheres;
	}

	void PolimiSoot_Analyzer::SetMinimumSectionAggregates(const int minimum_section_aggregates)
	{
		std::stringstream number;
		number << minimum_section_aggregates;
		bin_minimum_aggregates_section_ = minimum_section_aggregates;
		bin_minimum_aggregates_ = bin_label_ + number.str();
	}

	void PolimiSoot_Analyzer::SetMinimumSectionAggregates(const std::string bin_minimum_aggregates)
	{
		bin_minimum_aggregates_section_ = BIN_extract_section(bin_minimum_aggregates, "BIN");
		bin_minimum_aggregates_ = bin_minimum_aggregates;
	}

	void PolimiSoot_Analyzer::SetDensity(const int bin_index_zero, const int bin_index_final, const double bin_density_zero, const double bin_density_final)
	{
		if ( bin_index_zero <= 0 ||
			 bin_index_final <= bin_index_zero ||
			 bin_density_zero <= 500. || bin_density_zero >= 5000. ||
			 bin_density_final <= 500. || bin_density_final >= 5000. )
		OpenSMOKE::FatalErrorMessage("PolimiSoot_Analyzer: wrong density coefficients. Please check the values: BINstart DENSITYstart BINend DENSITYend (density in kg/m3)");

		bin_density_index_zero_ = bin_index_zero;
		bin_density_index_final_ = bin_index_final;
		bin_density_value_zero_ = bin_density_zero;
		bin_density_value_final_ = bin_density_final;
	}

	void PolimiSoot_Analyzer::SetPlanckAbsorptionCoefficient(const SootPlanckCoefficient soot_planck_coefficient)
	{
		soot_planck_coefficient_ = soot_planck_coefficient;
	}

	void PolimiSoot_Analyzer::SetPlanckAbsorptionCoefficient(const std::string label)
	{
		if (label == "Smooke")
			soot_planck_coefficient_ = SootPlanckCoefficient::SOOT_PLANCK_COEFFICIENT_SMOOKE;
		else if (label == "Kent")
			soot_planck_coefficient_ = SootPlanckCoefficient::SOOT_PLANCK_COEFFICIENT_KENT;
		else if (label == "Sazhin")
			soot_planck_coefficient_ = SootPlanckCoefficient::SOOT_PLANCK_COEFFICIENT_SAZHIN;
		else if (label == "Hubbard")
			soot_planck_coefficient_ = SootPlanckCoefficient::SOOT_PLANCK_COEFFICIENT_HUBBARD;
		else if (label == "none")
			soot_planck_coefficient_ = SootPlanckCoefficient::SOOT_PLANCK_COEFFICIENT_NONE;
		else
			OpenSMOKE::FatalErrorMessage("@PlanckCoefficient: available options: none | Smooke (default) | Kent | Sazhin | Hubbard");
	}

	void PolimiSoot_Analyzer::Setup()
	{
		iC_ = thermo_.IndexOfElement("C") - 1;
		iH_ = thermo_.IndexOfElement("H") - 1;
		iO_ = thermo_.IndexOfElement("O") - 1;

		nspecies_ = thermo_.NumberOfSpecies();

		// Check is soot sections are really available in the kinetic mechanism
		{
			unsigned int count = 0;
			for (unsigned int i = 0; i < nspecies_; i++)
			if (thermo_.NamesOfSpecies()[i].compare(0, bin_label_.size(), bin_label_) == 0)
				count++;

			if (count == 0)
				OpenSMOKE::FatalErrorMessage("No soot sections (" + bin_label_ + ") were found in the provided kinetic mechanism");
		}

		// BIN main properties
		for (unsigned int i = 0; i < nspecies_; i++)
		{
			if (thermo_.NamesOfSpecies()[i].compare(0, bin_label_.size(), bin_label_) == 0)
			{
				const double nc = thermo_.atomic_composition()(i, iC_);
				const double nh = thermo_.atomic_composition()(i, iH_);
				const double no = thermo_.atomic_composition()(i, iO_);

				const unsigned int section = static_cast<int>(std::round(std::log(nc / 24.) / std::log(2.) + 1.));

				const bool is_radical = BIN_is_radical(thermo_.NamesOfSpecies()[i], 'J');
				const std::string hydrogenation_level = BIN_hydrogenation_level(thermo_.NamesOfSpecies()[i], bin_label_, section);

				// Density
				{
					if (section <= bin_density_index_zero_)
					{
						bin_density_.push_back(bin_density_value_zero_);
					}
					else
					{
						const double m = (bin_density_value_final_ - bin_density_value_zero_) / (bin_density_index_final_ - bin_density_index_zero_);
						const double c = bin_density_value_zero_ - m * bin_density_index_zero_;
						bin_density_.push_back(c + m * section);
					}
				}

				bin_indices_.push_back(i);											// Index of bin in the gas phase kinetic scheme [-]
				bin_names_.push_back(thermo_.NamesOfSpecies()[i]);					// Name of bin in the gas phase kinetic scheme [-]
				bin_mw_.push_back(thermo_.MW(i));									// Molecular weight [kg/kmol]
				bin_m_.push_back(thermo_.MW(i) / PhysicalConstants::Nav_kmol);		// Mass of particle [kg]
				bin_section_.push_back(section);									// BIN section [-]
				bin_is_radical_.push_back(is_radical);								// BIN section [-]
				bin_hydrogenation_level_.push_back(hydrogenation_level);			// BIN hydrogenation level

				bin_ds_.push_back(std::pow(6. / PhysicalConstants::pi * thermo_.MW(i) / (bin_density_[bin_density_.size() - 1] / 1000.)
					/ (PhysicalConstants::Nav_mol), 1. / 3.) * 1.e-2);	// Diameter of particle [m]

				bin_V_.push_back(PhysicalConstants::pi / 6. * std::pow(bin_ds_[bin_ds_.size() - 1], 3.));	// Volume of particle [m3]
				bin_c_.push_back(nc);	// C
				bin_h_.push_back(nh);	// H
				bin_o_.push_back(no);	// O
				bin_h_over_c_.push_back(nh / nc);	// Ratio H/C
				bin_o_over_c_.push_back(no / nc);	// Ratio O/C

				if (nh > 0)	bin_o_over_h_.push_back(no / nh);	// Ratio O/H
				else    	bin_o_over_h_.push_back(0.);		// Ratio O/H

				// Thermophoretic effect
				if (section >= thermophoretic_effect_minimum_bin_)
					bin_indices_thermophoresis_.push_back(i);

				// Collect BINs in classes
				if (section >= bin_minimum_spheres_section_)
				{
					const int index = static_cast<int>(bin_indices_.size()) - 1;

					bin_indices_large_.push_back(index);
					bin_indices_large_global_.push_back(i);
					bin_density_large_.push_back(bin_density_[index]);
					bin_V_large_.push_back(bin_V_[index]);
				}
				else
				{
					const int index = static_cast<int>(bin_indices_.size() - 1);

					bin_indices_small_.push_back(index);
					bin_indices_small_global_.push_back(i);
					bin_density_small_.push_back(bin_density_[index]);
					bin_V_small_.push_back(bin_V_[index]);
				}
			}
		}

		// Identification of primary particle
		unsigned int bin_primary_particle_index = 0;
		for (unsigned int i = 0; i < bin_indices_.size(); i++)
			if (bin_section_[i] == (bin_minimum_aggregates_section_ - 1) && bin_is_radical_[i] == false && bin_hydrogenation_level_[i] == "A")
					bin_primary_particle_index = i;

		if (bin_primary_particle_index == 0)
			for (unsigned int i = 0; i < bin_indices_.size(); i++)
				if (bin_section_[i] == (bin_minimum_aggregates_section_ - 1) && bin_is_radical_[i] == true && bin_hydrogenation_level_[i] == "A")
					bin_primary_particle_index = i;

		if (bin_primary_particle_index == 0)
			OpenSMOKE::FatalErrorMessage("No primary particle was identified in the provided kinetic mechanism");

		std::cout << std::endl;
		std::cout << "Primary particle " << std::endl;
		std::cout << "------------------------------------------------" << std::endl;
		std::cout << " * Name:          " << bin_names_[bin_primary_particle_index] << std::endl;
		std::cout << " * Mass (ng):     " << bin_m_[bin_primary_particle_index]*1.e12 << std::endl;
		std::cout << " * Diameter (nm): " << bin_ds_[bin_primary_particle_index]*1.e9 << std::endl;
		std::cout << std::endl;


		// Collisional diameters
		bin_np_.resize(bin_indices_.size());	std::fill(bin_np_.begin(), bin_np_.end(), 0.);
		bin_dc_.resize(bin_indices_.size());	std::fill(bin_dc_.begin(), bin_dc_.end(), 0.);
		bin_d_.resize(bin_indices_.size());		std::fill(bin_d_.begin(), bin_d_.end(), 0.);
		for (unsigned int i = 0; i < bin_indices_.size(); i++)
		{
			if (bin_section_[i] >= bin_minimum_aggregates_section_)
			{
				bin_np_[i] = bin_m_[i] / bin_m_[bin_primary_particle_index];
				bin_dc_[i] = std::sqrt(5. / 3.) * bin_ds_[bin_primary_particle_index] * std::pow(bin_np_[i] / std::pow(1. + 2. / Df_, Df_ / 2.), 1. / Df_);
				bin_d_[i] = bin_dc_[i];
			}
			else
			{
				bin_d_[i] = bin_ds_[i];
			}
		}

		// Sub-classes: spherical particles and aggregates
		for (unsigned int i = 0; i < bin_indices_.size(); i++)
		{
			if (bin_section_[i] >= bin_minimum_spheres_section_)
			{
				if (bin_section_[i] >= bin_minimum_aggregates_section_)
				{
					bin_indices_large_aggregates_.push_back(i);
					bin_indices_large_aggregates_global_.push_back(bin_indices_[i]);
				}
				else
				{
					bin_indices_large_spherical_.push_back(i);
					bin_indices_large_spherical_global_.push_back(bin_indices_[i]);
				}
			}
		}

		// Physical diffusivities (correction factors)
		{
			int binReference = -1;
			double MWReference = 0.;
			double MWCut = 0.;
			for (unsigned int i=0;i<bin_section_.size();i++)
				if (bin_section_[i] == (physical_diffusion_bin_to_start_-1) && bin_hydrogenation_level_[i] == "A")
				{
					binReference = i;
					MWReference = bin_mw_[binReference];
					MWCut = MWReference;
					break;
				}

			if (binReference >= 0)
			{
				// Identify the molecular weight at which to cut the diffusion correction
				for (unsigned int i=0;i<bin_section_.size();i++)
					if (bin_section_[i] == (physical_diffusion_bin_to_cut_) && bin_hydrogenation_level_[i] == "A")
					{
						MWCut = bin_mw_[i];
						break;
					}

				// Index of reference BIN
				bin_physical_diffusion_reference_species_ = bin_indices_[binReference];

				// Correction coefficients
				bin_physical_diffusion_correction_factors_.resize(bin_indices_.size());
				std::fill(bin_physical_diffusion_correction_factors_.begin(), bin_physical_diffusion_correction_factors_.end(), 0.);
				for (unsigned int i = 0; i < bin_indices_.size(); i++)
				{
					if (bin_section_[i] >= physical_diffusion_bin_to_start_)
					{
						const double MWratio = std::min( bin_mw_[i] / MWReference, MWCut / MWReference );
						const double teta = std::pow(MWratio, physical_diffusion_exp_);
						bin_physical_diffusion_correction_factors_[i] = teta;
					}
				}

				std::cout << std::endl;
				std::cout << "Physical diffusion reference BIN " << std::endl;
				std::cout << "------------------------------------------------" << std::endl;
				std::cout << " * Name:          " << bin_names_[binReference] << std::endl;
				std::cout << " * Mass (ng):     " << bin_m_[binReference] * 1.e12 << std::endl;
				std::cout << " * Diameter (nm): " << bin_ds_[binReference] * 1.e9 << std::endl;
				std::cout << std::endl;
			}
		}

		// Extract relevant data about the BINs
		bin_hydrogenation_level_unique_values_ = extract_unique_values(bin_hydrogenation_level_);

		// Sub-classes: single sections (molecular+radical, any hydrogenation level)
		const unsigned int number_sections = *std::max_element(std::begin(bin_section_), std::end(bin_section_));
		bin_indices_sections_.resize(number_sections);
		bin_indices_sections_global_.resize(number_sections);
		for (unsigned int i = 0; i < bin_indices_.size(); i++)
		{
			bin_indices_sections_[bin_section_[i]-1].push_back(i);
			bin_indices_sections_global_[bin_section_[i]-1].push_back(bin_indices_[i]);
		}

		// Write species names for each section
		std::cout << std::endl;
		for (unsigned int i = 0; i < number_sections; i++)
		{
			std::cout << "Section " << i+1 << std::endl;
			std::cout << "------------------------------------------------" << std::endl;
			for (unsigned int j = 0; j < bin_indices_sections_[i].size(); j++)
				std::cout << bin_names_[bin_indices_sections_[i][j]] << std::endl;
		}

		// Write summary
		WriteBinData();

		// Final operations
		{
			// Memory allocation
			bin_omega_.resize(bin_indices_.size());
			bin_x_.resize(bin_indices_.size());
			bin_fv_.resize(bin_indices_.size());
			bin_rho_.resize(bin_indices_.size());
			bin_N_.resize(bin_indices_.size());

			
			{
				// Number of baskets (number of C atoms)
				bin_baskets_c_.resize(number_sections);
				for (unsigned int i = 0; i < bin_indices_.size(); i++)
					bin_baskets_c_[bin_section_[i]-1] = bin_c_[i];

				std::sort(bin_baskets_c_.begin(), bin_baskets_c_.end());
				bin_baskets_indices_.resize(bin_baskets_c_.size());


				for (unsigned int i = 0; i < bin_indices_.size(); i++)
					for (unsigned int k = 0; k < bin_baskets_c_.size(); k++)
					{
						if (bin_c_[i] == bin_baskets_c_[k])
							bin_baskets_indices_[k].push_back(i);
					}

				bin_baskets_d_.resize(bin_baskets_c_.size());
				bin_baskets_mw_.resize(bin_baskets_c_.size());
				bin_baskets_log10d_.resize(bin_baskets_c_.size());
				bin_baskets_dlog10d_.resize(bin_baskets_c_.size());
				dN_over_dlog10d_.resize(bin_baskets_c_.size());

				bin_baskets_V_.resize(bin_baskets_c_.size());
				bin_baskets_log10V_.resize(bin_baskets_c_.size());
				bin_baskets_dlog10V_.resize(bin_baskets_c_.size());
				dN_over_dlog10V_.resize(bin_baskets_c_.size());

				bin_baskets_m_.resize(bin_baskets_c_.size());
				bin_baskets_log10m_.resize(bin_baskets_c_.size());
				bin_baskets_dlog10m_.resize(bin_baskets_c_.size());
				dN_over_dlog10m_.resize(bin_baskets_c_.size());

				bin_baskets_N_.resize(bin_baskets_c_.size());
				bin_baskets_fv_.resize(bin_baskets_c_.size());
				bin_baskets_rho_.resize(bin_baskets_c_.size());
				bin_baskets_x_.resize(bin_baskets_c_.size());
				bin_baskets_omega_.resize(bin_baskets_c_.size());
			}

			for (unsigned int k = 0; k < bin_baskets_c_.size(); k++)
			{
				bin_baskets_mw_[k] = 0.;
				bin_baskets_d_[k] = 0.;
				bin_baskets_m_[k] = 0.;
				bin_baskets_V_[k] = 0.;

				for (unsigned int i = 0; i < bin_baskets_indices_[k].size(); i++)
				{
					const unsigned int j = bin_baskets_indices_[k][i];
					bin_baskets_d_[k] += bin_d_[j];
					bin_baskets_mw_[k] += bin_mw_[j];
					bin_baskets_V_[k] += bin_V_[j];
					bin_baskets_m_[k] += bin_m_[j];
				}

				bin_baskets_d_[k] /= double(bin_baskets_indices_[k].size());
				bin_baskets_mw_[k] /= double(bin_baskets_indices_[k].size());
				bin_baskets_V_[k] /= double(bin_baskets_indices_[k].size());
				bin_baskets_m_[k] /= double(bin_baskets_indices_[k].size());

				bin_baskets_log10d_[k] = std::log10(bin_baskets_d_[k]);
				bin_baskets_log10V_[k] = std::log10(bin_baskets_V_[k]);
				bin_baskets_log10m_[k] = std::log10(bin_baskets_m_[k]);
			}

			// Intervals
			bin_baskets_dlog10d_[0] = ((bin_baskets_log10d_[1] + bin_baskets_log10d_[0]) / 2. - bin_baskets_log10d_[0])*2.;
			bin_baskets_dlog10V_[0] = ((bin_baskets_log10V_[1] + bin_baskets_log10V_[0]) / 2. - bin_baskets_log10V_[0])*2.;
			bin_baskets_dlog10m_[0] = ((bin_baskets_log10m_[1] + bin_baskets_log10m_[0]) / 2. - bin_baskets_log10m_[0])*2.;

			for (unsigned int k = 1; k < bin_baskets_c_.size() - 1; k++)
			{
				bin_baskets_dlog10d_[k] = (bin_baskets_log10d_[k + 1] + bin_baskets_log10d_[k]) / 2. - (bin_baskets_log10d_[k] + bin_baskets_log10d_[k - 1]) / 2.;
				bin_baskets_dlog10V_[k] = (bin_baskets_log10V_[k + 1] + bin_baskets_log10V_[k]) / 2. - (bin_baskets_log10V_[k] + bin_baskets_log10V_[k - 1]) / 2.;
				bin_baskets_dlog10m_[k] = (bin_baskets_log10m_[k + 1] + bin_baskets_log10m_[k]) / 2. - (bin_baskets_log10m_[k] + bin_baskets_log10m_[k - 1]) / 2.;
			}

			const unsigned int k = static_cast<unsigned int>(bin_baskets_c_.size() - 1);
			bin_baskets_dlog10d_[k] = ((bin_baskets_log10d_[k] + bin_baskets_log10d_[k - 1]) / 2. - bin_baskets_log10d_[k - 1])*2.;
			bin_baskets_dlog10V_[k] = ((bin_baskets_log10V_[k] + bin_baskets_log10V_[k - 1]) / 2. - bin_baskets_log10V_[k - 1])*2.;
			bin_baskets_dlog10m_[k] = ((bin_baskets_log10m_[k] + bin_baskets_log10m_[k - 1]) / 2. - bin_baskets_log10m_[k - 1])*2.;
		}

		// Soot dimer
		for(unsigned int i=0;i<nspecies_;i++)
			if (thermo_.NamesOfSpecies()[i].compare(0, bin_minimum_spheres_.size(), bin_minimum_spheres_) == 0)
				soot_dimer_indices_global_.push_back(i);

		// PAH (340 nm, 1/2 aromatic rings)
		pah_1_2_rings_indices_global_.resize(0);
		if ( thermo_.IndexOfSpeciesWithoutError("C6H6") > 0 )		pah_1_2_rings_indices_global_.push_back( thermo_.IndexOfSpeciesWithoutError("C6H6")-1 );
		if ( thermo_.IndexOfSpeciesWithoutError("C7H8") > 0 )		pah_1_2_rings_indices_global_.push_back( thermo_.IndexOfSpeciesWithoutError("C7H8")-1 );
		if ( thermo_.IndexOfSpeciesWithoutError("INDENE") > 0 )		pah_1_2_rings_indices_global_.push_back( thermo_.IndexOfSpeciesWithoutError("INDENE")-1 );
		if ( thermo_.IndexOfSpeciesWithoutError("C10H8") > 0 )		pah_1_2_rings_indices_global_.push_back( thermo_.IndexOfSpeciesWithoutError("C10H8")-1 );
		if ( thermo_.IndexOfSpeciesWithoutError("C12H8") > 0 )		pah_1_2_rings_indices_global_.push_back( thermo_.IndexOfSpeciesWithoutError("C12H8")-1 );
		if ( thermo_.IndexOfSpeciesWithoutError("BIPHENYL") > 0 )	pah_1_2_rings_indices_global_.push_back( thermo_.IndexOfSpeciesWithoutError("BIPHENYL")-1 );
		if ( thermo_.IndexOfSpeciesWithoutError("FLUORENE") > 0 )	pah_1_2_rings_indices_global_.push_back( thermo_.IndexOfSpeciesWithoutError("FLUORENE")-1 );
		if ( thermo_.IndexOfSpeciesWithoutError("C6H5C2H") > 0 )	pah_1_2_rings_indices_global_.push_back( thermo_.IndexOfSpeciesWithoutError("C6H5C2H")-1 );
		if ( thermo_.IndexOfSpeciesWithoutError("C6H5C2H3") > 0 )	pah_1_2_rings_indices_global_.push_back( thermo_.IndexOfSpeciesWithoutError("C6H5C2H3")-1 );
		if ( thermo_.IndexOfSpeciesWithoutError("C6H5C2H5") > 0 )	pah_1_2_rings_indices_global_.push_back( thermo_.IndexOfSpeciesWithoutError("C6H5C2H5")-1 );
		if ( thermo_.IndexOfSpeciesWithoutError("C6H5CH2C6H5") > 0)	pah_1_2_rings_indices_global_.push_back( thermo_.IndexOfSpeciesWithoutError("C6H5CH2C6H5")-1 );
		if ( thermo_.IndexOfSpeciesWithoutError("C10H7CH3") > 0 )	pah_1_2_rings_indices_global_.push_back( thermo_.IndexOfSpeciesWithoutError("C10H7CH3")-1 );

		// From mechanisms > 2020
		if (thermo_.IndexOfSpeciesWithoutError("C6H5") > 0)		pah_1_2_rings_indices_global_.push_back(thermo_.IndexOfSpeciesWithoutError("C6H5") - 1);
		if (thermo_.IndexOfSpeciesWithoutError("C7H7") > 0)		pah_1_2_rings_indices_global_.push_back(thermo_.IndexOfSpeciesWithoutError("C7H7") - 1);
		if (thermo_.IndexOfSpeciesWithoutError("C10H7") > 0)	pah_1_2_rings_indices_global_.push_back(thermo_.IndexOfSpeciesWithoutError("C10H7") - 1);
		if (thermo_.IndexOfSpeciesWithoutError("C12H7") > 0)	pah_1_2_rings_indices_global_.push_back(thermo_.IndexOfSpeciesWithoutError("C12H7") - 1);
		if (thermo_.IndexOfSpeciesWithoutError("C12H9") > 0)	pah_1_2_rings_indices_global_.push_back(thermo_.IndexOfSpeciesWithoutError("C12H9") - 1);
		if (thermo_.IndexOfSpeciesWithoutError("CH3C6H4") > 0)	pah_1_2_rings_indices_global_.push_back(thermo_.IndexOfSpeciesWithoutError("CH3C6H4") - 1);
		if (thermo_.IndexOfSpeciesWithoutError("C6H4C2H") > 0)	pah_1_2_rings_indices_global_.push_back(thermo_.IndexOfSpeciesWithoutError("C6H4C2H") - 1);
		if (thermo_.IndexOfSpeciesWithoutError("C6H5C2H2") > 0)	pah_1_2_rings_indices_global_.push_back(thermo_.IndexOfSpeciesWithoutError("C6H5C2H2") - 1);
		if (thermo_.IndexOfSpeciesWithoutError("C10H7CH2") > 0)	pah_1_2_rings_indices_global_.push_back(thermo_.IndexOfSpeciesWithoutError("C10H7CH2") - 1);
		if (thermo_.IndexOfSpeciesWithoutError("C10H6CH3") > 0)	pah_1_2_rings_indices_global_.push_back(thermo_.IndexOfSpeciesWithoutError("C10H6CH3") - 1);
		if (thermo_.IndexOfSpeciesWithoutError("XYLENE") > 0)	pah_1_2_rings_indices_global_.push_back(thermo_.IndexOfSpeciesWithoutError("XYLENE") - 1);
		if (thermo_.IndexOfSpeciesWithoutError("RXYLENE") > 0)	pah_1_2_rings_indices_global_.push_back(thermo_.IndexOfSpeciesWithoutError("RXYLENE") - 1);
		if (thermo_.IndexOfSpeciesWithoutError("INDENYL") > 0)	pah_1_2_rings_indices_global_.push_back(thermo_.IndexOfSpeciesWithoutError("INDENYL") - 1);

		// PAH (340 nm, 3/4 aromatic rings)
		pah_3_4_rings_indices_global_.resize(0);
		if ( thermo_.IndexOfSpeciesWithoutError("C14H10") > 0 )		pah_3_4_rings_indices_global_.push_back( thermo_.IndexOfSpeciesWithoutError("C14H10")-1 );
		if ( thermo_.IndexOfSpeciesWithoutError("C16H10") > 0 )		pah_3_4_rings_indices_global_.push_back( thermo_.IndexOfSpeciesWithoutError("C16H10")-1 );

		// From mechanisms > 2020
		if (thermo_.IndexOfSpeciesWithoutError("C18H10") > 0)		pah_3_4_rings_indices_global_.push_back(thermo_.IndexOfSpeciesWithoutError("C18H10") - 1);
		if (thermo_.IndexOfSpeciesWithoutError("C18H14") > 0)		pah_3_4_rings_indices_global_.push_back(thermo_.IndexOfSpeciesWithoutError("C18H14") - 1);
		if (thermo_.IndexOfSpeciesWithoutError("C14H9") > 0)		pah_3_4_rings_indices_global_.push_back(thermo_.IndexOfSpeciesWithoutError("C14H9") - 1);
		if (thermo_.IndexOfSpeciesWithoutError("C16H9") > 0)		pah_3_4_rings_indices_global_.push_back(thermo_.IndexOfSpeciesWithoutError("C16H9") - 1);
		if (thermo_.IndexOfSpeciesWithoutError("C18H9") > 0)		pah_3_4_rings_indices_global_.push_back(thermo_.IndexOfSpeciesWithoutError("C18H9") - 1);
		if (thermo_.IndexOfSpeciesWithoutError("C6H5C2H4C6H5") > 0)	pah_3_4_rings_indices_global_.push_back(thermo_.IndexOfSpeciesWithoutError("C6H5C2H4C6H5") - 1);
		if (thermo_.IndexOfSpeciesWithoutError("C6H5CH2C6H5") > 0)	pah_3_4_rings_indices_global_.push_back(thermo_.IndexOfSpeciesWithoutError("C6H5CH2C6H5") - 1);

		// PAH (340 nm, more than 4 aromatic rings)
		// Actually, they are equivalent to small BINs
		pah_large_indices_global_.resize(0);
		pah_large_indices_global_ = bin_indices_small_global_;	

		// Summary PAH data
		std::cout << "PAHs with 1/2 rings (340 nm)" << std::endl;
		std::cout << "------------------------------------------------" << std::endl;
		for(unsigned int i=0;i<pah_1_2_rings_indices_global_.size();i++)
			std::cout << thermo_.NamesOfSpecies()[pah_1_2_rings_indices_global_[i]] << " " << pah_1_2_rings_indices_global_[i] << std::endl;	
		std::cout << std::endl;

		std::cout << "PAHs with 3/4 rings (400 nm)" << std::endl;
		std::cout << "------------------------------------------------" << std::endl;
		for(unsigned int i=0;i<pah_3_4_rings_indices_global_.size();i++)
			std::cout << thermo_.NamesOfSpecies()[pah_3_4_rings_indices_global_[i]] << " " << pah_3_4_rings_indices_global_[i] << std::endl;	
		std::cout << std::endl;

		std::cout << "PAHs with more than 4 rings (500 nm)" << std::endl;
		std::cout << "------------------------------------------------" << std::endl;
		for(unsigned int i=0;i<pah_large_indices_global_.size();i++)
			std::cout << thermo_.NamesOfSpecies()[pah_large_indices_global_[i]] << " " << pah_large_indices_global_[i] << std::endl;	
		std::cout << std::endl;

		std::cout << "Soot dimers" << std::endl;
		std::cout << "------------------------------------------------" << std::endl;
		for(unsigned int i=0;i<soot_dimer_indices_global_.size();i++)
			std::cout << thermo_.NamesOfSpecies()[soot_dimer_indices_global_[i]] << " " << soot_dimer_indices_global_[i] << std::endl;	
		std::cout << std::endl;

		std::cout << "Soot particles" << std::endl;
		std::cout << "------------------------------------------------" << std::endl;
		for(unsigned int i=0;i<bin_indices_large_global_.size();i++)
			std::cout << thermo_.NamesOfSpecies()[bin_indices_large_global_[i]] << " " << bin_indices_large_global_[i] << std::endl;	
		std::cout << std::endl;

		std::cout << "Soot spherical particles" << std::endl;
		std::cout << "------------------------------------------------" << std::endl;
		for (unsigned int i = 0; i<bin_indices_large_spherical_global_.size(); i++)
			std::cout << thermo_.NamesOfSpecies()[bin_indices_large_spherical_global_[i]] << " " << bin_indices_large_spherical_global_[i] << std::endl;
		std::cout << std::endl;

		std::cout << "Soot aggregates" << std::endl;
		std::cout << "------------------------------------------------" << std::endl;
		for (unsigned int i = 0; i<bin_indices_large_aggregates_global_.size(); i++)
			std::cout << thermo_.NamesOfSpecies()[bin_indices_large_aggregates_global_[i]] << " " << bin_indices_large_aggregates_global_[i] << std::endl;
		std::cout << std::endl;
	}

	double PolimiSoot_Analyzer::fv_large(const double rhoGas, const Eigen::VectorXd &omegaGas) const
	{
		double sum = 0.;
		for (unsigned int i = 0; i < bin_indices_large_global_.size(); i++)
		{
			int j = bin_indices_large_global_[i];
			sum += rhoGas*omegaGas(j) / bin_density_large_[i];
		}
		return sum;
	}

	double PolimiSoot_Analyzer::rho_large(const double rhoGas, const Eigen::VectorXd &omegaGas) const
	{
		double sum = 0.;
		for (unsigned int i = 0; i < bin_indices_large_global_.size(); i++)
		{
			int j = bin_indices_large_global_[i];
			sum += rhoGas*omegaGas(j);
		}
		return sum;
	}

	double PolimiSoot_Analyzer::fv_small(const double rhoGas, const Eigen::VectorXd &omegaGas) const
	{
		double sum = 0.;
		for (unsigned int i = 0; i < bin_indices_small_global_.size(); i++)
		{
			int j = bin_indices_small_global_[i];
			sum += rhoGas*omegaGas(j) / bin_density_small_[i];
		}
		return sum;
	}

	double PolimiSoot_Analyzer::fv(const double rhoGas, const Eigen::VectorXd &omegaGas) const
	{
		return (fv_small(rhoGas, omegaGas) + fv_large(rhoGas, omegaGas));
	}

	double PolimiSoot_Analyzer::omega_pah_1_2_rings(const Eigen::VectorXd &omegaGas) const
	{
		double sum = 0.;
		for (unsigned int i = 0; i < pah_1_2_rings_indices_global_.size(); i++)
		{
			const unsigned int j = pah_1_2_rings_indices_global_[i];
			sum += omegaGas(j);
		}
		return sum;
	}

	double PolimiSoot_Analyzer::omega_pah_3_4_rings(const Eigen::VectorXd &omegaGas) const
	{
		double sum = 0.;
		for (unsigned int i = 0; i < pah_3_4_rings_indices_global_.size(); i++)
		{
			const unsigned int j = pah_3_4_rings_indices_global_[i];
			sum += omegaGas(j);
		}
		return sum;
	}

	double PolimiSoot_Analyzer::omega_pah_large(const Eigen::VectorXd &omegaGas) const
	{
		double sum = 0.;
		for (unsigned int i = 0; i < pah_large_indices_global_.size(); i++)
		{
			const unsigned int j = pah_large_indices_global_[i];
			sum += omegaGas(j);
		}
		return sum;
	}

	double PolimiSoot_Analyzer::x_pah_1_2_rings(const Eigen::VectorXd& xGas) const
	{
		double sum = 0.;
		for (unsigned int i = 0; i < pah_1_2_rings_indices_global_.size(); i++)
		{
			const unsigned int j = pah_1_2_rings_indices_global_[i];
			sum += xGas(j);
		}
		return sum;
	}

	double PolimiSoot_Analyzer::x_pah_3_4_rings(const Eigen::VectorXd& xGas) const
	{
		double sum = 0.;
		for (unsigned int i = 0; i < pah_3_4_rings_indices_global_.size(); i++)
		{
			const unsigned int j = pah_3_4_rings_indices_global_[i];
			sum += xGas(j);
		}
		return sum;
	}

	double PolimiSoot_Analyzer::x_pah_large(const Eigen::VectorXd& xGas) const
	{
		double sum = 0.;
		for (unsigned int i = 0; i < pah_large_indices_global_.size(); i++)
		{
			const unsigned int j = pah_large_indices_global_[i];
			sum += xGas(j);
		}
		return sum;
	}

	double PolimiSoot_Analyzer::planck_coefficient(const double rhoGas, const double T, const Eigen::VectorXd &omegaGas) const
	{
		if (soot_planck_coefficient_ == SootPlanckCoefficient::SOOT_PLANCK_COEFFICIENT_SMOOKE)
			return (1307.*fv_large(rhoGas, omegaGas)*T);							// [1/m]	(Smooke et al. Combustion and Flame 2009)
		else if (soot_planck_coefficient_ == SootPlanckCoefficient::SOOT_PLANCK_COEFFICIENT_KENT)
			return (2262.*fv_large(rhoGas, omegaGas)*T);							// [1/m]	(Kent al. Combustion and Flame 1990)
		else if (soot_planck_coefficient_ == SootPlanckCoefficient::SOOT_PLANCK_COEFFICIENT_SAZHIN)
			return (1232.*rho_large(rhoGas, omegaGas)*(1. + 4.8e-4*(T - 2000.)));	// [1/m]	(Sazhin, Fluent 1994)
		else if (soot_planck_coefficient_ == SootPlanckCoefficient::SOOT_PLANCK_COEFFICIENT_HUBBARD)
			return (0.33*fv_large(rhoGas, omegaGas)*(-3.75e5 + 1735.*T));			// [1/m]	(Hubbard G.L., Tien C.L., Infrared Mean Absorption Coefficient of Luminous Flames and Smoke, Journal of Heat Transfer, 100, p.235, 1978)

		else
			return 0.;
	}

	void PolimiSoot_Analyzer::Analysis(const double T, const double P_Pa, const double rhoGas, const Eigen::VectorXd &omegaGas, const Eigen::VectorXd &xGas)
	{
		{
			const double small_eps = 1e-20;
			const double denominator_max = 1.e32;
			const double threshold = 1e-20;

			for (unsigned int i = 0; i < bin_indices_.size(); i++)
			{
				const unsigned int j = bin_indices_[i];
				bin_omega_[i] = omegaGas(j);				// mass fraction
				bin_x_[i] = xGas(j);						// mole fraction
				bin_rho_[i] = rhoGas*omegaGas(j);			// density [kg_soot/m3]
				bin_fv_[i] = bin_rho_[i] / bin_density_[i];	// volume fraction [m3_soot/m3]
				bin_N_[i] = bin_fv_[i] / bin_V_[i];			// [1/m3]
			}

			// Large PAHs (>4 rings)
			{
				fv_small_ = 0.;
				rho_small_ = 0.;
				N_small_ = 0.;
				omega_small_ = 0.;
				x_small_ = 0.;

				double atom_c = 0.;
				double atom_h = 0.;
				double atom_o = 0.;
				for (unsigned int i = 0; i < bin_indices_small_.size(); i++)
				{
					const unsigned int j = bin_indices_small_[i];
					fv_small_ += bin_fv_[j];
					rho_small_ += bin_rho_[j];
					N_small_ += bin_N_[j];
					omega_small_ += bin_omega_[j];
					x_small_ += bin_x_[j];
					atom_c += bin_x_[j] * bin_c_[j];
					atom_h += bin_x_[j] * bin_h_[j];
					atom_o += bin_x_[j] * bin_o_[j];
				}

				double denominator = 0.;
				if (fv_small_ < small_eps)	denominator = denominator_max;
				h_over_c_small_ = atom_h / (atom_c + denominator);
				o_over_c_small_ = atom_o / (atom_c + denominator);
				o_over_h_small_ = atom_o / (atom_h + denominator);

				d10_small_ = SootDiameters(1., 0., bin_indices_small_, threshold);
				d32_small_ = SootDiameters(3., 2., bin_indices_small_, threshold);
				d43_small_ = SootDiameters(4., 3., bin_indices_small_, threshold);
			}

			// Soot (total)
			{
				fv_large_ = 0.;
				rho_large_ = 0.;
				N_large_ = 0.;
				omega_large_ = 0.;
				x_large_ = 0.;

				double atom_c = 0.;
				double atom_h = 0.;
				double atom_o = 0.;
				for (unsigned int i = 0; i < bin_indices_large_.size(); i++)
				{
					const unsigned int j = bin_indices_large_[i];
					fv_large_ += bin_fv_[j];
					rho_large_ += bin_rho_[j];
					N_large_ += bin_N_[j];
					omega_large_ += bin_omega_[j];
					x_large_ += bin_x_[j];
					atom_c += bin_x_[j] * bin_c_[j];
					atom_h += bin_x_[j] * bin_h_[j];
					atom_o += bin_x_[j] * bin_o_[j];
				}

				double denominator = 0.;
				if (fv_large_ < small_eps)	denominator = denominator_max;
				h_over_c_large_ = atom_h / (atom_c + denominator);
				o_over_c_large_ = atom_o / (atom_c + denominator);
				o_over_h_large_ = atom_o / (atom_h + denominator);

				d10_large_ = SootDiameters(1., 0., bin_indices_large_, threshold);
				d32_large_ = SootDiameters(3., 2., bin_indices_large_, threshold);
				d43_large_ = SootDiameters(4., 3., bin_indices_large_, threshold);
			}

			// Soot: spherical particles
			{
				fv_large_spherical_ = 0.;
				rho_large_spherical_ = 0.;
				N_large_spherical_ = 0.;
				omega_large_spherical_ = 0.;
				x_large_spherical_ = 0.;
				
				double atom_c = 0.;
				double atom_h = 0.;
				double atom_o = 0.;
				for (unsigned int i = 0; i < bin_indices_large_spherical_.size(); i++)
				{
					const unsigned int j = bin_indices_large_spherical_[i];

					fv_large_spherical_ += bin_fv_[j];
					rho_large_spherical_ += bin_rho_[j];
					N_large_spherical_ += bin_N_[j];
					omega_large_spherical_ += bin_omega_[j];
					x_large_spherical_ += bin_x_[j];

					atom_c += bin_x_[j] * bin_c_[j];
					atom_h += bin_x_[j] * bin_h_[j];
					atom_o += bin_x_[j] * bin_o_[j];
				}

				double denominator = 0.;
				if (fv_large_spherical_ < small_eps)	denominator = denominator_max;
				h_over_c_large_spherical_ = atom_h / (atom_c + denominator);
				o_over_c_large_spherical_ = atom_o / (atom_c + denominator);
				o_over_h_large_spherical_ = atom_o / (atom_h + denominator);

				d10_large_spherical_ = SootDiameters(1., 0., bin_indices_large_spherical_, threshold);
				d32_large_spherical_ = SootDiameters(3., 2., bin_indices_large_spherical_, threshold);
				d43_large_spherical_ = SootDiameters(4., 3., bin_indices_large_spherical_, threshold);
			}

			// Soot: aggregates
			{
				fv_large_aggregates_ = 0.;
				rho_large_aggregates_ = 0.;
				N_large_aggregates_ = 0.;
				omega_large_aggregates_ = 0.;
				x_large_aggregates_ = 0.;

				double atom_c = 0.;
				double atom_h = 0.;
				double atom_o = 0.;
				for (unsigned int i = 0; i < bin_indices_large_aggregates_.size(); i++)
				{
					const unsigned int j = bin_indices_large_aggregates_[i];

					fv_large_aggregates_ += bin_fv_[j];
					rho_large_aggregates_ += bin_rho_[j];
					N_large_aggregates_ += bin_N_[j];
					omega_large_aggregates_ += bin_omega_[j];
					x_large_aggregates_ += bin_x_[j];

					atom_c += bin_x_[j] * bin_c_[j];
					atom_h += bin_x_[j] * bin_h_[j];
					atom_o += bin_x_[j] * bin_o_[j];
				}

				double denominator = 0.;
				if (fv_large_aggregates_ < small_eps)	denominator = denominator_max;
				h_over_c_large_aggregates_ = atom_h / (atom_c + denominator);
				o_over_c_large_aggregates_ = atom_o / (atom_c + denominator);
				o_over_h_large_aggregates_ = atom_o / (atom_h + denominator);

				d10_large_aggregates_ = SootDiameters(1., 0., bin_indices_large_aggregates_, threshold);
				d32_large_aggregates_ = SootDiameters(3., 2., bin_indices_large_aggregates_, threshold);
				d43_large_aggregates_ = SootDiameters(4., 3., bin_indices_large_aggregates_, threshold);
			}

			// Analysis of PAHs: mass fraction
			{
				omega_pah_1_2_rings_ = 0.;
				for (unsigned int i = 0; i < pah_1_2_rings_indices_global_.size(); i++)
				{
					const unsigned int j = pah_1_2_rings_indices_global_[i];
					omega_pah_1_2_rings_ += omegaGas(j);
				}

				omega_pah_3_4_rings_ = 0.;
				for (unsigned int i = 0; i < pah_3_4_rings_indices_global_.size(); i++)
				{
					const unsigned int j = pah_3_4_rings_indices_global_[i];
					omega_pah_3_4_rings_ += omegaGas(j);
				}

				omega_pah_large_ = 0.;
				for (unsigned int i = 0; i < pah_large_indices_global_.size(); i++)
				{
					const unsigned int j = pah_large_indices_global_[i];
					omega_pah_large_ += omegaGas(j);
				}
			}

			// Analysis of PAHs: mole fraction
			{
				x_pah_1_2_rings_ = 0.;
				for (unsigned int i = 0; i < pah_1_2_rings_indices_global_.size(); i++)
				{
					const unsigned int j = pah_1_2_rings_indices_global_[i];
					x_pah_1_2_rings_ += xGas(j);
				}

				x_pah_3_4_rings_ = 0.;
				for (unsigned int i = 0; i < pah_3_4_rings_indices_global_.size(); i++)
				{
					const unsigned int j = pah_3_4_rings_indices_global_[i];
					x_pah_3_4_rings_ +=xGas(j);
				}

				x_pah_large_ = 0.;
				for (unsigned int i = 0; i < pah_large_indices_global_.size(); i++)
				{
					const unsigned int j = pah_large_indices_global_[i];
					x_pah_large_ += omegaGas(j);
				}
			}

			// Analysis of PAHs: volume fraction
			{
				fv_pah_1_2_rings_ = omega_pah_1_2_rings_;
				fv_pah_3_4_rings_ = omega_pah_3_4_rings_;
				fv_pah_large_ = fv_small_;
			}

			// Analysis of PAHs: partial density (kg/m3)
			{
				rho_pah_1_2_rings_ = omega_pah_1_2_rings_ * rhoGas;
				rho_pah_3_4_rings_ = omega_pah_3_4_rings_ * rhoGas;
			}

			// Analysis of PAHs: number of molecules (#/m3)
			{
				const double Ctot = P_Pa / PhysicalConstants::R_J_kmol / T;
				N_pah_1_2_rings_ = Ctot*x_pah_1_2_rings_ * PhysicalConstants::Nav_kmol;
				N_pah_3_4_rings_ = Ctot*x_pah_3_4_rings_ * PhysicalConstants::Nav_kmol;
			}

		}
	}

	double PolimiSoot_Analyzer::SootDiameters(const double a, const double b, const std::vector<unsigned int>& bin_indices, const double threshold)
	{
		// (3,2) Surface mean diameter
		// (4,3) De Brouckere mean diameter

		double num = 0.;
		double den = 0.;
		double ntot = 0.;
		for (unsigned int i = 0; i < bin_indices.size(); i++)
		{
			const unsigned int j = bin_indices[i];
			const double n = bin_x_[j];
			const double d = bin_d_[j];

			ntot += n;
			num += n * std::pow(d, a);
			den += n * std::pow(d, b);
		}

		double dmean = 0.;
		if (ntot > threshold)
			dmean = std::pow(num / den, (1. / (a - b)));

		return dmean;
	}

	void PolimiSoot_Analyzer::Analysis(	OpenSMOKE::KineticsMap_CHEMKIN& kinetics, 
										const double T, const double P_Pa, const double rhoGas, const Eigen::VectorXd& omegaGas, const Eigen::VectorXd& xGas)
	{
		if (formation_rates_by_classes_ == true)
		{
			// Kinetic data
			Eigen::VectorXd r_complete(kinetics.NumberOfReactions());
			{
				thermo_.SetTemperature(T);
				thermo_.SetPressure(P_Pa);

				double MW;
				Eigen::VectorXd xGas_ = xGas;
				thermo_.MoleFractions_From_MassFractions(xGas_.data(), MW, omegaGas.data());
				const double cTot = P_Pa / (PhysicalConstants::R_J_kmol * T);
				const Eigen::VectorXd cGas = cTot * xGas;

				kinetics.SetTemperature(T);
				kinetics.SetPressure(P_Pa);
				kinetics.ReactionRates(cGas.data());
				kinetics.GiveMeReactionRates(r_complete.data());
			}

			Eigen::VectorXi mask_bin_indices_small_global(thermo_.NumberOfSpecies());
			mask_bin_indices_small_global.setZero();
			for (unsigned int i = 0; i < bin_indices_small_global_.size(); i++)
				mask_bin_indices_small_global(bin_indices_small_global_[i]) = 1;

			Eigen::VectorXi mask_bin_indices_large_global(thermo_.NumberOfSpecies());
			mask_bin_indices_large_global.setZero();
			for (unsigned int i = 0; i < bin_indices_large_global_.size(); i++)
				mask_bin_indices_large_global(bin_indices_large_global_[i]) = 1;

			Eigen::VectorXi mask_bin_indices_large_spherical_global(thermo_.NumberOfSpecies());
			mask_bin_indices_large_spherical_global.setZero();
			for (unsigned int i = 0; i < bin_indices_large_spherical_global_.size(); i++)
				mask_bin_indices_large_spherical_global(bin_indices_large_spherical_global_[i]) = 1;

			Eigen::VectorXi mask_bin_indices_large_aggregates_global(thermo_.NumberOfSpecies());
			mask_bin_indices_large_aggregates_global.setZero();
			for (unsigned int i = 0; i < bin_indices_large_aggregates_global_.size(); i++)
				mask_bin_indices_large_aggregates_global(bin_indices_large_aggregates_global_[i]) = 1;

			R_sum_small_.resize(soot_classes_->soot_class_names().size());
			R_sum_large_.resize(soot_classes_->soot_class_names().size());
			R_sum_large_spherical_.resize(soot_classes_->soot_class_names().size());
			R_sum_large_aggregates_.resize(soot_classes_->soot_class_names().size());
			Omega_sum_large_.resize(soot_classes_->soot_class_names().size());
			Omega_sum_small_.resize(soot_classes_->soot_class_names().size());
			Omega_sum_large_spherical_.resize(soot_classes_->soot_class_names().size());
			Omega_sum_large_aggregates_.resize(soot_classes_->soot_class_names().size());
			for (unsigned int k = 0; k < soot_classes_->soot_class_names().size(); k++)
			{
				Eigen::VectorXd r(kinetics.NumberOfReactions());
				r.setZero();
				for (unsigned int j = 0; j < soot_classes_->reaction_indices(k).size(); j++)
				{
					const unsigned int i = soot_classes_->reaction_indices(k)[j];
					r(i) = r_complete(i);
				}

				Eigen::VectorXd R_semi_complete(thermo_.NumberOfSpecies());
				kinetics.stoichiometry().FormationRatesFromReactionRates(R_semi_complete.data(),r.data());

				{
					Eigen::VectorXd R_small(thermo_.NumberOfSpecies());
					Eigen::VectorXd Omega_small(thermo_.NumberOfSpecies());
					for (unsigned int i = 0; i < thermo_.NumberOfSpecies(); i++)
					{
						R_small(i) = R_semi_complete(i) * mask_bin_indices_small_global[i];
						Omega_small(i) = R_small(i) * thermo_.MW(i);
					}
					R_sum_small_(k) = R_small.sum();
					Omega_sum_small_(k) = Omega_small.sum();
				}

				{
					Eigen::VectorXd R_large(thermo_.NumberOfSpecies());
					Eigen::VectorXd Omega_large(thermo_.NumberOfSpecies());
					for (unsigned int i = 0; i < thermo_.NumberOfSpecies(); i++)
					{
						R_large(i) = R_semi_complete(i) * mask_bin_indices_large_global[i];
						Omega_large(i) = R_large(i) * thermo_.MW(i);
					}
					R_sum_large_(k) = R_large.sum();
					Omega_sum_large_(k) = Omega_large.sum();
				}

				{
					Eigen::VectorXd R_large_spherical(thermo_.NumberOfSpecies());
					Eigen::VectorXd Omega_large_spherical(thermo_.NumberOfSpecies());
					for (unsigned int i = 0; i < thermo_.NumberOfSpecies(); i++)
					{
						R_large_spherical(i) = R_semi_complete(i) * mask_bin_indices_large_spherical_global[i];
						Omega_large_spherical(i) = R_large_spherical(i) * thermo_.MW(i);
					}
					R_sum_large_spherical_(k) = R_large_spherical.sum();
					Omega_sum_large_spherical_(k) = Omega_large_spherical.sum();
				}

				{
					Eigen::VectorXd R_large_aggregates(thermo_.NumberOfSpecies());
					Eigen::VectorXd Omega_large_aggregates(thermo_.NumberOfSpecies());
					for (unsigned int i = 0; i < thermo_.NumberOfSpecies(); i++)
					{
						R_large_aggregates(i) = R_semi_complete(i) * mask_bin_indices_large_aggregates_global[i];
						Omega_large_aggregates(i) = R_large_aggregates(i) * thermo_.MW(i);
					}
					R_sum_large_aggregates_(k) = R_large_aggregates.sum();
					Omega_sum_large_aggregates_(k) = Omega_large_aggregates.sum();
				}
			}
		}
	}

	void PolimiSoot_Analyzer::Distribution()
	{
		for (unsigned int k = 0; k < bin_baskets_c_.size(); k++)
		{
			bin_baskets_N_[k] = 0.;
			bin_baskets_fv_[k] = 0.;
			bin_baskets_rho_[k] = 0.;
			bin_baskets_omega_[k] = 0.;
			bin_baskets_x_[k] = 0.;

			for (unsigned int i = 0; i < bin_baskets_indices_[k].size(); i++)
			{
				const unsigned int j = bin_baskets_indices_[k][i];

				bin_baskets_N_[k] += bin_N_[j];
				bin_baskets_fv_[k] += bin_fv_[j];
				bin_baskets_rho_[k] += bin_rho_[j];
				bin_baskets_omega_[k] += bin_omega_[j];
				bin_baskets_x_[k] += bin_x_[j];
			}

			dN_over_dlog10d_[k] = bin_baskets_N_[k] / bin_baskets_dlog10d_[k];
			dN_over_dlog10V_[k] = bin_baskets_N_[k] / bin_baskets_dlog10V_[k];
			dN_over_dlog10m_[k] = bin_baskets_N_[k] / bin_baskets_dlog10m_[k];
		}

		{
			double small_eps = 1e-16;

			double m0 = 0.;
			double m1 = 0.;
			double m2 = 0.;
			double m3 = 0.;
			for (unsigned int k = 0; k < bin_baskets_c_.size(); k++)
			{
				m0 += bin_baskets_N_[k];
				m1 += bin_baskets_N_[k] * bin_baskets_d_[k];
				m2 += bin_baskets_N_[k] * bin_baskets_d_[k] * bin_baskets_d_[k];
				m3 += bin_baskets_N_[k] * bin_baskets_d_[k] * bin_baskets_d_[k] * bin_baskets_d_[k];
			}

			if (m0 < small_eps)	m0 = 1e32;
			if (m2 < small_eps)	m2 = 1e32;

			// Mean diameters [m]
			d10_N_large_ = m1 / m0;
			d32_N_large_ = m3 / m2;

			// Variance [m2]
			dvariance_N_ = 0.;
			for (unsigned int k = 0; k < bin_baskets_c_.size(); k++)
				dvariance_N_ += bin_baskets_N_[k] * std::pow(bin_baskets_d_[k] - d10_N_large_, 2.);
			dvariance_N_ /= m0;

			// Standard deviation [m]
			dstd_N_ = std::sqrt(dvariance_N_);
		}
	}

	void PolimiSoot_Analyzer::WriteDistribution(std::ofstream& fSootDistribution, const double t, const double x, const double y, const double z, const double T)
	{
		for (unsigned int k = 0; k < bin_baskets_c_.size(); k++)
		{
			fSootDistribution << std::scientific << std::setw(20) << std::left << t;
			fSootDistribution << std::scientific << std::setw(20) << std::left << x;
			fSootDistribution << std::scientific << std::setw(20) << std::left << y;
			fSootDistribution << std::scientific << std::setw(20) << std::left << z;
			fSootDistribution << std::scientific << std::setw(20) << std::left << T;

			fSootDistribution << std::scientific << std::setw(20) << std::left << bin_baskets_d_[k] * 1e9;		// [nm]
			fSootDistribution << std::scientific << std::setw(20) << std::left << bin_baskets_m_[k] * 1e12;		// [ng]
			fSootDistribution << std::scientific << std::setw(20) << std::left << bin_baskets_V_[k] * 1e27;		// [nm3]

			fSootDistribution << std::scientific << std::setw(20) << std::left << std::log10(bin_baskets_d_[k]);
			fSootDistribution << std::scientific << std::setw(20) << std::left << std::log10(bin_baskets_m_[k]);
			fSootDistribution << std::scientific << std::setw(20) << std::left << std::log10(bin_baskets_V_[k]);

			fSootDistribution << std::scientific << std::setw(20) << std::left << std::log10(bin_baskets_dlog10d_[k]);
			fSootDistribution << std::scientific << std::setw(20) << std::left << std::log10(bin_baskets_dlog10m_[k]);
			fSootDistribution << std::scientific << std::setw(20) << std::left << std::log10(bin_baskets_dlog10V_[k]);

			fSootDistribution << std::scientific << std::setw(20) << std::left << bin_baskets_fv_[k];
			fSootDistribution << std::scientific << std::setw(20) << std::left << bin_baskets_x_[k];
			fSootDistribution << std::scientific << std::setw(20) << std::left << bin_baskets_omega_[k];
			fSootDistribution << std::scientific << std::setw(20) << std::left << bin_baskets_rho_[k];
			fSootDistribution << std::scientific << std::setw(20) << std::left << bin_baskets_N_[k] / 1.e6; // [#/cm3]
			fSootDistribution << std::scientific << std::setw(20) << std::left << bin_baskets_N_[k] / bin_baskets_dlog10d_[k] / 1.e6; // [#/cm3]

			fSootDistribution << std::scientific << std::setw(20) << std::left << bin_baskets_fv_[k] / (fv_large_ + fv_small_);
			fSootDistribution << std::scientific << std::setw(20) << std::left << bin_baskets_x_[k] / (x_large_ + x_small_);
			fSootDistribution << std::scientific << std::setw(20) << std::left << bin_baskets_omega_[k] / (omega_large_ + omega_small_);
			fSootDistribution << std::scientific << std::setw(20) << std::left << bin_baskets_rho_[k] / (rho_large_ + rho_small_);
			fSootDistribution << std::scientific << std::setw(20) << std::left << bin_baskets_N_[k] / (N_large_ + N_small_);

			fSootDistribution << std::endl;
		}
	}

	void PolimiSoot_Analyzer::WriteDistributionLabel(std::ofstream &fSoot)
	{
		fSoot	<< std::setw(20) << std::left << "time[s](1)"
			<< std::setw(20) << std::left << "x[m](2)"
			<< std::setw(20) << std::left << "y[m](3)"
			<< std::setw(20) << std::left << "z[m](4)"
			<< std::setw(20) << std::left << "T[K](5)";

		fSoot	<< std::setw(20) << std::left << "d[nm](6)"
				<< std::setw(20) << std::left << "m[ng](7)"
				<< std::setw(20) << std::left << "V[nm3](8)"
				<< std::setw(20) << std::left << "log10d(9)"
				<< std::setw(20) << std::left << "log10m(10)"
				<< std::setw(20) << std::left << "log10V(11)"
				<< std::setw(20) << std::left << "dlog10d(12)"
				<< std::setw(20) << std::left << "dlog10m(13)"
				<< std::setw(20) << std::left << "dlog10V(14)";

		fSoot	<< std::setw(20) << std::left << "fv(15)"
				<< std::setw(20) << std::left << "x(16)"
				<< std::setw(20) << std::left << "y(17)"
				<< std::setw(20) << std::left << "rho[kg/m3](18)"
				<< std::setw(20) << std::left << "N[#/cm3](19)"
				<< std::setw(20) << std::left << "dN/dlgd10[cm-3](20)";

		fSoot	<< std::setw(20) << std::left << "fvN(21)"
				<< std::setw(20) << std::left << "xN(22)"
				<< std::setw(20) << std::left << "yN(23)"
				<< std::setw(20) << std::left << "rhoN[](24)"
				<< std::setw(20) << std::left << "NN[](25)";

		fSoot	<< std::endl << std::endl;
	}

	void PolimiSoot_Analyzer::WriteBinData()
	{
		std::ofstream fOutput;
		fOutput.open("PolimiSootModel.out", std::ios::out);
		fOutput.setf(std::ios::scientific);

		fOutput << std::setw(5) << std::left << "#";
		fOutput << std::setw(5) << std::left << "#";
		fOutput << std::setw(10) << std::left << "Name";
		fOutput << std::setw(10) << std::left << "Section";
		fOutput << std::setw(10) << std::left << "Radical";
		fOutput << std::setw(10) << std::left << "Hyd.Lev.";
		fOutput << std::setw(16) << std::fixed << std::left << std::setprecision(2) << "MW[kg/kmol]";	// [kg/kmol]
		fOutput << std::setw(12) << std::fixed << std::left << std::setprecision(5) << "d[nm]";			// [nm]
		fOutput << std::setw(12) << std::fixed << std::left << std::setprecision(5) << "dsph[nm]";		// [nm]
		fOutput << std::setw(12) << std::fixed << std::left << std::setprecision(5) << "dcol[nm]";		// [nm]
		fOutput << std::setw(16) << std::scientific << std::left << "m[ng]";							// [ng] 
		fOutput << std::setw(16) << std::scientific << std::left << "V[nm3]";							// [nm3]
		fOutput << std::setw(12) << std::scientific << std::left << "np[-]";							// [-]
		fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(0) << "C";	
		fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(0) << "H";
		fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(0) << "O";
		fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(4) << "H/C";
		fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(4) << "O/C";
		fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(4) << "O/H";
		fOutput << std::setw(18) << std::scientific << std::left << "Teta(corr.fact.)";					// [-]
		fOutput << std::endl;

		for (unsigned int i = 0; i < bin_indices_small_.size(); i++)
		{
			int j = bin_indices_small_[i];
			fOutput << std::setw(5) << std::left << j + 1;
			fOutput << std::setw(5) << std::left << i + 1;
			fOutput << std::setw(10) << std::left << bin_names_[j].c_str();
			fOutput << std::setw(10) << std::fixed << std::left << bin_section_[j];								// [-]
			fOutput << std::setw(10) << std::fixed << std::left << bin_is_radical_[j];							// [-]
			fOutput << std::setw(10) << std::fixed << std::left << bin_hydrogenation_level_[j];					// [-]
			fOutput << std::setw(16) << std::fixed << std::left << std::setprecision(2) << bin_mw_[j];			// [kg/kmol]
			fOutput << std::setw(12) << std::fixed << std::left << std::setprecision(5) << bin_d_[j]  * 1e9;	// [nm]
			fOutput << std::setw(12) << std::fixed << std::left << std::setprecision(5) << bin_ds_[j] * 1e9;	// [nm]
			fOutput << std::setw(12) << std::fixed << std::left << std::setprecision(5) << bin_dc_[j] * 1e9;	// [nm]
			fOutput << std::setw(16) << std::scientific << std::left << bin_m_[j] * 1e12;						// [ng] 
			fOutput << std::setw(16) << std::scientific << std::left << bin_V_[j] * 1e27;						// [nm3]
			fOutput << std::setw(12) << std::fixed << std::left << std::setprecision(2) << bin_np_[j];			// [-]
			fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(0) << bin_c_[j];
			fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(0) << bin_h_[j];
			fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(0) << bin_o_[j];
			fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(4) << bin_h_over_c_[j];
			fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(4) << bin_o_over_c_[j];
			fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(4) << bin_o_over_h_[j];
			if (bin_physical_diffusion_correction_factors_[j] > 0.)
				fOutput << std::setw(18) << std::scientific << std::left << bin_physical_diffusion_correction_factors_[j];
			else fOutput << std::setw(18) << std::scientific << std::left << "-";
			fOutput << std::endl;
		}
		fOutput << std::endl;

		for (unsigned int i = 0; i < bin_indices_large_.size(); i++)
		{
			int j = bin_indices_large_[i];
			fOutput << std::setw(5) << std::left << j + 1;
			fOutput << std::setw(5) << std::left << i + 1;
			fOutput << std::setw(10) << std::left << bin_names_[j].c_str();
			fOutput << std::setw(10) << std::fixed << std::left << bin_section_[j];								// [-]
			fOutput << std::setw(10) << std::fixed << std::left << bin_is_radical_[j];							// [-]
			fOutput << std::setw(10) << std::fixed << std::left << bin_hydrogenation_level_[j];					// [-]
			fOutput << std::setw(16) << std::fixed << std::left << std::setprecision(2) << bin_mw_[j];			// [kg/kmol]
			fOutput << std::setw(12) << std::fixed << std::left << std::setprecision(5) << bin_d_[j] * 1e9;		// [nm]
			fOutput << std::setw(12) << std::fixed << std::left << std::setprecision(5) << bin_ds_[j] * 1e9;	// [nm]
			fOutput << std::setw(12) << std::fixed << std::left << std::setprecision(5) << bin_dc_[j] * 1e9;	// [nm]
			fOutput << std::setw(16) << std::scientific << std::left << bin_m_[j] * 1e12;						// [ng] 
			fOutput << std::setw(16) << std::scientific << std::left << bin_V_[j] * 1e27;						// [nm3]
			fOutput << std::setw(12) << std::fixed << std::left << std::setprecision(2) << bin_np_[j];			// [-]
			fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(0) << bin_c_[j];
			fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(0) << bin_h_[j];
			fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(0) << bin_o_[j];
			fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(4) << bin_h_over_c_[j];
			fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(4) << bin_o_over_c_[j];
			fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(4) << bin_o_over_h_[j];
			if (bin_physical_diffusion_correction_factors_[j] > 0.)
				fOutput << std::setw(18) << std::scientific << std::left << bin_physical_diffusion_correction_factors_[j];
			else fOutput << std::setw(18) << std::scientific << std::left << "-";
			fOutput << std::endl;
		}

		fOutput.close();
	}

	void PolimiSoot_Analyzer::WriteLabelIntegralDataFile(std::ofstream& fOutput, unsigned int& count, const unsigned int width)
	{
		// Total soot (default: BIN5-BIN-Max)
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "fv(tot)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "x(tot)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "w(tot)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width+2, fOutput, "rho(tot)[kg/m3]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "N(tot)[#/cm3]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "H/C(tot)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "O/C(tot)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "O/H(tot)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "d10(tot)(nm)", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "d32(tot)(nm)", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "d43(tot)(nm)", count);

		// Spherical particles (default: BIN5-BIN12)
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "fv(sp)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "x(sp)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "w(sp)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width+2, fOutput, "rho(sp)[kg/m3]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "N(sp)[#/cm3]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "H/C(sp)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "O/C(sp)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "O/H(sp)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "d10(sp)(nm)", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "d32(sp)(nm)", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "d43(sp)(nm)", count);

		// Aggregates (default: BIN13-BIN-Max)
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "fv(agg)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "x(agg)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "w(agg)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width+2, fOutput, "rho(agg)[kg/m3]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "N(agg)[#/cm3]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "H/C(agg)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "O/C(agg)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "O/H(agg)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "d10(agg)(nm)", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "d32(agg)(nm)", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "d43(agg)(nm)", count);

		// Large soot precursors (PAHs 1-2rings)
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "fv(pah12)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "x(pah12)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "w(pah12)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width+2, fOutput, "rho(pah12)[kg/m3]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "N(pah12)[#/cm3]", count);

		// Large soot precursors (PAH>4rings)
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "fv(pah34)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "x(pah34)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "w(pah34)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width+2, fOutput, "rho(pah34)[kg/m3]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "N(pah34)[#/cm3]", count);

		// Large soot precursors (PAH>4rings)
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "fv(pah>4)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "x(pah>4)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "w(pah>4)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width+2, fOutput, "rho(pah>4)[kg/m3]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "N(pah>4)[#/cm3]", count);

		for (unsigned int k = 0; k < bin_baskets_indices_.size(); k++)
		{
			const std::string section = std::to_string(k + 1);
			OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "fv(" + section + ")[-]", count);
			OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "x(" + section + ")[-]", count);
			OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "w(" + section + ")[-]", count);
			OpenSMOKE::PrintTagOnASCIILabel(width+2, fOutput, "rho(" + section + ")[kg/m3]", count);
			OpenSMOKE::PrintTagOnASCIILabel(width,   fOutput, "N(" + section + ")[#/cm3]", count);
		}
	}

	void PolimiSoot_Analyzer::WriteIntegralDataFile(std::ofstream& fOutput, const unsigned int width)
	{
		// Total soot (default BIN5-BIN-Max)
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << fv_large();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << x_large();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << omega_large();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width+2) << rho_large();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << N_large() / 1.e6;
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << h_over_c_large();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << o_over_c_large();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << o_over_h_large();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << d10_large_ * 1.e9;
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << d32_large_ * 1.e9;
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << d43_large_ * 1.e9;

		// Spherical particles (default: BIN5-BIN12)
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << fv_large_spherical();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << x_large_spherical();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << omega_large_spherical();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width+2) << rho_large_spherical();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << N_large_spherical() / 1.e6;
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << h_over_c_large_spherical();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << o_over_c_large_spherical();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << o_over_h_large_spherical();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << d10_large_spherical_ * 1.e9;
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << d32_large_spherical_ * 1.e9;
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << d43_large_spherical_ * 1.e9;

		// Aggregates (default: BIN13-BIN-Max)
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << fv_large_aggregates();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << x_large_aggregates();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << omega_large_aggregates();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width+2) << rho_large_aggregates();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << N_large_aggregates() / 1.e6;
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << h_over_c_large_aggregates();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << o_over_c_large_aggregates();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << o_over_h_large_aggregates();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << d10_large_aggregates_ * 1.e9;
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << d32_large_aggregates_ * 1.e9;
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << d43_large_aggregates_ * 1.e9;

		// PAHs 1-2 rings
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << fv_pah_1_2_rings();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << x_pah_1_2_rings();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << omega_pah_1_2_rings();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width+2) << rho_pah_1_2_rings();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << N_pah_1_2_rings() / 1.e6;

		// PAHs 3-4 rings
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << fv_pah_3_4_rings();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << x_pah_3_4_rings();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << omega_pah_3_4_rings();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width+2) << rho_pah_3_4_rings();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << N_pah_3_4_rings() / 1.e6;

		// PAHs>4rings (large soot precursors)
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << fv_small();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << x_small();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << omega_small();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width+2) << rho_small();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << N_small() / 1.e6;

		// Individual sections
		for (unsigned int k = 0; k < bin_baskets_indices_.size(); k++)
		{
			double fv = 0.;
			double x = 0.;
			double w = 0.;
			double rho = 0.;
			double N = 0.;
			
			for (unsigned int i = 0; i < bin_baskets_indices_[k].size(); i++)
			{
				const unsigned int j = bin_baskets_indices_[k][i];

				fv += bin_fv_[j];
				x += bin_x_[j];
				w += bin_omega_[j];
				rho += bin_rho_[j];
				N += bin_N_[j];
			}

			fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << fv;
			fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << x;
			fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << w;
			fOutput << std::scientific << std::setprecision(9) << std::setw(width+2) << rho;
			fOutput << std::scientific << std::setprecision(9) << std::setw(width)   << N/1.e6;
		}
	}


	void PolimiSoot_Analyzer::WriteFormationSootClassesDataFile(std::ofstream& fOutput, const unsigned int width)
	{
		if (soot_classes_->is_active() == true)
		{
			fOutput << std::scientific << std::setprecision(9) << std::setw(width) << R_sum_large_.sum();
			fOutput << std::scientific << std::setprecision(9) << std::setw(width) << R_sum_large_spherical_.sum();
			fOutput << std::scientific << std::setprecision(9) << std::setw(width) << R_sum_large_aggregates_.sum();
			fOutput << std::scientific << std::setprecision(9) << std::setw(width) << R_sum_small_.sum();

			fOutput << std::scientific << std::setprecision(9) << std::setw(width) << Omega_sum_large_.sum();
			fOutput << std::scientific << std::setprecision(9) << std::setw(width) << Omega_sum_large_spherical_.sum();
			fOutput << std::scientific << std::setprecision(9) << std::setw(width) << Omega_sum_large_aggregates_.sum();
			fOutput << std::scientific << std::setprecision(9) << std::setw(width) << Omega_sum_small_.sum();

			for (unsigned int k = 0; k < soot_classes_->soot_class_names().size(); k++)
			{
				fOutput << std::scientific << std::setprecision(9) << std::setw(width) << R_sum_large_(k);
				fOutput << std::scientific << std::setprecision(9) << std::setw(width) << R_sum_large_spherical_(k);
				fOutput << std::scientific << std::setprecision(9) << std::setw(width) << R_sum_large_aggregates_(k);
				fOutput << std::scientific << std::setprecision(9) << std::setw(width) << R_sum_small_(k);
			}

			for (unsigned int k = 0; k < soot_classes_->soot_class_names().size(); k++)
			{
				fOutput << std::scientific << std::setprecision(9) << std::setw(width) << Omega_sum_large_(k);
				fOutput << std::scientific << std::setprecision(9) << std::setw(width) << Omega_sum_large_spherical_(k);
				fOutput << std::scientific << std::setprecision(9) << std::setw(width) << Omega_sum_large_aggregates_(k);
				fOutput << std::scientific << std::setprecision(9) << std::setw(width) << Omega_sum_small_(k);
			}
		}
	}
	
	void PolimiSoot_Analyzer::WriteLabelFormationSootClassesDataFile(std::ofstream& fOutput, unsigned int& count, const unsigned int width)
	{
		if (soot_classes_->is_active() == true)
		{
			OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "RTot(kmol/m3/s)", count);
			OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "RSp(kmol/m3/s)", count);
			OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "RAgg(kmol/m3/s)", count);
			OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "Rpah>4(kmol/m3/s)", count);

			OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "OmTot(kg/m3/s)", count);
			OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "OmSp(kg/m3/s)", count);
			OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "OmAgg(kg/m3/s)", count);
			OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "OmPah>4(kg/m3/s)", count);

			std::string label;
			for (unsigned int k = 0; k < soot_classes_->soot_class_names().size(); k++)
			{
				label = "RTot_" + soot_classes_->soot_class_names()[k];
				OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, label, count);

				label = "RSp_" + soot_classes_->soot_class_names()[k];
				OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, label, count);

				label = "RAgg_" + soot_classes_->soot_class_names()[k];
				OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, label, count);

				label = "Rpah>4_" + soot_classes_->soot_class_names()[k];
				OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, label, count);
			}

			for (unsigned int k = 0; k < soot_classes_->soot_class_names().size(); k++)
			{
				label = "OmTot_" + soot_classes_->soot_class_names()[k];
				OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, label, count);

				label = "OmSp_" + soot_classes_->soot_class_names()[k];
				OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, label, count);

				label = "OmAgg_" + soot_classes_->soot_class_names()[k];
				OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, label, count);

				label = "OmPah>4_" + soot_classes_->soot_class_names()[k];
				OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, label, count);
			}
		}
	}

	bool PolimiSoot_Analyzer::ClassesCorrectionCoefficients(std::vector<unsigned int>& indices, std::vector<double>& coefficients)
	{
		if (soot_classes_->is_active() == false)
		{
			return false;
		}
		else
		{
			indices.resize(0);
			coefficients.resize(0);
			for (unsigned int j = 0; j < list_class_corrections_names_.size(); j++)
			{
				const unsigned int k = soot_classes_->index_from_class_name(list_class_corrections_names_[j]);
				for (unsigned int i = 0; i < soot_classes_->reaction_indices(k).size(); i++)
				{
					indices.push_back(soot_classes_->reaction_indices(k)[i]);
					coefficients.push_back(list_class_corrections_coefficients_[j]);
				}
			}
			return true;
		}
	}
}
