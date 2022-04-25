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
|   Copyright(C) 2020 Alberto Cuoci                                       |
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

#include <queue>

// Default name of main dictionary
const Foam::word Foam::CVI::dictName("CVI");

// Constructor
Foam::CVI::CVI(const Foam::fvMesh& mesh, OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap)
:
	IOdictionary
	(
		IOobject
		(
			dictName,
			mesh.time().constant(),
			mesh,
			IOobject::MUST_READ_IF_MODIFIED,
			IOobject::NO_WRITE
		)
	),

	thermodynamicsMap_(thermodynamicsMap), 

	kineticsMap_(kineticsMap)
{

	
	std::vector< Eigen::VectorXi > disk_topology_indices;
	std::vector< Eigen::VectorXd > disk_source_terms;
	
	Foam::string folder_disks = this->lookup("folder");
	path_folder_disks_ = folder_disks;

	diskPreprocess_ = Switch(this->lookupOrDefault(word("preProcess"), word("on")));
	list_disks_ = readList<int>(this->lookup("disks"));
	exclude_negative_disk_source_terms_ = Switch(this->lookupOrDefault(word("excludeNegativeSourceTerms"), word("off")));
	time_target_ = this->lookup("targetTime");
	correction_coefficient_ = readScalar(this->lookup("correctionCoefficient"));
}

Foam::CVI::~CVI()
{}

void Foam::CVI::Setup(const Foam::fvMesh& mesh)
{
	Info << "CVI Setup" << endl;
}

void Foam::CVI::Analysis(const Foam::fvMesh& mesh)
{
	Info << "CVI Analysis" << endl;	

	if (diskPreprocess_ == true)
	{
		std::vector< Eigen::VectorXi > disk_topology_indices;
		std::vector< Eigen::VectorXd > disk_source_terms;

		disk_topology_indices.resize(list_disks_.size());
		disk_source_terms.resize(list_disks_.size());
		std::vector<double> disk_topology_volume(list_disks_.size());

		// Read topology
		{
			boost::property_tree::ptree ptree;
    			boost::property_tree::read_xml( (path_folder_disks_ / "DiskTopology.xml").string(), ptree );

			try
			{
				const int ndisks = ptree.get<int>("opensmoke.Disks");
			}
			catch(...)
			{
				FatalError << "Wrong number of disks in DiskTopology.xml file" << ::Foam::exit(FatalError);
			}

			try
			{
				for(unsigned int i=0;i<list_disks_.size();i++)
				{
					Info << "Reading topology for Disk " << list_disks_[i] << endl;

					std::string disk_name = "opensmoke.Disk" + std::to_string(list_disks_[i]) + ".";

					// Number of points
					const int np = ptree.get<int>(disk_name + "NumberOfPoints");
					disk_topology_indices[i].resize(np);

					// Volume
					disk_topology_volume[i] = ptree.get<double>(disk_name + "Volume");
					std::cout << " * Cells: " << np << " Volume: " << disk_topology_volume[i] << std::endl;

					// Cell indices
					std::stringstream data;
					data.str( ptree.get< std::string >(disk_name + "Cells") );  
					for(unsigned int j=0;j<np;j++)
					{
						std::string dummy;
						data >> dummy; disk_topology_indices[i](j) = std::stoi(dummy);
						data >> dummy;
						Info << j << " " << dummy << endl; 
					}
				}
			}
			catch(...)
			{
				FatalError << "Wrong input data in DiskTopology.xml file" << ::Foam::exit(FatalError);
			}
		}

		// Read disks
		for(unsigned int i=0;i<list_disks_.size();i++)
		{
			Info << "Reading source terms for Disk " << list_disks_[i] << endl;

			std::string disk_name = "Disk." + std::to_string(list_disks_[i]) + ".source." + time_target_ + ".xml";


			boost::property_tree::ptree ptree;
    			boost::property_tree::read_xml( (path_folder_disks_ / disk_name).string(), ptree );

			// Source terms
			disk_source_terms[i].resize(thermodynamicsMap_.NumberOfSpecies());
			disk_source_terms[i].setZero();

			try
			{
				const int number_of_species = ptree.get<int>("opensmoke.number-species");
				const double slice_volume = ptree.get<double>("opensmoke.slice-volume");

				std::stringstream data;
				data.str( ptree.get< std::string >("opensmoke.source-terms") );  
				for(unsigned int j=0;j<number_of_species;j++)
				{
					std::string name_species;
					data >> name_species;
					const int j_species = thermodynamicsMap_.IndexOfSpecies(name_species)-1;
					if (j_species < 0)
						Info << "Warning: " << name_species << " is not available in the current kinetic mechanism" << endl;

					std::string dummy;				
					data >> dummy; 
					data >> dummy; 				
					data >> dummy; 

					disk_source_terms[i](j_species) = correction_coefficient_*std::stod(dummy)*slice_volume/disk_topology_volume[i];	// [kg/s/m3]
					if (exclude_negative_disk_source_terms_ == true)
						if (disk_source_terms[i](j_species) < 0.)
							disk_source_terms[i](j_species) = 0.;
				}
			}
			catch(...)
			{
				FatalError << "Wrong input data in " << disk_name << " file"  << ::Foam::exit(FatalError);
			}
		}

		// Creation source terms
		{
			Info<< "Source terms from disks creation... " << endl;

			const unsigned int NC = thermodynamicsMap_.NumberOfSpecies();
			sourceFromDisk_.resize(NC);

			for (int i=0;i<NC;i++)
			{
				sourceFromDisk_.set
				(
					i,
					new volScalarField
					(
						IOobject
						(
							"SDisk:" + thermodynamicsMap_.NamesOfSpecies()[i],
							mesh.time().timeName(),
							mesh,
							IOobject::NO_READ,
							IOobject::AUTO_WRITE
						),
						mesh,
						dimensionedScalar("SDisk", dimensionSet(1, -3, -1, 0, 0), 0.0)
					)
				);
			}

			for(unsigned int i=0;i<list_disks_.size();i++)
				for(unsigned int j=0;j<disk_topology_indices[i].size();j++)
				{
					const int celli = disk_topology_indices[i](j);

					for(int k=0;k<NC;k++)
						sourceFromDisk_[k].ref()[celli] = disk_source_terms[i](k);
				}

			for(int k=0;k<NC;k++)
				sourceFromDisk_[k].write();

			FatalError << "Source terms from disks correctly written!"  << ::Foam::exit(FatalError);
		}
	}

	if (diskPreprocess_ == false)
	{
		Info<< "Source terms from disks creation... " << endl;

		const unsigned int NC = thermodynamicsMap_.NumberOfSpecies();
		sourceFromDisk_.resize(NC);

		for (int i=0;i<NC;i++)
		{
			sourceFromDisk_.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						"SDisk:" + thermodynamicsMap_.NamesOfSpecies()[i],
						mesh.time().timeName(),
						mesh,
						IOobject::MUST_READ,
						IOobject::AUTO_WRITE
					),
					mesh
				)
			);
		}
	}
}


