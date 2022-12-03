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
|   Copyright(C) 2020, 2021 Alberto Cuoci                                 |
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

// Include standard OpenFOAM classes
#include "fvCFD.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "pressureReference.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "fvcSmooth.H"
#include "interpolation.H"

// Include laminarSMOKE++ classes
#include "SparkModel.H"
#include "OpenSMOKEppReactingMixture.H"
#include "OpenSMOKEppReactingMixtureForRadiation.H"

// Lookup tables
#include "utilities/lookuptables/LookUpTable_ZY.h"

// Customized radiation model
#include "radiationModelOpenSMOKE++.H"

// Function for writing header lines in output files
std::string headerLabel(const std::string text, unsigned int& count)
{
	std::stringstream number; number << count;
	count++;
	return ( text + "(" + number.str() + ")" );
}

// Main code
int main(int argc, char *argv[])
{
	timeSelector::addOptions();

    	#include "addRegionOption.H"
    	#include "addDictOption.H"
    	#include "setRootCase.H"
    	#include "createTime.H"

    	instantList timeDirs = timeSelector::select0(runTime, args);

    	#include "createNamedMesh.H"	
	#include "createControl.H"

	// Read basic fields
	#include "createFields.H"
	#include "createFieldRefs.H"

	Info<< "Reading PostProcessing dictionary\n" << endl;
	IOdictionary postProcessingDictionary
	(
		IOobject
		(
			"PostProcessing",
			U.time().constant(),
			U.db(),
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	);

	// Default options
	Switch export_csv_file = false;
	Switch export_clustering_csv_file = false;
	Switch export_clustering_xml_file = false;
	Switch export_clustering_history = false;
	Switch export_opensmokeppxml_file = false;
	Switch calculate_scalar_dissipation_rate = false;
	Switch calculate_flame_curvature = false;
	Switch calculate_mixture_fraction_and_progress_variable = false;
	Switch export_thermo_and_transport_properties = false;
	Switch reconstruct_polimi_soot = false;
	Switch export_outlet = false;
	Switch export_disksxml = false;
	double tabulatedChemistryTemperatureThreshold = 300.;

	// Read options from PostProcessing dictionary
	export_csv_file = Switch(postProcessingDictionary.lookupOrDefault(word("exportCSV"), word("off")));
	export_clustering_csv_file = Switch(postProcessingDictionary.lookupOrDefault(word("exportClusteringCSV"), word("off")));
	export_clustering_xml_file = Switch(postProcessingDictionary.lookupOrDefault(word("exportClusteringXML"), word("off")));
	export_clustering_history = Switch(postProcessingDictionary.lookupOrDefault(word("exportClusteringHistory"), word("off")));
	export_opensmokeppxml_file = Switch(postProcessingDictionary.lookupOrDefault(word("exportOpenSMOKE++XML"), word("off")));
	export_thermo_and_transport_properties = Switch(postProcessingDictionary.lookupOrDefault(word("exportThermoAndTransportProperties"), word("off")));
	calculate_mixture_fraction_and_progress_variable = Switch(postProcessingDictionary.lookupOrDefault(word("tabulatedChemistry"), word("off")));
	calculate_scalar_dissipation_rate = Switch(postProcessingDictionary.lookupOrDefault(word("scalarDissipationRate"), word("off")));
	calculate_flame_curvature = Switch(postProcessingDictionary.lookupOrDefault(word("flameCurvature"), word("off")));
	reconstruct_polimi_soot = Switch(postProcessingDictionary.lookupOrDefault(word("reconstructPolimiSoot"), word("off")));
	export_outlet = Switch(postProcessingDictionary.lookupOrDefault(word("exportOutlet"), word("off")));
	export_disksxml = Switch(postProcessingDictionary.lookupOrDefault(word("exportDisksXML"), word("off")));

	// Forcing calculation of scalar dissipation rate and curvature
	if (export_csv_file == true || export_opensmokeppxml_file == true)
	{
		calculate_scalar_dissipation_rate = true;
		calculate_flame_curvature = true;
	}

	// Open global files
	autoPtr<std::ofstream> fClusteringHistory;
	autoPtr<std::ofstream> fClusteringPCAHistory;
	if (export_clustering_history == true)
	{
		fClusteringHistory.reset(new std::ofstream("Clustering.History", std::ios::out));
		fClusteringHistory().setf(std::ios::scientific);

		fClusteringHistory() << std::setw(20) << "time[s](1)";
		fClusteringHistory() << std::setw(20) << "nclusters(2)";
		fClusteringHistory() << std::setw(20) << "minSize(3)";
		fClusteringHistory() << std::setw(20) << "maxSize(4)";
		fClusteringHistory() << std::setw(20) << "meanSize(5)";
		fClusteringHistory() << std::setw(20) << "median(6)";
		fClusteringHistory() << std::endl;

		if (mixture.pca() == true)
		{
			fClusteringPCAHistory.reset(new std::ofstream("Clustering.PCA.History", std::ios::out));
			fClusteringPCAHistory().setf(std::ios::scientific);

			fClusteringPCAHistory() << std::setw(20) << "time[s](1)";
			fClusteringPCAHistory() << std::setw(20) << "exp-var[1](2)";	
			fClusteringPCAHistory() << std::setw(20) << "exp-var[2](3)";	
			fClusteringPCAHistory() << std::setw(20) << "exp-var[3](4)";	
			fClusteringPCAHistory() << std::setw(20) << "exp-var[4](5)";	
			
                        fClusteringPCAHistory() << std::setw(20) << "sum-var[1](6)";
                        fClusteringPCAHistory() << std::setw(20) << "sum-var[2](7)";
                        fClusteringPCAHistory() << std::setw(20) << "sum-var[3](8)";
                        fClusteringPCAHistory() << std::setw(20) << "sum-var[4](9)";

			fClusteringPCAHistory() << std::setw(20) << "exp-var[99%](10)";
			fClusteringPCAHistory() << std::setw(20) << "exp-var[98%](11)";
			fClusteringPCAHistory() << std::setw(20) << "exp-var[95%](12)";
                        fClusteringPCAHistory() << std::setw(20) << "exp-var[90%](13)";
                        fClusteringPCAHistory() << std::setw(20) << "exp-var[80%](14)";
                        fClusteringPCAHistory() << std::setw(20) << "exp-var[70%](15)";
                        fClusteringPCAHistory() << std::setw(20) << "exp-var[60%](16)";
                        fClusteringPCAHistory() << std::setw(20) << "exp-var[50%](17)";

			unsigned int count=18;
			fClusteringPCAHistory() << std::setw(20) << headerLabel("cos(W1,T)", count); 
			for (unsigned int i=0;i<mixture.pcaModel().enabled_species_indices().size();i++)
				fClusteringPCAHistory() << std::setw(20) << headerLabel("cos(W1," + mixture.thermodynamicsMap().NamesOfSpecies()[mixture.pcaModel().enabled_species_indices()[i]] + ")", count);
			fClusteringPCAHistory() << std::setw(20) << headerLabel("cos(W2,T)", count); 
			for (unsigned int i=0;i<mixture.pcaModel().enabled_species_indices().size();i++)
				fClusteringPCAHistory() << std::setw(20) << headerLabel("cos(W2," + mixture.thermodynamicsMap().NamesOfSpecies()[mixture.pcaModel().enabled_species_indices()[i]] + ")", count);  
			fClusteringPCAHistory() << std::endl;	
		}
	}
		
	// Options for lookup table tabulation
	OpenSMOKE::LookUpTable_ZY* lookup_table_zy;
	if (calculate_mixture_fraction_and_progress_variable == true)
	{
		Info << "Reading lookup table..." << endl;

		// Folder containing the tables
		const Foam::string table_folder = postProcessingDictionary.lookup("tabulatedChemistryFolder");
		const boost::filesystem::path table_folder_path = table_folder;
		const boost::filesystem::path path_to_xml = table_folder_path / "lookuptable.main.xml";

		// List of fields available in the tables
		std::vector<std::string> list_fields;
		List<word>  listFields(postProcessingDictionary.lookup("tabulatedChemistryFields"));
		for (unsigned int i=0;i<listFields.size();i++)
			list_fields.push_back(listFields[i]);	

		// Threshold temperature for scatter plots
		tabulatedChemistryTemperatureThreshold = postProcessingDictionary.lookupOrDefault<scalar>("tabulatedChemistryTemperatureThreshold", 300.);

		// Reading tables
		lookup_table_zy = new OpenSMOKE::LookUpTable_ZY(mixture.thermodynamicsMap());
		lookup_table_zy->Setup(path_to_xml, list_fields, true);
	}

	// Loop over times
	forAll(timeDirs, timeI)
    	{
       		runTime.setTime(timeDirs[timeI], timeI);
        	Info<< "Time = " << runTime.timeName() << endl;

        	// Handle geometry/topology changes
        	polyMesh::readUpdateState state = mesh.readUpdate();

		// Update main fields
		#include "updateFields.H"

		// Update thermodynamic and transport properties
		mixture.update_properties();

		// Update transport terms (i.e. fluxes)
		if (calculate_scalar_dissipation_rate == true)
			mixture.update_transport_terms(mesh, rho);

		#include "reconstructStoichiometricScalarDissipationRateAndStrainRate.H"

		#include "reconstructFlameCurvature.H"

		#include "reconstructMixtureFractionAndRPV.H"

		#include "exportCSV.H"

		#include "exportClusteringCSV.H"

		#include "exportClusteringXML.H"

		#include "exportClusteringHistory.H"

		#include "exportOpenSMOKE++XML.H"

		#include "exportThermoAndTransportProperties.H"

		#include "reconstructPolimiSoot.H"

		#include "exportOutlet.H"

		#include "exportDisksXML.H"
	}

	// Close global files
	if (export_clustering_history == true)
	{
		fClusteringHistory().close();
		if (mixture.pca() == true)
			fClusteringPCAHistory().close();
	}

	return 0;
}

