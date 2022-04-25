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

if (export_clustering_history == true)
{
	Info<< "Export Clustering History..." << endl;

	// Pointers to relevant fields
	const scalarField& TCells = T.internalField();

	//- Clustering of cells
	mixture.clusteringModel().Analysis(mesh);

	// Analysis of clusters
	const unsigned int nclusters = mixture.clusteringModel().nc();

	std::vector<unsigned int> ncells_per_cluster(nclusters);
	std::fill(ncells_per_cluster.begin(), ncells_per_cluster.end(), 0);
	
	forAll(TCells, celli)
	{
		ncells_per_cluster[static_cast<unsigned int>(mixture.clusteringModel().cluster_index()[celli])]++;
	}

	const unsigned int max_size = *std::max_element(ncells_per_cluster.begin(), ncells_per_cluster.end());
	const unsigned int min_size = *std::min_element(ncells_per_cluster.begin(), ncells_per_cluster.end());

        std::sort( ncells_per_cluster.begin(), ncells_per_cluster.end() );
        const unsigned int sum = std::accumulate( ncells_per_cluster.begin(), ncells_per_cluster.end(), 0 );
	const double mean_size = static_cast<double>(sum)/ncells_per_cluster.size();
        const unsigned int median = ( ncells_per_cluster[ncells_per_cluster.size()/2] + ncells_per_cluster[ (ncells_per_cluster.size()-1)/2 ] ) / 2.0 ;


	// Header line
	fClusteringHistory() << std::setw(20) << runTime.timeName();
	fClusteringHistory() << std::setw(20) << nclusters;
	fClusteringHistory() << std::setw(20) << min_size;
	fClusteringHistory() << std::setw(20) << max_size;
	fClusteringHistory() << std::setw(20) << mean_size;
	fClusteringHistory() << std::setw(20) << median;
	fClusteringHistory() << std::endl;


	// Additional info in case of PCA
	if (mixture.pca() == true)
	{
		Eigen::VectorXd explained;
		Eigen::MatrixXd cosines;

		mixture.pcaModel().Analysis(T,Y, explained, cosines);

		fClusteringPCAHistory() << std::setw(20) << runTime.timeName();
		fClusteringPCAHistory() << std::setw(20) << explained(0);
		fClusteringPCAHistory() << std::setw(20) << explained(1);
		fClusteringPCAHistory() << std::setw(20) << explained(2);
		fClusteringPCAHistory() << std::setw(20) << explained(3);

		for (unsigned int i=0;i<mixture.pcaModel().enabled_species_indices().size()+1;i++)
			fClusteringPCAHistory() << std::setw(20) << cosines(0,i);
		for (unsigned int i=0;i<mixture.pcaModel().enabled_species_indices().size()+1;i++)
			fClusteringPCAHistory() << std::setw(20) << cosines(1,i);

		fClusteringPCAHistory() << std::endl;
	}

}

