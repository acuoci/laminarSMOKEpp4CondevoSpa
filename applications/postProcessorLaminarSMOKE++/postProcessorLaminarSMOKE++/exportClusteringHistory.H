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
		Eigen::VectorXd explainedSum;
		Eigen::MatrixXd cosines;

		mixture.pcaModel().Analysis(T,Y, explained, cosines);

		explainedSum.resize(explained.size());
		explainedSum(0) = explained(0);
		for (unsigned int i=1;i<explainedSum.size();i++)
			explainedSum(i) = explainedSum(i-1)+explained(i);

		unsigned int index50 = 0;
		unsigned int index60 = 0;
		unsigned int index70 = 0;
		unsigned int index80 = 0;
		unsigned int index90 = 0;
		unsigned int index95 = 0;
		unsigned int index98 = 0;
		unsigned int index99 = 0;
		for (unsigned int i=0;i<explainedSum.size();i++)
		{
			if (index50 == 0 && explainedSum(i)>=0.50)	index50 = i+1;
                        if (index60 == 0 && explainedSum(i)>=0.60)      index60 = i+1;
                        if (index70 == 0 && explainedSum(i)>=0.70)      index70 = i+1;
                        if (index80 == 0 && explainedSum(i)>=0.80)      index80 = i+1;
                        if (index90 == 0 && explainedSum(i)>=0.90)      index90 = i+1;
                        if (index95 == 0 && explainedSum(i)>=0.95)      index95 = i+1;
                        if (index98 == 0 && explainedSum(i)>=0.98)      index98 = i+1;
                        if (index99 == 0 && explainedSum(i)>=0.99)      index99 = i+1;
		}

		fClusteringPCAHistory() << std::setw(20) << runTime.timeName();
		fClusteringPCAHistory() << std::setw(20) << explained(0);
		fClusteringPCAHistory() << std::setw(20) << explained(1);
		fClusteringPCAHistory() << std::setw(20) << explained(2);
		fClusteringPCAHistory() << std::setw(20) << explained(3);

                fClusteringPCAHistory() << std::setw(20) << explainedSum(0);
                fClusteringPCAHistory() << std::setw(20) << explainedSum(1);
                fClusteringPCAHistory() << std::setw(20) << explainedSum(2);
                fClusteringPCAHistory() << std::setw(20) << explainedSum(3);

                fClusteringPCAHistory() << std::setw(20) << index99;
               	fClusteringPCAHistory() << std::setw(20) << index98;
               	fClusteringPCAHistory() << std::setw(20) << index95;
               	fClusteringPCAHistory() << std::setw(20) << index90;
               	fClusteringPCAHistory() << std::setw(20) << index80;
               	fClusteringPCAHistory() << std::setw(20) << index70;
               	fClusteringPCAHistory() << std::setw(20) << index60;
               	fClusteringPCAHistory() << std::setw(20) << index50;

		for (unsigned int i=0;i<mixture.pcaModel().enabled_species_indices().size()+1;i++)
			fClusteringPCAHistory() << std::setw(20) << cosines(0,i);
		for (unsigned int i=0;i<mixture.pcaModel().enabled_species_indices().size()+1;i++)
			fClusteringPCAHistory() << std::setw(20) << cosines(1,i);

		fClusteringPCAHistory() << std::endl;
	}

}

