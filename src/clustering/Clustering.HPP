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
const Foam::word Foam::Clustering::dictName("Clustering");

// Constructor
Foam::Clustering::Clustering(const Foam::fvMesh& mesh)
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

	cluster_size_
	(
	    IOobject
	    (
		"cluster:size",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("cluster:size", dimensionSet(0, 0, 0, 0, 0), 0.0)
	),

	cluster_index_
	(
	    IOobject
	    (
		"cluster:index",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("cluster:index", dimensionSet(0, 0, 0, 0, 0), 0.0)
	)
{

	// Set key-fields (mandatory)
	{
		List<word>  listFields(this->lookup("fields"));
		SetKeyFields(listFields);
	}

	// Set accuracy coefficients (mandatory)
	epsilons_ = readList<double>(this->lookup("epsilon"));

	// Set accuracy coefficients type (mandatory)
	{
		List<word>  listEpsilonRel(this->lookup("epsilonTypes"));

		epsilons_rel_.resize(epsilons_.size());
		for(unsigned int i=0;i<epsilons_.size();i++)
		{
			if (listEpsilonRel[i] == "ABS" || listEpsilonRel[i] == "abs")
				epsilons_rel_[i] = false;
			else if (listEpsilonRel[i] == "REL" || listEpsilonRel[i] == "rel")
				epsilons_rel_[i] = true;
			else FatalError << "ClusteringModel: wrong epsilonTypes argument. Available options: ABS | REL | abs | rel" << ::Foam::exit(FatalError);
		}
	}

	// Set the algorithm to be used (mandatory)
	{
		word algorithm = lookup("algorithm");
		if (algorithm == "DMZ")
			algorithm_ = CLUSTERING_DMZ;
		else if (algorithm == "Binning")
			algorithm_ = CLUSTERING_BINNING;
		else FatalError << "ClusteringModel: wrong algorithm argument. Available options: DMZ | Binning" << ::Foam::exit(FatalError);
	}

	// Set the verbosity level to be used (optional, default: 1)
	verbosity_ = lookupOrDefault<label>("verbosity", 1);

	// Set the backward mapping policy (optional, default: correctedLinearInterpolation)
	{	
		word backward_mapping = lookupOrDefault<word>("speciesBackwardMapping", "correctedLinearInterpolation");
		if (backward_mapping == "correctedLinearInterpolation")
			species_backward_mapping_ = SPECIES_LINEAR_INTERPOLATION_WITH_CORRECTION;
		else if (backward_mapping == "linearInterpolation")
			species_backward_mapping_ = SPECIES_LINEAR_INTERPOLATION;
		else FatalError << "ClusteringModel: wrong speciesBackwardMapping argument. Available options: correctedLinearInterpolation | linearInterpolation" << ::Foam::exit(FatalError);
	}

	// Set the backward mapping policy for temperature (optional, default: linearInterpolation)
	{	
		word backward_mapping = lookupOrDefault<word>("temperatureBackwardMapping", "linearInterpolation");
		if (backward_mapping == "linearInterpolation")
			temperature_backward_mapping_ = TEMPERATURE_LINEAR_INTERPOLATION;
		else if (backward_mapping == "reconstructed")
			temperature_backward_mapping_ = TEMPERATURE_RECONSTRUCTED;
		else FatalError << "ClusteringModel: wrong temperatureBackwardMapping argument. Available options: linearInterpolation | reconstructed" << ::Foam::exit(FatalError);
	}
	
	// Analyze the list of provided input data
	{
		if (fields_.size() != epsilons_.size())
			FatalError << "ClusteringModel: the lists of key-fields and accuracy coefficients (epsilon) are not consistent" << ::Foam::exit(FatalError);
				
		for (unsigned int i=0;i<epsilons_.size();i++)
			if (epsilons_[i] <= 0.)
				FatalError << "ClusteringModel: the accuracy coefficients (epsilons) must be strictly positive" << ::Foam::exit(FatalError);
	}

	// Set the number of fields (i.e. features)
	D_ = fields_.size();
}

Foam::Clustering::~Clustering()
{}

void Foam::Clustering::Setup(const Foam::fvMesh& mesh)
{
	// Check for the existence of key-fields
	for (unsigned int i=0;i<D_;i++)
	{
		const volScalarField& ft = mesh.lookupObject<volScalarField>(fields_[i]);
	}	
}

void Foam::Clustering::Analysis(const Foam::fvMesh& mesh)
{
	if (algorithm_ == CLUSTERING_DMZ)
	{
		DMZAlgorithm dmz;
		if (Pstream::master())	dmz.SetVerbosityLevel(verbosity_);
		else			dmz.SetVerbosityLevel(0);
		dmz.SetNumberOfFeatures(D_);
		
		for(unsigned int i=0;i<D_;i++)
		{
			std::vector<double> psi(mesh.nCells());
			const volScalarField& field = mesh.lookupObject<volScalarField>(fields_[i]);
			forAll(field, celli)
				psi[celli] = field[celli];
			dmz.SetData(i, psi, epsilons_[i], epsilons_rel_[i]);
		}

		dmz.ClusterData();
		dmz.CalculateCentroids();
		dmz.CalculateDispersion();
		if (Pstream::master())	dmz.Print();

		dmz.GetData(nc_, gj_, lj_);
	}
	else if (algorithm_ == CLUSTERING_BINNING)
	{
		BinningAlgorithm binning;
		if (Pstream::master())	binning.SetVerbosityLevel(verbosity_);
		else			binning.SetVerbosityLevel(0);
		binning.SetNumberOfFeatures(D_);
		
		for(unsigned int i=0;i<D_;i++)
		{
			std::vector<double> psi(mesh.nCells());
			const volScalarField& field = mesh.lookupObject<volScalarField>(fields_[i]);
			forAll(field, celli)
				psi[celli] = field[celli];
			binning.SetData(i, psi, epsilons_[i], epsilons_rel_[i]);
		}

		binning.ClusterData();
		binning.CalculateCentroids();
		binning.CalculateDispersion();
		if (Pstream::master())	binning.Print();

		binning.GetData(nc_, gj_, lj_);
	}

	if (verbosity_ > 0)
	{
		// Set the cluster index
		{
			scalarField& clusterIndexCells = cluster_index_.ref();
			for (unsigned int i=0;i<gj_.size();i++)
				clusterIndexCells[i] = gj_[i];
		}

		// Set the size (i.e. the number of each cells in a cluster)
		{
			scalarField& clusterSizeCells = cluster_size_.ref();
			for (unsigned int i=0;i<gj_.size();i++)
				clusterSizeCells[i] = lj_[gj_[i]].size();
		}
	}		
}

const Foam::volScalarField& Foam::Clustering::cluster_index() const
{
	return cluster_index_;
}

Foam::volScalarField& Foam::Clustering::cluster_index()
{
	return cluster_index_;
}

const Foam::volScalarField& Foam::Clustering::cluster_size() const
{
	return cluster_size_;
}

Foam::volScalarField& Foam::Clustering::cluster_size()
{
	return cluster_size_;
}

