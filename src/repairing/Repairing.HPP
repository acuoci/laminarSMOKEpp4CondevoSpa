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
|   Copyright(C) 2022 Alberto Cuoci                                       |
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
const Foam::word Foam::Repairing::dictName("Repairing");

// Constructor
Foam::Repairing::Repairing(const Foam::fvMesh& mesh)
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
	)
{

	// Default values
	minYValue_ = 0.;
	maxYValue_ = 1.;
	patchThreshold_ = 0.98;

	// Set key-fields (mandatory)
	{
		List<word>  listMinFields(this->lookup("minFields"));
		List<word>  listMaxFields(this->lookup("maxFields"));
		minFields_ = listMinFields;
		maxFields_ = listMaxFields;

		minValues_ = readList<double>(this->lookup("minValues"));
		maxValues_ = readList<double>(this->lookup("maxValues"));
	}

	// Define patching values
	{	
		
		List<word>  listPatchNames(this->lookup("patchNames"));
		patchNames_ = listPatchNames;
		patchValues_ = readList<double>(this->lookup("patchValues"));
		iPatchNames_.resize(patchNames_.size());
	}

	// Check patching values
	double sum = 0.;
	for (unsigned int i=0; i<patchValues_.size(); i++)
		sum += patchValues_[i];
	if ( std::fabs(sum-1.)>1.e-9)
		FatalError << "RepairingModel: the sum of patch values must be equal to 1" << endl;
}

Foam::Repairing::~Repairing()
{}

void Foam::Repairing::Setup(const std::vector<std::string>& names)
{	
	for (unsigned int j=0; j<patchNames_.size(); j++)
	{
		iPatchNames_[j] = -1;	
		for (unsigned int k=0; k<names.size(); k++)
			if (patchNames_[j] == names[k])
			{
				iPatchNames_[j] = k;
				break;
			}
		if ( iPatchNames_[j] == -1)
			FatalError << "RepairingModel: the following patchNames is not available in the current kinetics: " << patchNames_[j] << endl;
	}
}

void Foam::Repairing::RepairMassFractions(const Foam::fvMesh& mesh, PtrList<volScalarField>& Y)
{
	Info << " * Repairing mass fractions..." << endl;

	const volScalarField& T = mesh.lookupObject<volScalarField>("T");
	const scalarField& TCells = T.internalField();

	// Internal cells
	forAll(TCells, celli)
	{
		for (unsigned int j=0; j<Y.size(); j++)
		{
			Y[j].ref()[celli] = std::max(Y[j].internalField()[celli], minYValue_);
			Y[j].ref()[celli] = std::min(Y[j].internalField()[celli], maxYValue_);
		}
	}

	// Boundaries
	{
		forAll(T.boundaryField(), patchi)
		{
			const fvPatchScalarField& pT = T.boundaryField()[patchi];

			forAll(pT, facei)
			{
				for (unsigned int j=0; j<Y.size(); j++)
				{
					Y[j].boundaryFieldRef()[patchi][facei] = std::max(Y[j].boundaryField()[patchi][facei], minYValue_);
					Y[j].boundaryFieldRef()[patchi][facei]= std::min(Y[j].boundaryField()[patchi][facei], maxYValue_);
				}
			}
		}
	}
}

void Foam::Repairing::RepairVelocity(const Foam::fvMesh& mesh)
{
	Info << " * Repairing velocity..." << endl;

	const volScalarField& T = mesh.lookupObject<volScalarField>("T");
	const scalarField& TCells = T.internalField();

	volVectorField& U = mesh.lookupObjectRef<volVectorField>("U");
	scalarField& uCells = U.component(0).ref();
	scalarField& vCells = U.component(1).ref();
	scalarField& wCells = U.component(2).ref();

	// Internal cells
	forAll(TCells, celli)
	{
		uCells[celli] = std::max(uCells[celli], -150.);
		uCells[celli] = std::min(uCells[celli],  150.);
		vCells[celli] = std::max(vCells[celli], -150.);
		vCells[celli] = std::min(vCells[celli],  150.);
		wCells[celli] = std::max(wCells[celli], -150.);
		wCells[celli] = std::min(wCells[celli],  150.);
	}
/*
	// Boundaries
	{
		forAll(T.boundaryField(), patchi)
		{
			const fvPatchScalarField& pT = T.boundaryField()[patchi];

			forAll(pT, facei)
			{
				for (unsigned int j=0; j<Y.size(); j++)
				{
					Y[j].boundaryFieldRef()[patchi][facei] = std::max(Y[j].boundaryField()[patchi][facei], minYValue_);
					Y[j].boundaryFieldRef()[patchi][facei]= std::min(Y[j].boundaryField()[patchi][facei], maxYValue_);
				}
			}
		}
	}
*/
}

void Foam::Repairing::RepairMinMaxFields(const Foam::fvMesh& mesh)
{	
	Info << " * Repairing min fields..." << endl;

	// Repair min values
	for (unsigned int j=0; j<minFields_.size(); j++)
	{
		Info << "   - " << minFields_[j] << endl;

		volScalarField& field = mesh.lookupObjectRef<volScalarField>(minFields_[j]);
		scalarField& fieldCells = field.ref();

		forAll(fieldCells, celli)
		{
			fieldCells[celli] = std::max(fieldCells[celli], minValues_[j]);
		}
	}

	// Boundaries
	{
		const volScalarField& T = mesh.lookupObject<volScalarField>("T");
		const scalarField& TCells = T.internalField();

		forAll(T.boundaryField(), patchi)
		{
			const fvPatchScalarField& pT = T.boundaryField()[patchi];

			for (unsigned int j=0; j<minFields_.size(); j++)
			{
				volScalarField& field = mesh.lookupObjectRef<volScalarField>(minFields_[j]);
				fvPatchScalarField& pfield = field.boundaryFieldRef()[patchi];

				forAll(pT, facei)
				{
					pfield[facei] = std::max(pfield[facei], minValues_[j]);
				}
			}
		}
	}

	Info << " * Repairing max fields..." << endl;

	// Repair max values
	for (unsigned int j=0; j<maxFields_.size(); j++)
	{
		Info << "   - " << maxFields_[j] << endl;

		volScalarField& field = mesh.lookupObjectRef<volScalarField>(maxFields_[j]);
		scalarField& fieldCells = field.ref();

		forAll(fieldCells, celli)
		{
			fieldCells[celli] = std::min(fieldCells[celli], maxValues_[j]);
		}
	}

	// Boundaries
	{
		const volScalarField& T = mesh.lookupObject<volScalarField>("T");
		const scalarField& TCells = T.internalField();

		forAll(T.boundaryField(), patchi)
		{
			const fvPatchScalarField& pT = T.boundaryField()[patchi];

			for (unsigned int j=0; j<maxFields_.size(); j++)
			{
				volScalarField& field = mesh.lookupObjectRef<volScalarField>(maxFields_[j]);
				fvPatchScalarField& pfield = field.boundaryFieldRef()[patchi];

				forAll(pT, facei)
				{
					pfield[facei] = std::min(pfield[facei], maxValues_[j]);
				}
			}
		}
	}
}

void Foam::Repairing::RepairPatch(const Foam::fvMesh& mesh, PtrList<volScalarField>& Y)
{
	Info << " * Patching mass fractions..." << endl;

	const volScalarField& T = mesh.lookupObject<volScalarField>("T");
	const scalarField& TCells = T.internalField();
	const unsigned int ns = Y.size();

	Eigen::VectorXd mass_fractions(ns);

	int countCells = 0;
	forAll(TCells, celli)
	{
		for(int i=0;i<ns;++i)
			mass_fractions(i) = Y[i].internalField()[celli];
		
		const double sum = mass_fractions.sum();
		if ( sum <= patchThreshold_)
		{
			countCells++;

			for(int i=0;i<ns;++i)
				Y[i].ref()[celli] = 0.;

			for (unsigned int i=0; i<patchNames_.size(); i++)
				Y[iPatchNames_[i]].ref()[celli] = patchValues_[i];
		}
		else
		{
			for(int i=0;i<ns;++i)
				Y[i].ref()[celli] = Y[i].internalField()[celli] / sum;
		}
	}

	Info << "   - " << "number of patched cells: " << countCells << endl;

	// Boundaries
	int countFaces = 0;
	{
		const volScalarField& T = mesh.lookupObject<volScalarField>("T");
		const scalarField& TCells = T.internalField();

		forAll(T.boundaryField(), patchi)
		{
			const fvPatchScalarField& pT = T.boundaryField()[patchi];

			forAll(pT, facei)
			{
				for(int i=0;i<ns;++i)
					mass_fractions(i) = Y[i].boundaryField()[patchi][facei];

				const double sum = mass_fractions.sum();
				if ( sum <= patchThreshold_)
				{
					countFaces++;

					for(int i=0;i<ns;++i)
						Y[i].boundaryFieldRef()[patchi][facei] = 0.;

					for (unsigned int i=0; i<patchNames_.size(); i++)
						Y[iPatchNames_[i]].boundaryFieldRef()[patchi][facei] = patchValues_[i];
				}
				else
				{
					for(int i=0;i<ns;++i)
						Y[i].boundaryFieldRef()[patchi][facei] = Y[i].boundaryField()[patchi][facei] / sum;
				}
			}
		}
	}

	Info << "   - " << "number of patched faces: " << countFaces << endl;
}



