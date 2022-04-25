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

if (export_clustering_csv_file == true)
{
	Info<< "Export Clustering CSV file..." << endl;

	// Recognize geometry
	const Vector<label> validDirs = (mesh.geometricD() + Vector<label>::one)/2;
     	Info<< "    Mesh has " << mesh.nGeometricD() << " geometric (non-empty/wedge) directions " << validDirs << endl;

	// Open file
	std::ofstream fCSV;
	fCSV.open( runTime.timeName() + "/status.clustering.csv", std::ios::out);
	fCSV.setf(std::ios::scientific);

	// Pointers to relevant fields
	const scalarField& TCells = T.internalField();
	const scalarField& pCells = p.internalField();
	const scalarField& csiCells = csi.internalField();
	const scalarField& tauCells = tau.internalField();

	// Cluster index
	volScalarField cIndex
	(
		IOobject
		(
			"cluster:index",
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		),
		mesh
	);

	// Cluster size
	volScalarField cSize
	(
		IOobject
		(
			"cluster:size",
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		),
		mesh
	);

	// Optional fields
	unsigned int npcs = 0;
	for (unsigned int i=1;i<=mixture.thermodynamicsMap().NumberOfSpecies();i++)
	{
		std::stringstream label; label<<i;
		IOobject header
        	(
        		"pc:" + label.str(),
        		mesh.time().timeName(),
        		mesh,
        		IOobject::NO_READ
    		);

		if (header.typeHeaderOk<volScalarField>(true))
			npcs++;
		else
			break;
	}
	
	// Principal Components
	PtrList<volScalarField> pcFields;
	pcFields.resize(npcs);
	for (unsigned int i=1;i<=npcs;i++)
	{
		std::stringstream label; label<<i;
		pcFields.set ( i-1, new volScalarField ( IOobject ( "pc:" + label.str(), runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE ), mesh ) );
	}

	// Header line
	fCSV << "t[s],";			// time (s)
	fCSV << "x[m],";			// x-coordinate (m)
	fCSV << "y[m],";			// y-coordinate (m)
	fCSV << "z[m],";			// z-coordinate (m)
	fCSV << "csi,";				// mixture fraction
	fCSV << "tau[s],";			// local residence time (s)
	fCSV << "T[K],";			// temperature (K)
	fCSV << "P[Pa],";			// pressure (Pa)
	fCSV << "Cindex,";			// cluster index
	fCSV << "Csize,";			// cluster size
	for (unsigned int i=1;i<=npcs;i++)
		fCSV << "PC" << i << ",";	// principal components
	fCSV << std::endl;

	// Loop over all the cells
	forAll(TCells, celli)
	{
		// Write on file
		fCSV << runTime.timeName() << ",";	// time (s)

		// Coordinates
		if (validDirs[0] != 0)	fCSV << mesh.C()[celli][0] << ",";	// x-coordinate (m)
		else			fCSV << 0 << ",";			// x-coordinate (m)
		if (validDirs[1] != 0)	fCSV << mesh.C()[celli][1] << ",";	// y-coordinate (m)
		else			fCSV << 0 << ",";			// y-coordinate (m)
		if (validDirs[2] != 0)	fCSV << mesh.C()[celli][2] << ",";	// z-coordinate (m)
		else			fCSV << 0 << ",";			// z-coordinate (m)

		fCSV << csiCells[celli] << ",";		// mixture fraction (-)
		fCSV << tauCells[celli] << ",";		// local residence time (s)
		fCSV << TCells[celli] << ",";		// temperature (K)
		fCSV << pCells[celli] << ",";		// pressure (Pa)
		fCSV << cIndex[celli] << ",";		// pressure (Pa)
		fCSV << cSize[celli] << ",";		// pressure (Pa)

		for (unsigned int i=0;i<npcs;i++)
			fCSV << pcFields[i][celli] << ",";	

		fCSV << std::endl;
	}

	fCSV.close();
}

