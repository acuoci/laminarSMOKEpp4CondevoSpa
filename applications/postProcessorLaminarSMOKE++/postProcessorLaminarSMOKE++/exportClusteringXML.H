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

if (export_clustering_xml_file == true)
{
	Info<< "Export Clustering XML file..." << endl;

	// Recognize geometry
	const Vector<label> validDirs = (mesh.geometricD() + Vector<label>::one)/2;
     	Info<< "    Mesh has " << mesh.nGeometricD() << " geometric (non-empty/wedge) directions " << validDirs << endl;

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

	// Open file
	std::ofstream fXML;
	fXML.open( runTime.timeName() + "/status.clustering.xml", std::ios::out);
	fXML.setf(std::ios::scientific);

	// Write file
	fXML << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
	fXML << "<opensmoke version=\"0.1a\">" << std::endl;
	
	fXML << "<classes>"<< static_cast<unsigned int>(max(cIndex).value()+1) << "</classes>" << std::endl;
	fXML << "<items>" << mesh.cells().size() << "</items>" << std::endl;
	fXML << "<original-components>" << mixture.thermodynamicsMap().NumberOfSpecies()+1 << "</original-components>" << std::endl;
	fXML << "<filtered-components>" << mixture.thermodynamicsMap().NumberOfSpecies()+1 << "</filtered-components>" << std::endl;
	fXML << "<removed-components>" << 0 << "</removed-components>" << std::endl;
	fXML << "<number-retained-species>" << mixture.thermodynamicsMap().NumberOfSpecies() << "</number-retained-species>" << std::endl;
	
	// Pointers to relevant fields
	const scalarField& TCells = T.internalField();
	const scalarField& csiCells = csi.internalField();
	const scalarField& cIndexCells = cIndex.internalField();

	fXML << "<data-original>" << std::endl;
	forAll(TCells, celli)
	{
		// Write on file
		fXML << static_cast<unsigned int>(cIndexCells[celli]+1) << " ";	// cluster index (from 1)
		fXML << csiCells[celli] << " ";					// mixture fraction (-)
		fXML << TCells[celli] << " ";					// temperature (K)
		
		// Sum of mass fractions
		double sum = 0.;
		for (unsigned int j = 0; j < mixture.thermodynamicsMap().NumberOfSpecies(); j++)
			sum += Y[j][celli];

		// Write mass fractions
		const double threshold = 1.e-64;
		for (unsigned int j = 0; j < mixture.thermodynamicsMap().NumberOfSpecies(); j++)
		{
			if (Y[j][celli] < 1.e-64)	fXML << 0. << " ";
			else				fXML << Y[j][celli]/sum << " ";
		}

		fXML << std::endl;
	}

	fXML << "</data-original>" << std::endl;
	fXML << "</opensmoke>" << std::endl;
	fXML.close();
}

