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

if (reconstruct_polimi_soot == true)
{
	OpenSMOKE::PolimiSoot_Analyzer sootAnalyzer(&mixture.thermodynamicsMap());
	sootAnalyzer.Initialize();
	sootAnalyzer.SetLabel("BIN");
	sootAnalyzer.SetFractalDiameter(1.8);
	sootAnalyzer.SetMinimumSectionSphericalParticles("BIN5");
	sootAnalyzer.SetMinimumSectionAggregates("BIN13");
	sootAnalyzer.SetDensity(10, 20, 1500., 1500.);
	sootAnalyzer.SetPhysicalDiffusionBinToStart(5);
	sootAnalyzer.SetPhysicalDiffusionBinToCut(10);
	sootAnalyzer.SetThermophoreticEffectMinimumBin(5);
	sootAnalyzer.Setup();

	const unsigned int ns = mixture.thermodynamicsMap().NumberOfSpecies();

	volScalarField sootPrecursorsY
	(
	    IOobject
	    (
		"sootPrecursorsY",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("dimensionless", dimensionSet(0, 0, 0, 0, 0), 0.)
	);

	volScalarField sootSphericalY
	(
	    IOobject
	    (
		"sootSphericalY",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("dimensionless", dimensionSet(0, 0, 0, 0, 0), 0.)
	);

 	volScalarField sootAggregatesY
	(
	    IOobject
	    (
		"sootAggregatesY",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("dimensionless", dimensionSet(0, 0, 0, 0, 0), 0.)
	);

        volScalarField sootFv
	(
	    IOobject
	    (
		"sootFv",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
    	    dimensionedScalar("dimensionless", dimensionSet(0, 0, 0, 0, 0), 0.)
	);

	volScalarField sootSphericalN
	(
	    IOobject
	    (
		"sootSphericalN",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("dimensionless", dimensionSet(0, -3, 0, 0, 0), 0.)
	);

	volScalarField sootAggregatesN
	(
	    IOobject
	    (
		"sootAggregatesN",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("dimensionless", dimensionSet(0, -3, 0, 0, 0), 0.)
	);

	volScalarField sootSphericalD32
	(
	    IOobject
	    (
		"sootSphericalD32",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("dimensionless", dimensionSet(0, 1, 0, 0, 0), 0.)
	);

	volScalarField sootAggregatesD32
	(
	    IOobject
	    (
		"sootAggregatesD32",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("dimensionless", dimensionSet(0, 1, 0, 0, 0), 0.)
	);

	volScalarField sootSphericalD43
	(
	    IOobject
	    (
		"sootSphericalD43",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("dimensionless", dimensionSet(0, 1, 0, 0, 0), 0.)
	);

	volScalarField sootAggregatesD43
	(
	    IOobject
	    (
		"sootAggregatesD43",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("dimensionless", dimensionSet(0, 1, 0, 0, 0), 0.)
	);

	volScalarField sootSphericalHC
	(
	    IOobject
	    (
		"sootSphericalHC",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("dimensionless", dimensionSet(0, 0, 0, 0, 0), 0.)
	);

	volScalarField sootAggregatesHC
	(
	    IOobject
	    (
		"sootAggregatesHC",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("dimensionless", dimensionSet(0, 0, 0, 0, 0), 0.)
	);


	Info << "Reconstructing Polimi Soot ..." << endl;
	{
		double mw;
		Eigen::VectorXd y(mixture.thermodynamicsMap().NumberOfSpecies());
		Eigen::VectorXd x(mixture.thermodynamicsMap().NumberOfSpecies());

		const scalarField& TCells = T.internalField();
		const scalarField& pCells = p.internalField();
		scalarField& sootPrecursorsYCells = sootPrecursorsY.ref(); 
		scalarField& sootSphericalYCells = sootSphericalY.ref();
		scalarField& sootAggregatesYCells = sootAggregatesY.ref();
		scalarField& sootFvCells = sootFv.ref();
		scalarField& sootSphericalNCells = sootSphericalN.ref();
		scalarField& sootSphericalD32Cells = sootSphericalD32.ref();
		scalarField& sootSphericalD43Cells = sootSphericalD43.ref();
		scalarField& sootAggregatesNCells = sootAggregatesN.ref();
		scalarField& sootAggregatesD32Cells = sootAggregatesD32.ref();
		scalarField& sootAggregatesD43Cells = sootAggregatesD43.ref();
		scalarField& sootSphericalHCCells = sootSphericalHC.ref();
		scalarField& sootAggregatesHCCells = sootAggregatesHC.ref();
		
		  
		// Internal cells
		forAll(TCells, celli)
		{	
			// Mass fractions
			for (unsigned int i=0;i<ns;i++)
				y(i) = Y[i].internalField()[celli];
	
			// Molecular weight and mole fractions
			mixture.thermodynamicsMap().MoleFractions_From_MassFractions(x.data(), mw, y.data());

			// Density
			const double rhoGas = pCells[celli]*mw/PhysicalConstants::R_J_kmol/TCells[celli];

			// Soot analysis
			sootAnalyzer.Analysis(TCells[celli], pCells[celli], rhoGas, y, x);
			

			// Transfer
			sootPrecursorsYCells[celli] = sootAnalyzer.omega_small();
			sootSphericalYCells[celli] = sootAnalyzer.omega_large_spherical();
			sootAggregatesYCells[celli] = sootAnalyzer.omega_large_aggregates();
			sootFvCells[celli] = sootAnalyzer.fv_large();
			sootSphericalNCells[celli] = sootAnalyzer.N_large_spherical();
			sootSphericalD32Cells[celli] = sootAnalyzer.SootDiameters(3., 2., sootAnalyzer.bin_indices_large_spherical(), 1.e-11);
			sootSphericalD43Cells[celli] = sootAnalyzer.SootDiameters(4., 3., sootAnalyzer.bin_indices_large_spherical(), 1.e-11);
			sootAggregatesNCells[celli] = sootAnalyzer.N_large_aggregates();
			sootAggregatesD32Cells[celli] = sootAnalyzer.SootDiameters(3., 2., sootAnalyzer.bin_indices_large_aggregates(), 1.e-11);
			sootAggregatesD43Cells[celli] = sootAnalyzer.SootDiameters(4., 3., sootAnalyzer.bin_indices_large_aggregates(), 1.e-11);
			sootSphericalHCCells[celli] = sootAnalyzer.h_over_c_large_spherical();
			sootAggregatesHCCells[celli] = sootAnalyzer.h_over_c_large_aggregates();
		}

		// Boundaries
		const volScalarField::Boundary& Tf = T.boundaryField();
		const volScalarField::Boundary& pf = p.boundaryField();
		volScalarField::Boundary& sootPrecursorsYf = sootPrecursorsY.boundaryFieldRef();
		volScalarField::Boundary& sootSphericalYf = sootSphericalY.boundaryFieldRef();
		volScalarField::Boundary& sootAggregatesYf = sootAggregatesY.boundaryFieldRef();
		volScalarField::Boundary& sootFvf = sootFv.boundaryFieldRef();
		volScalarField::Boundary& sootSphericalNf = sootSphericalN.boundaryFieldRef();
		volScalarField::Boundary& sootSphericalD32f = sootSphericalD32.boundaryFieldRef();
		volScalarField::Boundary& sootSphericalD43f = sootSphericalD43.boundaryFieldRef();
		volScalarField::Boundary& sootAggregatesNf = sootAggregatesN.boundaryFieldRef();
		volScalarField::Boundary& sootAggregatesD32f = sootAggregatesD32.boundaryFieldRef();
		volScalarField::Boundary& sootAggregatesD43f = sootAggregatesD43.boundaryFieldRef();
		volScalarField::Boundary& sootSphericalHCf = sootSphericalHC.boundaryFieldRef();
		volScalarField::Boundary& sootAggregatesHCf = sootAggregatesHC.boundaryFieldRef();

		forAll(T.boundaryField(), patchi)
		{
			const fvPatchScalarField& pT = Tf[patchi];
			const fvPatchScalarField& pp = pf[patchi];
			fvPatchScalarField& psootPrecursorsY = sootPrecursorsYf[patchi];
			fvPatchScalarField& psootSphericalY = sootSphericalYf[patchi];
			fvPatchScalarField& psootAggregatesY = sootAggregatesYf[patchi];
			fvPatchScalarField& psootFv = sootFvf[patchi];
			fvPatchScalarField& psootSphericalN = sootSphericalNf[patchi];
			fvPatchScalarField& psootSphericalD32 = sootSphericalD32f[patchi];
			fvPatchScalarField& psootSphericalD43 = sootSphericalD43f[patchi];
			fvPatchScalarField& psootAggregatesN = sootAggregatesNf[patchi];
			fvPatchScalarField& psootAggregatesD32 = sootAggregatesD32f[patchi];
			fvPatchScalarField& psootAggregatesD43 = sootAggregatesD43f[patchi];
			fvPatchScalarField& psootSphericalHC = sootSphericalHCf[patchi];
			fvPatchScalarField& psootAggregatesHC = sootAggregatesHCf[patchi];
	
			forAll(psootFv, facei)
			{
				// Mass fractions
				for (unsigned int i=0;i<ns;i++)
					y(i) = Y[i].boundaryField()[patchi][facei];

				// Molecular weight and mole fractions
				mixture.thermodynamicsMap().MoleFractions_From_MassFractions(x.data(), mw, y.data());

				// Density
				const double rhoGas = pp[facei]*mw/PhysicalConstants::R_J_kmol/pT[facei];

				// Soot analysis
				sootAnalyzer.Analysis(pT[facei], pp[facei], rhoGas, y, x);

				// Transfer
				psootPrecursorsY[facei] = sootAnalyzer.omega_small();
				psootSphericalY[facei] = sootAnalyzer.omega_large_spherical();
				psootAggregatesY[facei] = sootAnalyzer.omega_large_aggregates();
				psootFv[facei] = sootAnalyzer.fv_large();
				psootSphericalN[facei] = sootAnalyzer.N_large_spherical();
				psootSphericalD32[facei] = sootAnalyzer.SootDiameters(3., 2., sootAnalyzer.bin_indices_large_spherical(), 1.e-11);
				psootSphericalD43[facei] = sootAnalyzer.SootDiameters(4., 3., sootAnalyzer.bin_indices_large_spherical(), 1.e-11);
				psootAggregatesN[facei] = sootAnalyzer.N_large_aggregates();
				psootAggregatesD32[facei] = sootAnalyzer.SootDiameters(3., 2., sootAnalyzer.bin_indices_large_aggregates(), 1.e-11);
				psootAggregatesD43[facei] = sootAnalyzer.SootDiameters(4., 3., sootAnalyzer.bin_indices_large_aggregates(), 1.e-11);
				psootSphericalHC[facei] = sootAnalyzer.h_over_c_large_spherical();
				psootAggregatesHC[facei] = sootAnalyzer.h_over_c_large_aggregates();
			}
		}
	}

	Info << "Scatter plots ..." << endl;
	{
		// Solution
		{
			const scalarField& TCells = T.internalField();
			const scalarField& csiCells = csi.internalField();
			const scalarField& sootPrecursorsYCells = sootPrecursorsY.internalField();
			const scalarField& sootSphericalYCells = sootSphericalY.internalField();
			const scalarField& sootAggregatesYCells = sootAggregatesY.internalField();
			const scalarField& sootFvCells = sootFv.internalField();
			const scalarField& sootSphericalNCells = sootSphericalN.internalField();
			const scalarField& sootSphericalD32Cells = sootSphericalD32.internalField();
			const scalarField& sootSphericalD43Cells = sootSphericalD43.internalField();
			const scalarField& sootAggregatesNCells = sootAggregatesN.internalField();
			const scalarField& sootAggregatesD32Cells = sootAggregatesD32.internalField();
			const scalarField& sootAggregatesD43Cells = sootAggregatesD43.internalField();
			const scalarField& sootSphericalHCCells = sootSphericalHC.internalField();
			const scalarField& sootAggregatesHCCells = sootAggregatesHC.internalField();
						
			// Internal cells
			const boost::filesystem::path path_out = runTime.timeName() / "SootScatter.out";
			std::ofstream fOut(path_out.c_str(), std::ios::out);

			fOut << std::left << std::setw(16) << "csi(1)";
			fOut << std::left << std::setw(16) << "T[K](2)";
			fOut << std::left << std::setw(16) << "SPrecY(3)";
			fOut << std::left << std::setw(16) << "SSpheY(4)";
			fOut << std::left << std::setw(16) << "SAggrY(5)";
			fOut << std::left << std::setw(16) << "SootFv(6)";
			fOut << std::left << std::setw(16) << "SSpheN[1/m3](7)";
			fOut << std::left << std::setw(16) << "SSpheD32[m](8)";
			fOut << std::left << std::setw(16) << "SSpheD43[m](9)";
			fOut << std::left << std::setw(16) << "SSpheHC(10)";
			fOut << std::left << std::setw(16) << "SAggrN[1/m3](11)";
			fOut << std::left << std::setw(16) << "SAggrD32[m](12)";
			fOut << std::left << std::setw(16) << "SAggrD43[m](13)";
			fOut << std::left << std::setw(16) << "SAggrHC[1/m3](14)";
			for (unsigned int i=0;i<ns;i++)
			{
				std::stringstream count; count << (15+i); std::string label =  mixture.thermodynamicsMap().NamesOfSpecies()[i] + "(" + count.str() + ")";
				fOut << std::left << std::setw(16) << label;
			}
			fOut << std::endl;

			forAll(TCells, celli)
			{	
				if (TCells[celli] >= 310.)
				{
					fOut << std::left << std::setw(16) << std::scientific << csiCells[celli];
					fOut << std::left << std::setw(16) << std::scientific << TCells[celli];
					fOut << std::left << std::setw(16) << std::scientific << sootPrecursorsYCells[celli];
					fOut << std::left << std::setw(16) << std::scientific << sootSphericalYCells[celli];
					fOut << std::left << std::setw(16) << std::scientific << sootAggregatesYCells[celli];
					fOut << std::left << std::setw(16) << std::scientific << sootFvCells[celli];
					fOut << std::left << std::setw(16) << std::scientific << sootSphericalNCells[celli];
					fOut << std::left << std::setw(16) << std::scientific << sootSphericalD32Cells[celli];
					fOut << std::left << std::setw(16) << std::scientific << sootSphericalD43Cells[celli];
					fOut << std::left << std::setw(16) << std::scientific << sootSphericalHCCells[celli];
					fOut << std::left << std::setw(16) << std::scientific << sootAggregatesNCells[celli];
					fOut << std::left << std::setw(16) << std::scientific << sootAggregatesD32Cells[celli];
					fOut << std::left << std::setw(16) << std::scientific << sootAggregatesD43Cells[celli];
					fOut << std::left << std::setw(16) << std::scientific << sootAggregatesHCCells[celli];
					for (unsigned int i=0;i<ns;i++)
						fOut << std::left << std::setw(16) << std::scientific << Y[i].internalField()[celli];
					fOut << std::endl;
				}
			}
			fOut.close();
		}
	}

	Info << "PSDF ..." << endl;
	{
		const double sootFvThreshold = 1e-14;

		double mw;
		Eigen::VectorXd y(mixture.thermodynamicsMap().NumberOfSpecies());
		Eigen::VectorXd x(mixture.thermodynamicsMap().NumberOfSpecies());

		const scalarField& TCells = T.internalField();
		const scalarField& pCells = p.internalField();
		const scalarField& sootFvCells = sootFv.internalField();

		const std::vector<double> xc{5e-3, 10e-3, 15e-3, 20e-3, 25e-3};	// m
		const double yc = 17.5e-3;					// m
		const double zc = 0.;						// m

		for (unsigned int j=0;j<xc.size();j++)
		{
			// Identify cell			
			label index = mesh.findCell(point(xc[j], yc, zc));

			if (sootFvCells[index] > sootFvThreshold)
			{
				// Mass fractions
				for (unsigned int i=0;i<ns;i++)
					y(i) = Y[i].internalField()[index];
	
				// Molecular weight and mole fractions
				mixture.thermodynamicsMap().MoleFractions_From_MassFractions(x.data(), mw, y.data());

				// Density
				const double rhoGas = pCells[index]*mw/PhysicalConstants::R_J_kmol/TCells[index];

				// Soot analysis
				sootAnalyzer.Analysis(TCells[index], pCells[index], rhoGas, y, x);
				sootAnalyzer.Distribution();

				// Write on file
				const boost::filesystem::path path_out = runTime.timeName() / ("PSDF_" + std::to_string(j+1) + ".out");
				std::ofstream fOut(path_out.c_str(), std::ios::out);
				sootAnalyzer.WriteDistributionLabel(fOut);
				sootAnalyzer.WriteDistribution(fOut, std::stod(runTime.timeName()), xc[j], yc, zc, TCells[index]);
				fOut.close();
			}
		}

		// Location of maximum soot fv
		{
			// Find max
			scalar maxSootFv = -GREAT;
			label index = -1;

			forAll(sootFv, celli)
			{
   	 			if ( maxSootFv < sootFvCells[celli] )
    				{
       		 			maxSootFv = sootFvCells[celli];
        				index = celli;
   			 	}
			}

			Info << "Max soot volume fraction: " << sootFvCells[index] << " @Cell: " << index << " (" <<  mesh.C()[index].component(0) << ", " <<  mesh.C()[index].component(1) << ")" << endl;		

			if (sootFvCells[index] > sootFvThreshold)
			{
				// Mass fractions
				for (unsigned int i=0;i<ns;i++)
					y(i) = Y[i].internalField()[index];
	
				// Molecular weight and mole fractions
				mixture.thermodynamicsMap().MoleFractions_From_MassFractions(x.data(), mw, y.data());

				// Density
				const double rhoGas = pCells[index]*mw/PhysicalConstants::R_J_kmol/TCells[index];

				// Soot analysis
				sootAnalyzer.Analysis(TCells[index], pCells[index], rhoGas, y, x);
				sootAnalyzer.Distribution();

				// Write on file
				const boost::filesystem::path path_out = runTime.timeName() / ("PSDF_Max.out");
				std::ofstream fOut(path_out.c_str(), std::ios::out);
				sootAnalyzer.WriteDistributionLabel(fOut);
				sootAnalyzer.WriteDistribution(fOut, std::stod(runTime.timeName()),  mesh.C()[index].component(0),  mesh.C()[index].component(1),  mesh.C()[index].component(2), TCells[index]);
				fOut.close();
			}
		}
	}
	
	// Write
	Info<< " * Writing ..." << endl;
	sootPrecursorsY.write();
	sootSphericalY.write();
	sootAggregatesY.write();
	sootFv.write();
	sootSphericalN.write();
	sootSphericalD32.write();
	sootSphericalD43.write();
	sootAggregatesN.write();
	sootAggregatesD32.write();
	sootAggregatesD43.write();
	sootSphericalHC.write();
	sootAggregatesHC.write();	
}
