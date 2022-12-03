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

if (calculate_mixture_fraction_and_progress_variable == true)
{
	const unsigned int ns = mixture.thermodynamicsMap().NumberOfSpecies();
	const unsigned int nfields = lookup_table_zy->list_fields().size();

	volScalarField ProgressVariable_rec
	(
	    IOobject
	    (
		"rec:C",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("dimensionless", dimensionSet(0, 0, 0, 0, 0), 0.)
	);

	volScalarField MixtureFraction_rec
	(
	    IOobject
	    (
		"rec:Z",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
    	    dimensionedScalar("dimensionless", dimensionSet(0, 0, 0, 0, 0), 0.)
	);

        volScalarField YPAH12_rec
	(
	    IOobject
	    (
		"rec:YPAH12",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
    	    dimensionedScalar("dimensionless", dimensionSet(0, 0, 0, 0, 0), 0.)
	);

        volScalarField YPAH34_rec
	(
	    IOobject
	    (
		"rec:YPAH34",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
    	    dimensionedScalar("dimensionless", dimensionSet(0, 0, 0, 0, 0), 0.)
	);

        volScalarField YPAHLP_rec
	(
	    IOobject
	    (
		"rec:YPAHLP",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
    	    dimensionedScalar("dimensionless", dimensionSet(0, 0, 0, 0, 0), 0.)
	);

        volScalarField YSP_rec
	(
	    IOobject
	    (
		"rec:YSP",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
    	    dimensionedScalar("dimensionless", dimensionSet(0, 0, 0, 0, 0), 0.)
	);

        volScalarField YAGG_rec
	(
	    IOobject
	    (
		"rec:YAGG",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
    	    dimensionedScalar("dimensionless", dimensionSet(0, 0, 0, 0, 0), 0.)
	);

	PtrList<volScalarField> tabFields_;
	tabFields_.resize(nfields);

	forAll(tabFields_, i)
	{
		tabFields_.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					"tab:Y" + lookup_table_zy->list_fields()[i],
					runTime.timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
				mesh,
				dimensionedScalar("dimless", dimensionSet(0, 0, 0, 0, 0), 0.)
			)
		);
	}


	Info << "Reconstructing mixture fraction and progress variable ..." << endl;
	{
		std::vector<double> mass_fractions(mixture.thermodynamicsMap().NumberOfSpecies());

		const scalarField& TCells = T.internalField();
		scalarField& recCCells = ProgressVariable_rec.ref(); 
		scalarField& recZCells = MixtureFraction_rec.ref();
		scalarField& recYPAH12Cells = YPAH12_rec.ref();
		scalarField& recYPAH34Cells = YPAH34_rec.ref();
		scalarField& recYPAHLPCells = YPAHLP_rec.ref();   
		scalarField& recYSPCells = YSP_rec.ref(); 
		scalarField& recYAGGCells = YAGG_rec.ref(); 
		  
		// Internal cells
		forAll(TCells, celli)
		{	
			for (unsigned int i=0;i<ns;i++)
				mass_fractions[i] = Y[i].internalField()[celli];

			recCCells[celli] = lookup_table_zy->ReconstructProgressVariable(mass_fractions.data());
			recZCells[celli] = lookup_table_zy->ReconstructMixtureFraction(mass_fractions.data());

			lookup_table_zy->Interpolate(recZCells[celli], recCCells[celli]);
			for (unsigned int i=0;i<nfields;i++)
				tabFields_[i].ref()[celli] = lookup_table_zy->Y(i);

			// PAHs (if available)
			for (unsigned int i = 0; i < lookup_table_zy->indices_pah12().size(); i++)
				recYPAH12Cells[celli] += Y[lookup_table_zy->indices_pah12()[i]].internalField()[celli];
			for (unsigned int i = 0; i < lookup_table_zy->indices_pah34().size(); i++)
				recYPAH34Cells[celli] += Y[lookup_table_zy->indices_pah34()[i]].internalField()[celli];
			for (unsigned int i = 0; i < lookup_table_zy->indices_pahlp().size(); i++)
				recYPAHLPCells[celli] += Y[lookup_table_zy->indices_pahlp()[i]].internalField()[celli];

			// Spherical particles (if available)
			for (unsigned int i = 0; i < lookup_table_zy->indices_sp().size(); i++)
				recYSPCells[celli] += Y[lookup_table_zy->indices_sp()[i]].internalField()[celli];

			// Aggregates (if available)
			for (unsigned int i = 0; i < lookup_table_zy->indices_agg().size(); i++)
				recYAGGCells[celli] += Y[lookup_table_zy->indices_agg()[i]].internalField()[celli];
		}

		// Boundaries
		volScalarField::Boundary& recCf = ProgressVariable_rec.boundaryFieldRef();
		volScalarField::Boundary& recZf = MixtureFraction_rec.boundaryFieldRef();
		volScalarField::Boundary& recYPAH12f = YPAH12_rec.boundaryFieldRef();
		volScalarField::Boundary& recYPAH34f = YPAH34_rec.boundaryFieldRef();
		volScalarField::Boundary& recYPAHLPf = YPAHLP_rec.boundaryFieldRef();
		volScalarField::Boundary& recYSPf = YSP_rec.boundaryFieldRef();
		volScalarField::Boundary& recYAGGf = YAGG_rec.boundaryFieldRef();

		forAll(T.boundaryField(), patchi)
		{
			fvPatchScalarField& pC = recCf[patchi];
			fvPatchScalarField& pZ = recZf[patchi];
			fvPatchScalarField& pYPAH12 = recYPAH12f[patchi];
			fvPatchScalarField& pYPAH34 = recYPAH34f[patchi];
			fvPatchScalarField& pYPAHLP = recYPAHLPf[patchi];
			fvPatchScalarField& pYSP = recYSPf[patchi];
			fvPatchScalarField& pYAGG = recYAGGf[patchi];
	
			forAll(pZ, facei)
			{
				for (unsigned int i=0;i<ns;i++)
					mass_fractions[i] = Y[i].boundaryField()[patchi][facei];
				pC[facei] = lookup_table_zy->ReconstructProgressVariable(mass_fractions.data());
				pZ[facei] = lookup_table_zy->ReconstructMixtureFraction(mass_fractions.data());

				lookup_table_zy->Interpolate(pZ[facei], pC[facei]);
				for (unsigned int i=0;i<nfields;i++)
					tabFields_[i].boundaryFieldRef()[patchi][facei] = lookup_table_zy->Y(i); 

				// PAHs (if available)
				for (unsigned int i = 0; i < lookup_table_zy->indices_pah12().size(); i++)
					pYPAH12[facei] += Y[lookup_table_zy->indices_pah12()[i]].boundaryField()[patchi][facei];
				for (unsigned int i = 0; i < lookup_table_zy->indices_pah34().size(); i++)
					pYPAH34[facei] += Y[lookup_table_zy->indices_pah34()[i]].boundaryField()[patchi][facei];
				for (unsigned int i = 0; i < lookup_table_zy->indices_pahlp().size(); i++)
					pYPAHLP[facei] += Y[lookup_table_zy->indices_pahlp()[i]].boundaryField()[patchi][facei];

				// Spherical particles (if available)
				for (unsigned int i = 0; i < lookup_table_zy->indices_sp().size(); i++)
					pYSP[facei] += Y[lookup_table_zy->indices_sp()[i]].boundaryField()[patchi][facei];

				// Aggregates (if available)
				for (unsigned int i = 0; i < lookup_table_zy->indices_agg().size(); i++)
					pYAGG[facei] += Y[lookup_table_zy->indices_agg()[i]].boundaryField()[patchi][facei];
			}
		}
	}

	Info << "Scatter plots ..." << endl;
	{
		// Solution
		{
			const scalarField& TCells = T.internalField();
			const scalarField& recCCells = ProgressVariable_rec.internalField();
			const scalarField& recZCells = MixtureFraction_rec.internalField();
			const scalarField& recYPAH12Cells = YPAH12_rec.internalField();
			const scalarField& recYPAH34Cells = YPAH34_rec.internalField();
			const scalarField& recYPAHLPCells = YPAHLP_rec.internalField();
			const scalarField& recYSPCells = YSP_rec.internalField();
			const scalarField& recYAGGCells = YAGG_rec.internalField();			

			// Internal cells
			const boost::filesystem::path path_out = runTime.timeName() / "Scatter:Sol.out";
			std::ofstream fOut(path_out.c_str(), std::ios::out);

			fOut << std::left << std::setw(16) << "T[K](1)";
			fOut << std::left << std::setw(16) << "Z(2)";
			fOut << std::left << std::setw(16) << "C(3)";
			fOut << std::left << std::setw(16) << "YPAH12(4)";
			fOut << std::left << std::setw(16) << "YPAH34(5)";
			fOut << std::left << std::setw(16) << "YPAHLP(6)";
			fOut << std::left << std::setw(16) << "YSP(7)";
			fOut << std::left << std::setw(16) << "YAGG(8)";
			for (unsigned int i=0;i<ns;i++)
			{
				std::stringstream count; count << (9+i); std::string label =  mixture.thermodynamicsMap().NamesOfSpecies()[i] + "(" + count.str() + ")";
				fOut << std::left << std::setw(16) << label;
			}
			fOut << std::endl;

			forAll(TCells, celli)
			{	
				if (TCells[celli] >= tabulatedChemistryTemperatureThreshold)
				{
					fOut << std::left << std::setw(16) << std::scientific << TCells[celli];
					fOut << std::left << std::setw(16) << std::scientific << recZCells[celli];
					fOut << std::left << std::setw(16) << std::scientific << recCCells[celli];
					fOut << std::left << std::setw(16) << std::scientific << recYPAH12Cells[celli];
					fOut << std::left << std::setw(16) << std::scientific << recYPAH34Cells[celli];
					fOut << std::left << std::setw(16) << std::scientific << recYPAHLPCells[celli];
					fOut << std::left << std::setw(16) << std::scientific << recYSPCells[celli];
					fOut << std::left << std::setw(16) << std::scientific << recYAGGCells[celli];
					for (unsigned int i=0;i<ns;i++)
						fOut << std::left << std::setw(16) << std::scientific << Y[i].internalField()[celli];
					fOut << std::endl;
				}
			}
			fOut.close();
		}

		// Reconstruction
		{
			const scalarField& TCells = T.internalField();

			// Internal cells
			const boost::filesystem::path path_out = runTime.timeName() / "Scatter:Rec.out";
			std::ofstream fOut(path_out.c_str(), std::ios::out);

			for (unsigned int i=0;i<nfields;i++)
			{
				std::stringstream count; count << i+1; std::string label =  lookup_table_zy->list_fields()[i] + "(" + count.str() + ")";
				fOut << std::left << std::setw(16) << label;
			}
			fOut << std::endl;

			forAll(TCells, celli)
			{
				if (TCells[celli] >= tabulatedChemistryTemperatureThreshold)
				{
					for (unsigned int i=0;i<nfields;i++)
						fOut << std::left << std::setw(16) << std::scientific << tabFields_[i].ref()[celli];
					fOut << std::endl;
				}
			}
			fOut.close();
		}

		// Comparison
		{
			const scalarField& TCells = T.internalField();

			const scalarField& recYPAH12Cells = YPAH12_rec.internalField();
			const scalarField& recYPAH34Cells = YPAH34_rec.internalField();
			const scalarField& recYPAHLPCells = YPAHLP_rec.internalField();
			const scalarField& recYSPCells = YSP_rec.internalField();
			const scalarField& recYAGGCells = YAGG_rec.internalField();

			// Internal cells
			const boost::filesystem::path path_out = runTime.timeName() / "Scatter:Comparison.out";
			std::ofstream fOut(path_out.c_str(), std::ios::out);

			for (unsigned int i=0;i<nfields;i++)
			{
				// Solution
				{
					std::stringstream count; count << i+1; std::string label =  lookup_table_zy->list_fields()[i] + "[S]" + "(" + count.str() + ")";
					fOut << std::left << std::setw(16) << label;
				}
				// Reconstruction
				{
					std::stringstream count; count << i+1; std::string label =  lookup_table_zy->list_fields()[i] + "[R]" + "(" + count.str() + ")";
					fOut << std::left << std::setw(16) << label;
				}
			}
			fOut << std::endl;

			forAll(TCells, celli)
			{
				if (TCells[celli] >= tabulatedChemistryTemperatureThreshold)
				{
					for (unsigned int i=0;i<nfields;i++)
					{
						// Solution
						if (lookup_table_zy->list_fields()[i] == "PAH12")
							fOut << std::left << std::setw(16) << std::scientific << recYPAH12Cells[celli];
						else if (lookup_table_zy->list_fields()[i] == "PAH34")
							fOut << std::left << std::setw(16) << std::scientific << recYPAH34Cells[celli];
						else if (lookup_table_zy->list_fields()[i] == "PAHLP")
							fOut << std::left << std::setw(16) << std::scientific << recYPAHLPCells[celli];
						else if (lookup_table_zy->list_fields()[i] == "SP")
							fOut << std::left << std::setw(16) << std::scientific << recYSPCells[celli];
						else if (lookup_table_zy->list_fields()[i] == "AGG")
							fOut << std::left << std::setw(16) << std::scientific << recYAGGCells[celli];
						else if (lookup_table_zy->list_fields()[i] == "T")
							fOut << std::left << std::setw(16) << std::scientific << TCells[celli];
						else
						{
							const unsigned int j = mixture.thermodynamicsMap().IndexOfSpecies(lookup_table_zy->list_fields()[i]) - 1;
							fOut << std::left << std::setw(16) << std::scientific << Y[j].internalField()[celli];
						}

						// Reconstruction
						fOut << std::left << std::setw(16) << std::scientific << tabFields_[i].ref()[celli];
					}
					fOut << std::endl;
				}
			}
			fOut.close();
		}
	}
	
	// Write
	Info<< " * Writing ..." << endl;
	MixtureFraction_rec.write();
	ProgressVariable_rec.write();
	YPAH12_rec.write();
	YPAH34_rec.write();
	YPAHLP_rec.write();
	YSP_rec.write();
	YAGG_rec.write();
	for (unsigned int i=0;i<nfields;i++)
		tabFields_[i].write();
}
