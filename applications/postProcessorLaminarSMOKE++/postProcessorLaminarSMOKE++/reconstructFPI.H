if (reconstruct_fpi == true)
{
	// Target kinetics
	OpenSMOKE::ThermodynamicsMap_CHEMKIN* fpiThermodynamicsMap; 
	{
		Foam::string kinetics_folder = postProcessingDictionary.lookup("fpiKinetics");

		//- Reading homogeneous data
		Info<< " * Loading homogeneous kinetics...\n" << endl;
		boost::filesystem::path path_kinetics = kinetics_folder;
		boost::property_tree::ptree ptree;
    		boost::property_tree::read_xml( (path_kinetics / "kinetics.xml").c_str(), ptree );

		// Import map
		fpiThermodynamicsMap = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(ptree); 
	}


	// FPI Table
	Info << "Reading lookup table..." << endl;
	OpenSMOKE::LookUpTable_FPI* lookup_table_fpi;
	{
		// Folder containing the tables
		const Foam::string table_folder = postProcessingDictionary.lookup("fpiFolder");
		const boost::filesystem::path table_folder_path = table_folder;
		const boost::filesystem::path path_to_csv = table_folder_path / "fpi.main.csv";

		// Threshold temperature for scatter plots
		// tabulatedChemistryTemperatureThreshold = postProcessingDictionary.lookupOrDefault<scalar>("tabulatedChemistryTemperatureThreshold", 300.);

		// Reading tables
		lookup_table_fpi = new OpenSMOKE::LookUpTable_FPI(*fpiThermodynamicsMap);
		lookup_table_fpi->Setup(path_to_csv, true);
	}
			
	// Progress variable field
	volScalarField Yc
	(
	    IOobject
	    (
		"Yc",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("dimensionless", dimensionSet(0, 0, 0, 0, 0), 0.)
	);

	// Create reconstructed FPI mass fraction fields
	Info << "Create reconstructed FPI mass fraction fields..." << endl;
	PtrList<volScalarField> fpiY(fpiThermodynamicsMap->NumberOfSpecies());
	tmp<volScalarField> tYdefault;
	for (unsigned int i=0;i<fpiThermodynamicsMap->NumberOfSpecies();i++)
	{
		IOobject header
		(
			fpiThermodynamicsMap->NamesOfSpecies()[i],
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ
		);

		// Check if field exists and can be read
		if (header.typeHeaderOk<volScalarField>(true))
		{
			Info << "Reading " << fpiThermodynamicsMap->NamesOfSpecies()[i] << endl;

			fpiY.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						fpiThermodynamicsMap->NamesOfSpecies()[i],
						mesh.time().timeName(),
						mesh,
						IOobject::MUST_READ,
						IOobject::AUTO_WRITE
					),
					mesh
				)
			);
		}
		else
		{
			// Read Ydefault if not already read
			if (!tYdefault.valid())
			{
				word YdefaultName("Ydefault");

				IOobject timeIO ( YdefaultName, mesh.time().timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE );
				IOobject constantIO ( YdefaultName,mesh.time().constant(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE );
				IOobject time0IO ( YdefaultName,Time::timeName(0),mesh,IOobject::MUST_READ,IOobject::NO_WRITE );

				if (timeIO.typeHeaderOk<volScalarField>(true))
				{
					tYdefault = new volScalarField(timeIO, mesh);
				}
				else if (constantIO.typeHeaderOk<volScalarField>(true))
				{
					tYdefault = new volScalarField(constantIO, mesh);
				}
				else
				{
					tYdefault = new volScalarField(time0IO, mesh);
				}
			}

			fpiY.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						fpiThermodynamicsMap->NamesOfSpecies()[i],
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::AUTO_WRITE
					),
					tYdefault()
				)
			);
		}
	}

	Info << "Reconstructing FPI ..." << endl;
	{
		const scalarField& TCells = T.internalField();
		scalarField& ycCells = Yc.ref(); 
		  
		// Internal cells
		forAll(TCells, celli)
		{	
			// Reconstruct progress variable
			ycCells[celli] = 0.;
			for (unsigned int i=0;i<fpiWeights.size();i++)
				ycCells[celli] += Y[fpiIndices[i]].internalField()[celli] * fpiWeights[i];

			// Interpolation
			lookup_table_fpi->Interpolate(ycCells[celli]);
			for (unsigned int i=0;i<fpiThermodynamicsMap->NumberOfSpecies();i++)
				fpiY[i].ref()[celli] = std::max(lookup_table_fpi->Y(i),0.);
				
		}

		// Boundaries
		volScalarField::Boundary& ycf = Yc.boundaryFieldRef();

		forAll(T.boundaryField(), patchi)
		{
			fvPatchScalarField& pYc = ycf[patchi];
	
			forAll(pYc, facei)
			{
				pYc[facei] = 0.;
				for (unsigned int i=0;i<fpiWeights.size();i++)
					pYc[facei] += Y[fpiIndices[i]].boundaryField()[patchi][facei] * fpiWeights[i];

				// Interpolation
				lookup_table_fpi->Interpolate(pYc[facei]);
				for (unsigned int i=0;i<fpiThermodynamicsMap->NumberOfSpecies();i++)
					fpiY[i].boundaryFieldRef()[patchi][facei] = std::max(lookup_table_fpi->Y(i), 0.);
			}
		}
	}
	
	// Write
	Info<< " * Writing ..." << endl;
	Yc.write();
	for (unsigned int i=0;i<fpiThermodynamicsMap->NumberOfSpecies();i++)
		fpiY[i].write();
}
