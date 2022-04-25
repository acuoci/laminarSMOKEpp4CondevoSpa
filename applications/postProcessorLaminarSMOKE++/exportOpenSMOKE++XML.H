if (export_opensmokeppxml_file == true)
{
	Info<< "Export OpenSMOKE++XML file..." << endl;

	// Optional fields
	IOobject header
        (
        	"soot:fv",
        	mesh.time().timeName(),
        	mesh,
        	IOobject::NO_READ
    	);

	// Recognize geometry
	const Vector<label> validDirs = (mesh.geometricD() + Vector<label>::one)/2;
     	Info<< "    Mesh has " << mesh.nGeometricD() << " geometric (non-empty/wedge) directions " << validDirs << endl;

	// Open file
	std::ofstream fXML;
	fXML.open( runTime.timeName() + "/status.xml", std::ios::out);
	fXML.setf(std::ios::scientific);

	fXML << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
	fXML << "<opensmoke version=\"0.1a\">" << std::endl;
	
	fXML << "<Type>Flame2D</Type>" << std::endl;

	const unsigned int n_additional = 6 + 16 + (header.typeHeaderOk<volScalarField>(true) ? 5 : 0);

	unsigned int count = 2;
	fXML << "<additional>" << std::endl;

	fXML << n_additional << std::endl;

	fXML << "time [s] " << count++ << std::endl;
	fXML << "temperature [K] " << count++ << std::endl;
	fXML << "pressure [Pa] " << count++ << std::endl;
	fXML << "mol-weight [kg/kmol] " << count++ << std::endl;
	fXML << "density [kg/m3] " << count++ << std::endl;
	fXML << "heat-release [W/m3] " << count++ << std::endl;

	fXML << "x-coord [m] " << count++ << std::endl;
	fXML << "y-coord [m] " << count++ << std::endl;
	fXML << "z-coord [m] " << count++ << std::endl;
	fXML << "x-vel [m/s] " << count++ << std::endl;
	fXML << "y-vel [m/s] " << count++ << std::endl;
	fXML << "z-vel [m/s] " << count++ << std::endl;
	fXML << "volume [m3] " << count++ << std::endl;
	fXML << "mix-fract [-] " << count++ << std::endl;
	fXML << "tau [s] " << count++ << std::endl;
	fXML << "chi-st [1/s] " << count++ << std::endl;
	fXML << "strain-rate [1/s] " << count++ << std::endl;
	fXML << "curvature [1/m] " << count++ << std::endl;
	fXML << "enthalpy-sens [J/kg] " << count++ << std::endl;
	fXML << "cp [J/kg/K] " << count++ << std::endl;
	fXML << "therm-cond [W/m/K] " << count++ << std::endl;
	fXML << "viscosity [kg/m/s] " << count++ << std::endl;

	if (header.typeHeaderOk<volScalarField>(true))
	{
		fXML << "YPAH12 [-] " << count++ << std::endl;	// mass fraction of PAHs with 1-2 aromatic rings
		fXML << "YPAH34 [-] " << count++ << std::endl;	// mass fraction of PAHs with 3-4 aromatic rings
		fXML << "YPAHLP [-] " << count++ << std::endl;	// mass fraction of PAHs with more than 4 aromatic rings
		fXML << "YSoot [-] "  << count++ << std::endl;	// mass fraction of soot
		fXML << "fvSoot [-] " << count++ << std::endl;	// volume fraction of soot
	}

	fXML << "</additional>" << std::endl;

	fXML << "<t-p-mw>1 2 3</t-p-mw>" << std::endl;

	fXML << "<mass-fractions>" << std::endl;
	fXML << mixture.thermodynamicsMap().NumberOfSpecies() << std::endl;
	for (unsigned int j = 0; j < mixture.thermodynamicsMap().NumberOfSpecies(); j++)
		fXML << mixture.thermodynamicsMap().NamesOfSpecies()[j] << " " << mixture.thermodynamicsMap().MWs()[j] << " " << count++ << std::endl;
	fXML << "</mass-fractions>" << std::endl;

	fXML << "<profiles-size>" << std::endl;
	fXML << mesh.cells().size() << " " << n_additional + mixture.thermodynamicsMap().NumberOfSpecies() << std::endl;
	fXML << "</profiles-size>" << std::endl;


	// Pointers to relevant fields
	const scalarField& TCells = T.internalField();
	const scalarField& pCells = p.internalField();
	const scalarField& rhoCells = rho.internalField();
	const scalarField& hsCells = mixture.Hs().internalField();
	const scalarField& cpCells = Cp.internalField();
	const scalarField& QCells = mixture.Q().internalField();
	const scalarField& lambdaCells = lambda.internalField();
	const scalarField& muCells = mu.internalField();
	const scalarField& csiCells = csi.internalField();
	const scalarField& vCells = mesh.V();
	const vectorField& uCells = U.internalField();
	const scalarField& tauCells = tau.internalField();

	// Scalar dissipation rate and strain rate
	volScalarField chiSt
	(
		IOobject
		(
			"chiSt",
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		),
		mesh
	);

	volScalarField strainRate
	(
		IOobject
		(
			"strainRate",
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		),
		mesh
	);

	volScalarField curvature
	(
		IOobject
		(
			"curvature",
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		),
		mesh
	);

	const scalarField& chiStCells = chiSt.internalField();
	const scalarField& strainRateCells = strainRate.internalField();
	const scalarField& curvatureCells = curvature.internalField();

	PtrList<volScalarField> optionalFields;
	if (header.typeHeaderOk<volScalarField>(true))
    	{
		optionalFields.resize(5);

		optionalFields.set ( 0, new volScalarField ( IOobject ( "pah:Y12", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE ), mesh ) );
		optionalFields.set ( 1, new volScalarField ( IOobject ( "pah:Y34", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE ), mesh ) );
		optionalFields.set ( 2, new volScalarField ( IOobject ( "pah:YLP", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE ), mesh ) );
		optionalFields.set ( 3, new volScalarField ( IOobject ( "soot:Y",  runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE ), mesh ) );
		optionalFields.set ( 4, new volScalarField ( IOobject ( "soot:fv", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE ), mesh ) );	
    	}

	fXML << "<profiles>" << std::endl;

	// Loop over all the cells
	forAll(TCells, celli)
	{
		// Total concentration (kmol/m3)
		const double ctot = pCells[celli]/PhysicalConstants::R_J_kmol/TCells[celli];

		// Molecular weight (kg/kmol)
		const double mwtot = rhoCells[celli]/ctot;

		// Write on file
		fXML << runTime.timeName() << " ";	// time (s)
		fXML << TCells[celli] << " ";		// temperature (K)
		fXML << pCells[celli] << " ";		// pressure (Pa)
		fXML << mwtot << " ";			// molecular weight (kg/kmol)
		fXML << rhoCells[celli] << " ";		// density (kg/m3)
		fXML << QCells[celli] << " ";		// heat release rate (W/m3)

		// Coordinates
		if (validDirs[0] != 0)	fXML << mesh.C()[celli][0] << " ";	// x-coordinate (m)
		else			fXML << 0 << " ";			// x-coordinate (m)
		if (validDirs[1] != 0)	fXML << mesh.C()[celli][1] << " ";	// y-coordinate (m)
		else			fXML << 0 << " ";			// y-coordinate (m)
		if (validDirs[2] != 0)	fXML << mesh.C()[celli][2] << " ";	// z-coordinate (m)
		else			fXML << 0 << " ";			// z-coordinate (m)

		// Velocities
		if (validDirs[0] != 0)	fXML << uCells[celli][0] << " ";	// x-velocity (m/s)
		else			fXML << 0 << " ";			// x-velocity (m/s)
		if (validDirs[1] != 0)	fXML << uCells[celli][1] << " ";	// y-velocity (m/s)
		else			fXML << 0 << " ";			// y-velocity (m/s)
		if (validDirs[2] != 0)	fXML << uCells[celli][2] << " ";	// z-velocity (m/s)
		else			fXML << 0 << " ";			// z-velocity (m/s)

		fXML << vCells[celli] << " ";		// cell volume (m3)
		fXML << csiCells[celli] << " ";		// mixture fraction (-)
		fXML << tauCells[celli] << " ";		// local residence time (s)
		
		fXML << chiStCells[celli] << " ";	// scalar dissipation rate (1/s)
		fXML << strainRateCells[celli] << " ";	// strain rate (1/s)
		fXML << curvatureCells[celli] << " ";	// curvature (1/m)
		
		fXML << hsCells[celli] << " ";		// sensible enthalpy (J/kg)		
		fXML << cpCells[celli] << " ";		// constant-pressure specific heat (J/kg/K)
		fXML << lambdaCells[celli] << " ";	// thermal conductivity (W/m/K)
		fXML << muCells[celli] << " ";		// dynamic viscosity (kg/m/s)

		if (optionalFields.size() != 0)
		{
			fXML << optionalFields[0][celli] << " ";	// mass fraction of PAHs with 1-2 aromatic rings
			fXML << optionalFields[1][celli] << " ";	// mass fraction of PAHs with 3-4 aromatic rings
			fXML << optionalFields[2][celli] << " ";	// mass fraction of PAHs with more than 4 aromatic rings
			fXML << optionalFields[3][celli] << " ";	// mass fraction of soot
			fXML << optionalFields[4][celli] << " ";	// volume fraction of soot
		}

		// Sum of mass fractions
		double sum = 0.;
		for (unsigned int j = 0; j < mixture.thermodynamicsMap().NumberOfSpecies(); j++)
			sum += Y[j][celli];

		// Write mass fractions
		const double threshold = 1.e-32;
		for (unsigned int j = 0; j < mixture.thermodynamicsMap().NumberOfSpecies(); j++)
		{
			if (Y[j][celli] < 1e-32)	fXML << 0 << " ";
			else				fXML << Y[j][celli]/sum << " ";
		}

		fXML << std::endl;
	}
	fXML << "</profiles>" << std::endl;

	fXML << "</opensmoke>" << std::endl;
	fXML.close();
}

