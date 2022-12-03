if (export_csv_file == true)
{
	Info<< "Export CSV file..." << endl;

	// Recognize geometry
	const Vector<label> validDirs = (mesh.geometricD() + Vector<label>::one)/2;
     	Info<< "    Mesh has " << mesh.nGeometricD() << " geometric (non-empty/wedge) directions " << validDirs << endl;

	// Open file
	std::ofstream fCSV;
	fCSV.open( runTime.timeName() + "/status.csv", std::ios::out);
	fCSV.setf(std::ios::scientific);

	// Pointers to relevant fields
	const scalarField& TCells = T.internalField();
	const scalarField& pCells = p.internalField();
	const scalarField& rhoCells = rho.internalField();
	const scalarField& hsCells = mixture.Hs().internalField();
	const scalarField& cpCells = Cp.internalField();
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

	// Optional fields
	IOobject header
        (
        	"soot:fv",
        	mesh.time().timeName(),
        	mesh,
        	IOobject::NO_READ
    	);

	PtrList<volScalarField> optionalFields;
	if (header.typeHeaderOk<volScalarField>(true))
    	{
		optionalFields.resize(5);

		optionalFields.set ( 0, new volScalarField ( IOobject ( "pah:Y12", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE ), mesh ) );
		optionalFields.set ( 1, new volScalarField ( IOobject ( "pah:Y34", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE ), mesh ) );
		optionalFields.set ( 2, new volScalarField ( IOobject ( "pah:YLP", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE ), mesh ) );
		optionalFields.set ( 3, new volScalarField ( IOobject ( "soot:Y", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE ), mesh ) );
		optionalFields.set ( 4, new volScalarField ( IOobject ( "soot:fv", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE ), mesh ) );	
    	}

	// Header line
	fCSV << "t[s],";		// time (s)
	fCSV << "x[m],";		// x-coordinate (m)
	fCSV << "y[m],";		// y-coordinate (m)
	fCSV << "z[m],";		// z-coordinate (m)
	fCSV << "vx[m/s],";		// x-velocity (m/s)
	fCSV << "vy[m/s],";		// y-velocity (m/s)
	fCSV << "vz[m/s],";		// z-velocity (m/s)
	fCSV << "V[m3],";		// cell volume (m3)
	fCSV << "csi,";			// mixture fraction
	fCSV << "tau[s],";		// local residence time (s)
	fCSV << "T[K],";		// temperature (K)
	fCSV << "P[Pa],";		// pressure (Pa)
	fCSV << "chiSt[1/s],";		// scalar dissipation rate (1/s)
	fCSV << "strain-rate[1/s],";	// strain rate (1/s)
	fCSV << "curvature[1/m],";	// curvature (1/m)
	fCSV << "Hs[J/kg],";		// sensible enthalpy (J/kg)
	fCSV << "MW[kg/kmol],";		// molecular weight (kg/kmol)
	fCSV << "rho[kg/m3],";		// density (kg/m3)		
	fCSV << "cp[J/kg/K],";		// constant-pressure specific heat (J/kg/K)
	fCSV << "k[W/m/K],";		// thermal conductivity (W/m/K)
	fCSV << "mu[kg/m/s],";		// dynamic viscosity (kg/m/s)

	if (optionalFields.size() != 0)
	{
		fCSV << "YPAH12,";	// mass fraction of PAHs with 1-2 aromatic rings
		fCSV << "YPAH34,";	// mass fraction of PAHs with 3-4 aromatic rings
		fCSV << "YPAHLP,";	// mass fraction of PAHs with more than 4 aromatic rings
		fCSV << "YSoot,";	// mass fraction of soot
		fCSV << "fvSoot,";	// volume fraction of soot
	}

	for (unsigned int j = 0; j < mixture.thermodynamicsMap().NumberOfSpecies(); j++)
		fCSV << (mixture.thermodynamicsMap().NamesOfSpecies()[j] + ",");	// mass fractions of species

	fCSV << std::endl;

	// Loop over all the cells
	forAll(TCells, celli)
	{
		// Total concentration (kmol/m3)
		const double ctot = pCells[celli]/PhysicalConstants::R_J_kmol/TCells[celli];

		// Molecular weight (kg/kmol)
		const double mwtot = rhoCells[celli]/ctot;

		// Write on file
		fCSV << runTime.timeName() << ",";	// time (s)

		// Coordinates
		if (validDirs[0] != 0)	fCSV << mesh.C()[celli][0] << ",";	// x-coordinate (m)
		else			fCSV << 0 << ",";			// x-coordinate (m)
		if (validDirs[1] != 0)	fCSV << mesh.C()[celli][1] << ",";	// y-coordinate (m)
		else			fCSV << 0 << ",";			// y-coordinate (m)
		if (validDirs[2] != 0)	fCSV << mesh.C()[celli][2] << ",";	// z-coordinate (m)
		else			fCSV << 0 << ",";			// z-coordinate (m)

		// Velocities
		if (validDirs[0] != 0)	fCSV << uCells[celli][0] << ",";	// x-velocity (m/s)
		else			fCSV << 0 << ",";			// x-velocity (m/s)
		if (validDirs[1] != 0)	fCSV << uCells[celli][1] << ",";	// y-velocity (m/s)
		else			fCSV << 0 << ",";			// y-velocity (m/s)
		if (validDirs[2] != 0)	fCSV << uCells[celli][2] << ",";	// z-velocity (m/s)
		else			fCSV << 0 << ",";			// z-velocity (m/s)

		fCSV << vCells[celli] << ",";		// cell volume (m3)
		fCSV << csiCells[celli] << ",";		// mixture fraction (-)
		fCSV << tauCells[celli] << ",";		// local residence time (s)
		fCSV << TCells[celli] << ",";		// temperature (K)
		fCSV << pCells[celli] << ",";		// pressure (Pa)
		fCSV << chiStCells[celli] << ",";	// scalar dissipation rate (1/s)
		fCSV << strainRateCells[celli] << ",";	// strain rate (1/s)
		fCSV << curvatureCells[celli] << ",";	// curvature (1/m)
		fCSV << hsCells[celli] << ",";		// sensible enthalpy (J/kg)
		fCSV << mwtot << ",";			// molecular weight (kg/kmol)
		fCSV << rhoCells[celli] << ",";		// density (kg/m3)		
		fCSV << cpCells[celli] << ",";		// constant-pressure specific heat (J/kg/K)
		fCSV << lambdaCells[celli] << ",";	// thermal conductivity (W/m/K)
		fCSV << muCells[celli] << ",";		// dynamic viscosity (kg/m/s)

		if (optionalFields.size() != 0)
		{
			fCSV << optionalFields[0][celli] << ",";	// mass fraction of PAHs with 1-2 aromatic rings
			fCSV << optionalFields[1][celli] << ",";	// mass fraction of PAHs with 3-4 aromatic rings
			fCSV << optionalFields[2][celli] << ",";	// mass fraction of PAHs with more than 4 aromatic rings
			fCSV << optionalFields[3][celli] << ",";	// mass fraction of soot
			fCSV << optionalFields[4][celli] << ",";	// volume fraction of soot
		}

		// Sum of mass fractions
		double sum = 0.;
		for (unsigned int j = 0; j < mixture.thermodynamicsMap().NumberOfSpecies(); j++)
			sum += Y[j][celli];

		// Write mass fractions
		const double threshold = 1.e-32;
		for (unsigned int j = 0; j < mixture.thermodynamicsMap().NumberOfSpecies(); j++)
		{
			if (Y[j][celli] < 1e-32)	fCSV << 0 << ",";
			else				fCSV << Y[j][celli]/sum << ",";
		}

		fCSV << std::endl;
	}

	fCSV.close();
}

