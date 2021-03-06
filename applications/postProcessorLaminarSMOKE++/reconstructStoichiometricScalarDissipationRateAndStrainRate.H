if (calculate_scalar_dissipation_rate == true)
{
	Info<< "Reconstructing stoichiometric scalar dissipation rate and strain rate..." << endl;

	Info<< " * Reading csi field..." << endl;
	volScalarField csi
	(
		IOobject
		(
			"csi",
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		mesh
	);

	volScalarField chiSt
	(
	    IOobject
	    (
		"chiSt",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
    	    dimensionSet(0, 0, -1, 0, 0)
	);

	volScalarField strainRate
	(
		IOobject
		(
			"strainRate",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionSet(0, 0, -1, 0, 0)
	);

	dimensionedScalar epsilon_mag_grad_csi
	(
		"epsilon_mag_grad_csi",
		dimensionSet(0,-1,0,0,0,0,0),
		scalar(1.e-9)
	); 

	// Calculate
	Info<< " * Reconstructing ..." << endl;
	volVectorField grad_csi = fvc::grad(csi);
	volScalarField mag_grad_csi_uncorrected = mag(grad_csi); 
	volScalarField mag_grad_csi = mag(grad_csi) + epsilon_mag_grad_csi; 
	chiSt = 2.*alpha*mag(sqr(grad_csi));
	strainRate = -( sqr(grad_csi) && fvc::grad(U) )/sqr(mag_grad_csi)+ fvc::div(U); 

	// Clean
	Info<< " * Cleaning ..." << endl;
	{
		const scalarField& csiCells = csi.internalField();
		const scalarField& mag_grad_csi_uncorrectedCells = mag_grad_csi_uncorrected.internalField();

		scalarField& chiStCells = chiSt.ref(); 
		scalarField& strainRateCells = strainRate.ref();   

		forAll(csiCells, celli)
		{	
			if (csiCells[celli] <= 1e-3 || csiCells[celli] >= 0.999)
			{
				chiStCells[celli] = 0.;
				strainRateCells[celli] = 0.;
			}

			if (mag_grad_csi_uncorrectedCells[celli] < 1.e-9)
			{
				chiStCells[celli] = 0.;
				strainRateCells[celli] = 0.;
			}
		}
	}
	
	// Write
	Info<< " * Writing ..." << endl;
	chiSt.write();
	strainRate.write();
}
