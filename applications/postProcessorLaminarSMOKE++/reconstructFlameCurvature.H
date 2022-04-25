if (calculate_flame_curvature == true)
{
	Info<< "Reconstructing flame curvature..." << endl;

	volScalarField curvature
	(
		IOobject
		(
			"curvature",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionSet(0, -1, 0, 0, 0)
	);

	dimensionedScalar epsilon_mag_grad_T
	(
		"epsilon_mag_grad_T",
		dimensionSet(0,-1,0,1,0,0,0),
		scalar(1.e-6)
	); 

	// Calculate
	Info<< " * Calculating..." << endl;
	volVectorField grad_T = fvc::grad(T);
	volScalarField mag_grad_T_uncorrected = mag(grad_T); 
	volScalarField mag_grad_T = mag(grad_T) + epsilon_mag_grad_T;
	curvature = fvc::div(-grad_T/ mag_grad_T);

	// Cleaning
	const double minT = 320.;
	const double threshold = 1.e-6;
	Info<< " * Cleaning..." << endl;
	{
		const scalarField& TCells = T.internalField();
		const scalarField& mag_grad_T_uncorrectedCells = mag_grad_T_uncorrected.internalField();

		scalarField& curvatureCells = curvature.ref();   

		forAll(TCells, celli)
		{	
			if (TCells[celli] <= minT)
			{
				curvatureCells[celli] = 0.;
			}

			if (mag_grad_T_uncorrectedCells[celli] < threshold)
			{
				curvatureCells[celli] = 0.;
			}
		}
	}
	
	// Write
	Info<< " * Writing..." << endl;
	curvature.write();
}

