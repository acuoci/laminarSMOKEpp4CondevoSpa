{		
	// Const fields
	const scalarField& YcCells = Yc;
	const scalarField& pCells = p;

	// Non-const fields
	scalarField& TCells = T.primitiveFieldRef();
	scalarField& psiCells = psi.primitiveFieldRef();
	scalarField& cpCells = Cp.primitiveFieldRef();
	scalarField& muCells = mu.primitiveFieldRef();
	scalarField& lambdaCells = lambda.primitiveFieldRef();
	scalarField& alphaCells = alpha.primitiveFieldRef();
	scalarField& mwCells = MW.primitiveFieldRef();
	scalarField& cTotCells = cTot.primitiveFieldRef();

	// Loop over all the internal cells
	forAll(TCells, celli)
	{	
		// Interpolation
		lookup_table_fpi->Interpolate(YcCells[celli]);

		// Temperature (K)
		TCells[celli] = lookup_table_fpi->T();

		// Temperature (K)
		mwCells[celli] = lookup_table_fpi->MW();
			
		// Total concentration [kmol/m3]
		cTotCells[celli] = pCells[celli]/PhysicalConstants::R_J_kmol/TCells[celli];

		// Compressibility [s2/m2]
		psiCells[celli]  = cTotCells[celli]*mwCells[celli]/pCells[celli];

		// Constant pressure specific heat [J/kg/K]
		cpCells[celli]   = lookup_table_fpi->Cp();

		// Viscosity [kg/m/s]
		muCells[celli] = lookup_table_fpi->mu();
		
		// Thermal conductivity [W/m/K]
		lambdaCells[celli] = lookup_table_fpi->lambda();

		// Thermal diffusivity [m2/s]
		alphaCells[celli] = lambdaCells[celli]/cpCells[celli]/(pCells[celli]*psiCells[celli]);
	}

	// Boundaries
	{
		volScalarField::Boundary& ycBf 		= Yc.boundaryFieldRef();
		volScalarField::Boundary& pBf 		= p.boundaryFieldRef();
		volScalarField::Boundary& TBf 		= T.boundaryFieldRef();
		volScalarField::Boundary& psiBf 	= psi.boundaryFieldRef();
		volScalarField::Boundary& cpBf 		= Cp.boundaryFieldRef();

		volScalarField::Boundary& muBf 		= mu.boundaryFieldRef();
		volScalarField::Boundary& lambdaBf 	= lambda.boundaryFieldRef();
		volScalarField::Boundary& alphaBf 	= alpha.boundaryFieldRef();
		
		volScalarField::Boundary& mwBf 		= MW.boundaryFieldRef();
		volScalarField::Boundary& cTotBf 	= cTot.boundaryFieldRef();

		forAll(T.boundaryField(), patchi)
		{
			fvPatchScalarField& pyc = ycBf[patchi];
			fvPatchScalarField& pp = pBf[patchi];
			fvPatchScalarField& pT = TBf[patchi];
			fvPatchScalarField& ppsi = psiBf[patchi];
			fvPatchScalarField& pcp = cpBf[patchi];

			fvPatchScalarField& pmu = muBf[patchi];
			fvPatchScalarField& plambda = lambdaBf[patchi];
			fvPatchScalarField& palpha = alphaBf[patchi];
			fvPatchScalarField& pmw = mwBf[patchi];
			fvPatchScalarField& pcTot = cTotBf[patchi];
	
			forAll(pT, facei)
			{
				// Interpolation
				lookup_table_fpi->Interpolate(pyc[facei]);

				// Temperature (K)
				pT[facei] = lookup_table_fpi->T();

				// Temperature (K)
				pmw[facei] = lookup_table_fpi->MW();
			
				// Total concentration [kmol/m3]
				pcTot[facei] = pp[facei]/PhysicalConstants::R_J_kmol/pT[facei];

				// Compressibility [s2/m2]
				ppsi[facei]  = pcTot[facei]*pmw[facei]/pp[facei];

				// Constant pressure specific heat [J/kg/K]
				pcp[facei]   = lookup_table_fpi->Cp();

				// Viscosity [kg/m/s]
				pmu[facei] = lookup_table_fpi->mu();
		
				// Thermal conductivity [W/m/K]
				plambda[facei] = lookup_table_fpi->lambda();

				// Thermal diffusivity [m2/s]
				palpha[facei] = plambda[facei]/pcp[facei]/(pp[facei]*ppsi[facei]);
			}
		}
	}
