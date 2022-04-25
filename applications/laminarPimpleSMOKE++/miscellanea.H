if (export_opensmokepp_file == true)
{
	// Every point
	{
		Info<< "Exporting OpenSMOKE++ file ..." << endl;

		const unsigned int n_additional = 13;

		const std::string file_name = "output.xml";

		std::ofstream fXML;
		fXML.open(file_name.c_str(), std::ios::out);
		fXML.setf(std::ios::scientific);

		fXML << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
		fXML << "<opensmoke version=\"0.1a\">" << std::endl;

		fXML << "<Type>MultiDimFlame</Type>" << std::endl;

		fXML << "<additional>" << std::endl;
		fXML << n_additional << std::endl;
		fXML << "axial-coordinate [cm] 2" << std::endl;
		fXML << "temperature [K] 3" << std::endl;
		fXML << "pressure [Pa] 4" << std::endl;
		fXML << "mol-weight [kg/kmol] 5" << std::endl;
		fXML << "density [kg/m3] 6" << std::endl;
		fXML << "heat-release [W/m3] 7" << std::endl;
		fXML << "axial-velocity [m/s] 8" << std::endl;
		fXML << "vy [m/s] 9" << std::endl;
		fXML << "vz [m/s] 10" << std::endl;
		fXML << "volume [m3] 11" << std::endl;
		fXML << "csi [-] 12" << std::endl;
		fXML << "y [-] 13" << std::endl;
		fXML << "z [-] 14" << std::endl;
		fXML << "</additional>" << std::endl;

		fXML << "<t-p-mw>" << std::endl;
		fXML << "1 2 3" << std::endl;
		fXML << "</t-p-mw>" << std::endl;

		fXML << "<mass-fractions>" << std::endl;
		fXML << mixture.thermodynamicsMap().NumberOfSpecies() << std::endl;
		for (unsigned int i = 0; i < mixture.thermodynamicsMap().NumberOfSpecies(); i++)
			fXML << mixture.thermodynamicsMap().NamesOfSpecies()[i] << " " << mixture.thermodynamicsMap().MW(i) << " " << n_additional + (i+1) << std::endl;
		fXML << "</mass-fractions>" << std::endl;

	
		const scalarField& TCells = T.internalField();
		const scalarField& pCells = p.internalField();
		const scalarField& rhoCells = rho.internalField();
		const scalarField& csiCells = csi.internalField();
		const scalarField& vCells = mesh.V();
		const vectorField& uCells = U.internalField();

		// Select active cells
		int ntotcells = 0;
		forAll(TCells, celli) { ntotcells++; }
		std::vector<bool> activecell(ntotcells);
		std::fill(activecell.begin(), activecell.end(), false);

		int ncells = 0;
		forAll(TCells, celli)
		{
			if ( TCells[celli] > 320.)
			{
				activecell[celli] = true;
				ncells++;
			}
		}


		fXML << "<profiles-additional>" << std::endl;
		forAll(TCells, celli)
		{
			if (activecell[celli] == true)
			{
				const double ctot = pCells[celli]/8314./TCells[celli];
				const double mwtot = rhoCells[celli]/ctot;

				fXML << 1.e2*mesh.C()[celli][0] << " ";
				fXML << TCells[celli] << " ";
				fXML << pCells[celli] << " ";
				fXML << mwtot << " ";
				fXML << rhoCells[celli] << " ";
				fXML << 0 << " ";
				fXML << uCells[celli][0] << " ";
				fXML << uCells[celli][1] << " ";
				fXML << uCells[celli][2] << " ";
				fXML << vCells[celli] << " ";
				fXML << csiCells[celli] << " ";
				fXML << 1.e2*mesh.C()[celli][1] << " ";
				fXML << 1.e2*mesh.C()[celli][2] << " ";

				fXML << std::endl;
			}
		}
		fXML << "</profiles-additional>" << std::endl;

		for (unsigned int j = 0; j < mixture.thermodynamicsMap().NumberOfSpecies(); j++)
		{
			std::stringstream label; label << j+1;
			fXML << "<profiles-" << label.str() << ">" << std::endl;
			forAll(TCells, celli)
			{
				int index = 0;
				if (activecell[celli] == true)
				{
					if (Y[j][celli] < 1e-32)	fXML << 0 << " ";
					else				fXML << Y[j][celli] << " ";
				
					index++;
					if (index%100 == 0)		fXML << std::endl;
				}
			}
			fXML << std::endl;
			fXML << "</profiles-" << label.str() << ">" << std::endl;
		}

		fXML << "<profiles-size>" << std::endl;
		fXML << ncells << " " << mixture.thermodynamicsMap().NumberOfSpecies() + n_additional << std::endl;
		fXML << "</profiles-size>" << std::endl;
		fXML << "</opensmoke>" << std::endl;
	
		fXML.close();
	}


	// Axis
	{
		Info<< "Exporting OpenSMOKE++ file ..." << endl;

		const unsigned int n_additional = 13;

		const std::string file_name = "output.axis.xml";

		std::ofstream fXML;
		fXML.open(file_name.c_str(), std::ios::out);
		fXML.setf(std::ios::scientific);

		fXML << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
		fXML << "<opensmoke version=\"0.1a\">" << std::endl;

		fXML << "<Type>Axis</Type>" << std::endl;

		fXML << "<additional>" << std::endl;
		fXML << n_additional << std::endl;
		fXML << "axial-coordinate [cm] 2" << std::endl;
		fXML << "temperature [K] 3" << std::endl;
		fXML << "pressure [Pa] 4" << std::endl;
		fXML << "mol-weight [kg/kmol] 5" << std::endl;
		fXML << "density [kg/m3] 6" << std::endl;
		fXML << "heat-release [W/m3] 7" << std::endl;
		fXML << "axial-velocity [m/s] 8" << std::endl;
		fXML << "dummy [-] 9" << std::endl;
		fXML << "dummy [-] 10" << std::endl;
		fXML << "dummy [-] 11" << std::endl;
		fXML << "csi [-] 12" << std::endl;
		fXML << "x [-] 13" << std::endl;
		fXML << "y [-] 14" << std::endl;
		fXML << "</additional>" << std::endl;

		fXML << "<t-p-mw>" << std::endl;
		fXML << "1 2 3" << std::endl;
		fXML << "</t-p-mw>" << std::endl;

		fXML << "<mass-fractions>" << std::endl;
		fXML << mixture.thermodynamicsMap().NumberOfSpecies() << std::endl;
		for (unsigned int i = 0; i < mixture.thermodynamicsMap().NumberOfSpecies(); i++)
			fXML << mixture.thermodynamicsMap().NamesOfSpecies()[i] << " " << mixture.thermodynamicsMap().MW(i) << " " << n_additional + (i+1) << std::endl;
		fXML << "</mass-fractions>" << std::endl;

	
		const scalarField& TCells = T.internalField();
		const scalarField& pCells = p.internalField();
		const scalarField& rhoCells = rho.internalField();
		const scalarField& csiCells = csi.internalField();
		const scalarField& vCells = mesh.V();
		const vectorField& uCells = U.internalField();



		int ntotcells = 0;
		double minx = 1.e16;
		forAll(TCells, celli)
		{
			if (mesh.C()[celli][0] < minx)
				minx = mesh.C()[celli][0];
			ntotcells++;
		}

		// Select active cells
		std::vector<bool> activecell(ntotcells);
		std::fill(activecell.begin(), activecell.end(), false);
		
		int ncells = 0;
		forAll(TCells, celli)
		{
			if ( std::fabs( (mesh.C()[celli][0]-minx)/minx )<0.05 )
			{
				activecell[celli] = true;
				ncells++;
			}
		}


		fXML << "<profiles-additional>" << std::endl;
		forAll(TCells, celli)
		{
			if (activecell[celli] == true)
			{
				const double ctot = pCells[celli]/8314./TCells[celli];
				const double mwtot = rhoCells[celli]/ctot;

				fXML << 1.e2*mesh.C()[celli][2] << " ";
				fXML << TCells[celli] << " ";
				fXML << pCells[celli] << " ";
				fXML << mwtot << " ";
				fXML << rhoCells[celli] << " ";
				fXML << 0 << " ";
				fXML << uCells[celli][2] << " ";
				fXML << 0 << " ";
				fXML << 0 << " ";
				fXML << 0 << " ";
				fXML << csiCells[celli] << " ";
				fXML << 1.e2*mesh.C()[celli][0] << " ";
				fXML << 1.e2*mesh.C()[celli][1] << " ";

				fXML << std::endl;
			}
		}
		fXML << "</profiles-additional>" << std::endl;

		for (unsigned int j = 0; j < mixture.thermodynamicsMap().NumberOfSpecies(); j++)
		{
			std::stringstream label; label << j+1;
			fXML << "<profiles-" << label.str() << ">" << std::endl;
			forAll(TCells, celli)
			{
				int index = 0;
				if (activecell[celli] == true)
				{
					if (Y[j][celli] < 1e-32)	fXML << 0 << " ";
					else				fXML << Y[j][celli] << " ";
				
					index++;
					if (index%100 == 0)		fXML << std::endl;
				}
			}
			fXML << std::endl;
			fXML << "</profiles-" << label.str() << ">" << std::endl;
		}

		fXML << "<profiles-size>" << std::endl;
		fXML << ncells << " " << mixture.thermodynamicsMap().NumberOfSpecies() + n_additional << std::endl;
		fXML << "</profiles-size>" << std::endl;
		fXML << "</opensmoke>" << std::endl;
	
		fXML.close();
	}

	// To tecplot
	{
		Info<< "Soot volume fraction\n" << endl;
		volScalarField fv
		(
		    IOobject
		    (
			"soot:fv",
			mesh.time().timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
			),
			mesh
		);

		Info<< "Soot number density\n" << endl;
		volScalarField N
		(
		    IOobject
		    (
			"soot:N",
			mesh.time().timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
			),
			mesh
		);

		Info<< "Mixture fraction\n" << endl;
		volScalarField csi
		(
		    IOobject
		    (
			"csi",
			mesh.time().timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
			),
			mesh
		);

		std::string scheme = "cellPointFace";

		autoPtr< interpolation<scalar> > csiinterp = interpolation<scalar>::New(scheme, csi);
		autoPtr< interpolation<scalar> > rhointerp = interpolation<scalar>::New(scheme, rho);
		autoPtr< interpolation<scalar> > Pinterp = interpolation<scalar>::New(scheme, p);
		autoPtr< interpolation<scalar> > Tinterp = interpolation<scalar>::New(scheme, T);
		autoPtr< interpolation<scalar> > C2H4interp = interpolation<scalar>::New(scheme, Y[mixture.thermodynamicsMap().IndexOfSpecies("C2H4")-1]);
		autoPtr< interpolation<scalar> > O2interp = interpolation<scalar>::New(scheme, Y[mixture.thermodynamicsMap().IndexOfSpecies("O2")-1]);
		autoPtr< interpolation<scalar> > H2Ointerp = interpolation<scalar>::New(scheme, Y[mixture.thermodynamicsMap().IndexOfSpecies("H2O")-1]);
		autoPtr< interpolation<scalar> > CO2interp = interpolation<scalar>::New(scheme, Y[mixture.thermodynamicsMap().IndexOfSpecies("CO2")-1]);
		autoPtr< interpolation<scalar> > COinterp = interpolation<scalar>::New(scheme, Y[mixture.thermodynamicsMap().IndexOfSpecies("CO")-1]);
		autoPtr< interpolation<scalar> > H2interp = interpolation<scalar>::New(scheme, Y[mixture.thermodynamicsMap().IndexOfSpecies("H2")-1]);
		autoPtr< interpolation<scalar> > OHinterp = interpolation<scalar>::New(scheme, Y[mixture.thermodynamicsMap().IndexOfSpecies("OH")-1]);
		autoPtr< interpolation<scalar> > C2H2interp = interpolation<scalar>::New(scheme, Y[mixture.thermodynamicsMap().IndexOfSpecies("C2H2")-1]);
		autoPtr< interpolation<scalar> > C6H6interp = interpolation<scalar>::New(scheme, Y[mixture.thermodynamicsMap().IndexOfSpecies("C6H6")-1]);
		autoPtr< interpolation<scalar> > C10H8interp = interpolation<scalar>::New(scheme, Y[mixture.thermodynamicsMap().IndexOfSpecies("C10H8")-1]);
		autoPtr< interpolation<scalar> > N2interp = interpolation<scalar>::New(scheme, Y[mixture.thermodynamicsMap().IndexOfSpecies("N2")-1]);
		autoPtr< interpolation<scalar> > fvinterp = interpolation<scalar>::New(scheme, fv);
		autoPtr< interpolation<scalar> > Ninterp = interpolation<scalar>::New(scheme, N);

		scalar fvmax = max(fv).value();

		const double rmax = 0.01;
		const double zmax = 0.09;		// F3
//		const double zmax = 0.10;		// LF2
		const unsigned int nr = 51;
		const unsigned int nz = 301;
		
		const double dr = rmax/double(nr-1);
		const double dz = zmax/double(nz-1);

		const std::string file_name = "output.tec";
		std::ofstream fTecplot;
		fTecplot.open(file_name.c_str(), std::ios::out);
		fTecplot.setf(std::ios::scientific);

		fTecplot << "VARIABLES = " << "\"X\"" << "," << "\"Y\"";
		fTecplot << ", \"csi\"";
		fTecplot << ", \"rho\"";
		fTecplot << ", \"MW\"";
		fTecplot << ", \"P\"";
		fTecplot << ", \"T\"";
		fTecplot << ", \"C2H4\"";
		fTecplot << ", \"O2\"";
		fTecplot << ", \"H2O\"";
		fTecplot << ", \"CO2\"";
		fTecplot << ", \"CO\"";
		fTecplot << ", \"H2\"";
		fTecplot << ", \"OH\"";
		fTecplot << ", \"C2H2\"";
		fTecplot << ", \"C6H6\"";
		fTecplot << ", \"C10H8\"";
		fTecplot << ", \"N2\"";
		fTecplot << ", \"fv\"";
		fTecplot << ", \"fvnorm\"";
		fTecplot << ", \"N\"";
		fTecplot << std::endl;
		
		fTecplot << "ZONE I=" << nz << ", J=" << nr << ", DATAPACKING=POINT" << std::endl;

		for (unsigned int i=1;i<=nr;i++)
		{
			const double r = (i-1)*dr;
			for (unsigned int j=1;j<=nz;j++)
			{
				const double z = (j-1)*dz;
				vector pos(r,0.,z);

				label cellI = mesh.findCell(pos);

				scalar csii  = csiinterp->interpolate(pos, cellI);
				scalar Pi    = Pinterp->interpolate(pos, cellI);
				scalar rhoi  = rhointerp->interpolate(pos, cellI);
				scalar Ti    = Tinterp->interpolate(pos, cellI);
				scalar mwi = rhoi/(Pi/8314./Ti);
				scalar C2H4i = C2H4interp->interpolate(pos, cellI);
				scalar O2i = O2interp->interpolate(pos, cellI);
				scalar H2Oi = H2Ointerp->interpolate(pos, cellI);
				scalar CO2i = CO2interp->interpolate(pos, cellI);
				scalar COi = COinterp->interpolate(pos, cellI);
				scalar H2i = H2interp->interpolate(pos, cellI);
				scalar OHi = OHinterp->interpolate(pos, cellI);
				scalar C2H2i = C2H2interp->interpolate(pos, cellI);
				scalar C6H6i = C6H6interp->interpolate(pos, cellI);
				scalar C10H8i = C10H8interp->interpolate(pos, cellI);
				scalar N2i = N2interp->interpolate(pos, cellI);
				scalar fvi   = fvinterp->interpolate(pos, cellI);
				scalar Ni   = Ninterp->interpolate(pos, cellI);

	
				fTecplot << r << " " << z << " ";
				fTecplot << csii << " ";
				fTecplot << rhoi << " ";
				fTecplot << mwi << " ";
				fTecplot << Pi << " ";
				fTecplot << Ti << " ";
				fTecplot << C2H4i << " ";
				fTecplot << O2i << " ";
				fTecplot << H2Oi << " ";
				fTecplot << CO2i << " ";
				fTecplot << COi << " ";
				fTecplot << H2i << " ";
				fTecplot << OHi << " ";
				fTecplot << C2H2i << " ";
				fTecplot << C6H6i << " ";
				fTecplot << C10H8i << " ";
				fTecplot << N2i << " ";
				fTecplot << fvi*1e6 << " ";
				fTecplot << fvi/fvmax << " ";
				fTecplot << Ni << " ";
				fTecplot << std::endl;
			}
		}

		fTecplot.close();
	}


	// Pathline max soot
	{
		Info<< "Soot volume fraction\n" << endl;
		volScalarField fv
		(
		    IOobject
		    (
			"soot:fv",
			mesh.time().timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
			),
			mesh
		);

		Info<< "Soot number density\n" << endl;
		volScalarField N
		(
		    IOobject
		    (
			"soot:N",
			mesh.time().timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
			),
			mesh
		);

		Info<< "Mixture fraction\n" << endl;
		volScalarField csi
		(
		    IOobject
		    (
			"csi",
			mesh.time().timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
			),
			mesh
		);

		std::string scheme = "cellPointFace";
		autoPtr< interpolation<scalar> > fvinterp = interpolation<scalar>::New(scheme, fv);

		const double rmax = 0.020;
		const double zmax = 0.120;		
		const unsigned int nr = 1001;
		const unsigned int nz =  301;
		
		const double dr = rmax/double(nr-1);
		const double dz = zmax/double(nz-1);


		std::vector<double> maxr(nz);	
		std::vector<double> maxfv(nz);		
		for (unsigned int j=1;j<=nz;j++)
		{
			const double z = (j-1)*dz;
			maxr[j-1] = 0.;
			maxfv[j-1] = -1.e16;
			for (unsigned int i=1;i<=nr;i++)
			{
				const double r = (i-1)*dr;
				vector pos(r,0.,z);

				label cellI = mesh.findCell(pos);
				scalar fvi  = fvinterp->interpolate(pos, cellI);
				if (fvi > maxfv[j-1])
				{
					maxfv[j-1] = fvi;
					maxr[j-1] = r;
				}
			}
		}

		
		{
			Info<< "Exporting OpenSMOKE++ file ..." << endl;

			const unsigned int n_additional = 13;

			const std::string file_name = "output.maxfv.xml";

			std::ofstream fXML;
			fXML.open(file_name.c_str(), std::ios::out);
			fXML.setf(std::ios::scientific);

			fXML << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
			fXML << "<opensmoke version=\"0.1a\">" << std::endl;

			fXML << "<Type>Axis</Type>" << std::endl;

			fXML << "<additional>" << std::endl;
			fXML << n_additional << std::endl;
			fXML << "axial-coordinate [cm] 2" << std::endl;
			fXML << "temperature [K] 3" << std::endl;
			fXML << "pressure [Pa] 4" << std::endl;
			fXML << "mol-weight [kg/kmol] 5" << std::endl;
			fXML << "density [kg/m3] 6" << std::endl;
			fXML << "heat-release [W/m3] 7" << std::endl;
			fXML << "axial-velocity [m/s] 8" << std::endl;
			fXML << "dummy [-] 9" << std::endl;
			fXML << "dummy [-] 10" << std::endl;
			fXML << "dummy [-] 11" << std::endl;
			fXML << "csi [-] 12" << std::endl;
			fXML << "x [-] 13" << std::endl;
			fXML << "dummy [-] 14" << std::endl;
			fXML << "</additional>" << std::endl;

			fXML << "<t-p-mw>" << std::endl;
			fXML << "1 2 3" << std::endl;
			fXML << "</t-p-mw>" << std::endl;

			fXML << "<mass-fractions>" << std::endl;
			fXML << mixture.thermodynamicsMap().NumberOfSpecies() << std::endl;
			for (unsigned int i = 0; i < mixture.thermodynamicsMap().NumberOfSpecies(); i++)
				fXML << mixture.thermodynamicsMap().NamesOfSpecies()[i] << " " << mixture.thermodynamicsMap().MW(i) << " " << n_additional + (i+1) << std::endl;
			fXML << "</mass-fractions>" << std::endl;

			autoPtr< interpolation<scalar> > csiinterp = interpolation<scalar>::New(scheme, csi);
			autoPtr< interpolation<scalar> > rhointerp = interpolation<scalar>::New(scheme, rho);
			autoPtr< interpolation<scalar> > Pinterp = interpolation<scalar>::New(scheme, p);
			autoPtr< interpolation<scalar> > Tinterp = interpolation<scalar>::New(scheme, T);
			autoPtr< interpolation<vector> > Uinterp = interpolation<vector>::New(scheme, U);

			fXML << "<profiles-additional>" << std::endl;
			for (unsigned int j=1;j<=nz;j++)
			{
				const double z = (j-1)*dz;
				const double r = maxr[j-1];
				vector pos(r,0.,z);

				label cellI = mesh.findCell(pos);

				vector Ui    = Uinterp->interpolate(pos, cellI);
				scalar Pi    = Pinterp->interpolate(pos, cellI);
				scalar rhoi  = rhointerp->interpolate(pos, cellI);
				scalar Ti    = Tinterp->interpolate(pos, cellI);
				scalar csii  = csiinterp->interpolate(pos, cellI);
				scalar mwi   = rhoi/(Pi/8314./Ti);

				fXML << 1.e2*z << " ";
				fXML << Ti << " ";
				fXML << Pi << " ";
				fXML << mwi << " ";
				fXML << rhoi << " ";
				fXML << 0 << " ";
				fXML << Ui[2] << " ";
				fXML << 0 << " ";
				fXML << 0 << " ";
				fXML << 0 << " ";
				fXML << csii << " ";
				fXML << r << " ";
				fXML << 0 << " ";

				fXML << std::endl;
			}
			fXML << "</profiles-additional>" << std::endl;

			for (unsigned int j = 0; j < mixture.thermodynamicsMap().NumberOfSpecies(); j++)
			{
				autoPtr< interpolation<scalar> > Yinterp = interpolation<scalar>::New(scheme, Y[j]);

				std::stringstream label; label << j+1;
				fXML << "<profiles-" << label.str() << ">" << std::endl;
				for (unsigned int jj=1;jj<=nz;jj++)
				{
					const double z = (jj-1)*dz;
					const double r = maxr[jj-1];
					vector pos(r,0.,z);
					const int cellI = mesh.findCell(pos);
					scalar Yi    = Yinterp->interpolate(pos,  cellI);
					if (Yi < 1e-32)		fXML << 0  << " ";
					else			fXML << Yi << " ";
				
					if (jj%100 == 0)	fXML << std::endl;
				}
				fXML << std::endl;
				fXML << "</profiles-" << label.str() << ">" << std::endl;
			}

			fXML << "<profiles-size>" << std::endl;
			fXML << nz << " " << mixture.thermodynamicsMap().NumberOfSpecies() + n_additional << std::endl;
			fXML << "</profiles-size>" << std::endl;
			fXML << "</opensmoke>" << std::endl;
	
			fXML.close();
		}
	}


/*

	// Along the axis
	{
		Info<< "Exporting OpenSMOKE++ file ..." << endl;

		const unsigned int n_additional = 13;

		const std::string file_name = "output.axis.xml";

		std::ofstream fXML;
		fXML.open(file_name.c_str(), std::ios::out);
		fXML.setf(std::ios::scientific);

		fXML << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
		fXML << "<opensmoke version=\"0.1a\">" << std::endl;

		fXML << "<Type>Axis</Type>" << std::endl;

		fXML << "<additional>" << std::endl;
		fXML << n_additional << std::endl;
		fXML << "axial-coordinate [cm] 2" << std::endl;
		fXML << "temperature [K] 3" << std::endl;
		fXML << "pressure [Pa] 4" << std::endl;
		fXML << "mol-weight [kg/kmol] 5" << std::endl;
		fXML << "density [kg/m3] 6" << std::endl;
		fXML << "heat-release [W/m3] 7" << std::endl;
		fXML << "axial-velocity [m/s] 8" << std::endl;
		fXML << "dummy [-] 9" << std::endl;
		fXML << "dummy [-] 10" << std::endl;
		fXML << "dummy [-] 11" << std::endl;
		fXML << "csi [-] 12" << std::endl;
		fXML << "dummy [-] 13" << std::endl;
		fXML << "dummy [-] 14" << std::endl;
		fXML << "</additional>" << std::endl;

		fXML << "<t-p-mw>" << std::endl;
		fXML << "1 2 3" << std::endl;
		fXML << "</t-p-mw>" << std::endl;

		fXML << "<mass-fractions>" << std::endl;
		fXML << mixture.thermodynamicsMap().NumberOfSpecies() << std::endl;
		for (unsigned int i = 0; i < mixture.thermodynamicsMap().NumberOfSpecies(); i++)
			fXML << mixture.thermodynamicsMap().NamesOfSpecies()[i] << " " << mixture.thermodynamicsMap().MW(i) << " " << n_additional + (i+1) << std::endl;
		fXML << "</mass-fractions>" << std::endl;


		label patchi = mesh.boundaryMesh().findPatchID("axis");
	
		volScalarField::Boundary& pBf 		= p.boundaryFieldRef();
		volScalarField::Boundary& TBf 		= T.boundaryFieldRef();
		volScalarField::Boundary& rhoBf 	= rho.boundaryFieldRef();
		volScalarField::Boundary& csiBf 	= csi.boundaryFieldRef();
		volVectorField::Boundary& uBf 		= U.boundaryFieldRef();

		fvPatchScalarField& pp = pBf[patchi];
		fvPatchScalarField& pT = TBf[patchi];
		fvPatchScalarField& prho = rhoBf[patchi];
		fvPatchScalarField& pcsi = csiBf[patchi];
		fvPatchVectorField& pu = uBf[patchi];

		// Select active cells
		int ntotfaces = 0;
		forAll(pT, facei) { ntotfaces++; }
		std::vector<bool> activeface(ntotfaces);
		std::fill(activeface.begin(), activeface.end(), false);

		int nfaces = 0;
		forAll(pT, facei)
		{
			if ( pT[facei] > 0.)
			{
				activeface[facei] = true;
				nfaces++;
			}
		}

		fXML << "<profiles-additional>" << std::endl;
		forAll(pT, facei)
		{
			if (activeface[facei] == true)
			{
				const double ctot = pp[facei]/8314./pT[facei];
				const double mwtot = prho[facei]/ctot;

				fXML << 1.e2*mesh.boundaryMesh()[patchi].faceCentres()[facei].z() << " ";
				fXML << pT[facei] << " ";
				fXML << pp[facei] << " ";
				fXML << mwtot << " ";
				fXML << prho[facei] << " ";
				fXML << 0 << " ";
				fXML << pu[facei][2] << " ";
				fXML << 0 << " ";
				fXML << 0 << " ";
				fXML << 0 << " ";
				fXML << pcsi[facei] << " ";
				fXML << 0 << " ";
				fXML << 0 << " ";

				fXML << std::endl;
			}
		}
		fXML << "</profiles-additional>" << std::endl;

		for (unsigned int j = 0; j < mixture.thermodynamicsMap().NumberOfSpecies(); j++)
		{
			std::stringstream label; label << j+1;
			fXML << "<profiles-" << label.str() << ">" << std::endl;
			forAll(pT, facei)
			{
				int index = 0;
				if (activeface[facei] == true)
				{
					if (Y[j].boundaryField()[patchi][facei] < 1e-24)	fXML << 0 << " ";
					else							fXML << Y[j].boundaryField()[patchi][facei] << " ";
				
					index++;
					if (index%100 == 0)					fXML << std::endl;
				}
			}
			fXML << std::endl;
			fXML << "</profiles-" << label.str() << ">" << std::endl;
		}

		fXML << "<profiles-size>" << std::endl;
		fXML << nfaces << " " << mixture.thermodynamicsMap().NumberOfSpecies() + n_additional << std::endl;
		fXML << "</profiles-size>" << std::endl;
		fXML << "</opensmoke>" << std::endl;
	
		fXML.close();
	}

*/
}
