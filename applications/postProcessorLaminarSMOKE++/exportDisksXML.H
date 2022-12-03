if (export_disksxml == true)
{
	const dictionary& disksDictionary = postProcessingDictionary.subDict("DisksXML");

	const Foam::string north_name_u = disksDictionary.lookup("north");
	const Foam::string south_name_u = disksDictionary.lookup("south");
	const Foam::string east_name_u = disksDictionary.lookup("east");
	const Foam::string west_name_u = disksDictionary.lookup("west");

	const scalar Ri   = readScalar(disksDictionary.lookup("Ri")) / 1000.;
	const scalar Re   = readScalar(disksDictionary.lookup("Re")) / 1000.;
	const scalar Htot = readScalar(disksDictionary.lookup("H")) / 1000.;


	bool disk_exist = false;
	label patchID[4];
		
	patchID[0] = mesh.boundaryMesh().findPatchID(north_name_u);
	patchID[1] = mesh.boundaryMesh().findPatchID(south_name_u);
	patchID[2] = mesh.boundaryMesh().findPatchID(west_name_u);
	patchID[3] = mesh.boundaryMesh().findPatchID(east_name_u);
	
	if (patchID[0]!=-1 && patchID[1]!=-1 && patchID[3]!=-1)
		disk_exist = true;

	if (disk_exist == true)
	{
		label ns = Y.size();

		std::string filename = "Disk.xml";
		std::ofstream fDisk(filename.c_str(), std::ios::out);
		fDisk.setf(std::ios::scientific);

		fDisk << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
		fDisk << "<opensmoke version=\"0.1a\">" << std::endl;

		fDisk << "<AxialCoordinate>" << std::endl;
		fDisk << "z" << std::endl;
		fDisk << "</AxialCoordinate>" << std::endl;
		fDisk << "<RadialCoordinate>" << std::endl;
		fDisk << "x" << std::endl;
		fDisk << "</RadialCoordinate>" << std::endl;

		fDisk << "<NumberOfSpecies>" << std::endl;
		fDisk << ns << std::endl;
		fDisk << "</NumberOfSpecies>" << std::endl;
		fDisk << "<NamesOfSpecies>" << std::endl;
		for(unsigned int i=0;i<ns;i++)
			fDisk << mixture.thermodynamicsMap().NamesOfSpecies()[i] << std::endl;
		fDisk << "</NamesOfSpecies>" << std::endl;	

		
		for(int k=0;k<4;k++)
		{
			label patchi = patchID[k];

			std::string side_string = "";
			if (k==0) side_string = "North";
			if (k==1) side_string = "South";
			if (k==2) side_string = "West";
			if (k==3) side_string = "East";

			if (patchi != -1)
			{
				fDisk << "<Data side=\"" << side_string << "\">" << std::endl;				

				const int number_points = 10;
				
				fDisk << "<NumberOfPoints>" << std::endl;
				fDisk << number_points << std::endl;
				fDisk << "</NumberOfPoints>" << std::endl;

				fDisk << "<Points>" << std::endl;

				if (side_string == "North")				
				{
					const double delta = (Re-Ri)/static_cast<double>(number_points);
					for (unsigned int i=0;i<number_points;i++)
						fDisk 	<< (Ri+delta/2.) + delta*i 	<< " " 
							<< 0. 				<< " "
							<< Htot/2. 			<< " "
							<< std::endl;
				}
				else if (side_string == "South")				
				{
					const double delta = (Re-Ri)/static_cast<double>(number_points);
					for (unsigned int i=0;i<number_points;i++)
						fDisk 	<< (Ri+delta/2.) + delta*i 	<< " " 
							<< 0. 				<< " "
							<< -Htot/2. 			<< " "
							<< std::endl;
				}

				else if (side_string == "East")				
				{
					const double delta = Htot/static_cast<double>(number_points);
					for (unsigned int i=0;i<number_points;i++)
						fDisk 	<< Re 					<< " " 
							<< 0. 					<< " "
							<< (-Htot/2.+delta/2.) + delta*i 	<< " "
							<< std::endl;
				}

				else if (side_string == "West")				
				{
					const double delta = Htot/static_cast<double>(number_points);
					for (unsigned int i=0;i<number_points;i++)
						fDisk 	<< Ri 					<< " " 
							<< 0. 					<< " "
							<< (-Htot/2.+delta/2.) + delta*i 	<< " "
							<< std::endl;
				}

				fDisk << "</Points>" << std::endl;

				Eigen::VectorXd omegaAvg(ns); omegaAvg.setZero();
				double TAvg = 0.;
				double areaTot = 0.;
				forAll(mesh.boundaryMesh()[patchi].faceCentres(), facei)
				{
					const double area = mesh.magSf().boundaryField()[patchi][facei]; 
					areaTot += area;
						
					TAvg += T.boundaryField()[patchi][facei] * area;				
					for(unsigned int i=0;i<ns;i++)
						omegaAvg(i) += Y[i].boundaryField()[patchi][facei] * area;	
				}
				omegaAvg /= areaTot;
				TAvg /= areaTot;

				fDisk << "<MassFractions>" << std::endl;
				for (unsigned int i=0;i<number_points;i++)
				{
					for(unsigned int i=0;i<ns;i++)
						fDisk << omegaAvg(i)/omegaAvg.sum() << " ";
					fDisk << std::endl;
				}
				fDisk << "</MassFractions>" << std::endl;

				fDisk << "<Temperature>" << std::endl;
				for (unsigned int i=0;i<number_points;i++)
					fDisk << TAvg << std::endl;
				fDisk << "</Temperature>" << std::endl;

				fDisk << "</Data>" << std::endl;
			}

			else
			{
				fDisk << "<Data side=\"" << side_string << "\">" << std::endl;				

				int number_points = 0;			
				fDisk << "<NumberOfPoints>" << std::endl;
				fDisk << number_points << std::endl;
				fDisk << "</NumberOfPoints>" << std::endl;
				fDisk << "</Data>" << std::endl;
			}
		}

		fDisk << "</opensmoke>" << std::endl;
		fDisk.close();

	}

}


