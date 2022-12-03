if (export_outlet == true)
{
	Info << "Reconstructing outlet..." << endl;
        {
                std::vector<double> mass_fractions(mixture.thermodynamicsMap().NumberOfSpecies());

                const scalarField& TCells = T.internalField();

		double YO2in_massflowrate = 0.;
		double YO2in_mass_avg = 0.;
                double YO2out_massflowrate = 0.;
                double YO2out_mass_avg = 0.;

                double XO2in_massflowrate = 0.;
                double XO2in_mass_avg = 0.;
                double XO2out_massflowrate = 0.;
                double XO2out_mass_avg = 0.;
		double massFlowRate_inlet = 0.;

		const unsigned int iO2 = mixture.thermodynamicsMap().IndexOfSpecies("O2")-1;
		
		forAll(T.boundaryField(), patchi)
                {
			if ((T.boundaryField()[patchi].patch().name()).substr(0,5) == "inlet")
                        {
                                const fvPatchScalarField& pT = T.boundaryField()[patchi];

				massFlowRate_inlet += -sum(phi.boundaryField()[patchi]);
                                const scalar area = gSum(mesh.magSf().boundaryField()[patchi]);
				
				YO2in_massflowrate += -sum(phi.boundaryField()[patchi]*Y[iO2].boundaryField()[patchi]);
				XO2in_massflowrate += -sum(phi.boundaryField()[patchi]*X[iO2].boundaryField()[patchi]);
			}

			YO2in_mass_avg = YO2in_massflowrate / massFlowRate_inlet;
                        XO2in_mass_avg = XO2in_massflowrate / massFlowRate_inlet;
		}

		forAll(T.boundaryField(), patchi)
                {
			if (T.boundaryField()[patchi].patch().name() == "outlet")
			{
				const fvPatchScalarField& pT = T.boundaryField()[patchi];
				const scalar massFlowRate = sum(phi.boundaryField()[patchi]);
				const scalar area = gSum(mesh.magSf().boundaryField()[patchi]);

				const scalar rho_mass_avg = sum(phi.boundaryField()[patchi]*rho.boundaryField()[patchi]) / massFlowRate;
				const scalar rho_area_avg = sum(mesh.magSf().boundaryField()[patchi]*rho.boundaryField()[patchi]) / area;

                                const scalar MW_mass_avg = sum(phi.boundaryField()[patchi]*MW.boundaryField()[patchi]) / massFlowRate;
                                const scalar MW_area_avg = sum(mesh.magSf().boundaryField()[patchi]*MW.boundaryField()[patchi]) / area;

                                const scalar U_mass_avg = massFlowRate/rho_mass_avg/area;
                                const scalar U_area_avg = massFlowRate/rho_area_avg/area;

                                const scalar T_mass_avg = sum(phi.boundaryField()[patchi]*T.boundaryField()[patchi]) / massFlowRate;
                                const scalar T_area_avg = sum(mesh.magSf().boundaryField()[patchi]*T.boundaryField()[patchi]) / area;

                                YO2out_massflowrate = sum(phi.boundaryField()[patchi]*Y[iO2].boundaryField()[patchi]);
                                YO2out_mass_avg = YO2out_massflowrate / massFlowRate;

                                XO2out_massflowrate = sum(phi.boundaryField()[patchi]*X[iO2].boundaryField()[patchi]);
                                XO2out_mass_avg = XO2out_massflowrate / massFlowRate;

				std::vector<double> Y_mass_avg(Y.size());
				std::vector<double> Y_area_avg(Y.size());
                                std::vector<double> X_mass_avg(Y.size());
                                std::vector<double> X_area_avg(Y.size());
				
				double Y_mass_dry = 1.;
				double X_mass_dry = 1.;
				for (unsigned int i=0;i<Y.size();i++)
                                {
                                        Y_mass_avg[i] = sum(phi.boundaryField()[patchi]*Y[i].boundaryField()[patchi]) / massFlowRate;
                                        Y_area_avg[i] = sum(mesh.magSf().boundaryField()[patchi]*Y[i].boundaryField()[patchi]) / area;
                                        X_mass_avg[i] = sum(phi.boundaryField()[patchi]*X[i].boundaryField()[patchi]) / massFlowRate;
                                        X_area_avg[i] = sum(mesh.magSf().boundaryField()[patchi]*X[i].boundaryField()[patchi]) / area;					
                                	if (mixture.thermodynamicsMap().NamesOfSpecies()[i] == "H2O" || mixture.thermodynamicsMap().NamesOfSpecies()[i] == "h2o" )
					{
						Y_mass_dry = 1.-Y_mass_avg[i];
						X_mass_dry = 1.-X_mass_avg[i];
					}
				//        Info << mixture.thermodynamicsMap().NamesOfSpecies()[i] << "\t" << Y_area_avg << "\t" << Y_mass_avg << endl;
                                }

				double YNOX_mass_avg = 0.;
				double YNOX_area_avg = 0.;
				double XNOX_mass_avg = 0.;
				double XNOX_area_avg = 0.;
				if (mixture.noxPostProcessor() == true)
				{
					YNOX_mass_avg = sum(phi.boundaryField()[patchi]*mixture.noxPostProcessorModel().YNO().boundaryField()[patchi]) / massFlowRate;
                                        YNOX_area_avg = sum(mesh.magSf().boundaryField()[patchi]*mixture.noxPostProcessorModel().YNO().boundaryField()[patchi]) / area;
                                        XNOX_mass_avg = (YNOX_mass_avg*MW_mass_avg/28.);
                                        XNOX_area_avg = (YNOX_area_avg*MW_area_avg/28.);
				}

				Info << "Outlet mass flow rate:       " << massFlowRate << " kg/s" << endl;
				Info << "Outlet volumetric flow rate: " << massFlowRate/rho_mass_avg << " m3/s" << endl; 
				Info << "Outlet surface area:         " << area*1.e6 << " mm2" << endl;

				Info << "Average values (mass-flow-rate)                                " << endl;
				Info << "---------------------------------------------------------------------" << endl;
				Info << "Density (kg/m3):            " << rho_mass_avg << endl;
				Info << "Molecular weight (kg/kmol): " << MW_mass_avg << endl;
				Info << "Velocity (m/s):             " << U_mass_avg << endl;
				Info << "Temperature (C):            " << T_mass_avg-273.15 << endl;
				Info << endl;
				
				std::cout << "Outlet (wet):     Y                     X" << std::endl;
				for (unsigned int i=0;i<Y.size();i++)
				{
					std::cout << std::setw(18) << std::left << mixture.thermodynamicsMap().NamesOfSpecies()[i];
					std::cout << std::setw(21) << std::right << std::scientific << Y_mass_avg[i];
					std::cout << std::setw(21) << std::right << std::scientific << X_mass_avg[i];
					std::cout << std::endl;
				}
				if (mixture.noxPostProcessor() == true)
				{
                                        std::cout << std::setw(18) << std::left << "NOX";
                                        std::cout << std::setw(21) << std::right << std::scientific << YNOX_mass_avg;
                                        std::cout << std::setw(21) << std::right << std::scientific << XNOX_mass_avg;
                                        std::cout << std::endl;
				}
				std::cout << std::endl;

                                std::cout << "Outlet (dry):     Y                     X" << std::endl;
                                for (unsigned int i=0;i<Y.size();i++)
                                {
                                        std::cout << std::setw(18) << std::left << mixture.thermodynamicsMap().NamesOfSpecies()[i];
                                        std::cout << std::setw(21) << std::right << std::scientific << Y_mass_avg[i]/Y_mass_dry;
                                        std::cout << std::setw(21) << std::right << std::scientific << X_mass_avg[i]/X_mass_dry;
                                        std::cout << std::endl;;
                                }
                                if (mixture.noxPostProcessor() == true)
                                {
                                        std::cout << std::setw(18) << std::left << "NOX";
                                        std::cout << std::setw(21) << std::right << std::scientific << YNOX_mass_avg/Y_mass_dry;
                                        std::cout << std::setw(21) << std::right << std::scientific << XNOX_mass_avg/X_mass_dry;
                                        std::cout << std::endl;
                                }                 
		                std::cout << std::endl;

                                Info << "---------------------------------------------------------------------" << endl;
                        }
		}
		
		std::cout << "In/Out comparison" << std::endl;
			
		std::cout << std::left << "O2 Mass fraction:";
                                std::cout << std::setw(21) << std::right << YO2in_mass_avg;
                                std::cout << std::setw(21) << std::right << YO2out_mass_avg;
                                std::cout << std::endl;

                                std::cout << std::left << "O2 Mole fraction:";
                                std::cout << std::setw(21) << std::right << XO2in_mass_avg;
                                std::cout << std::setw(21) << std::right << XO2out_mass_avg;
                                std::cout << std::endl;			
	
				std::cout << std::left << "O2 Mass flow rate (kg/s):";
				std::cout << std::setw(21) << std::right << YO2in_massflowrate;
                                std::cout << std::setw(21) << std::right << YO2out_massflowrate;
				std::cout << std::setw(21) << std::right << (YO2in_massflowrate-YO2out_massflowrate)/YO2in_massflowrate;
				std::cout << std::endl;
	
				Info << "---------------------------------------------------------------------" << endl;	
		
		
	}				
}
