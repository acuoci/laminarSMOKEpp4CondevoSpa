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

void Foam::OpenSMOKEppReactingMixture::chemistry_direct_integration(const double t0, const double tf, const fvMesh& mesh, const Foam::volScalarField& rho)
{
	if (hmom_ == false && drg_ == false && clustering_ == false && noxPostProcessor_ == false)
	{
		//- Initial conditions
		const double deltat = tf - t0;
		scalarField& TCells = T_.ref();
		scalarField& pCells = p_.ref();
		scalarField& QCells = Q_.ref();
		scalarField& cpuChemistryCells = cpuChemistry_.ref();

		const scalarField& rhoCells = rho.internalField();
		const scalarField& cpCells = Cp_.internalField();
		const scalarField& vCells = mesh.V();

		// Umbalances
		double maxUmbalance = 0.;
		double sumUmbalance = 0.;

		// Integration via OpenSMOKE++ solver
		if (homogeneousReactions_ == true && odeParameterBatchReactorHomogeneous_.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
		{
			const unsigned int neqs = ns_+1;
		
			// Min and max values
			Eigen::VectorXd yMin(neqs); yMin.setConstant(0.);  yMin(ns_) = Tmin_;
			Eigen::VectorXd yMax(neqs); yMax.setConstant(1.);  yMax(ns_) = Tmax_;
			Eigen::VectorXd y0(neqs);
			Eigen::VectorXd yf(neqs);

			double cpuTotal = 0.;

			Info <<" * Solving homogeneous chemistry (OpenSMOKE++ solver)... "<<endl;
			{			
				unsigned int counter = 0;
			
				double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

				forAll(TCells, celli)
				{
					double tStartLocal = OpenSMOKE::OpenSMOKEGetCpuTime();

					//- Solving for celli:	
					if (TCells[celli] > direct_integration_minimum_temperature_for_chemistry_)
					{
						for(unsigned int i=0;i<ns_;i++)
							y0(i) = Y_[i].internalField()[celli];
						y0(ns_) = TCells[celli];

						// Check and normalize the composition
						{
							double sum = 0.;
							for(unsigned int i=0;i<ns_;i++)
							{
								if (y0(i) < 0.)	y0(i) = 0.;
								sum += y0(i);
							}
							for(unsigned int i=0;i<ns_;i++)
								y0(i) /= sum;
						}

						//if (constPressureBatchReactor == true)
						{
							// Set reactor
							batchReactorHomogeneousConstantPressure_->SetReactor(pCells[celli]);
							batchReactorHomogeneousConstantPressure_->SetEnergyEquation(solveForEnergyEquation_);

							// Set initial conditions
							odeSolverConstantPressure_->SetInitialConditions(t0, y0);

							// Additional ODE solver options
							if (celli == 0)
							{
								// Set linear algebra options
								odeSolverConstantPressure_->SetLinearAlgebraSolver(odeParameterBatchReactorHomogeneous_.linear_algebra());
								odeSolverConstantPressure_->SetFullPivoting(odeParameterBatchReactorHomogeneous_.full_pivoting());

								// Set relative and absolute toleransces
								odeSolverConstantPressure_->SetAbsoluteTolerances(odeParameterBatchReactorHomogeneous_.absolute_tolerance());
								odeSolverConstantPressure_->SetRelativeTolerances(odeParameterBatchReactorHomogeneous_.relative_tolerance());

								// Set minimum and maximum values
								odeSolverConstantPressure_->SetMinimumValues(yMin);
								odeSolverConstantPressure_->SetMaximumValues(yMax);
							}

							// Solve
							OdeSMOKE::OdeStatus status = odeSolverConstantPressure_->Solve(tf);

							// Transfer solution
							odeSolverConstantPressure_->Solution(yf);
						
							// Fill heat release field
							QCells[celli] = batchReactorHomogeneousConstantPressure_->QR();
						}
					}		
					else
					{
						for(unsigned int i=0;i<ns_;i++)
							yf(i) = Y_[i].internalField()[celli];
						yf(ns_) = TCells[celli];
					}

					// Check mass fractions and umbalance
					const double umbalance = normalizeMassFractions(yf, celli);
					maxUmbalance = std::max(umbalance, maxUmbalance);
					sumUmbalance += umbalance;

				//	if (strangAlgorithm != STRANG_COMPACT)
					{
						// Assign mass fractions
						for(int i=0;i<ns_;i++)
							Y_[i].ref()[celli] = yf(i);

						//- Allocating final values: temperature
						if (solveForEnergyEquation_ == true)
							TCells[celli] = yf(ns_);
					}
				
					double tEndLocal = OpenSMOKE::OpenSMOKEGetCpuTime();
					cpuChemistryCells[celli] = (tEndLocal-tStartLocal)*1000.;

					cpuTotal += (tEndLocal-tStartLocal);

					if (counter%(int(0.20*mesh.nCells())+1) == 0)
						Info <<"   Accomplished: " << counter << "/" << mesh.nCells() << " Total(ms): " << cpuTotal*1000 << endl;

					counter++;				
				}

				double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			
				Info << "   Max umbalance: " << maxUmbalance << " - Mean umbalance: " << sumUmbalance/double(mesh.nCells()) << endl;
				Info << "   Homogeneous chemistry solved in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per reactor)" << endl;
			}
		}
	
	    	Info<< " * T gas min/max (after chemistry) = " << min(T_).value() << ", " << max(T_).value() << endl;

	}

	else if (hmom_ == false && drg_ == true && clustering_ == false && noxPostProcessor_ == false)
	{
		//- Initial conditions
		const double deltat = tf - t0;
		scalarField& TCells = T_.ref();
		scalarField& pCells = p_.ref();
		scalarField& QCells = Q_.ref();
		scalarField& cpuChemistryCells = cpuChemistry_.ref();

		const scalarField& rhoCells = rho.internalField();
		const scalarField& cpCells = Cp_.internalField();
		const scalarField& vCells = mesh.V();

		scalarField& drgSpeciesCells = drgModel_().important_species().ref();
		scalarField& drgReactionsCells = drgModel_().important_reactions().ref();

		// Umbalances
		double maxUmbalance = 0.;
		double sumUmbalance = 0.;

		// Integration via OpenSMOKE++ solver
		if (homogeneousReactions_ == true && odeParameterBatchReactorHomogeneous_.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
		{
			// Min and max values
			Eigen::VectorXd yMin;
			Eigen::VectorXd yMax;
	
			// Initial and final values
			Eigen::VectorXd y0;
			Eigen::VectorXd yf;

			// Mass fractions, mole fractions, and concentrations
			Eigen::VectorXd omega(ns_);
			Eigen::VectorXd x(ns_);
			Eigen::VectorXd c(ns_);

			double cpuDRG = 0.;
			double cpuSolution = 0.;
			double cpuTotal = 0.;
			int speciesDRG = 0;
			int reactionsDRG = 0;

			Info <<" * Solving homogeneous chemistry (OpenSMOKE++ solver) + DRG... "<<endl;
			{			
				unsigned int counter = 0;
			
				double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

				forAll(TCells, celli)
				{
					double tStartLocal = OpenSMOKE::OpenSMOKEGetCpuTime();

					//- Solving for celli:	
					if (TCells[celli] > direct_integration_minimum_temperature_for_chemistry_)
					{
						// DRG Analysis
						double tStartDRG = OpenSMOKE::OpenSMOKEGetCpuTime();
						{
							// Mass fractions
							for(unsigned int i=0;i<ns_;i++)
								omega(i) = Y_[i].internalField()[celli];

							// Concentrations
							double mw;
							thermodynamicsMap_->MoleFractions_From_MassFractions(x.data(), mw, omega.data());
							double cTot = pCells[celli]/PhysicalConstants::R_J_kmol/TCells[celli];
							c = cTot*x;

							// Analysis
							drgModel_().Analysis(TCells[celli], pCells[celli], c);

							// Fill the relevant DRG fields
							drgSpeciesCells[celli]   = drgModel_().number_important_species();
							drgReactionsCells[celli] = drgModel_().number_important_reactions();
						}
						double tEndDRG = OpenSMOKE::OpenSMOKEGetCpuTime();

						// Recover the relevant variables
						double tStartSolution = OpenSMOKE::OpenSMOKEGetCpuTime();
						{
							const unsigned int neqs = drgModel_().number_important_species()+1;

							y0.resize(neqs);
							yf.resize(neqs);
							yMin.resize(neqs);	yMin.setConstant(0.); yMin(neqs-1) = Tmin_;
							yMax.resize(neqs);	yMax.setConstant(1.); yMax(neqs-1) = Tmax_;

							for (unsigned int i=0;i<drgModel_().number_important_species();++i)	
							{
								const unsigned int j = drgModel_().indices_important_species()[i];
								y0(i) = omega(j);
							}
							y0(neqs-1) = TCells[celli];
						}

						//if (constPressureBatchReactor == true)
						{
							// Set reactor
							batchReactorHomogeneousConstantPressure_DRG_->SetReactor(pCells[celli]);
							batchReactorHomogeneousConstantPressure_DRG_->SetEnergyEquation(solveForEnergyEquation_);

							// Set DRG
							batchReactorHomogeneousConstantPressure_DRG_->SetDRGOptions( drgModel_().indices_important_species(), drgModel_().indices_unimportant_reactions(), drgModel_().important_reactions_steps() );
							batchReactorHomogeneousConstantPressure_DRG_->SetMassFractions(omega);

							// Set initial conditions
							odeSolverConstantPressure_DRG_->SetInitialConditions(t0, y0);

							// Additional ODE solver options
							{
								// Set linear algebra options
								odeSolverConstantPressure_DRG_->SetLinearAlgebraSolver(odeParameterBatchReactorHomogeneous_.linear_algebra());
								odeSolverConstantPressure_DRG_->SetFullPivoting(odeParameterBatchReactorHomogeneous_.full_pivoting());

								// Set relative and absolute toleransces
								odeSolverConstantPressure_DRG_->SetAbsoluteTolerances(odeParameterBatchReactorHomogeneous_.absolute_tolerance());
								odeSolverConstantPressure_DRG_->SetRelativeTolerances(odeParameterBatchReactorHomogeneous_.relative_tolerance());

								// Set minimum and maximum values
								odeSolverConstantPressure_DRG_->SetMinimumValues(yMin);
								odeSolverConstantPressure_DRG_->SetMaximumValues(yMax);
							}

							// Solve
							OdeSMOKE::OdeStatus status = odeSolverConstantPressure_DRG_->Solve(tf);

							// Transfer solution
							odeSolverConstantPressure_DRG_->Solution(yf);
						
							// Fill heat release field
							QCells[celli] = batchReactorHomogeneousConstantPressure_DRG_->QR();
						}

						// Distribute the solution
						{
							// Check the solution (no side effects on the OF fields)
							{
								Eigen::VectorXd yff(ns_+1);
								yff.head(ns_) = omega;
								yff(ns_) = yf(yf.size()-1);

								// Check mass fractions and umbalance
								const double umbalance = normalizeMassFractions(yff, celli);
								maxUmbalance = std::max(umbalance, maxUmbalance);
								sumUmbalance += umbalance;
							}
						
							for (unsigned int i=0;i<drgModel_().number_important_species();++i)	
							{
								const unsigned int j = drgModel_().indices_important_species()[i];
								Y_[j].ref()[celli] = yf(i);
							}

							//- Allocating final values: temperature
							if (solveForEnergyEquation_ == true)
								TCells[celli] = yf(yf.size()-1);
						}

						double tEndSolution = OpenSMOKE::OpenSMOKEGetCpuTime();

						cpuDRG += (tEndDRG - tStartDRG);
						cpuSolution += (tEndSolution - tStartSolution);
					}		
					else
					{
						// Fill the relevant DRG fields
						drgSpeciesCells[celli]   = drgModel_().number_key_species();
						drgReactionsCells[celli] = 0;
					}

					speciesDRG += drgSpeciesCells[celli];
					reactionsDRG += drgReactionsCells[celli];

					double tEndLocal = OpenSMOKE::OpenSMOKEGetCpuTime();
					cpuChemistryCells[celli] = (tEndLocal-tStartLocal)*1000.;
	
					cpuTotal += (tEndLocal-tStartLocal);

					if (counter%(int(0.20*mesh.nCells())+1) == 0)
						Info 	<< "   Accomplished: " << counter << "/" << mesh.nCells() 
							<< " ns=" << speciesDRG/double(mesh.nCells()) << " nr=" << reactionsDRG/double(mesh.nCells()) 
							<< " DRG(%)=" << cpuDRG/cpuTotal*100. << " Sol(%)=" << cpuSolution/cpuTotal*100. << " Tot(ms)=" << cpuTotal*1000 
							<< endl;

					counter++;
				}

				double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			
				Info << "   Max umbalance: " << maxUmbalance << " - Mean umbalance: " << sumUmbalance/double(mesh.nCells()) << endl;
				Info << "   Homogeneous chemistry solved in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per reactor)" << endl;
			}
		}
	
	    	Info<< " * T gas min/max (after chemistry) = " << min(T_).value() << ", " << max(T_).value() << endl;

	}

	else if (hmom_ == true && drg_ == false && clustering_ == false && noxPostProcessor_ == false)
	{
		//- Initial conditions
		const double deltat = tf - t0;
		scalarField& TCells = T_.ref();
		scalarField& pCells = p_.ref();
		scalarField& QCells = Q_.ref();
		scalarField& cpuChemistryCells = cpuChemistry_.ref();

		const scalarField& rhoCells = rho.internalField();
		const scalarField& cpCells = Cp_.internalField();
		const scalarField& vCells = mesh.V();

		const scalarField& fvCells = hmomModel_().fv().internalField(); 
		scalarField& Norm_M00Cells = hmomModel_().Norm_M00().ref();
		scalarField& Norm_M10Cells = hmomModel_().Norm_M10().ref();
		scalarField& Norm_M01Cells = hmomModel_().Norm_M01().ref();
		scalarField& Norm_N0Cells  = hmomModel_().Norm_N0().ref();
		
		// Umbalances
		double maxUmbalance = 0.;
		double sumUmbalance = 0.;

		// Integration via OpenSMOKE++ solver
		if (homogeneousReactions_ == true && odeParameterBatchReactorHomogeneous_.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
		{
			const unsigned int nhmom  = 4;
			const unsigned int neqs = ns_ + 1 + nhmom;
		
			// Min and max values
			Eigen::VectorXd yMin(neqs); yMin.setConstant(0.);  yMin(ns_) = Tmin_;	for(unsigned int i=1;i<=nhmom;i++) yMin(ns_+i) = 0.;
			Eigen::VectorXd yMax(neqs); yMax.setConstant(1.);  yMax(ns_) = Tmax_;	for(unsigned int i=1;i<=nhmom;i++) yMax(ns_+i) = 1.e16;

			Eigen::VectorXd y0(neqs);
			Eigen::VectorXd yf(neqs);

			Info <<" * Solving homogeneous chemistry (OpenSMOKE++ solver) + HMOM... "<<endl;
			{			
				unsigned int counter = 0;
			
				double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

				forAll(TCells, celli)
				{
					double tStartLocal = OpenSMOKE::OpenSMOKEGetCpuTime();

					//- Solving for celli:	
					if (TCells[celli] > direct_integration_minimum_temperature_for_chemistry_)
					{
						for(unsigned int i=0;i<ns_;i++)
							y0(i) = Y_[i].internalField()[celli];
						y0(ns_)   = TCells[celli];
						y0(ns_+1) = Norm_M00Cells[celli];
						y0(ns_+2) = Norm_M10Cells[celli];
						y0(ns_+3) = Norm_M01Cells[celli];
						y0(ns_+4) = Norm_N0Cells[celli];

						// Check and normalize the composition
						{
							double sum = 0.;
							for(unsigned int i=0;i<ns_;i++)
							{
								if (y0(i) < 0.)	y0(i) = 0.;
								sum += y0(i);
							}
							for(unsigned int i=0;i<ns_;i++)
								y0(i) /= sum;
						}

						//if (constPressureBatchReactor == true)
						{
							// Set reactor
							batchReactorHomogeneousConstantPressure_HMOM_->SetReactor(pCells[celli]);
							batchReactorHomogeneousConstantPressure_HMOM_->SetEnergyEquation(solveForEnergyEquation_);

							if ( fvCells[celli] < hmomModel_().minimum_fv_for_integration() )
								batchReactorHomogeneousConstantPressure_HMOM_->IntegrateHMOM(false);
							else
								batchReactorHomogeneousConstantPressure_HMOM_->IntegrateHMOM(true);
							
							// Set variables
							if (celli == 0)
							{
								batchReactorHomogeneousConstantPressure_HMOM_->SetHMOMOptions(	hmomModel_().index_H(), hmomModel_().index_OH(), hmomModel_().index_H2(), hmomModel_().index_H2O(), 
																hmomModel_().index_C2H2(), hmomModel_().index_O2(), hmomModel_().pah_species_indices() );
								batchReactorHomogeneousConstantPressure_HMOM_->SetPAHGasConsumption(hmomModel_().pah_gas_consumption());
							}

							// Set initial conditions
							odeSolverConstantPressure_HMOM_->SetInitialConditions(t0, y0);

							// Additional ODE solver options
							if (celli == 0)
							{
								// Set linear algebra options
								odeSolverConstantPressure_HMOM_->SetLinearAlgebraSolver(odeParameterBatchReactorHomogeneous_.linear_algebra());
								odeSolverConstantPressure_HMOM_->SetFullPivoting(odeParameterBatchReactorHomogeneous_.full_pivoting());

								// Set relative and absolute toleransces
								odeSolverConstantPressure_HMOM_->SetAbsoluteTolerances(odeParameterBatchReactorHomogeneous_.absolute_tolerance());
								odeSolverConstantPressure_HMOM_->SetRelativeTolerances(odeParameterBatchReactorHomogeneous_.relative_tolerance());

								// Set minimum and maximum values
								odeSolverConstantPressure_HMOM_->SetMinimumValues(yMin);
								odeSolverConstantPressure_HMOM_->SetMaximumValues(yMax);
							}

							// Solve
							OdeSMOKE::OdeStatus status = odeSolverConstantPressure_HMOM_->Solve(tf);

							// Transfer solution
							odeSolverConstantPressure_HMOM_->Solution(yf);
					
							// Fill heat release field
							QCells[celli] = batchReactorHomogeneousConstantPressure_HMOM_->QR();
						}
					}		
					else
					{
						for(unsigned int i=0;i<ns_;i++)
							yf(i) = Y_[i].internalField()[celli];
						yf(ns_) = TCells[celli];
						yf(ns_+1) = Norm_M00Cells[celli];
						yf(ns_+2) = Norm_M10Cells[celli];
						yf(ns_+3) = Norm_M01Cells[celli];
						yf(ns_+4) = Norm_N0Cells[celli];
					}

					// Check mass fractions and umbalance
					const double umbalance = normalizeMassFractions(yf, celli);
					maxUmbalance = std::max(umbalance, maxUmbalance);
					sumUmbalance += umbalance;

				//	if (strangAlgorithm != STRANG_COMPACT)
					{
						// Assign mass fractions
						for(int i=0;i<ns_;i++)
							Y_[i].ref()[celli] = yf(i);

						//- Allocating final values: temperature
						if (solveForEnergyEquation_ == true)
							TCells[celli] = yf(ns_);

						// Moments
						Norm_M00Cells[celli] = yf(ns_+1);
						Norm_M10Cells[celli] = yf(ns_+2);
						Norm_M01Cells[celli] = yf(ns_+3);
						Norm_N0Cells[celli]  = yf(ns_+4);
					}
				
					double tEndLocal = OpenSMOKE::OpenSMOKEGetCpuTime();
					cpuChemistryCells[celli] = (tEndLocal-tStartLocal)*1000.;

					if (counter%(int(0.20*mesh.nCells())+1) == 0)
						Info <<"   Accomplished: " << counter << "/" << mesh.nCells() << endl;

					counter++;
				}

				double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			
				Info << "   Max umbalance: " << maxUmbalance << " - Mean umbalance: " << sumUmbalance/double(mesh.nCells()) << endl;
				Info << "   Homogeneous chemistry solved in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per reactor)" << endl;
			}
		}
	
	    	Info<< " * T gas min/max (after chemistry) = " << min(T_).value() << ", " << max(T_).value() << endl;
	}

	if (hmom_ == false && drg_ == false && clustering_ == true && noxPostProcessor_ == false)
	{
		//- Clustering of cells
		clusteringModel_().Analysis(mesh);		
	
		//- Initial conditions
		const double deltat = tf - t0;
		scalarField& TCells = T_.ref();
		scalarField& pCells = p_.ref();
		scalarField& QCells = Q_.ref();
		scalarField& cpuChemistryCells = cpuChemistry_.ref();

		const scalarField& rhoCells = rho.internalField();
		const scalarField& cpCells = Cp_.internalField();
		const scalarField& cvCells = Cpv_.internalField();
		const scalarField& vCells = mesh.V();

		// Umbalances
		double maxUmbalance = 0.;
		double sumUmbalance = 0.;

		// Integration via OpenSMOKE++ solver
		if (homogeneousReactions_ == true && odeParameterBatchReactorHomogeneous_.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
		{
			const unsigned int neqs = ns_+1;
			const unsigned int nclusters = clusteringModel_().nc();
		
			// Min and max values
			Eigen::VectorXd yMin(neqs); yMin.setConstant(0.);  yMin(ns_) = Tmin_;
			Eigen::VectorXd yMax(neqs); yMax.setConstant(1.);  yMax(ns_) = Tmax_;
			Eigen::VectorXd y0(neqs);
			Eigen::VectorXd yf(neqs);

			double cpuTotal = 0.;

			Info <<" * Solving homogeneous chemistry with clustering of cells (OpenSMOKE++ solver)... "<<endl;
			{			
				unsigned int counter = 0;
			
				double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

				for (unsigned int cluster=0;cluster<nclusters;cluster++)
				{
					double tStartLocal = OpenSMOKE::OpenSMOKEGetCpuTime();

					bool chemistry_active = false;
					if (TCells[clusteringModel_().lj()[cluster][0]] > direct_integration_minimum_temperature_for_chemistry_)
						chemistry_active = true;

					//- Solving for single cluster	
					if (chemistry_active == true)
					{
						// Calculate the average value
						double pc = 0.;
						double mc = 0.;
						double vc = 0.;
						{
							y0.setZero();
							double Cvmc = 0.;
							for (unsigned int cell=0;cell<clusteringModel_().lj()[cluster].size();cell++)
							{
								const unsigned int j=clusteringModel_().lj()[cluster][cell];
								const double mj = rhoCells[j]*vCells[j];
								const double Cvj_mj = cvCells[j]*rhoCells[j]*vCells[j];

								pc += pCells[j]*vCells[j];
								for(unsigned int i=0;i<ns_;i++)
									y0(i) += Y_[i].internalField()[j]*mj;
								y0(ns_) += TCells[j]*Cvj_mj;

								mc += mj;
								Cvmc += Cvj_mj;
								vc += vCells[j];
							}
							y0 /= mc;
							y0(ns_) *= mc/Cvmc;
							pc /= vc;
						}
						const double rhoc = mc/vc;

						// Check and normalize the composition
						{
							double sum = 0.;
							for(unsigned int i=0;i<ns_;i++)
							{
								if (y0(i) < 0.)	y0(i) = 0.;
								sum += y0(i);
							}
							for(unsigned int i=0;i<ns_;i++)
								y0(i) /= sum;
						}

						//if (constPressureBatchReactor == true)
						{
							// Set reactor
							batchReactorHomogeneousConstantPressure_->SetReactor(pc);
							batchReactorHomogeneousConstantPressure_->SetEnergyEquation(solveForEnergyEquation_);

							// Set initial conditions
							odeSolverConstantPressure_->SetInitialConditions(t0, y0);

							// Additional ODE solver options
							if (cluster == 0)
							{
								// Set linear algebra options
								odeSolverConstantPressure_->SetLinearAlgebraSolver(odeParameterBatchReactorHomogeneous_.linear_algebra());
								odeSolverConstantPressure_->SetFullPivoting(odeParameterBatchReactorHomogeneous_.full_pivoting());

								// Set relative and absolute toleransces
								odeSolverConstantPressure_->SetAbsoluteTolerances(odeParameterBatchReactorHomogeneous_.absolute_tolerance());
								odeSolverConstantPressure_->SetRelativeTolerances(odeParameterBatchReactorHomogeneous_.relative_tolerance());

								// Set minimum and maximum values
								odeSolverConstantPressure_->SetMinimumValues(yMin);
								odeSolverConstantPressure_->SetMaximumValues(yMax);
							}

							// Solve
							OdeSMOKE::OdeStatus status = odeSolverConstantPressure_->Solve(tf);

							// Transfer solution
							odeSolverConstantPressure_->Solution(yf);
						
							// Fill heat release field
							for (unsigned int cell=0;cell<clusteringModel_().lj()[cluster].size();cell++)
							{
								const unsigned int j=clusteringModel_().lj()[cluster][cell];
								QCells[j] = batchReactorHomogeneousConstantPressure_->QR();
							}
						}

						// Backward-mapping
						{
							// Store mass fractions at previous time step
							std::vector<Eigen::VectorXd> yold;
							if (	solveForEnergyEquation_ == true && 
								clusteringModel_().temperature_backward_mapping() == Clustering::TEMPERATURE_RECONSTRUCTED)
							{
								yold.resize(clusteringModel_().lj()[cluster].size());
								for (unsigned int cell=0;cell<clusteringModel_().lj()[cluster].size();cell++)
								{
									const unsigned int j=clusteringModel_().lj()[cluster][cell];
									
									yold[cell].resize(ns_);
									for(int i=0;i<ns_;i++)
										yold[cell][i] = Y_[i].ref()[j];
								}
							}

							// Difference between new and old values
							Eigen::VectorXd delta = yf-y0;

							// Species backward mapping
							if (clusteringModel_().species_backward_mapping() == Clustering::SPECIES_LINEAR_INTERPOLATION_WITH_CORRECTION)
							{
								const double eps = 1e-16;
								double mtot = 0.;

								Eigen::VectorXd rho(ns_);
								
								for (unsigned int cell=0;cell<clusteringModel_().lj()[cluster].size();cell++)
								{
									const unsigned int j=clusteringModel_().lj()[cluster][cell];

									for(int i=0;i<ns_;i++)
									{
										if (delta(i)>=0.)
										{
											rho(i) = rhoCells[j]*(Y_[i].ref()[j] + delta(i));
										}
										else
										{
											rho(i) = rhoCells[j]*(Y_[i].ref()[j] + delta(i)*((Y_[i].ref()[j]+eps)/(y0(i)+eps) ) );
										}
									}

									const double rhotot = rho.sum();
									mtot += rhotot*vCells[j];
									
									// Assign mass fractions
									for(int i=0;i<ns_;i++)
										Y_[i].ref()[j] = rho(i)/rhotot;
								}	
							}
							else if (clusteringModel_().species_backward_mapping() == Clustering::SPECIES_LINEAR_INTERPOLATION)
							{
								for (unsigned int cell=0;cell<clusteringModel_().lj()[cluster].size();cell++)
								{
									const unsigned int j=clusteringModel_().lj()[cluster][cell];
						
									// Assign mass fractions
									for(int i=0;i<ns_;i++)
										Y_[i].ref()[j] += delta(i);
								}
							}

							//- Temperature backward mapping
							if (	solveForEnergyEquation_ == true && 
								clusteringModel_().temperature_backward_mapping() == Clustering::TEMPERATURE_LINEAR_INTERPOLATION )
							{
								for (unsigned int cell=0;cell<clusteringModel_().lj()[cluster].size();cell++)
								{
									const unsigned int j=clusteringModel_().lj()[cluster][cell];
									TCells[j] += delta(ns_);
								}
							}
							else if (	solveForEnergyEquation_ == true && 
									clusteringModel_().temperature_backward_mapping() == Clustering::TEMPERATURE_RECONSTRUCTED)
							{
								Eigen::VectorXd x(ns_);
								Eigen::VectorXd ynew(ns_);
								Eigen::VectorXd h_species(ns_);

								for (unsigned int cell=0;cell<clusteringModel_().lj()[cluster].size();cell++)
								{
									const unsigned int j=clusteringModel_().lj()[cluster][cell];

									// Cell mass fractions
									for(unsigned int i=0;i<ns_;++i)
										ynew(i) = Y_[i].ref()[j];

									// Set thermodynamic map
									thermodynamicsMap_->SetPressure(pCells[j]);
									thermodynamicsMap_->SetTemperature(TCells[j]);

									// Calculate the mole fractions from the mass fractions (and the molecular weight in kg/kmol)
									double mw;
									thermodynamicsMap_->MoleFractions_From_MassFractions(x.data(), mw, yold[cell].data());

									// Mixture enthalpy at previous time level [J/kg]
									const double H0 = thermodynamicsMap_->hMolar_Mixture_From_MoleFractions(x.data()) / mw;

									// Enthalpies of species at previous time level
									thermodynamicsMap_->hMolar_Species(h_species.data());	// total ethalpies of species [J/kmol]
									for(unsigned int i=0;i<ns_;++i)
										h_species(i) /= thermodynamicsMap_->MW(i);	// total enthalpies of species [J/kg]

									// Change of enthalpy (linearized)
									double deltah = 0.;
									for(unsigned int i=0;i<ns_;++i)
										deltah += h_species(i)*(ynew(i)-yold[cell](i)); // total change of enthalpy [J/kg]

									// New enthalpy [J/kg]
									const double H = H0 + deltah;

									// Calculate the mole fractions from the mass fractions (and the molecular weight in kg/kmol)
									thermodynamicsMap_->MoleFractions_From_MassFractions(x.data(), mw, ynew.data());

									// Get new temperature from enthalpy
									const double temp = thermodynamicsMap_->GetTemperatureFromEnthalpyAndMoleFractions(H*mw, pCells[j], x.data(), TCells[j]);
									if (temp != 0.) 
									{
										TCells[j] = temp;
										TCells[j] = std::max(TCells[j], Tmin_);
										TCells[j] = std::min(TCells[j], Tmax_);
									}
								} // loop over cells in a cluster

							} // temperature backward-mapping

						} // backward-mapping

					} // loop over clusters		

					// Check mass fractions and umbalance
				//	const double umbalance = normalizeMassFractions(yf, celli);
				//	maxUmbalance = std::max(umbalance, maxUmbalance);
				//	sumUmbalance += umbalance;
				
					double tEndLocal = OpenSMOKE::OpenSMOKEGetCpuTime();

					// Fill cpuChemistry field
					for (unsigned int cell=0;cell<clusteringModel_().lj()[cluster].size();cell++)
					{
						const unsigned int j=clusteringModel_().lj()[cluster][cell];
						cpuChemistryCells[j] = (tEndLocal-tStartLocal)*1000.;
					}

					cpuTotal += (tEndLocal-tStartLocal);

					if (counter%(int(0.20*nclusters)+1) == 0)
						Info <<"   Accomplished: " << counter << "/" << nclusters << " Total(ms): " << cpuTotal*1000 << endl;

					counter++;		
				}

				double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			
			//	Info << "   Max umbalance: " << maxUmbalance << " - Mean umbalance: " << sumUmbalance/double(mesh.nCells()) << endl;
				Info << "   Homogeneous chemistry solved in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per cell)" << endl;
				Info << "   Homogeneous chemistry solved in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(nclusters)*1000. << " ms per cluster)" << endl;
			}
		}
	
	    	Info<< " * T gas min/max (after chemistry) = " << min(T_).value() << ", " << max(T_).value() << endl;

	}

	if (hmom_ == false && drg_ == true && clustering_ == true && noxPostProcessor_ == false)
	{
		//- Clustering of cells
		clusteringModel_().Analysis(mesh);		
	
		//- Initial conditions
		const double deltat = tf - t0;
		scalarField& TCells = T_.ref();
		scalarField& pCells = p_.ref();
		scalarField& QCells = Q_.ref();
		scalarField& cpuChemistryCells = cpuChemistry_.ref();

		const scalarField& rhoCells = rho.internalField();
		const scalarField& cpCells = Cp_.internalField();
		const scalarField& cvCells = Cpv_.internalField();
		const scalarField& vCells = mesh.V();

		scalarField& drgSpeciesCells = drgModel_().important_species().ref();
		scalarField& drgReactionsCells = drgModel_().important_reactions().ref();

		// Umbalances
		double maxUmbalance = 0.;
		double sumUmbalance = 0.;

		// Integration via OpenSMOKE++ solver
		if (homogeneousReactions_ == true && odeParameterBatchReactorHomogeneous_.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
		{
			// Total number of clusters
			const unsigned int nclusters = clusteringModel_().nc();
		
			// Min and max values
			Eigen::VectorXd yMin; 
			Eigen::VectorXd yMax;
			Eigen::VectorXd y0;
			Eigen::VectorXd yf;

			// Mass fractions, mole fractions, and concentrations
			Eigen::VectorXd omega(ns_);
			Eigen::VectorXd x(ns_);
			Eigen::VectorXd c(ns_);

			double cpuDRG = 0.;
			double cpuSolution = 0.;
			double cpuTotal = 0.;
			int speciesDRG = 0;
			int reactionsDRG = 0;

			Info <<" * Solving homogeneous chemistry with clustering of cells (OpenSMOKE++ solver) + DRG... "<<endl;
			{			
				unsigned int counter = 0;
			
				double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

				for (unsigned int cluster=0;cluster<nclusters;cluster++)
				{
					double tStartLocal = OpenSMOKE::OpenSMOKEGetCpuTime();

					bool chemistry_active = false;
					if (TCells[clusteringModel_().lj()[cluster][0]] > direct_integration_minimum_temperature_for_chemistry_)
						chemistry_active = true;

					//- Solving for single cluster	
					if (chemistry_active == true)
					{
						// Calculate the average value and apply the DRG
						const double tStartDRG = OpenSMOKE::OpenSMOKEGetCpuTime();

						double pc = 0.;
						double mc = 0.;
						double vc = 0.;
						double Tc = 0.;
						omega.setZero();
						{
							double Cvmc = 0.;
							for (unsigned int cell=0;cell<clusteringModel_().lj()[cluster].size();cell++)
							{
								const unsigned int j=clusteringModel_().lj()[cluster][cell];
								const double mj = rhoCells[j]*vCells[j];
								const double Cvj_mj = cvCells[j]*rhoCells[j]*vCells[j];

								pc += pCells[j]*vCells[j];
								for(unsigned int i=0;i<ns_;i++)
									omega(i) += Y_[i].internalField()[j]*mj;
								Tc += TCells[j]*Cvj_mj;

								mc += mj;
								Cvmc += Cvj_mj;
								vc += vCells[j];
							}
							omega /= mc;
							Tc /= Cvmc;
							pc /= vc;

							// Concentrations
							double mw;
							thermodynamicsMap_->MoleFractions_From_MassFractions(x.data(), mw, omega.data());
							double cTot = pc/PhysicalConstants::R_J_kmol/Tc;
							c = cTot*x;

							// Analysis
							drgModel_().Analysis(Tc, pc, c);
							//Info << "Done" << endl;
							// Fill the relevant DRG fields
							for (unsigned int cell=0;cell<clusteringModel_().lj()[cluster].size();cell++)
							{
								const unsigned int j=clusteringModel_().lj()[cluster][cell];
								drgSpeciesCells[j]   = drgModel_().number_important_species();
								drgReactionsCells[j] = drgModel_().number_important_reactions();
							}
						}

						const double tEndDRG = OpenSMOKE::OpenSMOKEGetCpuTime();

						// Select the relevant species only
						const double tStartSolution = OpenSMOKE::OpenSMOKEGetCpuTime();
						{
							const unsigned int neqs = drgModel_().number_important_species()+1;

							y0.resize(neqs);
							yf.resize(neqs);
							yMin.resize(neqs);	yMin.setConstant(0.); yMin(neqs-1) = Tmin_;
							yMax.resize(neqs);	yMax.setConstant(1.); yMax(neqs-1) = Tmax_;

							for (unsigned int i=0;i<drgModel_().number_important_species();++i)	
							{
								const unsigned int j = drgModel_().indices_important_species()[i];
								y0(i) = omega(j);
							}
							y0(neqs-1) = Tc;							
						}

						//if (constPressureBatchReactor == true)
						{
							// Set reactor
							batchReactorHomogeneousConstantPressure_DRG_->SetReactor(pc);
							batchReactorHomogeneousConstantPressure_DRG_->SetEnergyEquation(solveForEnergyEquation_);

							// Set DRG
							batchReactorHomogeneousConstantPressure_DRG_->SetDRGOptions( drgModel_().indices_important_species(), drgModel_().indices_unimportant_reactions(), drgModel_().important_reactions_steps() );
							batchReactorHomogeneousConstantPressure_DRG_->SetMassFractions(omega);

							// Set initial conditions
							odeSolverConstantPressure_DRG_->SetInitialConditions(t0, y0);

							// Additional ODE solver options
							{
								// Set linear algebra options
								odeSolverConstantPressure_DRG_->SetLinearAlgebraSolver(odeParameterBatchReactorHomogeneous_.linear_algebra());
								odeSolverConstantPressure_DRG_->SetFullPivoting(odeParameterBatchReactorHomogeneous_.full_pivoting());

								// Set relative and absolute toleransces
								odeSolverConstantPressure_DRG_->SetAbsoluteTolerances(odeParameterBatchReactorHomogeneous_.absolute_tolerance());
								odeSolverConstantPressure_DRG_->SetRelativeTolerances(odeParameterBatchReactorHomogeneous_.relative_tolerance());

								// Set minimum and maximum values
								odeSolverConstantPressure_DRG_->SetMinimumValues(yMin);
								odeSolverConstantPressure_DRG_->SetMaximumValues(yMax);
							}

							// Solve
							OdeSMOKE::OdeStatus status = odeSolverConstantPressure_DRG_->Solve(tf);

							// Transfer solution
							odeSolverConstantPressure_DRG_->Solution(yf);

							// Fill heat release field
							for (unsigned int cell=0;cell<clusteringModel_().lj()[cluster].size();cell++)
							{
								const unsigned int j=clusteringModel_().lj()[cluster][cell];
								QCells[j] = batchReactorHomogeneousConstantPressure_DRG_->QR();
							}
						}

						// Complete vector (all the species + temperature) before ODE integration
						Eigen::VectorXd y00(ns_+1);
						y00.head(ns_) = omega;
						y00(ns_) = Tc;

						// Complete vector (all the species + temperature) after ODE integration
						Eigen::VectorXd yff(ns_+1);
						yff.head(ns_) = omega;
						yff(ns_) = yf(yf.size()-1);
						for (unsigned int i=0;i<drgModel_().number_important_species();++i)	
						{
							const unsigned int j = drgModel_().indices_important_species()[i];
							yff(j) = yf(i);
						}

						// Backward-mapping
						{
							// Store mass fractions at previous time step
							std::vector<Eigen::VectorXd> yold;
							if (	solveForEnergyEquation_ == true && 
								clusteringModel_().temperature_backward_mapping() == Clustering::TEMPERATURE_RECONSTRUCTED)
							{
								yold.resize(clusteringModel_().lj()[cluster].size());
								for (unsigned int cell=0;cell<clusteringModel_().lj()[cluster].size();cell++)
								{
									const unsigned int j=clusteringModel_().lj()[cluster][cell];
									
									yold[cell].resize(ns_);
									for(int i=0;i<ns_;i++)
										yold[cell][i] = Y_[i].ref()[j];
								}
							}

							// Difference between new and old values
							Eigen::VectorXd delta = yff-y00;

							// Species backward mapping
							if (clusteringModel_().species_backward_mapping() == Clustering::SPECIES_LINEAR_INTERPOLATION_WITH_CORRECTION)
							{
								const double eps = 1e-16;
								double mtot = 0.;

								Eigen::VectorXd rho(ns_);
								
								for (unsigned int cell=0;cell<clusteringModel_().lj()[cluster].size();cell++)
								{
									const unsigned int j=clusteringModel_().lj()[cluster][cell];

									for(int i=0;i<ns_;i++)
									{
										if (delta(i)>=0.)
										{
											rho(i) = rhoCells[j]*(Y_[i].ref()[j] + delta(i));
										}
										else
										{
											rho(i) = rhoCells[j]*(Y_[i].ref()[j] + delta(i)*((Y_[i].ref()[j]+eps)/(y00(i)+eps) ) );
										}
									}

									const double rhotot = rho.sum();
									mtot += rhotot*vCells[j];
									
									// Assign mass fractions
									for(int i=0;i<ns_;i++)
										Y_[i].ref()[j] = rho(i)/rhotot;
								}	
							}
							else if (clusteringModel_().species_backward_mapping() == Clustering::SPECIES_LINEAR_INTERPOLATION)
							{
								for (unsigned int cell=0;cell<clusteringModel_().lj()[cluster].size();cell++)
								{
									const unsigned int j=clusteringModel_().lj()[cluster][cell];
						
									// Assign mass fractions
									for(int i=0;i<ns_;i++)
										Y_[i].ref()[j] += delta(i);
								}
							}

							//- Temperature backward mapping
							if (	solveForEnergyEquation_ == true && 
								clusteringModel_().temperature_backward_mapping() == Clustering::TEMPERATURE_LINEAR_INTERPOLATION )
							{
								for (unsigned int cell=0;cell<clusteringModel_().lj()[cluster].size();cell++)
								{
									const unsigned int j=clusteringModel_().lj()[cluster][cell];
									TCells[j] += delta(ns_);
								}
							}
							else if (	solveForEnergyEquation_ == true && 
									clusteringModel_().temperature_backward_mapping() == Clustering::TEMPERATURE_RECONSTRUCTED)
							{
								Eigen::VectorXd x(ns_);
								Eigen::VectorXd ynew(ns_);
								Eigen::VectorXd h_species(ns_);

								for (unsigned int cell=0;cell<clusteringModel_().lj()[cluster].size();cell++)
								{
									const unsigned int j=clusteringModel_().lj()[cluster][cell];

									// Cell mass fractions
									for(unsigned int i=0;i<ns_;++i)
										ynew(i) = Y_[i].ref()[j];

									// Set thermodynamic map
									thermodynamicsMap_->SetPressure(pCells[j]);
									thermodynamicsMap_->SetTemperature(TCells[j]);

									// Calculate the mole fractions from the mass fractions (and the molecular weight in kg/kmol)
									double mw;
									thermodynamicsMap_->MoleFractions_From_MassFractions(x.data(), mw, yold[cell].data());

									// Mixture enthalpy at previous time level [J/kg]
									const double H0 = thermodynamicsMap_->hMolar_Mixture_From_MoleFractions(x.data()) / mw;

									// Enthalpies of species at previous time level
									thermodynamicsMap_->hMolar_Species(h_species.data());	// total ethalpies of species [J/kmol]
									for(unsigned int i=0;i<ns_;++i)
										h_species(i) /= thermodynamicsMap_->MW(i);	// total enthalpies of species [J/kg]

									// Change of enthalpy (linearized)
									double deltah = 0.;
									for(unsigned int i=0;i<ns_;++i)
										deltah += h_species(i)*(ynew(i)-yold[cell](i)); // total change of enthalpy [J/kg]

									// New enthalpy [J/kg]
									const double H = H0 + deltah;

									// Calculate the mole fractions from the mass fractions (and the molecular weight in kg/kmol)
									thermodynamicsMap_->MoleFractions_From_MassFractions(x.data(), mw, ynew.data());

									// Get new temperature from enthalpy
									const double temp = thermodynamicsMap_->GetTemperatureFromEnthalpyAndMoleFractions(H*mw, pCells[j], x.data(), TCells[j]);
									if (temp != 0.) 
									{
										TCells[j] = temp;
										TCells[j] = std::max(TCells[j], Tmin_);
										TCells[j] = std::min(TCells[j], Tmax_);
									}
								} // loop over cells in a cluster

							} // temperature backward-mapping

						} // backward-mapping

						const double tEndSolution = OpenSMOKE::OpenSMOKEGetCpuTime();

						cpuDRG += (tEndDRG-tStartDRG);
						cpuSolution += (tEndSolution-tStartSolution);

					} // loop over clusters		
					else
					{
						// Fill the relevant DRG fields
						for (unsigned int cell=0;cell<clusteringModel_().lj()[cluster].size();cell++)
						{
							const unsigned int j=clusteringModel_().lj()[cluster][cell];
							drgSpeciesCells[j]   = drgModel_().number_key_species();
							drgReactionsCells[j] = 0;
						}
					}

					speciesDRG += drgSpeciesCells[clusteringModel_().lj()[cluster][0]];
					reactionsDRG += drgReactionsCells[clusteringModel_().lj()[cluster][0]];

					double tEndLocal = OpenSMOKE::OpenSMOKEGetCpuTime();

					// Fill cpuChemistry field
					for (unsigned int cell=0;cell<clusteringModel_().lj()[cluster].size();cell++)
					{
						const unsigned int j=clusteringModel_().lj()[cluster][cell];
						cpuChemistryCells[j] = (tEndLocal-tStartLocal)*1000.;
					}

					cpuTotal += (tEndLocal-tStartLocal);

					if (counter%(int(0.20*nclusters)+1) == 0)
						Info 	<< "   Accomplished: " << counter << "/" << nclusters 
							<< " ns=" << speciesDRG/double(nclusters) << " nr=" << reactionsDRG/double(nclusters) 
							<< " DRG(%)=" << cpuDRG/cpuTotal*100. << " Sol(%)=" << cpuSolution/cpuTotal*100. << " Tot(ms)=" << cpuTotal*1000 
							<< endl;
					counter++;	
				}

                                double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

                        //      Info << "   Max umbalance: " << maxUmbalance << " - Mean umbalance: " << sumUmbalance/double(mesh.nCells()) << endl;
                                Info << "   Homogeneous chemistry solved in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per cell)" << endl;
                                Info << "   Homogeneous chemistry solved in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(nclusters)*1000. << " ms per cluster)" << endl;
                        }
                }

		Info<< " * T gas min/max (after chemistry) = " << min(T_).value() << ", " << max(T_).value() << endl;
	}

	else if (noxPostProcessor_ == true)
	{
		//- Initial conditions
		const double deltat = tf - t0;
		
		scalarField& cpuChemistryCells = cpuChemistry_.ref();

		const scalarField& TCells = T_.internalField();
		const scalarField& pCells = p_.internalField();
		const scalarField& rhoCells = rho.internalField();

		scalarField& YNOCells = noxPostProcessorModel_().YNO().ref();
	//	scalarField& YHCNCells = noxPostProcessorModel_().YHCN().ref();
	//	scalarField& YNH3Cells = noxPostProcessorModel_().YNH3().ref();
	//	scalarField& YN2OCells = noxPostProcessorModel_().YN2O().ref();

		// Integration via OpenSMOKE++ solver
		if (homogeneousReactions_ == true && odeParameterBatchReactorHomogeneous_.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
		{
			const unsigned int neqs = 1;
		
			// Min values
			Eigen::VectorXd yMin(neqs); 
			yMin.setConstant(0.);

			// Max values
			Eigen::VectorXd yMax(neqs); 
			yMax.setConstant(1.); 

			// Mass fractions of rlevant species
			Eigen::VectorXd omega(6);
			omega.setZero();

			// Solution vectors
			Eigen::VectorXd y0(neqs);
			Eigen::VectorXd yf(neqs);

			Info <<" * Solving NOx (OpenSMOKE++ solver)... "<<endl;
			{			
				unsigned int counter = 0;
			
				double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

				forAll(TCells, celli)
				{
					double tStartLocal = OpenSMOKE::OpenSMOKEGetCpuTime();

					//- Solving for celli:	
					if (TCells[celli] > direct_integration_minimum_temperature_for_chemistry_)
					{
						// Mass fractions of relevant species
						{
							omega(0) = Y_[noxPostProcessorModel_().index_N2()].internalField()[celli];
							omega(1) = Y_[noxPostProcessorModel_().index_O2()].internalField()[celli];
							omega(2) = Y_[noxPostProcessorModel_().index_H2O()].internalField()[celli];
							omega(3) = Y_[noxPostProcessorModel_().index_Fuel()].internalField()[celli];

							if (noxPostProcessorModel_().index_O() >= 0)
								omega(4) = Y_[noxPostProcessorModel_().index_O()].internalField()[celli];

							if (noxPostProcessorModel_().index_OH() >= 0)
								omega(5) = Y_[noxPostProcessorModel_().index_OH()].internalField()[celli];
						}
				
						// Unknowns
						y0(0) = std::max(YNOCells[celli],0.);
					//	y0(1) = std::max(YHCNCells[celli],0.);
					//	y0(2) = std::max(YNH3Cells[celli],0.);
					//	y0(3) = std::max(YN2OCells[celli],0.);

						//if (constPressureBatchReactor == true)
						{
						//	batchReactorHomogeneousConstantPressure_NOxPP_->SetDebug(-1);
						//	if (TCells[celli] > 1600.)
						//		batchReactorHomogeneousConstantPressure_NOxPP_->SetDebug(celli);
							// Set reactor
							batchReactorHomogeneousConstantPressure_NOxPP_->SetReactor(TCells[celli], pCells[celli], rhoCells[celli]);
							batchReactorHomogeneousConstantPressure_NOxPP_->SetMassFractions(omega);

							// Turn on integration of chemical step
							batchReactorHomogeneousConstantPressure_NOxPP_->IntegrateNOxPP(true);
							
							// Set variables
							if (celli == 0)
								batchReactorHomogeneousConstantPressure_NOxPP_->SetNOxPPOptions( noxPostProcessorModel_().index_OH(), noxPostProcessorModel_().index_O(), 
																 noxPostProcessorModel_().index_CH(), noxPostProcessorModel_().index_Fuel() );

							// Set initial conditions
							odeSolverConstantPressure_NOxPP_->SetInitialConditions(t0, y0);

							// Additional ODE solver options
							if (celli == 0)
							{
								// Set linear algebra options
								odeSolverConstantPressure_NOxPP_->SetLinearAlgebraSolver(odeParameterBatchReactorHomogeneous_.linear_algebra());
								odeSolverConstantPressure_NOxPP_->SetFullPivoting(odeParameterBatchReactorHomogeneous_.full_pivoting());

								// Set relative and absolute toleransces
								odeSolverConstantPressure_NOxPP_->SetAbsoluteTolerances(odeParameterBatchReactorHomogeneous_.absolute_tolerance());
								odeSolverConstantPressure_NOxPP_->SetRelativeTolerances(odeParameterBatchReactorHomogeneous_.relative_tolerance());

								// Set minimum and maximum values
								odeSolverConstantPressure_NOxPP_->SetMinimumValues(yMin);
								odeSolverConstantPressure_NOxPP_->SetMaximumValues(yMax);
							}

							// Solve
							OdeSMOKE::OdeStatus status = odeSolverConstantPressure_NOxPP_->Solve(tf);

							// Transfer solution
							odeSolverConstantPressure_NOxPP_->Solution(yf);
						}
					}		
					else
					{
						yf(0) = YNOCells[celli];
					//	yf(1) = YHCNCells[celli];
					//	yf(2) = YNH3Cells[celli];
					//	yf(3) = YN2OCells[celli];

					}

				//	if (strangAlgorithm != STRANG_COMPACT)
					{
						YNOCells[celli] = yf(0);
					//	YHCNCells[celli] = yf(1);
					//	YNH3Cells[celli] = yf(2);
					//	YN2OCells[celli] = yf(3);
					}
				
					double tEndLocal = OpenSMOKE::OpenSMOKEGetCpuTime();
					cpuChemistryCells[celli] = (tEndLocal-tStartLocal)*1000.;

					if (counter%(int(0.20*mesh.nCells())+1) == 0)
						Info <<"   Accomplished: " << counter << "/" << mesh.nCells() << endl;

					counter++;
				}

				double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			
				Info << "   NOx chemistry solved in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per reactor)" << endl;
			}
		}

//				double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			
			//	Info << "   Max umbalance: " << maxUmbalance << " - Mean umbalance: " << sumUmbalance/double(mesh.nCells()) << endl;
//				Info << "   Homogeneous chemistry solved in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per cell)" << endl;
//				Info << "   Homogeneous chemistry solved in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(nclusters)*1000. << " ms per cluster)" << endl;
//			}
//		}
	
	    	Info<< " * T gas min/max (after chemistry) = " << min(T_).value() << ", " << max(T_).value() << endl;

	}
}
