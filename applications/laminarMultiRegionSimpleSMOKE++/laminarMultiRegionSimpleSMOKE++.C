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
|   Copyright(C) 2022 Alberto Cuoci                                       |
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

// This is a steady state solver
#define STEADYSTATE 1

// This is a multi-region solver
#define MULTIREGIONSOLVER 1

// Include standard OpenFOAM classes
#include "fvCFD.H"
#include "multivariateScheme.H"
#include "simpleControl.H"

#include "fvModels.H"
#include "fvConstraints.H"
#include "coordinateSystem.H"
#include "pressureReference.H"

#include "interpolation.H"
#include "fvcSmooth.H"
#include "OFstream.H"

// Multiregion classes
#include "regionProperties.H"
#include "solidThermo.H"

// Include laminarSMOKE++ classes
#include "SparkModel.H"
#include "OpenSMOKEppReactingMixture.H"
#include "OpenSMOKEppReactingMixtureForRadiation.H"

// Customized radiation model
#include "radiationModelOpenSMOKE++.H"

// Main code
int main(int argc, char *argv[])
{
	#include "postProcess.H"

	#include "setRootCaseLists.H"
	#include "createTime.H"

	// Create multiple meshes
	#include "createMeshes.H"
	fvMesh& mesh = fluidRegions[0];
	fvMesh& solidMesh = solidRegions[0];

	// Create fields
	#include "createFields.H"
	#include "createSolidFields.H"
	#include "createFieldRefs.H"

	// Create numerical controls
	#include "createControl.H"
	#include "createTimeControls.H"
	#include "initContinuityErrs.H"

    	// Create simple controls
	pressureReference pressureReference(p, rho, simple.dict());
        /* 
	PtrList<simpleControl> simpleControlSolid(solidRegions.size());
	forAll(solidRegions, i)
	{
	    	simpleControlSolid.set( i, new simpleControl( solidRegions[i] ) );
	}
	*/
	simpleControl simpleControlSolid(solidMesh);
	
	Info<< nl << "Starting time loop" << nl << endl;
	if (mixture.noxPostProcessor() == true)
	{
		// Update mixture properties
		mixture.update_properties_for_steady_state();

		// Update fluxes and transport terms in species/energy equations
		mixture.update_transport_terms(mesh, rho);

		while (simple.loop(runTime))
		{
			Info<< "Time = " << runTime.timeName() << nl << endl;

			// Additional equations (HMOM, etc.)
			mixture.solve_additional_equations(mesh, rho, phi);																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																				 

			// On-the-fly post processing
			#include "ontheflyPostProcessing.H"

			// Write fields
			runTime.write();

			// Write CPU time on the screen
			Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl;
		}
	}
	else
	{
	while (simple.loop(runTime))
	{
		// Check for regions to be solved
		const bool is_fluid_active = runTime.controlDict().lookupOrDefault<Switch>("fluid", true);
		const bool is_solid_active = runTime.controlDict().lookupOrDefault<Switch>("solid", true);

		Info<< "Time = " << runTime.timeName() << nl << endl;
		
	 	const double t0 = runTime.value() - runTime.deltaT().value();
	 	const double tf = runTime.value();

		// Pressure-velocity SIMPLE corrector
         	{
			// Update SIMPLE control for solid phase
			/*
			forAll(simpleControlSolid, i)
			{
				simpleControlSolid[i].fluidSolutionControl::read() && simpleControlSolid[i].readResidualControls();

				// Store previous field: solid
				if (!simpleControlSolid[i].endIfConverged(runTime))
				{
					simpleControlSolid[i].storePrevIterFields();
				}
			}
			*/
                        {
                                simpleControlSolid.fluidSolutionControl::read() && simpleControlSolid.readResidualControls();

                                // Store previous field: solid
                                if (!simpleControlSolid.endIfConverged(runTime))
                                {
                                	simpleControlSolid.storePrevIterFields();
                                }
                        }

			// Fluid region equations
			if (is_fluid_active == true)
			{
				if (mixture.solveForMomentumEquation() == true)
				{
		    			#include "UEqn.H"

		    			// Update mixture properties
		    			mixture.update_properties_for_steady_state();

		    			// Update local Jacobian matrices
		    			mixture.update_chemical_source_terms();

		    			// Update fluxes and transport terms in species/energy equations
		    			mixture.update_transport_terms(mesh, rho);

		    			// Species equations
		    			#include "YEqn.H"

		    			// Energy equation
		    			#include "EEqn.H" 

                    			// Pressure equation
		    			{
		        			#include "pEqn.H"
		    			}
				}
				else
				{
			    		// Update mixture properties
			    		mixture.update_properties_for_steady_state();

		    			// Update local Jacobian matrices
					mixture.update_chemical_source_terms();

		    			// Update fluxes and transport terms in species/energy equations
		   	 		mixture.update_transport_terms(mesh, rho);

		    			// Species equations
		    			#include "YEqn.H"

		    			// Energy equation
		    			#include "EEqn.H"  
				}
	    
				// Passive scalars
	        		#include "csiEqn.H"
				#include "tauEqn.H"
			}

			// Solid region equations
			if (is_solid_active == true)
			{
				// Solve solid regions
				// forAll(solidRegions, i)
				{
					//Info << "\nSolving for solid region " << solidRegions[i].name() << endl;
					//#include "setRegionSolidFields.H"
					#include "solidEEqn.H"
				}
			}
         	}																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																								 

		// On-the-fly post processing
		#include "ontheflyPostProcessing.H"

		// Write fields
		runTime.write();

		// Write CPU time on the screen
		Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl;
	}
	}

	Info<< "End" << endl;

	return 0;
}

