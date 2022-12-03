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

// This is a multi-region solver
#define MULTIREGIONSOLVER 1

// Include standard OpenFOAM classes
#include "fvCFD.H"
#include "multivariateScheme.H"
#include "pimpleMultiRegionControl.H"
#include "pressureControl.H"
#include "pressureReference.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "coordinateSystem.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "OFstream.H"
#include "compressibleCourantNo.H"

// Multiregion classes
#include "regionProperties.H"
#include "solidRegionDiffNo.H"
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

	// To be able to use postProcess
	//#define CREATE_MESH createMeshes.H
	//#include "postProcess.H"

	//#include "addCheckCaseOptions.H"
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
	#include "createRhoUfIfPresent.H"

	// Create numerical controls
    	#include "initContinuityErrs.H"

    	pimpleMultiRegionControl pimples(fluidRegions, solidRegions);
    	#include "createFluidPressureControls.H"

	// The fluid region is unique
	pimpleNoLoopControl& pimple = pimples.pimple(0);
	pressureControl& pressureControl = pressureControlFluid[0];
	pressureReference pressureReference(p, pimple.dict(), false);

    	#include "createTimeControls.H"
   	#include "readSolidTimeControls.H"
    	#include "compressibleMultiRegionCourantNo.H"
    	#include "solidRegionDiffusionNo.H"
    	#include "setInitialMultiRegionDeltaT.H"

	Info<< nl << "Starting time loop" << nl << endl;

	while (pimples.run(runTime))
    	{
        	#include "readTimeControls.H"
        	#include "readSolidTimeControls.H"

       	 	#include "compressibleMultiRegionCourantNo.H"
        	#include "solidRegionDiffusionNo.H"
        	#include "setMultiRegionDeltaT.H"

        	runTime++;

        	Info<< "Time = " << runTime.timeName() << nl << endl;

		// Check for regions to be solved
		const bool is_fluid_active = runTime.controlDict().lookupOrDefault<Switch>("fluid", true);
		const bool is_solid_active = runTime.controlDict().lookupOrDefault<Switch>("solid", true);

		// Solve the continuity equation for density
		#include "rhoEqn.H"

		// PIMPLE loop
        	while (pimples.loop())
        	{
			if (is_fluid_active == true)
			{

				// Momentum equations
				#include "UEqn.H"

				// Update fluxes and transport terms in species/energy equations
				mixture.update_transport_terms(mesh, rho);

				// Species equations 
				#include "YEqn.H"

				// Additional equations (HMOM, etc.)
				//mixture.solve_additional_equations(mesh, rho, phi);

				// Energy equation
				#include "EEqn.H"

				// Chemical step
				{
					const double t0 = runTime.value() - runTime.deltaT().value();
					const double tf = runTime.value();
					mixture.chemistry_direct_integration(t0, tf, mesh, rho);
				}

				// Update mixture properties
				mixture.update_properties();

				// Pressure corrector loop
				while (pimple.correct())
				{
					#include "pEqn.H"
				}

				// Update the density from the new pressure field
				rho = mixture.rho();

				// Passive scalars
				if (is_fluid_active == true)
				{
            				#include "csiEqn.H"
					#include "tauEqn.H"
				}
			}

			
			if (is_solid_active == true)
			{
			//	forAll(solidRegions, i)
            			{
                	//		Info<< "\nSolving for solid region " << solidRegions[i].name() << endl;

					// To solve the solid energy equation
                	//		#include "setRegionSolidFields.H"
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

	Info<< "End" << endl;

	return 0;
}

