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
|   Copyright(C) 2020 Alberto Cuoci                                       |
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

// Include standard OpenFOAM classes
#include "fvCFD.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "pressureControl.H"
#include "fvOptions.H"
#include "fvcSmooth.H"

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
	#include "createMesh.H"
	#include "createControl.H"
	#include "createTimeControls.H"
	#include "initContinuityErrs.H"
	#include "createFields.H"
	#include "createFieldRefs.H"
	#include "compressibleCourantNo.H"
	#include "setInitialDeltaT.H"

	
	Info<< nl << "Starting time loop" << nl << endl;

	while (pimple.run(runTime))
	{
		// Set the time step 
		#include "readTimeControls.H"
		#include "compressibleCourantNo.H"
		#include "setDeltaT.H"

		// Advance in time
		runTime++;
		Info<< "Time = " << runTime.timeName() << nl << endl;

		// Solve the continuity equation for density
		#include "rhoEqn.H"

		// PIMPLE loop
		while (pimple.loop())
		{
			// Momentum equations
			#include "UEqn.H"

			// Update fluxes and transport terms in species/energy equations
			mixture.update_transport_terms(mesh, rho);

			// Species equations 
			#include "YEqn.H"

			// Additional equations (HMOM, etc.)
			mixture.solve_additional_equations(mesh, rho, phi);

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
				if (pimple.consistent())
				{
					#include "pcEqn.H"
				}
				else
				{
					#include "pEqn.H"
				}
			}
		}

		// Update the density from the new pressure field
		rho = mixture.rho();

		// Passive scalars
            	#include "csiEqn.H"
		#include "tauEqn.H"

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

