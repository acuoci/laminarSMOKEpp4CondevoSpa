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

// This is a steady state solver
#define STEADYSTATE 1

// Include standard OpenFOAM classes
#include "fvCFD.H"
#include "multivariateScheme.H"
#include "simpleControl.H"
#include "pressureControl.H"
#include "interpolation.H"
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
	
	Info<< nl << "Starting time loop" << nl << endl;

	while (simple.loop(runTime))
	{
		Info<< "Time = " << runTime.timeName() << nl << endl;
		
	 	const double t0 = runTime.value() - runTime.deltaT().value();
	 	const double tf = runTime.value();

		// Pressure-velocity SIMPLE corrector
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
		    		if (simple.consistent())
		    		{
		    			#include "pcEqn.H"
		    		}
		    		else
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

