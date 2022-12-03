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

if (!pimple.flow())
{
    if (pimple.models())
    {
        fvModels.correct();
    }

    // Thermophysics
    {
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
    }
}
else
{
    if (!mesh.steady() && pimple.firstPimpleIter())
    {
	// Solve the continuity equation for density
	#include "rhoEqn.H"
    }

    if (pimple.models())
    {
        fvModels.correct();
    }

    // Momentum equations
    #include "UEqn.H"

    // Thermophysics
    {
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
    }

    // --- PISO loop
    while (pimple.correct())
    {
        #include "../laminarPimpleSMOKE++/pEqn.H"
    }

    if (!mesh.steady() && pimple.finalPimpleIter())
    {
        rho = mixture.rho();
    }
}
