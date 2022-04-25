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

{
	// Heat dissipation due to viscous forces
	if (mixture.viscous_heating_source_term() == true)
	{
		volSymmTensorField rhoTau
		(
		    IOobject
		    (
			"mix:rhoTau",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		    ),
		    mu*dev(twoSymm(fvc::grad(U)))
		);

		Qdissipation = rhoTau && fvc::grad(U);
	}

	// Radiative heat transfer
	radiation->correct();

	// Energy equation in terms of temperature
	if (mixture.energyEquation() == ENERGY_EQ_TEMPERATURE)
	{
		fvScalarMatrix TEqn
		(
			  Cp*fvm::ddt(rho,T)
			+ Cp*fvm::div(phi, T)
                        - dpdt		    
			==
			  fvm::laplacian(lambda,T) 
			- ( mixture.Cps_Js() & fvc::grad(T) )
			+ radiation->divq(T)
			+ Qdissipation
			+ Cp*fvOptions(rho, T)
		);

		TEqn.relax();

		fvOptions.constrain(TEqn);

		TEqn.solve("T");

		fvOptions.correct(T);

		Info<< " * T gas min/max (after transport)   = " << min(T).value() << ", " << max(T).value() << endl;

	}

	// Energy equation in terms of enthalpy
	else if (mixture.energyEquation() == ENERGY_EQ_ENTHALPY)
	{
		fvScalarMatrix HEqn
		(
			  fvm::ddt(rho, hs)
			+ mvConvection->fvmDiv(phi, hs) 
			- dpdt
			+ ( 
				mixture.conduction_via_temperature_gradient() == true
				? ( -fvc::laplacian(lambda,T) ) 
				: ( -fvc::laplacian(lambda/Cp,hs) + fvc::div(mixture.HsConduction(), "div(mix:hs_Js)") )
			  )
			+ fvc::div(mixture.hs_Js())
			==
			//  radiation->divq(T) +	// TODO
                          Qdissipation
			+ fvOptions(rho, hs)
		);

		Info << "Energy equation (enthalpy form): need to fix radiative heat transfer" << endl;

		HEqn.relax();

		fvOptions.constrain(HEqn);

		HEqn.solve("H");

		fvOptions.correct(hs);
	}

	//Spark
	spark.SetSpark(runTime.value(), T);
}

