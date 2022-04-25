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

// Discretization scheme: convection
tmp<fv::convectionScheme<scalar> > mvConvection
(
	fv::convectionScheme<scalar>::New ( mesh, fields, phi, mesh.divScheme("div(phi,Yi_hs)") )
);

// Discretization scheme: convection via correction velocity
tmp<fv::convectionScheme<scalar> > mvConvectionCorrection
(
	fv::convectionScheme<scalar>::New ( mesh, fields, phic, mesh.divScheme("div(phic,Yi)") )
);

// Governing equation in case the diffusion fluxes are written as j=-rho*Dmixi*nabla(Yi)
if (mixture.diffusion_via_mole_fraction_gradient() == false)
{
	volScalarField Yt(0.0*Y[0]);

	forAll(Y, i)
	{
		if (mixture.active(i))
		{
			volScalarField& Yi = Y[i];
			volScalarField& Dmixi = Dmix[i];

			fvScalarMatrix YiEqn
			(
				fvm::ddt(rho, Yi)
				+ mvConvection->fvmDiv(phi, Yi)
				+ mvConvectionCorrection->fvmDiv(phic, Yi)
				- fvm::laplacian(rho*Dmixi, Yi)
				==
				fvOptions(rho, Yi)
			);

			// Add Soret effect
			if (mixture.soret_effect() == true)
			{ 
				if (mixture.soret_list(i) == true)
					YiEqn -= fvc::laplacian(rho*Dsoret[mixture.soret_index(i)]/T, T, "laplacian(teta,Yi)");
			}

			YiEqn.relax();

			fvOptions.constrain(YiEqn);

			YiEqn.solve("Yi");

			fvOptions.correct(Yi);

			Yi.max(0.0);
			Yt += Yi;
		}
	}

	if (mixture.closure_on_inert() == true)
	{
		Y[mixture.inertIndex()] = scalar(1.0) - Yt;
		Y[mixture.inertIndex()].max(0.0);
	}
	else
	{
		eps = 1./Yt-1.;
		forAll(Y, i)
		{
			volScalarField& Yi = Y[i];
			Yi *= (1.+eps);
		}

		Info<< " * Min/Max unbalance   = " << -max(-eps).value() << ", " << max(eps).value() << endl;
	}
}

// Governing equation in case the diffusion fluxes are written as j=-rho*(Wi/W)*Dmixi*nabla(Xi)
else
{
	volScalarField Yt(0.0*Y[0]);

	forAll(Y, i)
	{
		if (mixture.active(i))
		{
			volScalarField& Yi = Y[i];
			volScalarField& Dmixi = Dmix[i];
			dimensionedScalar MWi = mixture.mw_species(i);

			fvScalarMatrix YiEqn
			(
				fvm::ddt(rho, Yi)
				+ mvConvection->fvmDiv(phi, Yi)
				+ mvConvectionCorrection->fvmDiv(phic, Yi)
				- fvm::laplacian(rho*Dmixi-rho*Dmixi*MW*Yi/MWi, Yi)
				==
				- mwDiffusionCorrection[i]
				+ fvOptions(rho, Yi)
			);

			// Add Soret effect
			if (mixture.soret_effect() == true)
			{ 
				if (mixture.soret_list(i) == true)
					YiEqn -= fvc::laplacian(rho*Dsoret[mixture.soret_index(i)]/T, T, "laplacian(teta,Yi)");
			}

			YiEqn.relax();

			fvOptions.constrain(YiEqn);

			YiEqn.solve("Yi");

			fvOptions.correct(Yi);

			Yi.max(0.0);
			Yt += Yi;
		}
	}

	if (mixture.closure_on_inert() == true)
	{
		Y[mixture.inertIndex()] = scalar(1.0) - Yt;
		Y[mixture.inertIndex()].max(0.0);
	}
	else
	{
		eps = 1./Yt-1.;
		forAll(Y, i)
		{
			volScalarField& Yi = Y[i];
			Yi *= (1.+eps);
		}

		Info<< " * Min/Max unbalance   = " << -max(-eps).value() << ", " << max(eps).value() << endl;
	}
}


