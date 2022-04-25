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

// Create const references

// Compressibility [s2/m2]
const volScalarField& psi = mixture.psi();

// Dynamic viscosity [kg/m/s]
const volScalarField& mu = mixture.mu();

// Constant pressure specific heat [J/kg/K]
const volScalarField& Cp = mixture.Cp();

// Constant volume specific heat [J/kg/K]
const volScalarField& Cpv = mixture.Cpv();

// Thermal conductivity [W/m/K]
const volScalarField& lambda = mixture.lambda();

// Thermal diffusivity [m2/s]
const volScalarField& alpha = mixture.alpha();

// Molecular weight [kg/kmol]
const volScalarField& MW = mixture.MW();

// Mass diffusion correction flux [kg/m2/s]
volVectorField& Jc = mixture.Jc();

// Correction on the diffusion fluxes associated to the molecular weight
const PtrList<volScalarField>& mwDiffusionCorrection = mixture.mwDiffusionCorrection();

// Mass diffusion correction flow [kg/s]
surfaceScalarField& phic = mixture.phic();

// Create non-const references

// Temperature [K]
volScalarField& T = mixture.T();


// Reading/Writing options
IOobject::readOption  tauReadOption = IOobject::NO_READ;
IOobject::writeOption tauWriteOption = IOobject::NO_WRITE;


// Residence time (options)
if (mixture.solveForTauEquation() == true)
{
	tauReadOption = IOobject::MUST_READ;
	tauWriteOption = IOobject::AUTO_WRITE;
}

//Residence time field
Info<< "Reading tau (residence time) field\n" << endl;
volScalarField tau
(
	IOobject
	(
		"tau",
		runTime.timeName(),
		mesh,
		tauReadOption,
		tauWriteOption
	),
	mesh
);

// Passive scalar (options)
IOobject::readOption  csiReadOption = IOobject::NO_READ;
IOobject::writeOption csiWriteOption = IOobject::NO_WRITE;

if (mixture.solveForCsiEquation() == true)
{
	csiReadOption = IOobject::MUST_READ;
	csiWriteOption = IOobject::AUTO_WRITE;
}


// Passive scalar field
Info<< "Reading passive scalar csi field\n" << endl;
volScalarField csi
(
	IOobject
	(
		"csi",
		runTime.timeName(),
		mesh,
		csiReadOption,
		csiWriteOption
	),
	mesh
);

// Create the reacting mixture object as interface with the radiation model
#include "createOpenSMOKEppReactingMixtureForRadiation.H"

// Create the radiation model
#include "createRadiationModel.H"

// Radiative heat loss [W/m3]
volScalarField& Qrad = mixture.Qrad();


