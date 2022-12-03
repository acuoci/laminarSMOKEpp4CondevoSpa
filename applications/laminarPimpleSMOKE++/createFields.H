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

//- Detect spark
SparkModel spark(mesh);

// Create the reacting multicomponent mixture
Foam::autoPtr<Foam::OpenSMOKEppReactingMixture> pMixture = OpenSMOKEppReactingMixture::New(mesh);
OpenSMOKEppReactingMixture& mixture = pMixture();

// Set references
PtrList<volScalarField>& Y = mixture.Y();
PtrList<volScalarField>& X = mixture.X();
PtrList<volScalarField>& Dmix = mixture.Dmix();
PtrList<volScalarField>& Dsoret = mixture.Dsoret();
PtrList<volVectorField>& J = mixture.J();
PtrList<volScalarField>& CpSpecies = mixture.CpSpecies();
PtrList<volScalarField>& hsSpecies = mixture.hsSpecies();
volScalarField& p = mixture.p();
volScalarField& hs = mixture.Hs();


// Creating density field as the product between p and psi
Info<< " * Creating density field rho ..." << endl;
#if STEADYSTATE == 1
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mixture.rho()
);
#else
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh
    ),
    mixture.rho()
);
#endif

Info << "   RhoMin = " << min(rho).value() << " RhoMax" << max(rho).value() << endl;


// Reading the velocity field
Info<< " * Reading velocity field U ..." << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "   Umin: " << min(U).value() << " Umax: " << max(U).value() << endl;

Info<< "Creating mass umbalances\n" << endl;
volScalarField eps
(
    IOobject
    (
        "mix:eps",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("mix:eps", dimensionSet(0, 0, 0, 0, 0), 0.)
);

volScalarField Qdissipation
(
    IOobject
    (
	"mix:Qdiss",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("mix:Qdiss", dimensionSet(1, -1, -3, 0, 0), 0. )
);

// Create mass flux field
#include "compressibleCreatePhi.H"

// Set flux for mesh
mesh.setFluxRequired(p.name());

#if STEADYSTATE != 1
Info<< " * Creating field dpdt ..." << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0.0)
);
#endif

// Creating the kinetic energy field
Info<< " * Creating kinetic energy field K ..." << endl;
volScalarField K("K", 0.5*magSqr(U));

// Set interpolation schemes for species mass fractions
multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;
forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(hs);

// Spark
spark.Setup(mesh);

// Calculate the initial mass in the domain
dimensionedScalar initialMass = fvc::domainIntegrate(rho);

// Include multi-reference frame
#include "createMRF.H"

// Include additional optional and constraint fields
#include "createFvModels.H"
#include "createFvConstraints.H"

// Check for the gravity
#include "readGravitationalAcceleration.H"
if ( (mag(g)).value() != 0.)
{
	Info << endl;
	Info << "------------------------------------------------------------------" << endl;
	Info << " ERROR:   gravity is turned on                                    " << endl;
	Info << "          Please consider that this solver is not conceived       " << endl;
	Info << "          for bouyancy driven flows. If you think gravity         " << endl;
	Info << "          plays a non-negligible role, consider the possibility   " << endl;
	Info << "          to use the laminarBouyantPimpleSMOKE solver (for        " << endl;
	Info << "          unsteady flows) or the laminarBouyantSimpleSMOKE solver " << endl;
	Info << "          (for steady-state flows)                                " << endl;
	Info << "------------------------------------------------------------------" << endl;
	Info << endl;
	abort();
}

