/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "radiationModelOpenSMOKE++.H"
#include "absorptionEmissionModelOpenSMOKE++.H"
#include "scatterModelOpenSMOKE++.H"
#include "sootModelOpenSMOKE++.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(radiationModelOpenSMOKEpp, 0);
    defineRunTimeSelectionTable(radiationModelOpenSMOKEpp, T);
    defineRunTimeSelectionTable(radiationModelOpenSMOKEpp, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::IOobject Foam::radiationModelOpenSMOKEpp::createIOobject(const fvMesh& mesh) const
{
    IOobject io
    (
        "radiationProperties",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.typeHeaderOk<IOdictionary>(true))
    {
        io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
        return io;
    }
    else
    {
        io.readOpt() = IOobject::NO_READ;
        return io;
    }
}


void Foam::radiationModelOpenSMOKEpp::initialise()
{
    solverFreq_ = max(1, lookupOrDefault<label>("solverFreq", 1));

    absorptionEmission_.reset
    (
        radiationModels::absorptionEmissionModelOpenSMOKEpp::New(*this, mesh_).ptr()
    );

    scatter_.reset(radiationModels::scatterModelOpenSMOKEpp::New(*this, mesh_).ptr());

    soot_.reset(radiationModels::sootModelOpenSMOKEpp::New(*this, mesh_).ptr());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModelOpenSMOKEpp::radiationModelOpenSMOKEpp(const volScalarField& T)
:
    IOdictionary
    (
        IOobject
        (
            "radiationProperties",
            T.time().constant(),
            T.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(T.mesh()),
    time_(T.time()),
    T_(T),
    coeffs_(dictionary::null),
    solverFreq_(0),
    firstIter_(true),
    absorptionEmission_(nullptr),
    scatter_(nullptr),
    soot_(nullptr)
{}


Foam::radiationModelOpenSMOKEpp::radiationModelOpenSMOKEpp(const word& type, const volScalarField& T)
:
    IOdictionary(createIOobject(T.mesh())),
    mesh_(T.mesh()),
    time_(T.time()),
    T_(T),
    coeffs_(subOrEmptyDict(type + "Coeffs")),
    solverFreq_(1),
    firstIter_(true),
    absorptionEmission_(nullptr),
    scatter_(nullptr),
    soot_(nullptr)
{
    initialise();
}


Foam::radiationModelOpenSMOKEpp::radiationModelOpenSMOKEpp
(
    const word& type,
    const dictionary& dict,
    const volScalarField& T
)
:
    IOdictionary
    (
        IOobject
        (
            "radiationProperties",
            T.time().constant(),
            T.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    mesh_(T.mesh()),
    time_(T.time()),
    T_(T),
    coeffs_(subOrEmptyDict(type + "Coeffs")),
    solverFreq_(1),
    firstIter_(true),
    absorptionEmission_(nullptr),
    scatter_(nullptr),
    soot_(nullptr)
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::radiationModelOpenSMOKEpp::~radiationModelOpenSMOKEpp()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiationModelOpenSMOKEpp::correct()
{
    if (firstIter_ || (time_.timeIndex() % solverFreq_ == 0))
    {
        calculate();
        firstIter_ = false;
    }

    if (!soot_.empty())
    {
        soot_->correct();
    }
}


bool Foam::radiationModelOpenSMOKEpp::read()
{
    if (regIOobject::read())
    {
        coeffs_ = subOrEmptyDict(type() + "Coeffs");

        solverFreq_ = lookupOrDefault<label>("solverFreq", 1);
        solverFreq_ = max(1, solverFreq_);

        return true;
    }
    else
    {
        return false;
    }
}


Foam::tmp<Foam::fvScalarMatrix> Foam::radiationModelOpenSMOKEpp::Sh
(
    const OpenSMOKEppReactingMixtureForRadiation& thermo,
    const volScalarField& h
) const
{
    const volScalarField Cpv(thermo.Cpv());
    const volScalarField T3(pow3(T_));

    return
    (
        Ru()
      - fvm::Sp(4.0*Rp()*T3/Cpv, h)
      - Rp()*T3*(T_ - 4.0*h/Cpv)
    );
}


Foam::tmp<Foam::fvScalarMatrix> Foam::radiationModelOpenSMOKEpp::ST
(
    const dimensionedScalar& rhoCp,
    volScalarField& T
) const
{
    return
    (
        Ru()/rhoCp
      - fvm::Sp(Rp()*pow3(T)/rhoCp, T)
    );
}

Foam::tmp<Foam::fvScalarMatrix> Foam::radiationModelOpenSMOKEpp::divq( volScalarField& T ) const
{
    return
    (
        Ru() - fvm::Sp(Rp()*pow3(T), T)
    );
}

void Foam::radiationModelOpenSMOKEpp::Qloss( volScalarField& T, volScalarField& Qloss)
{
	// TODO: Ru is an external source
	//       Thus it can be neglected in the expression below, unless the user
	//	 explicitly defined it
        Qloss = Rp()*pow4(T);	// - Ru();
}


const Foam::radiationModels::absorptionEmissionModelOpenSMOKEpp&
Foam::radiationModelOpenSMOKEpp::absorptionEmission() const
{
    if (!absorptionEmission_.valid())
    {
        FatalErrorInFunction
            << "Requested radiation absorptionEmission model, but model is "
            << "not active" << abort(FatalError);
    }

    return absorptionEmission_();
}


const Foam::radiationModels::sootModelOpenSMOKEpp& Foam::radiationModelOpenSMOKEpp::soot() const
{
    if (!soot_.valid())
    {
        FatalErrorInFunction
            << "Requested radiation sootModel model, but model is "
            << "not active" << abort(FatalError);
    }

    return soot_();
}


// ************************************************************************* //
