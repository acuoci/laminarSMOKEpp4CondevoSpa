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

#include "backwardDiffusionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::backwardDiffusionFvPatchScalarField::backwardDiffusionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedUserDefinedFvPatchScalarField(p, iF)
{
	alfa()    = 0.0;
	beta()    = 0.0;
        eta()     = 0.0;
	omega0()  = 0.0;
	rho0()    = 0.0;
        epsilon() = 0.0;

	nameInternal_ = internalField().name();
}


Foam::backwardDiffusionFvPatchScalarField::backwardDiffusionFvPatchScalarField
(
    const backwardDiffusionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedUserDefinedFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::backwardDiffusionFvPatchScalarField::backwardDiffusionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedUserDefinedFvPatchScalarField(p, iF)
{
    // Name of field
    nameInternal_ = internalField().name();
    
    // Set the nominal value
    omega0() = scalarField("omega0", dict, p.size());
    rho0()   = scalarField("rho0", dict, p.size());

    // Fixed value condition is forced
    alfa() = 1000.;
    eta()  = 1.;

    // Calculating epsilon
    epsilon() = 0;

    // Calculating beta
    const double Dmix = 1e-10;
    beta() = Dmix*this->patch().deltaCoeffs();

    // Read value if available
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        evaluate();
    }
}


Foam::backwardDiffusionFvPatchScalarField::backwardDiffusionFvPatchScalarField
(
    const backwardDiffusionFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedUserDefinedFvPatchScalarField(tppsf, iF)
{}


void Foam::backwardDiffusionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

    // Calculating alfa
    const volVectorField& U = db().lookupObject<volVectorField>("U");
    tmp<vectorField> n = patch().nf();
    alfa() = -(n & U.boundaryField()[patchi]);
    
    // Calculating eta
    const volScalarField& rho = db().lookupObject<volScalarField>("rho");
    eta() = rho0() / rho.boundaryField()[patchi];

    nameInternal_ = internalField().name();

    bool soretEffect  = db().foundObject<volScalarField>("mix:Dsoret:" + nameInternal_);

    // Calculating epsilon
    if (soretEffect == true)
    {
        const volScalarField& Dsoret = db().lookupObject<volScalarField>("mix:Dsoret:" + nameInternal_);
    	const volScalarField& T = db().lookupObject<volScalarField>("T");
    	epsilon() = -T.boundaryField()[patchi].snGrad() / T.boundaryField()[patchi] * Dsoret.boundaryField()[patchi];
    }
    else
    { 
	epsilon() = 0.;
    }

    // Calculating beta
    nameInternal_ = internalField().name();

    const volScalarField& Dmix = db().lookupObject<volScalarField>("mix:Dmix:" + nameInternal_);
    beta() = Dmix.boundaryField()[patchi]*this->patch().deltaCoeffs();

    if (debug)
    {
    }

    mixedUserDefinedFvPatchScalarField::updateCoeffs();
}

void Foam::backwardDiffusionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "omega0", omega0());
    writeEntry(os, "rho0", rho0());
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        backwardDiffusionFvPatchScalarField
    );
}

// ************************************************************************* //
