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

#include "gradientEnthalpyFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "OpenSMOKEppReactingMixture.H"

Foam::gradientEnthalpyFvPatchScalarField::gradientEnthalpyFvPatchScalarField
(
	const fvPatch& p, 
	const DimensionedField<scalar, volMesh>& iF
)
:
fixedGradientFvPatchScalarField(p, iF)
{}

Foam::gradientEnthalpyFvPatchScalarField::gradientEnthalpyFvPatchScalarField
(
	const gradientEnthalpyFvPatchScalarField& ptf,
	const fvPatch& p,
	const DimensionedField<scalar, volMesh>& iF,
	const fvPatchFieldMapper& mapper
)
:
fixedGradientFvPatchScalarField(ptf, p, iF, mapper)
{}

Foam::gradientEnthalpyFvPatchScalarField::
gradientEnthalpyFvPatchScalarField
(
	const fvPatch& p,
	const DimensionedField<scalar, volMesh>& iF,
	const dictionary& dict
)
:
fixedGradientFvPatchScalarField(p, iF, dict)
{}

Foam::gradientEnthalpyFvPatchScalarField::
gradientEnthalpyFvPatchScalarField
(
	const gradientEnthalpyFvPatchScalarField& tppsf,
	const DimensionedField<scalar, volMesh>& iF
)
:
fixedGradientFvPatchScalarField(tppsf, iF)
{}


void Foam::gradientEnthalpyFvPatchScalarField::updateCoeffs()
{
	if (updated())
	{
		return;
	}

	const OpenSMOKEppReactingMixture& mixture = this->db().objectRegistry::lookupObject<OpenSMOKEppReactingMixture> ("SimulationOptions");
	const label patchi = patch().index();

	const volScalarField& p = this->db().objectRegistry::lookupObject<volScalarField> ("p");
	const volScalarField& T = this->db().objectRegistry::lookupObject<volScalarField> ("T");

	const fvPatchScalarField& pw = p.boundaryField()[patchi];
	fvPatchScalarField& Tw = const_cast<fvPatchScalarField&>(T.boundaryField()[patchi]);

	Tw.evaluate();
	
	gradient() = 	mixture.Cpv(pw, Tw, patchi)*Tw.snGrad() +
			patch().deltaCoeffs()*
			(
			    mixture.hs(pw, Tw, patchi)
			  - mixture.hs(pw, Tw, patch().faceCells())
			);
	
	fixedGradientFvPatchScalarField::updateCoeffs();
}

namespace Foam
{
	makePatchTypeField
	(
		fvPatchScalarField,
		gradientEnthalpyFvPatchScalarField
	);
}

