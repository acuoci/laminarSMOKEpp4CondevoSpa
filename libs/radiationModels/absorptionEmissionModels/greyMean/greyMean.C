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

#include "greyMean.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "extrapolatedCalculatedFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{
    defineTypeNameAndDebug(greyMean, 0);

    addToRunTimeSelectionTable
    (
        absorptionEmissionModelOpenSMOKEpp,
        greyMean,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::greyMean::greyMean
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& modelName
)
:
    absorptionEmissionModelOpenSMOKEpp(dict, mesh),
    coeffsDict_(dict.subDict(modelName + "Coeffs")),
    speciesNames_(0),
    specieIndex_(label(0)),
    lookUpTablePtr_(),
    thermo_(mesh.lookupObject<OpenSMOKEppReactingMixtureForRadiation>("OpenSMOKEppReactingMixtureForRadiation")),
    Yj_(nSpecies_)
{
    label nFunc = 0;
    forAllConstIter(dictionary, coeffsDict_, iter)
    {
        // safety:
        if (!iter().isDict())
        {
            continue;
        }
        const word& key = iter().keyword();
        speciesNames_.insert(key, nFunc);
        const dictionary& dict = iter().dict();
        coeffs_[nFunc].initialise(dict);
        nFunc++;
    }

    if (coeffsDict_.found("lookUpTableFileName"))
    {
        const word name = coeffsDict_.lookup("lookUpTableFileName");
        if (name != "none")
        {
            lookUpTablePtr_.set
            (
                new interpolationLookUpTable<scalar>
                (
                    fileName(coeffsDict_.lookup("lookUpTableFileName")),
                    mesh.time().constant(),
                    mesh
                )
            );

            if (!mesh.foundObject<volScalarField>("ft"))
            {
                FatalErrorInFunction
                    << "specie ft is not present to use with "
                    << "lookUpTableFileName " << nl
                    << exit(FatalError);
            }
        }
    }

    // Check that all the species on the dictionary are present in the
    // look-up table and save the corresponding indices of the look-up table

    label j = 0;
    forAllConstIter(HashTable<label>, speciesNames_, iter)
    {
        if (!lookUpTablePtr_.empty())
        {
            if (lookUpTablePtr_().found(iter.key()))
            {
                label index = lookUpTablePtr_().findFieldIndex(iter.key());

                Info<< "specie: " << iter.key() << " found on look-up table "
                    << " with index: " << index << endl;

                specieIndex_[iter()] = index;
            }
            else if (mesh.foundObject<volScalarField>(iter.key()))
            {
                Yj_.set(j, &mesh.lookupObjectRef<volScalarField>(iter.key()));
                specieIndex_[iter()] = 0;
                j++;
                Info<< "specie: " << iter.key() << " is being solved" << endl;
            }
            else
            {
                FatalErrorInFunction
                    << "specie: " << iter.key()
                    << " is neither in look-up table: "
                    << lookUpTablePtr_().tableName()
                    << " nor is being solved" << nl
                    << exit(FatalError);
            }
        }
        else if (mesh.foundObject<volScalarField>(iter.key()))
        {
            Yj_.set(j, &mesh.lookupObjectRef<volScalarField>(iter.key()));
            specieIndex_[iter()] = 0;
            j++;
        }
        else
        {
            FatalErrorInFunction
                << " there is not lookup table and the specie" << nl
                << iter.key() << nl
                << " is not found " << nl
                << exit(FatalError);

        }
    }

	gas_correction_coefficient_  = coeffsDict_.lookupOrDefault<scalar>(word("gasCorrectionCoefficient"),  scalar(1.));
	Info << "Gas correction coefficient: " << gas_correction_coefficient_ << endl;

    	word soot_radiation(coeffsDict_.lookup("sootModel"));

	if (soot_radiation == "none")
		soot_planck_coefficient_ = SOOT_RADIATION_PLANCK_COEFFICIENT_NONE;
	else if (soot_radiation == "Smooke")
		soot_planck_coefficient_ = SOOT_RADIATION_PLANCK_COEFFICIENT_SMOOKE;
	else if (soot_radiation == "Kent")
		soot_planck_coefficient_ = SOOT_RADIATION_PLANCK_COEFFICIENT_KENT;
	else if (soot_radiation == "Sazhin")
		soot_planck_coefficient_ = SOOT_RADIATION_PLANCK_COEFFICIENT_SAZHIN;
	else
	{
		FatalErrorIn( "Foam::radiation::greyMeanAbsorptionEmission::greyMeanAbsorptionEmission")
	    		<< "Wrong sootModel. Available models: none | Smooke | Kent | Sazhin " << abort(FatalError);
	}

	soot_correction_coefficient_ = coeffsDict_.lookupOrDefault<scalar>(word("sootCorrectionCoefficient"), scalar(1.));
	Info << "Soot radiation model: " << soot_radiation << endl;
	Info << " - Correction coeff.: " << soot_correction_coefficient_ << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::greyMean::~greyMean()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::greyMean::aCont
(
    const label bandI
) const
{
    const OpenSMOKEppReactingMixtureForRadiation& mixture = dynamic_cast<const OpenSMOKEppReactingMixtureForRadiation&>(thermo_);

    const volScalarField& T = thermo_.T();
    const volScalarField& p = thermo_.p();


    tmp<volScalarField> ta
    (
        volScalarField::New
        (
            "aCont" + name(bandI),
            mesh(),
            dimensionedScalar(dimless/dimLength, 0),
            extrapolatedCalculatedFvPatchVectorField::typeName
        )
    );

    scalarField& a = ta.ref().primitiveFieldRef();

    forAll(a, celli)
    {
        forAllConstIter(HashTable<label>, speciesNames_, iter)
        {
            label n = iter();
            scalar Xipi = 0.0;
            if (specieIndex_[n] != 0)
            {
                // Specie found in the lookUpTable.
                const volScalarField& ft =
                    mesh_.lookupObject<volScalarField>("ft");

                const List<scalar>& Ynft = lookUpTablePtr_().lookUp(ft[celli]);
                // moles x pressure [atm]
                Xipi = Ynft[specieIndex_[n]]*paToAtm(p[celli]);
            }
            else
            {
                scalar invWt = 0.0;
                forAll(mixture.Y(), s)
                {
                    invWt += mixture.Y(s)[celli]/mixture.Wi(s);
                }

                label index = mixture.species(iter.key());
                scalar Xk = mixture.Y(index)[celli]/(mixture.Wi(index)*invWt);

                Xipi = Xk*paToAtm(p[celli]);
            }

            const absorptionCoeffs::coeffArray& b = coeffs_[n].coeffs(T[celli]);

            scalar Ti = T[celli];
            // negative temperature exponents
            if (coeffs_[n].invTemp())
            {
                Ti = 1.0/T[celli];
            }
            a[celli] +=
                Xipi
               *(
                    ((((b[5]*Ti + b[4])*Ti + b[3])*Ti + b[2])*Ti + b[1])*Ti
                  + b[0]
                );

	    a[celli] *= gas_correction_coefficient_;
        }
    }

    // Soot contribution
    if (soot_planck_coefficient_ != SOOT_RADIATION_PLANCK_COEFFICIENT_NONE)
    {
	const scalarField& fvsoot = T.db().objectRegistry::lookupObject<volScalarField>("soot:fv").internalField();

	forAll(a, celli)
    	{
		const scalar fvi = fvsoot[celli];
		const scalar Ti = T[celli];

		double aSoot = 0.;
															
		if (soot_planck_coefficient_ == SOOT_RADIATION_PLANCK_COEFFICIENT_SMOOKE)
			aSoot = 1307.*fvi*Ti;						// [1/m]	(Smooke et al. Combustion and Flame 2009)
		else if (soot_planck_coefficient_ == SOOT_RADIATION_PLANCK_COEFFICIENT_KENT)
			aSoot = 2262.*fvi*Ti;						// [1/m]	(Kent al. Combustion and Flame 1990)
		else if (soot_planck_coefficient_ == SOOT_RADIATION_PLANCK_COEFFICIENT_SAZHIN)
			aSoot= 1232.*fvi*(1. + 4.8e-4*(Ti - 2000.));			// [1/m]	(Sazhin, Fluent 1994)	

		a[celli] += aSoot*soot_correction_coefficient_;
        }
    }

    ta.ref().correctBoundaryConditions();
    return ta;
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::greyMean::eCont
(
    const label bandI
) const
{
    return aCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::greyMean::ECont
(
    const label bandI
) const
{
    return absorptionEmissionModelOpenSMOKEpp::ECont(bandI);
}


// ************************************************************************* //
