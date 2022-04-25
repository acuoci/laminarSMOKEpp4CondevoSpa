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

Class
    Foam::radiationModels::opticallyThin

Description
    Optically thin model

SourceFiles
    opticallyThin.C

\*---------------------------------------------------------------------------*/

#ifndef opticallyThin_H
#define opticallyThin_H

#include "radiationModelOpenSMOKE++.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{

/*---------------------------------------------------------------------------*\
                         Class opticallyThin Declaration
\*---------------------------------------------------------------------------*/

class opticallyThin
:
    public radiationModelOpenSMOKEpp
{
public:

    //- Runtime type information
    TypeName("opticallyThin");


    // Constructors

        //- Construct from components
        opticallyThin(const volScalarField& T);

        //- Construct from components
        opticallyThin(const dictionary& dict, const volScalarField& T);

        //- Disallow default bitwise copy construction
        opticallyThin(const opticallyThin&) = delete;


    //- Destructor
    virtual ~opticallyThin();


    // Member Functions

        // Edit

            //- Solve radiation equation(s)
            void calculate();

            //- Read radiationProperties dictionary
            bool read();

            //- Source term component (for power of T^4)
            tmp<volScalarField> Rp() const;

            //- Source term component (constant)
            tmp<volScalarField::Internal> Ru() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const opticallyThin&) = delete;

private:

	//- Total radiative heat flux [W/m2]
        volScalarField Qr_;

        //- Absorption coefficient
        volScalarField a_;

        //- Emission coefficient
        volScalarField e_;

        //- Emission contribution
        volScalarField E_;

	//- Ambient temperature
	scalar ambientTemperature_;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
