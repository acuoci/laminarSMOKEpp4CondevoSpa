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
    Foam::radiationModels::scatterModelOpenSMOKEpp

Description
    Base class for radiation scattering

\*---------------------------------------------------------------------------*/

#ifndef scatterModelOpenSMOKEpp_H
#define scatterModelOpenSMOKEpp_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{

/*---------------------------------------------------------------------------*\
                        Class scatterModelOpenSMOKEpp Declaration
\*---------------------------------------------------------------------------*/

class scatterModelOpenSMOKEpp
{

protected:

    // Protected data

        //- Reference to the fvMesh
        const fvMesh& mesh_;

public:

    //- Runtime type information
    TypeName("scatterModelOpenSMOKEpp");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            scatterModelOpenSMOKEpp,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& mesh
            ),
            (dict, mesh)
        );


    // Constructors

        //- Construct from components
        scatterModelOpenSMOKEpp(const dictionary& dict, const fvMesh& mesh);


    // Selector

        static autoPtr<scatterModelOpenSMOKEpp> New
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~scatterModelOpenSMOKEpp();


    // Member Functions

        //- Return scatter coefficient
        virtual tmp<volScalarField> sigmaEff() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
