/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
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
    Foam::radiationModels::sootModelOpenSMOKEpp

Description
    Base class for soot models

\*---------------------------------------------------------------------------*/

#ifndef sootModelOpenSMOKEpp_H
#define sootModelOpenSMOKEpp_H

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
                            Class sootModelOpenSMOKEpp Declaration
\*---------------------------------------------------------------------------*/

class sootModelOpenSMOKEpp
{
protected:

    // Protected data

        //- Radiation model dictionary
        const dictionary dict_;

        //- Reference to the fvMesh
        const fvMesh& mesh_;


public:

    //- Runtime type information
    TypeName("sootModelOpenSMOKEpp");

    //- Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            sootModelOpenSMOKEpp,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& mesh,
                const word& modelType
            ),
            (dict, mesh, modelType)
        );


    // Constructors

        //- Construct from components
        sootModelOpenSMOKEpp
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& modelType
        );


    //- Selector
    static autoPtr<sootModelOpenSMOKEpp> New
    (
        const dictionary& dict,
        const fvMesh& mesh
    );


    //- Destructor
    virtual ~sootModelOpenSMOKEpp();


    // Member Functions

        // Access

            //- Reference to the mesh
            inline const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Reference to the dictionary
            inline const dictionary& dict() const
            {
                return dict_;
            }


        // Edit

            //- Main update/correction routine
            virtual void correct() = 0;


        // Access

            //- Return const reference to soot
            virtual const volScalarField& soot() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
