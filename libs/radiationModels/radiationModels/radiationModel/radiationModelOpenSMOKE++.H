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
    Foam::radiationModelOpenSMOKEpp

Description
    Top level model for radiation modelling

SourceFiles
    radiationModelOpenSMOKE++.C
    radiationModelNewOpenSMOKE++.C

\*---------------------------------------------------------------------------*/

#ifndef radiationModelOpenSMOKEpp_H
#define radiationModelOpenSMOKEpp_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "fvMatricesFwd.H"
#include "Switch.H"
#include "OpenSMOKEppReactingMixtureForRadiation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class OpenSMOKEppReactingMixtureForRadiation;
class fvMesh;

namespace radiationModels
{
class absorptionEmissionModelOpenSMOKEpp;
class scatterModelOpenSMOKEpp;
class sootModelOpenSMOKEpp;
}

/*---------------------------------------------------------------------------*\
                       Class radiationModelOpenSMOKEpp Declaration
\*---------------------------------------------------------------------------*/

class radiationModelOpenSMOKEpp
:
    public IOdictionary
{
protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Reference to the time database
        const Time& time_;

        //- Reference to the temperature field
        const volScalarField& T_;

        //- Radiation model dictionary
        dictionary coeffs_;

        //- Radiation solver frequency - number flow solver iterations per
        //  radiation solver iteration
        label solverFreq_;

        //- Flag to enable radiation model to be evaluated on first iteration
        bool firstIter_;


        // References to the radiation sub-models

            //- Absorption/emission model
            autoPtr<radiationModels::absorptionEmissionModelOpenSMOKEpp>
                absorptionEmission_;

            //- Scatter model
            autoPtr<radiationModels::scatterModelOpenSMOKEpp> scatter_;

            //- Soot model
            autoPtr<radiationModels::sootModelOpenSMOKEpp> soot_;


private:

    // Private Member Functions

        //- Create IO object if dictionary is present
        IOobject createIOobject(const fvMesh& mesh) const;

        //- Initialise
        void initialise();


public:

    //- Runtime type information
    TypeName("radiationModelOpenSMOKEpp");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            radiationModelOpenSMOKEpp,
            T,
            (
                const volScalarField& T
            ),
            (T)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            radiationModelOpenSMOKEpp,
            dictionary,
            (
                const dictionary& dict,
                const volScalarField& T
            ),
            (dict, T)
        );


    // Constructors

        //- Null constructor
        radiationModelOpenSMOKEpp(const volScalarField& T);

        //- Construct from components
        radiationModelOpenSMOKEpp(const word& type, const volScalarField& T);

        //- Construct from components
        radiationModelOpenSMOKEpp
        (
            const word& type,
            const dictionary& dict,
            const volScalarField& T
        );

        //- Disallow default bitwise copy construction
        radiationModelOpenSMOKEpp(const radiationModelOpenSMOKEpp&) = delete;


    // Selectors

        //- Return a reference to the selected radiation model
        static autoPtr<radiationModelOpenSMOKEpp> New(const volScalarField& T);

        //- Return a reference to the selected radiation model
        static autoPtr<radiationModelOpenSMOKEpp> New
        (
            const dictionary& dict,
            const volScalarField& T
        );


    //- Destructor
    virtual ~radiationModelOpenSMOKEpp();


    // Member Functions

        // Edit

            //- Main update/correction routine
            virtual void correct();

            //- Solve radiation equation(s)
            virtual void calculate() = 0;

            //- Read radiationProperties dictionary
            virtual bool read() = 0;


        // Access

            //- Source term component (for power of T^4)
            virtual tmp<volScalarField> Rp() const = 0;

            //- Source term component (constant)
            virtual tmp<volScalarField::Internal> Ru() const = 0;

            //- Energy source term
            virtual tmp<fvScalarMatrix> Sh
            (
                const OpenSMOKEppReactingMixtureForRadiation& thermo,
                const volScalarField& h
            ) const;

            //- Temperature source term
            virtual tmp<fvScalarMatrix> ST
            (
                const dimensionedScalar& rhoCp,
                volScalarField& T
            ) const;

            //- Access to absorptionEmission model
            const radiationModels::absorptionEmissionModelOpenSMOKEpp&
                absorptionEmission() const;

            //- Access to soot Model
            const radiationModels::sootModelOpenSMOKEpp& soot() const;

	    //- Temperature source term Q=div(qrad)
            virtual tmp<fvScalarMatrix> divq(volScalarField& T) const;

	    //- Qloss [W/m3]
	    virtual void Qloss( volScalarField& T, volScalarField& Qloss);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const radiationModelOpenSMOKEpp&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define addToRadiationRunTimeSelectionTables(model)                            \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        radiationModelOpenSMOKEpp,                                             \
        model,                                                                 \
        dictionary                                                             \
    );                                                                         \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        radiationModelOpenSMOKEpp,                                             \
        model,                                                                 \
        T                                                                      \
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
