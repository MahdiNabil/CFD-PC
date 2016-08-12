/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 Alex Rattner
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

Namespace
    Foam::conductivityModels

Description
    A namespace for various incompressible conductivityModel implementations.
    Basically a RAS-made clone of viscosityModel

Class
    Foam::conductivityModel

Description
    An abstract base class for incompressible conductivityModels.

    The strain rate is defined by:

        mag(symm(grad(U)))

    but will probably never be used


SourceFiles
    conductivityModel.C
    conductivityModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef conductivityModel_H
#define conductivityModel_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dimensionedScalar.H"
#include "tmp.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class conductivityModel Declaration
\*---------------------------------------------------------------------------*/

class conductivityModel
{

protected:

    // Protected data

        word name_;
        dictionary conductivityProperties_;

        const volVectorField& U_;
        const surfaceScalarField& phi_;


    // Private Member Functions

        //- Disallow copy construct
        conductivityModel(const conductivityModel&);

        //- Disallow default bitwise assignment
        void operator=(const conductivityModel&);


public:

    //- Runtime type information
    TypeName("conductivityModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            conductivityModel,
            dictionary,
            (
                const word& name,
                const dictionary& conductivityProperties,
                const volVectorField& U,
                const surfaceScalarField& phi
            ),
            (name, conductivityProperties, U, phi)
        );


    // Selectors

        //- Return a reference to the selected conductivity model
        static autoPtr<conductivityModel> New
        (
            const word& name,
            const dictionary& conductivityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Constructors

        //- Construct from components
        conductivityModel
        (
            const word& name,
            const dictionary& conductivityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~conductivityModel()
    {}


    // Member Functions

        //- Return the phase transport properties dictionary
        const dictionary& conductivityProperties() const
        {
            return conductivityProperties_;
        }

        //- Return the strain rate
        tmp<volScalarField> strainRate() const;

        //- Return the laminar conductivity
        virtual tmp<volScalarField> lambda() const = 0;

        //- Correct the laminar conductivity
        virtual void correct() = 0;

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& conductivity) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
