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

Class
    Foam::conductivityModels::Fourier

Description
    A Fourier thermal conductivity model by RAS.

SourceFiles
    Fourier.C

\*---------------------------------------------------------------------------*/

#ifndef Fourier_H
#define Fourier_H

#include "conductivityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace conductivityModels
{

/*---------------------------------------------------------------------------*\
                           Class Newtonian Declaration
\*---------------------------------------------------------------------------*/

class Fourier
:
    public conductivityModel
{
    // Private data

        dimensionedScalar lambda0_;

        volScalarField lambda_;


public:

    //- Runtime type information
    TypeName("Fourier");


    // Constructors

        //- construct from components
        Fourier
        (
            const word& name,
            const dictionary& conductivityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    ~Fourier()
    {}


    // Member Functions

        //- Return the conductivity
        tmp<volScalarField> lambda() const
        {
            return lambda_;
        }

        //- Correct the conductivity (not appropriate, conductivity constant)
        void correct()
        {}

        //- Read transportProperties dictionary
        bool read(const dictionary& conductivityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace conductivityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
