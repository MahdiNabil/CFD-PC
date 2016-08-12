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

    RAS analogy to the Newtonian viscosity model

\*---------------------------------------------------------------------------*/

#include "Fourier.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace conductivityModels
{
    defineTypeNameAndDebug(Fourier, 0);
    addToRunTimeSelectionTable(conductivityModel, Fourier, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conductivityModels::Fourier::Fourier
(
    const word& name,
    const dictionary& conductivityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    conductivityModel(name, conductivityProperties, U, phi),
    lambda0_(conductivityProperties_.lookup("lambda")),
    lambda_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U_.mesh(),
        lambda0_
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::conductivityModels::Fourier::read
(
    const dictionary& conductivityProperties
)
{
    conductivityModel::read(conductivityProperties);

    conductivityProperties_.lookup("lambda") >> lambda0_;
    lambda_ = lambda0_;

    return true;
}


// ************************************************************************* //
