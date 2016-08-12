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

    Added by RAS for thermal transport, basically a clone of the viscosity model

\*---------------------------------------------------------------------------*/

#include "conductivityModel.H"
#include "volFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(conductivityModel, 0);
    defineRunTimeSelectionTable(conductivityModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conductivityModel::conductivityModel
(
    const word& name,
    const dictionary& conductivityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    name_(name),
    conductivityProperties_(conductivityProperties),
    U_(U),
    phi_(phi)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

//- I don't think this will ever be used, but W/E
Foam::tmp<Foam::volScalarField> Foam::conductivityModel::strainRate() const
{
    return sqrt(2.0)*mag(symm(fvc::grad(U_)));
}


bool Foam::conductivityModel::read(const dictionary& conductivityProperties)
{
    conductivityProperties_ = conductivityProperties;

    return true;
}


// ************************************************************************* //
