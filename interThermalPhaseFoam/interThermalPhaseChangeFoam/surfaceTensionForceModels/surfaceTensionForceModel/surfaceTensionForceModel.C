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

\*---------------------------------------------------------------------------*/

#include "surfaceTensionForceModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfaceTensionForceModel, 0);
    defineRunTimeSelectionTable(surfaceTensionForceModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceTensionForceModel::surfaceTensionForceModel
(
    const word& name,
    const dictionary& surfaceTensionForceProperties,
    const interfaceProperties& interface,
    const volScalarField& alpha1
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            alpha1.time().constant(),
            alpha1.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    name_(name),
    surfaceTensionForceProperties_(surfaceTensionForceProperties),
    interface_(interface),
    alpha1_(alpha1)

{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::surfaceTensionForceModel::pcap() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "DummyPC",
                alpha1_.time().timeName(),
                alpha1_.mesh()
            ),
            alpha1_.mesh(),
            dimensionedScalar("0", dimensionSet(1, -1, -2, 0, 0), 0)
        )
    );
}


Foam::tmp<Foam::surfaceScalarField> 
Foam::surfaceTensionForceModel::phi_c(const surfaceScalarField& rAUf_) const
{
    return 
    tmp<surfaceScalarField>( this->Fstff() * rAUf_ * alpha1_.mesh().magSf() );
}


bool Foam::surfaceTensionForceModel::read
(
    const dictionary& surfaceTensionForceProperties
)
{
    surfaceTensionForceProperties_ = surfaceTensionForceProperties;
    
    return true;
}

// ************************************************************************* //
