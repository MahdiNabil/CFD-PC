/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 Alex Rattner
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

#include "SST.H"
#include "addToRunTimeSelectionTable.H"
#include "EvalSSF.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceTensionForceModels
{
    defineTypeNameAndDebug(SST, 0);
    addToRunTimeSelectionTable(surfaceTensionForceModel, SST, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceTensionForceModels::SST::SST
(
	const word& name,
	const dictionary& surfaceTensionForceProperties,
	const interfaceProperties& interface,
	const volScalarField& alpha1
)
:
    surfaceTensionForceModel(name, surfaceTensionForceProperties, interface, alpha1),
	mesh_(alpha1.mesh()),
	Fstffv
    (
        IOobject
        (
            "SurfaceTensionForce",
            alpha1_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
		mesh_,
		dimensionedScalar( "dummy", dimensionSet(1,-2,-2,0,0,0,0), 0 )
    );

	//Sharp surface force (capillary forces on cell faces)
    volVectorField fc
    (
        IOobject
        (
            "fc",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("fc0", dimMass/(dimLength*dimLength*dimTime*dimTime), vector(0,0,0))
    );

	//Sharp surface force (capillary forces on cell faces)
    surfaceScalarField fcf
    (
        IOobject
        (
            "fcf",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("fc0", dimMass/(dimLength*dimLength*dimTime*dimTime), 0)
    );

	//Sharp surface force (capillary forces on cell faces)
    surfaceScalarField fcf_filter
    (
        IOobject
        (
            "fcf_filter",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("fc0", dimMass/(dimLength*dimLength*dimTime*dimTime), 0)
    );

{
	correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::surfaceTensionForceModels::SST::correct()
{
	Fstffv = fvc::interpolate(interface_.sigmaK())*fvc::snGrad(alpha1_);
}

bool Foam::surfaceTensionForceModels::SST::read(const dictionary& surfaceTensionForceProperties)
{
	surfaceTensionForceModel::read(surfaceTensionForceProperties);

	return true;
}


// ************************************************************************* //
