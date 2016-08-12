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

#include "thermalPhaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::thermalPhaseChangeModel>
Foam::thermalPhaseChangeModel::New
(
    const word& name,
    const dictionary& thermalPhaseChangeProperties,
    const twoPhaseThermalMixture& twoPhaseProperties,
    const volScalarField& T,
    const volScalarField& alpha1
)
{
    const word modelType(thermalPhaseChangeProperties.lookup("model"));

    Info<< "Selecting phase change model: " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "thermalPhaseChangeModel::New"
        )   << "Unknown thermalPhaseChangeModel type "
            << modelType << endl << endl
            << "Valid  thermalPhaseChangeModels are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<thermalPhaseChangeModel>
    (
        cstrIter()
        (
            name,
            thermalPhaseChangeProperties,
            twoPhaseProperties,
            T,
            alpha1
        )
    );
}


// ************************************************************************* //
