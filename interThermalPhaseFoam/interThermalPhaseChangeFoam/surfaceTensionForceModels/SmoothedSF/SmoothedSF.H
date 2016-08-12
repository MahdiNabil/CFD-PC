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
    Foam::surfaceTensionForceModel::SmoothedSF

Description
    Default surface tension force model, from Brackbill (1992)
    Adding some smoothing

SourceFiles
    SmoothedSF.C

\*--------------------------------------------------------------------*/

#ifndef SmoothedSF_H
#define SmoothedSF_H

#include "surfaceTensionForceModel.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceTensionForceModels
{

/*--------------------------------------------------------------------*\
                            Class SmoothedSF
\*--------------------------------------------------------------------*/

class SmoothedSF
:
    public surfaceTensionForceModel
{

// Private Member Functions

    //- This phase change model needs access to the mesh:
    const fvMesh& mesh_;

    //- Surface tension forces on cell faces
    surfaceScalarField Fstffv;

public:

    //- Runtime type information
    TypeName("SmoothedSF");


    // Constructors

    //- construct from components
    SmoothedSF
    (
        const word& name,
        const dictionary& surfaceTensionForceProperties,
        const interfaceProperties& interface,
        const volScalarField& alpha1
    );



    //- Destructor
    ~SmoothedSF()
    {}


    // Member Functions

        //- Correct the surface tension model
        void correct();
        
        //- Returns surface tension force
        tmp<surfaceScalarField> Fstff() const
        {
            return Fstffv;
        }


        //- Read the transportProperties dictionary and update
        bool read(const dictionary& surfaceTensionForceProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceTensionForceModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
