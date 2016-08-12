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
    Foam::thermalPhaseChangeModel::noPhaseChange

Description
    Dummy model for no phase change. Just set none as the model type.

SourceFiles
    noPhaseChange.C

\*--------------------------------------------------------------------*/

#ifndef noPhaseChange_H
#define noPhaseChange_H

#include "thermalPhaseChangeModel.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalPhaseChangeModels
{

/*--------------------------------------------------------------------*\
                         Class noPhaseChange
\*--------------------------------------------------------------------*/

class noPhaseChange
:
    public thermalPhaseChangeModel
{

// Private Member Functions

    //- This phase change model needs access to the mesh:
    const fvMesh& mesh_;

    //- Phase change heat
    volScalarField Q_pc_;

public:

    //- Runtime type information
    TypeName("none");


    // Constructors

    //- construct from components
    noPhaseChange
    (
        const word& name,
        const dictionary& thermalPhaseChangeProperties,
        const twoPhaseThermalMixture& twoPhaseProperties,
        const volScalarField& T,
        const volScalarField& alpha1
    );



    //- Destructor
    ~noPhaseChange()
    {}


    // Member Functions

        //- Correct the thermalPhaseChange model
        void correct()
        {

        }
        
        //- Return phase change heat
        tmp<volScalarField> Q_pc() const
        {
            return Q_pc_;
        }


        //- Read the transportProperties dictionary and update
        bool read(const dictionary& thermalPhaseChangeProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermalPhaseChangeModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
