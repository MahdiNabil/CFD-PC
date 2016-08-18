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
    Foam::thermalPhaseChangeModel::InterfaceEquilibrium_SplitDilatation

Description
    This model forces wall and interface cells to the equilibrium state,
    causing evaporation or condensation as appropriate. The phase change
    heat term is limited to not produce/sink more than the volume of a 
    cell in a single time step(like the courantAlpha limit).

    This model uses different interface thresholds for identifying
    evaporating and condensing cell pairs.

    This is similar to the previous HiLoFacePair code, but the phase change
    rate is numerically relaxed by some user-defined value (should be <= 1,
    where 0 is no phase change, and 1 is no relaxation). This seems to help
    with overshooting where the temperature in a cell can go from superheated
    to subcooled in one time step.

    This version splits the volumetric generation/sinking and applies it away
    from the interface so as to not spread out the interface or shoot liquid 
    through it


SourceFiles
    InterfaceEquilibrium_SplitDilatation.C

\*--------------------------------------------------------------------*/

#ifndef InterfaceEquilibrium_SplitDilatation_H
#define InterfaceEquilibrium_SplitDilatation_H

#include "thermalPhaseChangeModel.H"
#include "MeshGraph.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalPhaseChangeModels
{

/*--------------------------------------------------------------------*\
                Class InterfaceEquilibrium_SplitDilatation
\*--------------------------------------------------------------------*/

class InterfaceEquilibrium_SplitDilatation
:
    public thermalPhaseChangeModel
{

// Private Member Functions

    //- Calculate the phase change heat rate
    void calcQ_pc();

    //- Calculate the volume generation/sinking field
    void calcPCV();

    //- This phase change model needs access to the mesh:
    const fvMesh& mesh_;

    //- Phase change heat
    volScalarField Q_pc_;

    //- Graph of mesh grid, used for finding interface cells
    MeshGraph InterfaceMeshGraph;

    //- Field for tracking Interface cells
    volScalarField InterfaceField_;

    //- Field for blocking evaporation on wall cells
    volScalarField WallField;

    //- Volume change field (from phase change)
    volScalarField PCVField;

    // Scalar parameters for identifying condensation and evaporation int.
    // cell pairs:
    scalar CondThresh;
    scalar EvapThresh;
    scalar RelaxFac; //The relaxation factor (0-1)

public:

    //- Runtime type information
    TypeName("InterfaceEquilibrium_SplitDilatation");


    // Constructors

    //- construct from components
    InterfaceEquilibrium_SplitDilatation
    (
        const word& name,
        const dictionary& thermalPhaseChangeProperties,
        const twoPhaseThermalMixture& twoPhaseProperties,
        const volScalarField& T,
        const volScalarField& alpha1
    );


    //- Destructor
    ~InterfaceEquilibrium_SplitDilatation()
    {}


    // Member Functions

        //- Return phase change heat
        tmp<volScalarField> Q_pc() const
        {
            return Q_pc_;
        }

        //- Correct the thermalPhaseChange model
        void correct()
        {
            calcQ_pc();
            calcPCV();
        }

        //- Return interface cells marker field:
        tmp<volScalarField> InterfaceField() const
        {
            return InterfaceField_;
        }
        
        //- Gets volume generation (split and applied slightly away from 
        //  interface)
        tmp<volScalarField> PCV() const
        {
            return PCVField;
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
