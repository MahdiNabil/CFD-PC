/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 Alex Rattner, Sanjay Adhikari
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
    Foam::thermalPhaseChangeModel::interfacialResistance

Description
    This model evaluates the phase change rate based on the interfacial
    resistance applied over the diffuse interface

SourceFiles
    interfacialResistance.C

\*--------------------------------------------------------------------*/

#ifndef interfacialResistance_H
#define interfacialResistance_H

#include "thermalPhaseChangeModel.H"
#include "MeshGraph.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalPhaseChangeModels
{

/*--------------------------------------------------------------------*\
                              Class interfacialResistance
\*--------------------------------------------------------------------*/

class interfacialResistance
:
    public thermalPhaseChangeModel
{

// Private Member Functions

    //- Calculate the phase change heat rate
    void calcQ_pc();

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

    //- magnitude of alphaGradient field
    volScalarField interfaceArea;

    //Some constants for interfacial heat transfer
    dimensionedScalar R_g;       // specific ideal gas constant 
    dimensionedScalar sigmaHat;  // accomodation coefficient
    double v_lv;                 // change in specific volume 
    double hi;                   // interfacial resistance
    

public:

    //- Runtime type information
    TypeName("interfacialResistance");


    // Constructors

    //- construct from components
    interfacialResistance
    (
        const word& name,
        const dictionary& thermalPhaseChangeProperties,
        const twoPhaseThermalMixture& twoPhaseProperties,
        const volScalarField& T,
        const volScalarField& alpha1
    );


    //- Destructor
    ~interfacialResistance()
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
        }

        //- Return interface cells marker field:
        tmp<volScalarField> InterfaceField() const
        {
            return InterfaceField_;
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
