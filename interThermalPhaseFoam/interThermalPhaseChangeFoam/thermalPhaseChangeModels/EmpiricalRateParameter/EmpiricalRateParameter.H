/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 Alex Rattner and Mahdi Nabil
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
    Foam::thermalPhaseChangeModel::EmpiricalRateParameter

Description
    This model implements the empirical rate parameter phase change
    formulation of Yang, Z., Peng, X. F., & Ye, P. (2008). Numerical and
    experimental investigation of two phase flow during boiling in a coiled
    tube. International Journal of Heat and Mass Transfer, 51(5), 1003-1016.


SourceFiles
    EmpiricalRateParameter.C

\*--------------------------------------------------------------------*/

#ifndef EmpiricalRateParameter_H
#define EmpiricalRateParameter_H

#include "thermalPhaseChangeModel.H"
#include "MeshGraph.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalPhaseChangeModels
{

/*--------------------------------------------------------------------*\
                              Class EmpiricalRateParameter
\*--------------------------------------------------------------------*/

class EmpiricalRateParameter
:
    public thermalPhaseChangeModel
{

// Private Member Functions

    //- Calculate and return the phase change heat rate
    void calcQ_pc();

    //- Phase change heat
    volScalarField Q_pc_;

    // declaring rl and rv
    dimensionedScalar rl;
    dimensionedScalar rv;

public:

    //- Runtime type information
    TypeName("EmpiricalRateParameter");


    // Constructors

    //- construct from components
    EmpiricalRateParameter
    (
        const word& name,
        const dictionary& thermalPhaseChangeProperties,
        const twoPhaseThermalMixture& twoPhaseProperties,
        const volScalarField& T,
        const volScalarField& alpha1
    );


    //- Destructor
    ~EmpiricalRateParameter()
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

        //- Read the transportProperties dictionary and update
        bool read(const dictionary& thermalPhaseChangeProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermalPhaseChangeModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
