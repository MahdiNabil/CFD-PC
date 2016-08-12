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
    Foam::thermalPhaseChangeModel

Description
    Parent class for thermal phase change models.

SourceFiles
    thermalPhaseChangeModel.C
    newThermalPhaseChangeModel.C

\*---------------------------------------------------------------------------*/

#ifndef thermalPhaseChangeModel_H
#define thermalPhaseChangeModel_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "tmp.H"
#include "autoPtr.H"
#include "twoPhaseThermalMixture.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class thermalPhaseChangeModel
\*---------------------------------------------------------------------------*/

class thermalPhaseChangeModel : public IOdictionary
{

protected:

    // Protected data
    word name_;
    dictionary thermalPhaseChangeProperties_;
    // Thermal mixture props
    const twoPhaseThermalMixture& twoPhaseProperties_;
    // Pointers to important fields
    const volScalarField& T_;
    const volScalarField& alpha1_;
    // Get these properties from the thermalPhaseChange dictionary       
    dimensionedScalar T_sat_;
    dimensionedScalar h_lv_;
    // Switch to turn off and on the alpha1 and v_pc source terms
    Switch sw_PCV;
    Switch sw_alpha1Gen;


    // Private Member Functions

    //- Disallow copy construct
    thermalPhaseChangeModel(const thermalPhaseChangeModel&);

    //- Disallow default bitwise assignment
    void operator=(const thermalPhaseChangeModel&);


public:

    //- Runtime type information
    TypeName("thermalPhaseChangeModel");


    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        thermalPhaseChangeModel,
        dictionary,
        (
            const word& name,
            const dictionary& thermalPhaseChangeProperties,
            const twoPhaseThermalMixture& twoPhaseProperties,
            const volScalarField& T,
            const volScalarField& alpha1
        ),
        (name, thermalPhaseChangeProperties, twoPhaseProperties, T, alpha1)
    );


    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<thermalPhaseChangeModel> New
        (
            const word& name,
            const dictionary& thermalPhaseChangeProperties,
            const twoPhaseThermalMixture& twoPhaseProperties,
            const volScalarField& T,
            const volScalarField& alpha1
        );


    // Constructors

        //- Construct from components
        thermalPhaseChangeModel
        (
            const word& name,
            const dictionary& thermalPhaseChangeProperties,
            const twoPhaseThermalMixture& twoPhaseProperties,
            const volScalarField& T,
            const volScalarField& alpha1
        );


    //- Destructor
    virtual ~thermalPhaseChangeModel()
    {}


    // Member Functions
        
        //- Return the volumetric rate of phase change heat - i.e. the heat
        //  source/sink rate per unit volume due to phase change. Evaporation
        //  is positive, and condensation is negative, units (W/m^3)
        virtual tmp<volScalarField> Q_pc() const = 0;

        //- Return Volumetric generation rate of volume due to phase change
        //  i.e. how much fluid volume is generated due to evaporation
        //  or destroyed due to condensation, units (1/s)
        virtual tmp<volScalarField> PCV() const;

        //- Return volumetric production rate of the alpha1 field due to phase
        //  change. This is positive in the case of condensation, or negative
        //  for evaporation, units (1/s)
        virtual tmp<volScalarField> alpha1Gen() const;

        //- Correct the thermalPhaseChange model
        virtual void correct() = 0;

        //- Read the thermalPhaseChange dictionary and update
        virtual bool read(const dictionary& thermalPhaseChangeProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
