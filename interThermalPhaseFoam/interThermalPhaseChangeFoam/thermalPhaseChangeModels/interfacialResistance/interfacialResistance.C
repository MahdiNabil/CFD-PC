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

\*---------------------------------------------------------------------------*/

#include "interfacialResistance.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalPhaseChangeModels
{
    defineTypeNameAndDebug(interfacialResistance, 0);
    addToRunTimeSelectionTable
    (
        thermalPhaseChangeModel, 
        interfacialResistance,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalPhaseChangeModels::interfacialResistance::interfacialResistance
(
        const word& name,
        const dictionary& thermalPhaseChangeProperties,
        const twoPhaseThermalMixture& twoPhaseProperties,
        const volScalarField& T,
        const volScalarField& alpha1
)
:
    thermalPhaseChangeModel
    (
        name,
        thermalPhaseChangeProperties,
        twoPhaseProperties,
        T,
        alpha1
    ),
    mesh_(T.mesh()),
    Q_pc_
    (
        IOobject
        (
            "PhaseChangeHeat",
            T_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar( "dummy", dimensionSet(1,-1,-3,0,0,0,0), 0 )
    ),
    InterfaceMeshGraph( mesh_, alpha1 ),
    InterfaceField_
    (
        IOobject
        (
            "InterfaceField",
            T_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        scalar(0)
    ),
    WallField
    (
        IOobject
        (
            "WallField",
            T_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        scalar(0)
    ),
    interfaceArea //Is initialized to zero
    (
        IOobject
        (
            "interfaceArea",
            T_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar( "dummy", dimensionSet(0,2,0,0,0,0,0), 0 )
    ),
    R_g( thermalPhaseChangeProperties_.lookup("R_g") ),
    sigmaHat( thermalPhaseChangeProperties_.lookup("sigmaHat") ),
    v_lv
    (
          (1.0/twoPhaseProperties_.rho2().value())
        - (1.0/twoPhaseProperties_.rho1().value())
    ),
    hi
    (
         (2.0*sigmaHat.value()/(2.0 - sigmaHat.value()))
        *(h_lv_.value()*h_lv_.value()/(T_sat_.value()*v_lv))
        *pow(1.0/(2.0*3.1416*R_g.value()*T_sat_.value()),0.5)
    )

{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::thermalPhaseChangeModels::interfacialResistance::calcQ_pc()
{
    // Now add wall cells to the interfaceField:
    labelList WallCells;
    forAll( mesh_.boundary(), pI )
    {
        if( isA<wallFvPatch>( mesh_.boundary()[pI] ) )    
        {
            WallCells.append( mesh_.boundary()[pI].faceCells() );
        }
    }
    WallField = 0;
    forAll( WallCells, cI )
    {   
        WallField[WallCells[cI]] = 1;
    }

    // Reset all Q_pc to 0
    Q_pc_ = dimensionedScalar( "dummy", dimensionSet(1,-1,-3,0,0,0,0), 0 );

    // Compute some helpful props:
    // For some reason dT is dimensionless
    const dimensionedScalar& dT = 
         alpha1_.time().deltaTValue()
        *dimensionedScalar( "dummy", dimensionSet(0,0,1,0,0,0,0), 1.0 );
    const dimensionedScalar& rho1 = twoPhaseProperties_.rho1();
    const dimensionedScalar& rho2 = twoPhaseProperties_.rho2();

    interfaceArea.internalField() = mag(fvc::grad(alpha1_))*mesh_.V();

    // decaying Phase Change Heat per unit volume
    Q_pc_.internalField() =
         twoPhaseProperties_.rho()*twoPhaseProperties_.cp()
        *((1.0 - exp( -hi*interfaceArea*dT.value()
                      /(  mesh_.V()*twoPhaseProperties_.rho()
                         *twoPhaseProperties_.cp() )
                       ))*(T_ - T_sat_)/dT.value() );


    // Fluid availability limits
    // Get cond/evap limits
    volScalarField LimCond = (1.0 - alpha1_)*(rho2*h_lv_/dT);
    // No evaporation on wall cells!
    volScalarField LimEvap = (1.0 - WallField)*alpha1_*rho1*h_lv_/dT;

    // Apply fluid limiting
    volScalarField Q_pc_fluid = 
        neg(Q_pc_)*max(Q_pc_, -LimCond) + pos(Q_pc_)*min(Q_pc_, LimEvap) ;

    // Volume-based limiting (i.e. relative phase change rate can't exceed |1|
    // per time step
    volScalarField PCV_fac = 
        dT*(Q_pc_ / h_lv_)*(   (scalar(1.0)/twoPhaseProperties_.rho2()) 
                             - (scalar(1.0)/twoPhaseProperties_.rho1()) );

    //Again, don't allow evap on wall   
    volScalarField Q_pc_vol = 
        Q_pc_ * mag( min( max(1.0/(PCV_fac+SMALL), -1.0), (1.0-WallField) ) );

    //Composite limit
    Q_pc_ =
          neg(Q_pc_)*max( max( Q_pc_, Q_pc_fluid ), Q_pc_vol)
        + pos(Q_pc_)*min( min( Q_pc_, Q_pc_fluid ), Q_pc_vol);
}


bool Foam::thermalPhaseChangeModels::interfacialResistance::
read(const dictionary& thermalPhaseChangeProperties)
{
    thermalPhaseChangeModel::read(thermalPhaseChangeProperties);
    return true;
}


// ************************************************************************* //
