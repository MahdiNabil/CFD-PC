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

#include "InterfaceEquilibrium.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalPhaseChangeModels
{
    defineTypeNameAndDebug(InterfaceEquilibrium, 0);
    addToRunTimeSelectionTable
    (
        thermalPhaseChangeModel, 
        InterfaceEquilibrium, 
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalPhaseChangeModels::InterfaceEquilibrium::InterfaceEquilibrium
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
        T, alpha1
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
    )
{
    // Read in the cond/evap int. thresholds
    thermalPhaseChangeProperties_.lookup("CondThresh") >> CondThresh;
    thermalPhaseChangeProperties_.lookup("EvapThresh") >> EvapThresh;
    thermalPhaseChangeProperties_.lookup("RelaxFac") >> RelaxFac;   

    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::thermalPhaseChangeModels::InterfaceEquilibrium::calcQ_pc()
{
    // Get the sets of interface cell face pairs for evaporation/condensation
    std::vector<MeshGraph::CellFacePair> CondIntCellFacePairs;
    std::vector<MeshGraph::CellFacePair> EvapIntCellFacePairs;

    // Find internal interface cell pairs using graph traversal
    InterfaceMeshGraph.Reset();
    InterfaceMeshGraph.GetInterfaceCellFacePairs
    (
        CondIntCellFacePairs,
        CondThresh
    );
    
    InterfaceMeshGraph.Reset();
    InterfaceMeshGraph.GetInterfaceCellFacePairs
    (
        EvapIntCellFacePairs,
        EvapThresh
    );

    // Compute the interpolated T field to see which face pairs are actually
    // evaporating/condensing:
    surfaceScalarField Tf = fvc::interpolate(T_);

    // Reset interface field, then interpolate
    InterfaceField_ = 0;

    // Loop through cond cells:
    for
    (
        std::vector<MeshGraph::CellFacePair>::iterator it =
            CondIntCellFacePairs.begin();
        it != CondIntCellFacePairs.end(); 
        it++
    )
    {
        // Check that temp is below T_sat for condensation
        if ( Tf[(*it).f] <= T_sat_.value() )
        {
            InterfaceField_[(*it).c1] = 1;  
            InterfaceField_[(*it).c2] = 1;
        }
    }

    // Loop through evap cells:
    for
    (
        std::vector<MeshGraph::CellFacePair>::iterator it =
            EvapIntCellFacePairs.begin();
        it != EvapIntCellFacePairs.end(); 
        it++
    )
    {
        // Check that temp is above T_sat for evaporation
        if ( Tf[(*it).f] >= T_sat_.value() )
        {
            InterfaceField_[(*it).c1] = 1;
            InterfaceField_[(*it).c2] = 1;
        }
    }


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
        InterfaceField_[WallCells[cI]] = 1;
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

    //Unlimited phase change heat
    Q_pc_ = 
         InterfaceField_
        *twoPhaseProperties_.rho()
        *twoPhaseProperties_.cp()
        *((T_ - T_sat_)/dT);

    // Fluid availability limits
    // Get cond/evap limits
    volScalarField LimCond = (1.0 - alpha1_)*(rho2*h_lv_/dT);
    // No evaporation on wall cells!
    volScalarField LimEvap = (1.0 - WallField)*alpha1_*rho1*h_lv_/dT;


    // Apply fluid limiting
    volScalarField Q_pc_fluid = 
        neg(Q_pc_)*max(Q_pc_, -LimCond) + pos(Q_pc_)*min(Q_pc_, LimEvap) ;

    // Volume-based limiting (i.e. relative phase change rate can't exceed 
    // |1| per time step
    volScalarField PCV_fac = 
        dT*(Q_pc_ / h_lv_)*(   (scalar(1.0)/twoPhaseProperties_.rho2()) 
                             - (scalar(1.0)/twoPhaseProperties_.rho1()) );

    // Volume generation/sink based limited
    // Again, allow regular evap on wall 
    volScalarField Q_pc_vol =
        Q_pc_ * mag( min( max(1.0/(PCV_fac+SMALL), -1.0), (1.0-WallField) ) );

    //Composite limit
    Q_pc_ = 
           neg(Q_pc_)*max( max( Q_pc_, Q_pc_fluid ), Q_pc_vol) 
         + pos(Q_pc_)*min( min( Q_pc_, Q_pc_fluid ), Q_pc_vol);

    //Under relax phase change rate per user specification
    Q_pc_ = RelaxFac * Q_pc_;
}


bool Foam::thermalPhaseChangeModels::InterfaceEquilibrium::
read(const dictionary& thermalPhaseChangeProperties)
{
    thermalPhaseChangeModel::read(thermalPhaseChangeProperties);

    //Read in the cond/evap int. thresholds
    thermalPhaseChangeProperties_.lookup("CondThresh") >> CondThresh;
    thermalPhaseChangeProperties_.lookup("EvapThresh") >> EvapThresh;
    thermalPhaseChangeProperties_.lookup("RelaxFac") >> RelaxFac;

    return true;
}


// ************************************************************************* //
