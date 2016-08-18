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

#include "InterfaceEquilibrium_SplitDilatation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalPhaseChangeModels
{
    defineTypeNameAndDebug(InterfaceEquilibrium_SplitDilatation, 0);
    addToRunTimeSelectionTable
    (
        thermalPhaseChangeModel,
        InterfaceEquilibrium_SplitDilatation,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalPhaseChangeModels::InterfaceEquilibrium_SplitDilatation::
InterfaceEquilibrium_SplitDilatation
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
    PCVField
    (
        IOobject
        (
            "PhaseChangeVolume",
            T_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar( "dummy", dimensionSet(0,0,-1,0,0,0,0), 0 )
    )
{
    // Read in the cond/evap int. thresholds
    thermalPhaseChangeProperties_.lookup("CondThresh") >> CondThresh;
    thermalPhaseChangeProperties_.lookup("EvapThresh") >> EvapThresh;
    thermalPhaseChangeProperties_.lookup("RelaxFac") >> RelaxFac;   

    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::thermalPhaseChangeModels::InterfaceEquilibrium_SplitDilatation::calcQ_pc()
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

    // Unlimited phase change heat
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
          neg(Q_pc_)*max(Q_pc_, -LimCond)
        + pos(Q_pc_)*min(Q_pc_, LimEvap) ;

    // Volume-based limiting (i.e. relative phase change rate can't exceed |1|
    // per time step
    volScalarField PCV_fac = 
        dT*(Q_pc_/h_lv_)*(   (scalar(1.0)/twoPhaseProperties_.rho2())
                           - (scalar(1.0)/twoPhaseProperties_.rho1()) );

    // Again, don't allow evap on wall   
    volScalarField Q_pc_vol =
           Q_pc_ 
         * mag( min( max(1.0/(PCV_fac + SMALL), -1.0), (1.0 - WallField) ) );

    // Composite limit
    Q_pc_ =
          neg(Q_pc_)*max( max( Q_pc_, Q_pc_fluid ), Q_pc_vol) 
        + pos(Q_pc_)*min( min( Q_pc_, Q_pc_fluid ), Q_pc_vol);

    // Under relax phase change rate per user specification
    Q_pc_ = RelaxFac * Q_pc_;

}


//- Gets volume generation (split and applied slightly away from interface)
void Foam::thermalPhaseChangeModels::InterfaceEquilibrium_SplitDilatation::calcPCV()
{
    // Get some local references to helpful fields
    const volScalarField& Q_pc_ = this->Q_pc();
    // Direction of interface in each cell
    const volVectorField gradAlpha = fvc::grad( alpha1_ );
    dimensionedScalar epsInvLength
    (
        "dummy",
        dimensionSet(0,-1,0,0,0,0,0),
        SMALL
    );
    const volVectorField nAlpha = gradAlpha/( mag(gradAlpha) + epsInvLength );
    
    
    // Init the PCV field:
    volScalarField LiquidVolGen =
        (Q_pc_/h_lv_)*(scalar(-1.0)/twoPhaseProperties_.rho1());
    volScalarField VaporVolGen = 
        (Q_pc_/h_lv_)*(scalar(1.0)/twoPhaseProperties_.rho2());

    // Scale by mesh volume:
    LiquidVolGen.internalField() = LiquidVolGen.internalField() * mesh_.V();
    VaporVolGen.internalField() = VaporVolGen.internalField() * mesh_.V();
    
    // Now transport the volumetric generation away from the interface
    // Scan through all cells and apply rule

    // Some helpful constants
    label nFaceNeighbors = mesh_.faceNeighbour().size();
    dimensionedScalar epsQ_pc( "dummy", dimensionSet(0,-1,0,0,0,0,0), SMALL );
    const labelList& FaceOwners = mesh_.faceOwner();
    const labelList& FaceNeighbours = mesh_.faceNeighbour();
    const surfaceVectorField& FaceVectors = mesh_.Sf();

    // Try two passes to move away PCV:
    for (int k = 0; k<4; k++)
    {

        forAll(InterfaceField_, cI)
        {
            // First, skip if this cell is not on the interface
            if ( InterfaceField_[cI] != 1 )
            {
                continue;
            }

            // Also skip if no phase change is happening here:
            if
            (
                   (mag( LiquidVolGen[cI] ) <= SMALL)
                && (mag( VaporVolGen[cI] ) <= SMALL)
            )
            {
                continue;
            }


            // Get list of faces on the cell (that point to other cells)
            const cell& curCell = mesh_.cells()[cI];
            labelList curCellFaces;
            forAll(curCell, fI)
            {
                // Only consider faces that point to other cells (> size of
                // faceNeighbor list)
                // This face is shared by other cells, so possibly transport vol
                // Generation across it
                if (curCell[fI] < nFaceNeighbors) 
                {
                    curCellFaces.append( curCell[fI] );
                }
            }

            // Get cell neighbors and face area vectors for each shared face
            List<vector> curFaceVectors;
            labelList curFaceNeighbours;
            forAll(curCellFaces, fI)
            {
                label curFaceNeighbour = 
                       (cI == FaceOwners[curCellFaces[fI]])
                     ? FaceNeighbours[curCellFaces[fI]]
                     : FaceOwners[curCellFaces[fI]];
                     
                curFaceNeighbours.append( curFaceNeighbour );

                if ( cI < curFaceNeighbour ) // Points away from current cell
                {
                    curFaceVectors.append( FaceVectors[curCellFaces[fI]] );
                }
                else // Points into current cell
                {
                    curFaceVectors.append( -1.0*FaceVectors[curCellFaces[fI]] );
                }


            }


            // Now get relative surface area vectors from each face (i.e. those
            // pointing in each direction sum to 1)
            // Fraction of volume generation to move
            scalarList LiquidFaceFractions, VaporFaceFractions;    
            // Where to move the volume generation to
            labelList LiquidFaceTargets, VaporFaceTargets; 
            scalar TotalLiquidFacing = 0;
            scalar TotalVaporFacing = 0;

            forAll(curCellFaces, fI)
            {
                // Get relative facing to the interface normal (what fraction of
                // volgen should go through this face:
                scalar IntFacing = ( curFaceVectors[fI] & gradAlpha[cI] );
                if      ( IntFacing > 0 ) //Facing the liquid side
                {
                    LiquidFaceFractions.append( IntFacing );
                    LiquidFaceTargets.append( curFaceNeighbours[fI] );
                    TotalLiquidFacing += IntFacing;
                }
                else if ( IntFacing < 0 ) // Facing the vapor side
                {
                    VaporFaceFractions.append( -1.0*IntFacing );
                    VaporFaceTargets.append( curFaceNeighbours[fI] );
                    TotalVaporFacing -= IntFacing;
                }
            }

            // Normalize both face fraction lists
            forAll( LiquidFaceFractions, fI)
            {
                LiquidFaceFractions[fI] /= TotalLiquidFacing;
            }
            forAll( VaporFaceFractions, fI)
            {
                VaporFaceFractions[fI] /= TotalVaporFacing;
            }

            // Finally correct the PCV field - move volumetric generation away
            // from interface cells:
            dimensionedScalar curLiquidVolGen = LiquidVolGen[cI];
            dimensionedScalar curVaporVolGen = VaporVolGen[cI];
            //Apply liquid change:
            forAll( LiquidFaceFractions, fI)
            {
                LiquidVolGen[cI] -=
                    curLiquidVolGen.value()*LiquidFaceFractions[fI];
                LiquidVolGen[LiquidFaceTargets[fI]] +=
                    curLiquidVolGen.value()*LiquidFaceFractions[fI];
                
            }

            forAll( VaporFaceFractions, fI)
            {
                VaporVolGen[cI] -=
                    curVaporVolGen.value()*VaporFaceFractions[fI];
                VaporVolGen[VaporFaceTargets[fI]] +=
                    curVaporVolGen.value()*VaporFaceFractions[fI];
            }

        }

    //End try two passes
    }

    //Combine two parts of generation:
    PCVField = LiquidVolGen + VaporVolGen;

    //Renormalize by cell volume:
    PCVField.internalField() = PCVField.internalField() / mesh_.V();

}



bool Foam::thermalPhaseChangeModels::InterfaceEquilibrium_SplitDilatation::
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
