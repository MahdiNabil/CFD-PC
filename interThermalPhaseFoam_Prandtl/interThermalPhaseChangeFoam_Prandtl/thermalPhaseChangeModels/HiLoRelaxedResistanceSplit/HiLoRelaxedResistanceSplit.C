/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 Alex Rattner
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

#include "HiLoRelaxedResistanceSplit.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalPhaseChangeModels
{
    defineTypeNameAndDebug(HiLoRelaxedResistanceSplit, 0);
    addToRunTimeSelectionTable(thermalPhaseChangeModel, HiLoRelaxedResistanceSplit, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalPhaseChangeModels::HiLoRelaxedResistanceSplit::HiLoRelaxedResistanceSplit
(
		const word& name,
		const dictionary& thermalPhaseChangeProperties,
		const twoPhaseThermalMixture& twoPhaseProperties,
		const volScalarField& T,
		const volScalarField& alpha1
)
:
    thermalPhaseChangeModel(name, thermalPhaseChangeProperties, twoPhaseProperties, T, alpha1),
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
	//Read in the cond/evap int. thresholds
	thermalPhaseChangeProperties_.lookup("CondThresh") >> CondThresh;
	thermalPhaseChangeProperties_.lookup("EvapThresh") >> EvapThresh;
	thermalPhaseChangeProperties_.lookup("RelaxFac") >> RelaxFac;	

	correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::thermalPhaseChangeModels::HiLoRelaxedResistanceSplit::calcQ_pc()
{
	//Get the sets of interface cell face pairs for evaporation/condensation
	std::vector<MeshGraph::CellFacePair> CondIntCellFacePairs, EvapIntCellFacePairs;


	//Find internal interface cell pairs using graph traversal
	InterfaceMeshGraph.Reset();
	InterfaceMeshGraph.GetInterfaceCellFacePairs(CondIntCellFacePairs, CondThresh);
	InterfaceMeshGraph.Reset();
	InterfaceMeshGraph.GetInterfaceCellFacePairs(EvapIntCellFacePairs, EvapThresh);

	//Compute the interpolated T field to see which face pairs are actually evaporating/condensing:
	surfaceScalarField Tf = fvc::interpolate(T_);

	//Reset interface field, then interpolate
	InterfaceField_ = 0;

	//Loop through cond cells:
	for (std::vector<MeshGraph::CellFacePair>::iterator it = CondIntCellFacePairs.begin(); it != CondIntCellFacePairs.end(); it++)
	{
		//Check that temp is below T_sat for condensation
		if ( Tf[(*it).f] <= T_sat_.value() )
		{   InterfaceField_[(*it).c1] = 1;  InterfaceField_[(*it).c2] = 1;  }
	}

	//Loop through evap cells:
	for (std::vector<MeshGraph::CellFacePair>::iterator it = EvapIntCellFacePairs.begin(); it != EvapIntCellFacePairs.end(); it++)
	{
		//Check that temp is above T_sat for evaporation
		if ( Tf[(*it).f] >= T_sat_.value() )
		{   InterfaceField_[(*it).c1] = 1;  InterfaceField_[(*it).c2] = 1;  }
	}

	//Now add wall cells to the interfaceField:
	labelList WallCells;
	forAll( mesh_.boundary(), pI )
	{
		if( isA<wallFvPatch>( mesh_.boundary()[pI] ) )    
		{  WallCells.append( mesh_.boundary()[pI].faceCells() );  }
	}
	WallField = 0;
	forAll( WallCells, cI )
	{   
		WallField[WallCells[cI]] = 1;
		InterfaceField_[WallCells[cI]] = 1;
	}

	//List total int. cells
	//Info<< "Total interface cells: " << gSum(InterfaceField_) << endl;

	//Reset all Q_pc to 0
	Q_pc_ = dimensionedScalar( "dummy", dimensionSet(1,-1,-3,0,0,0,0), 0 );

	//*******************************************************************************************
	//*******************************************************************************************
	//Alternate approximate thermal-resistance based phase change model
	const dimensionedScalar eps_gAlpha( "dummy", dimensionSet(0,-1,0,0,0,0,0), SMALL );
	const dimensionedScalar eps_A("dummy", dimensionSet(0,2,0,0,0,0,0), SMALL );
	const surfaceVectorField g_alpha1f = fvc::interpolate( fvc::grad(alpha1_) );
	const surfaceVectorField n_alpha1f = g_alpha1f / ( mag(g_alpha1f) + eps_gAlpha );
	volScalarField cVols
    (
        IOobject
        (
            "CellVolume",
            T_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
		mesh_,
		dimensionedScalar( "dummy", dimensionSet(0,3,0,0,0,0,0), 1 )
    );
	cVols.internalField() *= mesh_.V();

	//Approximate interface area per cell:
	const volScalarField A_int = fvc::surfaceSum( mag( mesh_.Sf() & n_alpha1f ) )/2;

	//Approximate heat transfer length per cell:
	const volScalarField L_HT = cVols / ( fvc::surfaceSum( mag( mesh_.Sf() & n_alpha1f ) ) + eps_A );

	//Phase change heating rate
	Q_pc_ = InterfaceField_ * twoPhaseProperties_.lambda() * A_int * (T_ - T_sat_)/(L_HT * cVols);


//Comparison
//Info<< "Q_pc original: " << gSum( Q_pc_.internalField() * mesh_.V() ) << endl;
//Info<< "Q_pc new: " << gSum( Q_pcB_.internalField() * mesh_.V() ) << endl;

	//*******************************************************************************************
	//*******************************************************************************************
}


//- Gets volume generation (split and applied slightly away from interface)
void Foam::thermalPhaseChangeModels::HiLoRelaxedResistanceSplit::calcPCV()
{
	//Get some local references	to helpful fields
	const volScalarField& Q_pc_ = this->Q_pc();
	//Direction of interface in each cell
	const volVectorField gradAlpha = fvc::grad( alpha1_ );
	dimensionedScalar epsInvLength( "dummy", dimensionSet(0,-1,0,0,0,0,0), SMALL );
	const volVectorField nAlpha = gradAlpha/( mag(gradAlpha) + epsInvLength );
	//const surfaceVectorField gradAlphaf = fvc:sngrad( alpha1 );
	
	//Init the PCV field:
	volScalarField LiquidVolGen = (Q_pc_ / h_lv_)*( scalar(-1.0)/twoPhaseProperties_.rho1() );
	volScalarField VaporVolGen = (Q_pc_ / h_lv_)*( scalar(1.0)/twoPhaseProperties_.rho2() );

	//Scale by mesh volume:
	LiquidVolGen.internalField() = LiquidVolGen.internalField() * mesh_.V();
	VaporVolGen.internalField() = VaporVolGen.internalField() * mesh_.V();
	
	//Now transport the volumetric generation away from the interface
	//Scan through all cells and apply rule

	//Some helpful constants
	label nFaceNeighbors = mesh_.faceNeighbour().size();
	dimensionedScalar epsQ_pc( "dummy", dimensionSet(0,-1,0,0,0,0,0), SMALL );
	const labelList& FaceOwners = mesh_.faceOwner();
	const labelList& FaceNeighbours = mesh_.faceNeighbour();
	const surfaceVectorField& FaceVectors = mesh_.Sf();

//Try two passes to move away PCV:
for (int k = 0; k<4; k++)
{

	forAll(InterfaceField_, cI)
	{
		//First, skip if this cell is not on the interface
		if ( InterfaceField_[cI] != 1 )
		{  continue;  }

		//Also skip if no phase change is happening here:
		if ( (mag( LiquidVolGen[cI] ) <= SMALL) && (mag( VaporVolGen[cI] ) <= SMALL) )
		{  continue;  }

//Check nAlpha
//Info<< "CurCell: " << mesh_.C()[cI] << ",   nAlpha: " << nAlpha[cI] << endl;

		//Get list of faces on the cell (that point to other cells)
		const cell& curCell = mesh_.cells()[cI];
		labelList curCellFaces;
		forAll(curCell, fI)
		{
			//Only consider faces that point to other cells (> size of faceNeighbor list)
			if (curCell[fI] < nFaceNeighbors) //This face is shared by other cells, so possibly transport vol Generation across it
			{  curCellFaces.append( curCell[fI] );  }
		}

		//Get cell neighbors and face area vectors for each shared face
		List<vector> curFaceVectors;
		labelList curFaceNeighbours;
		forAll(curCellFaces, fI)
		{
			label curFaceNeighbour = (cI == FaceOwners[curCellFaces[fI]]) ? FaceNeighbours[curCellFaces[fI]] : FaceOwners[curCellFaces[fI]];
			curFaceNeighbours.append( curFaceNeighbour );
			if ( cI < curFaceNeighbour ) //Points away from current cell
			{
				curFaceVectors.append( FaceVectors[curCellFaces[fI]] );
				//Check face vector facing:
				//Info<< "CurCell: " << mesh_.C()[cI] << ",   CurFaceNeighbour: " << mesh_.C()[curFaceNeighbour] << ",  vector: " << FaceVectors[curCellFaces[fI]] << endl;
			}
			//else if ( cI == FaceNeighbours[curCellFaces[fI]] ) //We are the neighbor of this face (flip it)
			else //Points into current cell
			{
				curFaceVectors.append( -1.0*FaceVectors[curCellFaces[fI]] );
				//Check face vector facing:
				//Info<< "CurCell: " << mesh_.C()[cI] << ",   CurFaceNeighbour: " << mesh_.C()[curFaceNeighbour] << ",  vector: " << -1.0*FaceVectors[curCellFaces[fI]] << endl;
			}


		}


//Info<< cI << ": " << curCellFaces.size() << endl;
		//Now get relative surface area vectors from each face (i.e. those pointing in each direction sum to 1)
		scalarList LiquidFaceFractions, VaporFaceFractions;    //Fraction of volume generation to move
		labelList LiquidFaceTargets, VaporFaceTargets; //Where to move the volume generation to
		scalar TotalLiquidFacing = 0;
		scalar TotalVaporFacing = 0;

		forAll(curCellFaces, fI)
		{
			//Get relative facing to the interface normal (what fraction of volgen should go through this face:
			scalar IntFacing = ( curFaceVectors[fI] & gradAlpha[cI] );
			if      ( IntFacing > 0 ) //Facing the liquid side
			{
				LiquidFaceFractions.append( IntFacing );
				LiquidFaceTargets.append( curFaceNeighbours[fI] );
				TotalLiquidFacing += IntFacing;
			}
			else if ( IntFacing < 0 ) //Facing the vapor side
			{
				VaporFaceFractions.append( -1.0*IntFacing );
				VaporFaceTargets.append( curFaceNeighbours[fI] );
				TotalVaporFacing -= IntFacing;
			}
		}

		//Normalize both face fraction lists
		forAll( LiquidFaceFractions, fI)
		{   LiquidFaceFractions[fI] /= TotalLiquidFacing;  }
		forAll( VaporFaceFractions, fI)
		{   VaporFaceFractions[fI] /= TotalVaporFacing;  }

		//Finally correct the PCV field - move volumetric generation away from interface cells:
		dimensionedScalar curLiquidVolGen = LiquidVolGen[cI];
		dimensionedScalar curVaporVolGen = VaporVolGen[cI];
		//Apply liquid change:
		forAll( LiquidFaceFractions, fI)
		{
			//PCVField[cI]                    -= curLiquidVolGen.value()*LiquidFaceFractions[fI];
			//PCVField[LiquidFaceTargets[fI]] += curLiquidVolGen.value()*LiquidFaceFractions[fI];
			LiquidVolGen[cI]                    -= curLiquidVolGen.value()*LiquidFaceFractions[fI];
			LiquidVolGen[LiquidFaceTargets[fI]] += curLiquidVolGen.value()*LiquidFaceFractions[fI];
			
		}

		forAll( VaporFaceFractions, fI)
		{
			//PCVField[cI]                    -= curVaporVolGen.value()*VaporFaceFractions[fI];
			//PCVField[VaporFaceTargets[fI]]  += curVaporVolGen.value()*VaporFaceFractions[fI];
			VaporVolGen[cI]                   -= curVaporVolGen.value()*VaporFaceFractions[fI];
			VaporVolGen[VaporFaceTargets[fI]] += curVaporVolGen.value()*VaporFaceFractions[fI];
		}

	}

//End try two passes
}

	//Combine two parts of generation:
	PCVField = LiquidVolGen + VaporVolGen;

	//Renormalize by cell volume:
	PCVField.internalField() = PCVField.internalField() / mesh_.V();

    //PCVField = ( (Q_pc_ / h_lv_)*( (scalar(1.0)/twoPhaseProperties_.rho2()) - (scalar(1.0)/twoPhaseProperties_.rho1()) ) );
}



bool Foam::thermalPhaseChangeModels::HiLoRelaxedResistanceSplit::read(const dictionary& thermalPhaseChangeProperties)
{
	thermalPhaseChangeModel::read(thermalPhaseChangeProperties);

	//Read in the cond/evap int. thresholds
	thermalPhaseChangeProperties_.lookup("CondThresh") >> CondThresh;
	thermalPhaseChangeProperties_.lookup("EvapThresh") >> EvapThresh;
	thermalPhaseChangeProperties_.lookup("RelaxFac") >> RelaxFac;

	return true;
}


// ************************************************************************* //
