/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 Alex Rattner
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

#include "LimitedFacePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalPhaseChangeModels
{
    defineTypeNameAndDebug(LimitedFacePair, 0);
    addToRunTimeSelectionTable(thermalPhaseChangeModel, LimitedFacePair, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalPhaseChangeModels::LimitedFacePair::LimitedFacePair
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
    InterfaceField
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
	correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::thermalPhaseChangeModels::LimitedFacePair::calcQ_pc()
{
//- Find the set of cells on the interface (or on a wall)	
	//First get the interface cells:
	labelList InterfaceCells;
	labelList WallCells;

	//Find internal interface cells using graph traversal
	InterfaceMeshGraph.Reset();
	//InterfaceMeshGraph.GetInterfaceCells(InterfaceCells, 0.5);
	InterfaceMeshGraph.GetDoubleInterfaceCells(InterfaceCells, 0.5);

/*
	//Try deadbanding to find other interface cells:
	scalar dbHigh = 0.99;
	scalar dbLow  = 0.01;
	forAll (alpha1_, j)
	{
		if ( (alpha1_[j] >= dbLow) && (alpha1_[j] < dbHigh) )
		{  InterfaceCells.append( j );  }

	}
*/


	//Spit out internal interface cells count
	//We can't access collective operations directly, have to go through OpenFoam templates (like gSum)
	labelList myIntCellsCount;
	myIntCellsCount.append( InterfaceCells.size() );
	Info<< "Internal interface cells: " << gSum( myIntCellsCount ) << endl;

	//Also add cells touching the wall - these can also be phase change sites
	//This probably only needs to be done once, at the beginning of the simulation,
	// but is more convenient to do at every time step (as long as it isn't too expensive)

	forAll( mesh_.boundary(), pI )
	{
		if( isA<wallFvPatch>( mesh_.boundary()[pI] ) )    
		{
			InterfaceCells.append( mesh_.boundary()[pI].faceCells() );
			WallCells.append( mesh_.boundary()[pI].faceCells() );
		}
	}
	WallField = 0;
	forAll( WallCells, cI )
	{   WallField[WallCells[cI]] = 1.0;  }

	//Now mark the InterfaceField - mostly for convenient reading of output
	InterfaceField = 0;
	forAll( InterfaceCells, cI )
	{   InterfaceField[InterfaceCells[cI]] = 1.0;   }


	//Spit out total interface cells count
	myIntCellsCount.clear();
	myIntCellsCount.append( InterfaceCells.size() );
	Info<< "Total interface cells: " << gSum( myIntCellsCount ) << endl;

//- This next section evaluates the phase change heat transfer on interface/wall cells
	//Reset all Q_pc to 0
	//Compute some helpful props:
	//For some reason dT is dimensionless
	const dimensionedScalar& dT = alpha1_.time().deltaTValue() * dimensionedScalar( "dummy", dimensionSet(0,0,1,0,0,0,0), 1.0 );
	const dimensionedScalar& rho1 = twoPhaseProperties_.rho1();
	const dimensionedScalar& rho2 = twoPhaseProperties_.rho2();

	//Unlimited phase change heat
	Q_pc_ = InterfaceField*twoPhaseProperties_.rho()*twoPhaseProperties_.cp()*((T_-T_sat_)/dT);

volScalarField Q_pc_org = Q_pc_;

	//This limiting scheme only works for condensation!!!
	//Get cond/evap limits
	volScalarField LimCond = (1.0-alpha1_)*( rho2*h_lv_ / dT );
	//No evaporation on wall cells!
	volScalarField LimEvap = (1.0-WallField)*alpha1_*rho1*h_lv_ / dT;
	//Apply fluid limiting
	volScalarField Q_pc_fluid = neg(Q_pc_org)*max(Q_pc_org, -LimCond) + pos(Q_pc_org)*min(Q_pc_org, LimEvap) ;
	
	//Apply total limiting
	//Q_pc_ = neg(Q_pc_)*max(Q_pc_, -LimCond) + pos(Q_pc_)*min(Q_pc_, LimEvap) ;

	//This limiting scheme only works for evaporation!!!
	//Try time-based limiting (i.e. relative phase change rate can't exceed |1| per time step
	volScalarField PCV_fac = dT*(Q_pc_org / h_lv_)*( (scalar(1.0)/twoPhaseProperties_.rho2()) - (scalar(1.0)/twoPhaseProperties_.rho1()) );
	//Volume generation/sink based limited	
	volScalarField Q_pc_vol = Q_pc_org * mag( min( max(PCV_fac, -1.0), (1.0-WallField) ) );
	//Available fluid limiting

	//Composite limit
	Q_pc_ = neg(Q_pc_org)*max( max( Q_pc_org, Q_pc_fluid ), Q_pc_vol) + pos(Q_pc_org)*min( min( Q_pc_org, Q_pc_fluid ), Q_pc_vol);


//Print out Comparison of limits
/*
forAll (InterfaceField, j)
{
	if (InterfaceField[j] == 0)
	{  continue;  }
	else
	{
		Info<< "Cell " << j << ": Q_pc: " << Q_pc_[j] << ", Vol limit: " << Q_pc_vol[j] << ", alpha = " << alpha1_[j] << ", T = " << T_[j] << ", Fluid limit: " << Q_pc_fluid[j] << endl;
	}
}
*/

}


bool Foam::thermalPhaseChangeModels::LimitedFacePair::read(const dictionary& thermalPhaseChangeProperties)
{
	thermalPhaseChangeModel::read(thermalPhaseChangeProperties);

	return true;
}


// ************************************************************************* //
