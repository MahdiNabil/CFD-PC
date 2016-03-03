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

#include "DropwiseSGS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalPhaseChangeModels
{
    defineTypeNameAndDebug(DropwiseSGS, 0);
    addToRunTimeSelectionTable(thermalPhaseChangeModel, DropwiseSGS, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalPhaseChangeModels::DropwiseSGS::DropwiseSGS
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
   	Q_pc_sgs_
    (
        IOobject
        (
            "sgsPhaseChangeHeat",
            T_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh_,
		dimensionedScalar( "dummy", dimensionSet(1,-1,-3,0,0,0,0), 0 )
    ),
	qFlux_sgs_ // Declare the flux from integral as a volumeScalarField
    (
        IOobject
        (
            "sgsPhaseChangeHeatFlux",
            T_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh_,
		dimensionedScalar( "dummy", dimensionSet(1,0,-3,0,0,0,0), 0 )
    ),
	wet
    (
        IOobject
        (
            "wet",
            T_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar( "dummy", dimensionSet(0,0,0,0,0,0,0), 0 )
	),
	faceTime
    (
        IOobject
        (
            "faceTime",
            T_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar( "dummy", dimensionSet(0,0,1,0,0,0,0), 0 )
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
	//Read in the cond/evap int. thresholds and other thermal properties
	thermalPhaseChangeProperties_.lookup("CondThresh") >> CondThresh;
	thermalPhaseChangeProperties_.lookup("EvapThresh") >> EvapThresh;
	thermalPhaseChangeProperties_.lookup("RelaxFac") >> RelaxFac;
	thermalPhaseChangeProperties_.lookup("Gamma") >> Gamma;
	thermalPhaseChangeProperties_.lookup("C_1") >> C_1;
	thermalPhaseChangeProperties_.lookup("C_2") >> C_2;
	thermalPhaseChangeProperties_.lookup("R_g") >> R_g;
	thermalPhaseChangeProperties_.lookup("C_3") >> C_3;

	//Set other constant fluid properties
	const IOdictionary& transportProperties = mesh_.lookupObject<IOdictionary>("transportProperties");
	const dictionary& phase1Properties(transportProperties.subDict("phase1"));
	sigma = dimensionedScalar(transportProperties.lookup("sigma")).value();
	k_l = dimensionedScalar(phase1Properties.lookup("lambda")).value();
	
		
Info<< sigma << endl;
		
	correct();
	//GSLIntegral();
	decayHt();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::thermalPhaseChangeModels::DropwiseSGS::calcQ_pc()
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

	//Reset all Q_pc to 0
	Q_pc_ = dimensionedScalar( "dummy", dimensionSet(1,-1,-3,0,0,0,0), 0 );
	Q_pc_sgs_ = dimensionedScalar( "dummy", dimensionSet(1,-1,-3,0,0,0,0), 0 );

	//Compute some helpful props:
	//For some reason dT is dimensionless
	const dimensionedScalar& dT = alpha1_.time().deltaTValue() * dimensionedScalar( "dummy", dimensionSet(0,0,1,0,0,0,0), 1.0 );
	const dimensionedScalar& rho1 = twoPhaseProperties_.rho1();
	const dimensionedScalar& rho2 = twoPhaseProperties_.rho2();

	//Unlimited phase change heat
	Q_pc_ = InterfaceField_*twoPhaseProperties_.rho()*twoPhaseProperties_.cp()*((T_-T_sat_)/dT);

	//Fluid availability limits
	//Get cond/evap limits
	volScalarField LimCond = (1.0-alpha1_)*( rho2*h_lv_ / dT );
	//No evaporation on wall cells!
	volScalarField LimEvap = (1.0-WallField)*alpha1_*rho1*h_lv_ / dT;

	//Apply fluid limiting
	volScalarField Q_pc_fluid = neg(Q_pc_)*max(Q_pc_, -LimCond) + pos(Q_pc_)*min(Q_pc_, LimEvap) ;

	//Volume-based limiting (i.e. relative phase change rate can't exceed |1| per time step
	volScalarField PCV_fac = dT*(Q_pc_ / h_lv_)*( (scalar(1.0)/twoPhaseProperties_.rho2()) - (scalar(1.0)/twoPhaseProperties_.rho1()) );

	//Again, allow regular evap on wall	
	volScalarField Q_pc_vol = Q_pc_ * mag( min( max(1.0/(PCV_fac+SMALL), -1.0), (1.0-WallField) ) );

	//Composite limit
	Q_pc_ = neg(Q_pc_)*max( max( Q_pc_, Q_pc_fluid ), Q_pc_vol) + pos(Q_pc_)*min( min( Q_pc_, Q_pc_fluid ), Q_pc_vol);

	//Under relax phase change rate per user specification
	Q_pc_ = RelaxFac * Q_pc_;

	//Now apply subgridscale model on the wall patches
	surfaceScalarField alpha1f = fvc::interpolate(alpha1_); //Alpha1 on faces

	//- Applying the subGrid Model
	forAll(mesh_.boundary(),pI)
	{
		if(isA<wallFvPatch>(mesh_.boundary()[pI]))
		{
			scalarField& faceTimePatch = faceTime.boundaryField()[pI];
			faceTimePatch += dT.value();
			scalarField& wetPatch = wet.boundaryField()[pI];  
			scalarField& alphaPatch = alpha1f.boundaryField()[pI];
			
			forAll(wetPatch, fI) 
			{
				wetPatch[fI] = ( (alphaPatch[fI] > 0.9) || ((wetPatch[fI] == 1.0) && (alphaPatch[fI] > 0.1)) ) ? 1.0 : 0.0;
			}

			faceTimePatch = (1.0-wetPatch)*faceTimePatch;

			scalarField& qFlux_sgsPatch = qFlux_sgs_.boundaryField()[pI];
			qFlux_sgsPatch = max(0.0, (-1.0/faceTimePatch + C_4*(faceTimePatch/exp(faceTimePatch)) + C_5)*(C_6/C_7)); 	
		}	
	}
	//Calculate volumetric SGS phase change rate in the cell	
	Q_pc_sgs_ = fvc::surfaceIntegrate( (1-alpha1f)*qFlux_sgs_*mesh_.magSf() );
}

double Foam::thermalPhaseChangeModels::DropwiseSGS::GSLFunction(double r, void *params)
{
	
	struct GSLFunction_params * p = static_cast<struct GSLFunction_params *>(params);
	double Gamma1 = (p->Gamma);      // Cp/Cv
	double h_lv1 = (p->h_lv);        // latent heat of vaporization
	double k_l1 = (p->k_l);          // thermal conductivity of liquid phase
	double rho_l1 = (p->rho_l);      // density of liquid phase
	double rho_v1 = (p->rho_v);      // density of vapor phase
	double R_g1 = (p->R_g);          // specific ideal gas constant
	double sigma1 = (p->sigma);      // surface tension
	double T_sat1 = (p->T_sat);      // temperature of the saturated vapor phase
	double T_w1 = (p->T_w);          // temperature of the wall
	double C_11 = (p->C_1);          // constant 1 from Rose 1998
	double C_21 = (p->C_2);          // constant 2 multiplied by the correction factor from Rose 1998

//Info<< "T_f: " << T_w1;
	double GSLFunction = pow(r, -2.0/3.0) * ((T_sat1-T_w1) - (2*sigma1*T_sat1/(r*rho_l1*h_lv1))) / 
	(C_11*r/k_l1 + C_21*T_sat1*(Gamma1+1)/(h_lv1*h_lv1*rho_v1*(Gamma1-1))*pow(R_g1*T_sat1/(2*M_PI),0.5));
	
	return  GSLFunction;	
}

void Foam::thermalPhaseChangeModels::DropwiseSGS::GSLIntegral()
{
	qFlux_sgs_ = dimensionedScalar( "dummy", dimensionSet(1,0,-3,0,0,0,0), 0 ); // set the flux to 0 initialy
	
	// some interpolated variables that would be used in the calculation
	const dimensionedScalar& rho_l = twoPhaseProperties_.rho1();  // Density of the liquid face
	const dimensionedScalar& rho_v = twoPhaseProperties_.rho2();  // Density of the vapor phase
	const surfaceScalarField Tf = fvc::interpolate(T_);	    // Temperature interpolated to cell faces
	
	forAll( mesh_.boundary(), pI )
	{
		if( isA<wallFvPatch>( mesh_.boundary()[pI] ) )    
		{
			const fvPatch& fPatch = mesh_.boundary()[pI]; 
			const scalarField WallFaceAreas = fPatch.magSf();
			const scalarField& Wall_T1 = Tf.boundaryField()[pI];
			const scalarField rMax = C_3*Foam::sqrt( WallFaceAreas ); // Maximum SGS drop size radius
			scalarField rMin = (2.0*sigma*T_sat_.value()) / (  max( (T_sat_.value()-Wall_T1), SMALL)   *rho_l.value()*h_lv_.value()); // Minimum SGS drop size radius
			rMin = min(rMin, rMax);
			
			gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);  // A sufficiently large number	
			forAll(fPatch, fI)  // Now loop over all the cell Faces of the patch and calculate the heat flux 
			{
				double result, error;
				struct GSLFunction_params params = { Gamma, h_lv_.value(), k_l, rho_l.value(), rho_v.value(), R_g, sigma, T_sat_.value(), Wall_T1[fI], C_1, C_2 };
				gsl_function F;
				F.function = &Foam::thermalPhaseChangeModels::DropwiseSGS::GSLFunction;
				F.params = &params;	
				gsl_integration_qags (&F, rMin[fI], rMax[fI], 0, 1e-7, 1000, w, &result, &error);
				qFlux_sgs_.boundaryField()[pI][fI] = -result * pow(rMax[fI], -1.0/3.0) / 3;
			}
			gsl_integration_workspace_free (w);
		}
	}
}

void Foam::thermalPhaseChangeModels::DropwiseSGS::decayHt()
{
	qFlux_sgs_ = dimensionedScalar( "dummy", dimensionSet(1,0,-3,0,0,0,0), 0 ); // set the flux to 0 initialy

	// some interpolated variables that would be used in the calculation
	const surfaceScalarField Tf = fvc::interpolate(T_);	    // Temperature interpolated to cell faces
	
	forAll( mesh_.boundary(), pI )
	{
		if( isA<wallFvPatch>( mesh_.boundary()[pI] ) )    
		{
			const fvPatch& fPatch = mesh_.boundary()[pI]; 
			const scalarField WallFaceAreas = fPatch.magSf();
			const scalarField& Wall_T1 = Tf.boundaryField()[pI];
			scalarField decayHtConst = -(T_sat_.value()-Wall_T1) * 5e5; // we have chosen 5e5 arbitrarily
			
			forAll(fPatch, fI)  // Now loop over all the cell Faces of the patch and calculate the heat flux 
			{
				qFlux_sgs_.boundaryField()[pI][fI] = decayHtConst[fI];
			}
		}
	}
}

bool Foam::thermalPhaseChangeModels::DropwiseSGS::read(const dictionary& thermalPhaseChangeProperties)
{
	thermalPhaseChangeModel::read(thermalPhaseChangeProperties);

	//Read in the cond/evap int. thresholds
	thermalPhaseChangeProperties_.lookup("CondThresh") >> CondThresh;
	thermalPhaseChangeProperties_.lookup("EvapThresh") >> EvapThresh;
	thermalPhaseChangeProperties_.lookup("RelaxFac") >> RelaxFac;

	return true;
}


// ************************************************************************* //
