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

#include "SST.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceTensionForceModels
{
    defineTypeNameAndDebug(SST, 0);
    addToRunTimeSelectionTable(surfaceTensionForceModel, SST, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceTensionForceModels::SST::SST
(
    const word& name,
    const dictionary& surfaceTensionForceProperties,
    const interfaceProperties& interface,
    const volScalarField& alpha1
)
:
    surfaceTensionForceModel
    (
        name,
        surfaceTensionForceProperties,
        interface,
        alpha1
    ),
    mesh_(alpha1.mesh()),
    Cpc( surfaceTensionForceProperties.lookupOrDefault<scalar>("Cpc", 0.5) ),
    phi_c_thresholdFactor
    (
        surfaceTensionForceProperties.lookupOrDefault<scalar>
        ("ThresholdFactor", 0.01) 
    ),
    filterFactor
    (
        surfaceTensionForceProperties.lookupOrDefault<scalar>
        ("NormalFilterFactor", 0.9)
    ),
    pc
    (
        IOobject
        (
            "pc",
            alpha1_.time().timeName(),
            alpha1.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1.mesh()
    ),
    Fstffv
    (
        IOobject
        (
            "SurfaceTensionForce",
            alpha1_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar( "dummy", dimensionSet(1,-2,-2,0,0,0,0), 0 )
    ),
    //Sharp surface force (capillary forces on cell faces)
    fc
    (
        IOobject
        (
            "fc",
            alpha1_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector
        ("fc0", dimMass/(dimLength*dimLength*dimTime*dimTime), vector(0,0,0))
    ),
    //Sharp surface force (capillary forces on cell faces)
    fcf
    (
        IOobject
        (
            "fcf",
            alpha1_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar
        ("fc0", dimMass/(dimLength*dimLength*dimTime*dimTime), 0)
    ),
    //Sharp surface force (capillary forces on cell faces)
    fcf_filter
    (
        IOobject
        (
            "fcf_filter",
            alpha1_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar
        ("fc0", dimMass/(dimLength*dimLength*dimTime*dimTime), 0)
    )

{

    //Set reference pressure stuff: 
    setRefCell
    (
        pc,
        pc,
        mesh_.solutionDict().subDict("PIMPLE"),
        pcRefCell,
        pcRefValue
    );
    pcRefValue = 0;

    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::surfaceTensionForceModels::SST::correct()
{
    //Step 1: smoothing the phase fraction field (2 passes)
    scalar CSK = 0.5;
    //pass 1
    volScalarField alpha1s =
          CSK*(fvc::average(fvc::interpolate(alpha1_))) 
        + (1.0 - CSK) * alpha1_;
        
    //pass 2
    alpha1s = 
        CSK*(fvc::average(fvc::interpolate(alpha1s))) + (1.0 - CSK)*alpha1s;

    //Step 2: initialize interface curvature Kappa
    const volVectorField gradAlpha = fvc::grad(alpha1s);
    
    const dimensionedScalar 
    deltaN("deltaN", dimensionSet(0,-1,0,0,0,0,0), 1E-16);
    
    const volVectorField ns(gradAlpha/(mag(gradAlpha) + deltaN));
    volScalarField K = fvc::div(ns);

    //Step 3: smooth curvature field (2 passes)
    volScalarField w = Foam::sqrt( mag( alpha1_*(1.0 - alpha1_) ) + 1.0E-6);
    volScalarField factor = 2.0*Foam::sqrt( mag( alpha1_*(1.0 - alpha1_) ) );
    //pass 1
    volScalarField Ks_star = 
        fvc::average(fvc::interpolate(K*w))/fvc::average(fvc::interpolate(w));
    volScalarField Ks = factor*K + (1.0 - factor)*Ks_star;

    //pass 2
    Ks_star = 
        fvc::average(fvc::interpolate(Ks*w))/fvc::average(fvc::interpolate(w));
    Ks = factor*K + (1.0 - factor)*Ks_star;

    //Step 4: compute smoothed curvature on faces
    surfaceScalarField Kf = fvc::interpolate(w*Ks)/fvc::interpolate(w); 

    //Step 5: compute interface delta function from sharpened interface field
    volScalarField alpha1_pc =
        1.0/(1.0-Cpc)*(min( max(alpha1_,Cpc/2.0), (1.0-Cpc/2.0) ) - Cpc/2.0);

    surfaceScalarField deltasf = fvc::snGrad(alpha1_pc);

    //Step 6: compute surface tension force on faces
    const dimensionedScalar dummyA("DummyA", dimensionSet(0,2,0,0,0,0,0), 1.0);

    fcf = -interface_.sigma()*Kf*deltasf;

    //Step 7: filter surface tension forces to only be normal to interfaces
    fcf_filter =
          (deltasf/(mag(deltasf)+deltaN)) 
         *(
              filterFactor*fcf_filter 
            + (1.0-filterFactor)
            *(   fvc::interpolate( fvc::grad(pc)
               - (fvc::grad(pc) & ns)*ns ) & (mesh_.Sf()/mesh_.magSf()) )
          );
        
        
        
        
    fcf = fcf - fcf_filter;

    //Step 8: produce fc on cell centers
    fc = fvc::average(fcf*mesh_.Sf()/mesh_.magSf());

    //Step 9: solve for capillary pressure field:
    fvScalarMatrix pcEqn
    (
        fvm::laplacian(pc) == fvc::div(fcf*mesh_.magSf() )
    );

    //Get reference to pc
    //Reference point should be located outside of the film thickness 
    //(condensate) in order to have stable pressure solution runtime
    if (pc.needReference())
    {
        pcEqn.setReference(pcRefCell, getRefCellValue(pc, pcRefCell));
    }

    pcEqn.solve();

    //Fstffv = fcf - fvc::snGrad(pc) ;
    Fstffv = fvc::snGrad(pc);
}


Foam::tmp<Foam::surfaceScalarField> 
Foam::surfaceTensionForceModels::SST::
phi_c(const surfaceScalarField& rAUf_) const
{
    const surfaceScalarField phi_c_i( Fstffv * rAUf_ * mesh_.magSf() );

    //Apply limiting (filtering)
    const dimensionedScalar
    dummyFlux("dummyFlux", dimensionSet(0,3,-1,0,0,0,0), 1.0);

    dimensionedScalar phi_c_thresh
    ( 
          phi_c_thresholdFactor*dummyFlux*gSum( mag(phi_c_i.field()) ) 
        / gSum( SMALL + pos( mag(phi_c_i.field()) - SMALL ) ) 
    );

    //Return filtered phi_c
    return tmp<surfaceScalarField>
    ( phi_c_i - max( min(phi_c_i, phi_c_thresh), -phi_c_thresh ) );
}

bool Foam::surfaceTensionForceModels::SST::
read(const dictionary& surfaceTensionForceProperties)
{
    surfaceTensionForceModel::read(surfaceTensionForceProperties);

    return true;
}


// ************************************************************************* //
