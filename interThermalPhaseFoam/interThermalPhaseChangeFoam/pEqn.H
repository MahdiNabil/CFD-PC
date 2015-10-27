
Switch CombinedHSPressure("no");
runTime.controlDict().readIfPresent("CombinedHSPressure", CombinedHSPressure);
const scalar CHSP = CombinedHSPressure ? 1.0 : 0;

Switch ModifiedST("no");
runTime.controlDict().readIfPresent("ModifiedST", ModifiedST);
const scalar MST = ModifiedST ? 1.0 : 0;

//Info<< "CHSP: " << CHSP << endl;
//Info<< "MST : " << MST  << endl;

{
    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rAUf(fvc::interpolate(rAU));

    U = rAU*UEqn.H();

    surfaceScalarField phiU
    (
        "phiU",
        (fvc::interpolate(U) & mesh.Sf())
        + fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi)
    );
   fvOptions.makeRelative(phiU);

    adjustPhi(phiU, U, p_rgh);

// Added Condition for Modified Surface Tension

		//surface tension flux (capillary)
		//surfaceScalarField phi_c( ( (fvc::interpolate(fc) & mesh.Sf()) - (fvc::snGrad(pc) * mesh.magSf()) )*rAUf );

/*
	if ( ModifiedST )
	{
		phiU += phi_c - (1.0-CHSP)*ghf*fvc::snGrad(rho)*rAUf*mesh.magSf();;
	}
	else
	{
		phi = phiU +
		(
		    fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
		  - (1.0-CHSP)*ghf*fvc::snGrad(rho) 
		)*rAUf*mesh.magSf();
	}
*/

	//Calculate the flux from surface tension:
	surfaceScalarField phi_c
    (
        IOobject
        (
            "phi_c",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("flux", dimLength*dimLength*dimLength/dimTime, 0)
    );


	if (ModifiedST) //Using SSF model
	{
		phi_c =  ( (fcf * mesh.magSf()) - (fvc::snGrad(pc) * mesh.magSf()) )*rAUf;

		//Apply limiting (filtering)
		const dimensionedScalar dummyFlux("dummyFlux", dimensionSet(0,3,-1,0,0,0,0), 1.0);
		dimensionedScalar phi_c_thresh( 0.01* gAverage( mag(phi_c.field()) ) * dummyFlux );
		phi_c = phi_c - max( min(phi_c, phi_c_thresh), -phi_c_thresh );
	}
	else //Using default surface tension model
	{
		phi_c = fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1) * rAUf*mesh.magSf();
	}



	//Calculate the flux from gravitational forces
	surfaceScalarField phi_rgh
    (
        IOobject
        (
            "phi_rgh",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("flux", dimLength*dimLength*dimLength/dimTime, 0)
    );


	if (!CombinedHSPressure)
	{
		phi_rgh =  -ghf*fvc::snGrad(rho) * rAUf*mesh.magSf();
	}


	phi = phiU + phi_c + phi_rgh;


   // Update the fixedFluxPressure BCs to ensure flux consistency
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_rgh.boundaryField(),
        (
            phiU.boundaryField()
          - fvOptions.relative(mesh.Sf().boundaryField() & U.boundaryField())
        )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
    );

    while (pimple.correctNonOrthogonal())
    {

	//No need to correct p_rgh for phase change - already incorporated into div(phi) 
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phi) - phaseChangeModel->PCV()
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

	//This is -= (unlike phaseChangeFoam) because the sign of the governing equations is changed.
        if (pimple.finalNonOrthogonalIter())
        {
            phi -= p_rghEqn.flux();
        }

    }

// Added Condition for Modified Surface Tension Case 
	U += rAU*fvc::reconstruct((phi - phiU)/rAUf);

    U.correctBoundaryConditions();


    #include "continuityErrs.H"
	if ( ModifiedST )
	{
		p == p_rgh + (1.0-CHSP)*rho*gh + pc;    //Add in surface tension

		if (p_rgh.needReference())
		{
		    p += dimensionedScalar
		    (
		       "p",
		        p.dimensions(),
		        pRefValue - getRefCellValue(p, pRefCell)
		    );
		    p_rgh = p - (1.0-CHSP)*rho*gh - pc; //Add in surface tension
		}
	}
else
	{
		p == p_rgh + (1.0-CHSP)*rho*gh;    

		if (p_rgh.needReference())
		{
		    p += dimensionedScalar
		    (
		       "p",
		        p.dimensions(),
		        pRefValue - getRefCellValue(p, pRefCell)
		    );
		    p_rgh = p - (1.0-CHSP)*rho*gh; 
		}
	}
}






