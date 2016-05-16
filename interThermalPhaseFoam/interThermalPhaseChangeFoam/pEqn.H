
Switch CombinedHSPressure("no");
runTime.controlDict().readIfPresent("CombinedHSPressure", CombinedHSPressure);
const scalar CHSP = CombinedHSPressure ? 1.0 : 0;

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


	phi_c   = stfModel->phi_c( rAUf );		

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

	phi = phiU + stfModel->phi_c( rAUf ) + phi_rgh;


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
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phi) - PCV_accum
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
	p == p_rgh + (1.0-CHSP)*rho*gh + stfModel->pcap();    //Add in surface tension

	if (p_rgh.needReference())
	{
	    p += dimensionedScalar
	    (
	       "p",
	        p.dimensions(),
	        pRefValue - getRefCellValue(p, pRefCell)
	    );
	    p_rgh = p - (1.0-CHSP)*rho*gh - stfModel->pcap(); //Add in surface tension
	}

}





