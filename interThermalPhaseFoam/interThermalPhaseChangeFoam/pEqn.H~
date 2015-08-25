Info<<"Test 1"<< endl;

{
    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rAUf(fvc::interpolate(rAU));
    
Info<<"Test 2"<< endl;

    U = rAU*UEqn.H();

Info<<"Test 3"<< endl;

    surfaceScalarField phiU
    (
        "phiU",
        (fvc::interpolate(U) & mesh.Sf())
        + fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi)
    );

Info<<"Test 4"<< endl;

    adjustPhi(phiU, U, p_rgh);

Info<<"Test 5"<< endl;

    phi = phiU +
    (
        fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
      - ghf*fvc::snGrad(rho) 
    )*rAUf*mesh.magSf();

    Info<<"Test 6"<< endl;

    while (pimple.correctNonOrthogonal())
    {
        Info<<"Test 7"<< endl;

	//No need to correct p_rgh for phase change - already incorporated into div(phi) 
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phi) - phaseChangeModel->PCV()
        );

        Info<<"Test 8"<< endl;

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        Info<<"Test 9"<< endl;

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        Info<<"Test 10"<< endl;

	//This is -= (unlike phaseChangeFoam) because the sign of the governing equations is changed.
        if (pimple.finalNonOrthogonalIter())
        {
            phi -= p_rghEqn.flux();
        }
        Info<<"Test 11"<< endl;
    }

    U += rAU*fvc::reconstruct((phi - phiU)/rAUf);
    U.correctBoundaryConditions();

    Info<<"Test 12"<< endl;

    #include "continuityErrs.H"

    p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }
}
Info<<"Test 12"<< endl;
