{
    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rAUf(fvc::interpolate(rAU));

    U = rAU*UEqn.H();
    surfaceScalarField phiU
    (
        "phiU",
        (fvc::interpolate(U) & mesh.Sf())
        + fvc::ddtCorr(U, phi)
    );

    adjustPhi(phiU, U, p_rgh);



    phi = phiU +
    (
        fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
      - ghf*fvc::snGrad(rho) 
    )*rAUf*mesh.magSf();

    Info<<"Test 1"<< endl;
    while (pimple.correctNonOrthogonal())
    {
        Info<<"Test 2"<< endl;

	//No need to correct p_rgh for phase change - already incorporated into div(phi) 
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phi) - phaseChangeModel->PCV()
        );

        Info<<"Test 3"<< endl;

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        Info<<"Test 4"<< endl;

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        Info<<"Test 5"<< endl;

	//This is -= (unlike phaseChangeFoam) because the sign of the governing equations is changed.
        if (pimple.finalNonOrthogonalIter())
        {
            phi -= p_rghEqn.flux();
        }
        Info<<"Test 6"<< endl;
    }

    U += rAU*fvc::reconstruct((phi - phiU)/rAUf);
    U.correctBoundaryConditions();

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
