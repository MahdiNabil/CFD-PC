surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(1, 0, -1, 0, 0), 0)
);

{
    label nAlphaCorr(readLabel(pimple.dict().lookup("nAlphaCorr")));

    label nAlphaSubCycles(readLabel(pimple.dict().lookup("nAlphaSubCycles")));

    // Compute limited phi
    surfaceScalarField phic(mag(phi/mesh.magSf()));
    
    phic = min(interface.cAlpha()*phic, max(phic));

    // Don't include PCV in the div term here, 
    // it is already there from the updated phi from the previous time step:
    volScalarField divU(fvc::div(phi) - phaseChangeModel->PCV());
    
    //This term is needed for an explicit correction to the mass flux (rhoPhi)
    dimensionedScalar totalDeltaT = runTime.deltaT();

    if (nAlphaSubCycles > 1)
    {
        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            #include "alphaEqn.H"
        }
    }
    else
    {
        #include "alphaEqn.H"
    }

    // The original interFoam code did not update the interface
    // curvature/forces if there were multiple outer PIMPLE correctors
    // We should do it even if that is the case:
    interface.correct();

    rho == alpha1*rho1 + (scalar(1) - alpha1)*rho2;
}
