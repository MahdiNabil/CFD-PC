{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");


    surfaceScalarField phir(phic*interface.nHatf());

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        surfaceScalarField phiAlpha
        (
            fvc::flux
            (
                phi,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, scalar(1) - alpha1, alpharScheme),
                alpha1,
                alpharScheme
            )
        );

        //Linear source term
        volScalarField Sp
        (
            IOobject
            (
                "Sp",
                runTime.timeName(),
                mesh
            ),
            divU
        );

        //alpha1 source term due to condensation and divergence of velocity
        volScalarField Su
        (
            IOobject
            (
                "Su",
                runTime.timeName(),
                mesh
            ),
            phaseChangeModel->alpha1Gen()
        );

        MULES::explicitSolve
        (
            geometricOneField(),
            alpha1,
            phi,
            phiAlpha,
            Sp,
            Su,
            1,
            0
        );

        rhoPhi +=
            (runTime.deltaT()/totalDeltaT)
           *(phiAlpha*(rho1 - rho2) + phi*rho2); 
    }

    Info<< "Liquid phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value()
        << endl;
}
