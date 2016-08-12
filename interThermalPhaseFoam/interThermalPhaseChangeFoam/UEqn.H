// Combined pressure model
Switch CombinedHSPressure("no");
runTime.controlDict().readIfPresent("CombinedHSPressure", CombinedHSPressure);
const scalar CHSP = CombinedHSPressure ? 1.0 : 0;

// Improved transport models
Switch ImprovedTransportBlending("no");
transportProperties.readIfPresent("ImprovedTransportBlending", ImprovedTransportBlending);
const scalar ImpTB = ImprovedTransportBlending ? 1.0 : 0.0;

fvVectorMatrix UEqn
(
      fvm::ddt(rho, U)
    + fvm::div(rhoPhi, U)
    + turbulence->divDevRhoReff(rho, U)
    //ASR - the default divDevRhoReff uses the viscosity values on cell centers
    // better to use the face values with this improved blending model
    + ImpTB*(   -fvm::laplacian(twoPhaseProperties.muf(), U)
               + fvm::laplacian(twoPhaseProperties.mu(), U) )
    - CHSP*rho*g
);

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                  stfModel->Fstff()                
                  - (1.0-CHSP)*ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );
    }




