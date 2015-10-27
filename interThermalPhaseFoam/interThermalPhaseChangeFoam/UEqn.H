
Switch CombinedHSPressure("no");
runTime.controlDict().readIfPresent("CombinedHSPressure", CombinedHSPressure);
const scalar CHSP = CombinedHSPressure ? 1.0 : 0;

Switch ModifiedST("no");
runTime.controlDict().readIfPresent("ModifiedST", ModifiedST);
const scalar MST = ModifiedST ? 1.0 : 0;

fvVectorMatrix UEqn
(
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      + turbulence->divDevRhoReff(rho, U)
      - CHSP*rho*g
);

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        if ( ModifiedST )
		{
        solve
        (
            UEqn
         ==
	    (MST*fc) + 
            fvc::reconstruct
            (
                (
				  - MST*fvc::snGrad(pc)
                  - (1.0-CHSP)*ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );
		}
		else
		{
		solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
                  - (1.0-CHSP)*ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );
		}
    }




