
    //Update kEff
    kEff = 
          twoPhaseProperties.lambdaf()
        + fvc::interpolate( rho*twoPhaseProperties.cp()*turbulence->nut() );



    surfaceScalarField alphaEffRho
    (
         "alphaEff",
        fvc::interpolate(rho)*(   twoPhaseProperties.alphaf()
                                + fvc::interpolate( turbulence->nut() ) )
    ); 

    //Use limited alpha1:
    const volScalarField limAlpha1( min(max(alpha1, scalar(0)), scalar(1)) );

    //Latent heat accounted for in source term
    T.correctBoundaryConditions();
    H = ( (T - T_0)*(limAlpha1*rho1*cp1 + (1.0 - limAlpha1)*rho2*cp2) )/rho;

    //Scalar factor for artificial diffusion in the energy equation
    const scalar RelaxFac = 1.0;

    //Energy-temperature loop
    label nEnergyLoops(readLabel(pimple.dict().lookup("nEnergyLoops")));
    for (int EEqnCount=0; EEqnCount < nEnergyLoops; EEqnCount++)
    {

        //Form and solve the energy equation
        fvScalarMatrix EEqn
        (
            fvm::ddt(rho, H)
            + fvm::div(rhoPhi, H)
            - fvc::laplacian(kEff, T)
            - RelaxFac*(   fvm::laplacian(alphaEffRho, H)
                           - fvc::laplacian(alphaEffRho, H) )
            + phaseChangeModel->Q_pc()
        );

        EEqn.solve();
        //Now reevaluate T for the updated enthalpy fields
        T = T_0 + rho*H/(limAlpha1*rho1*cp1 + (1.0 - limAlpha1)*rho2*cp2);

     }


