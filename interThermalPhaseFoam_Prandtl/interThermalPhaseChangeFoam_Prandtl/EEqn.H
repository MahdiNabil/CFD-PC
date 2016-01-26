
	//First get the thermal conductivities from the dictionary
	//const dictionary& phase2Properties(transportProperties.subDict("phase2"));
	//lambda0_ = dimensionedScalar(phase2Properties.lookup("lambda")).value();
	transportProperties.lookup("lambda0") >> lambda0_;
	transportProperties.lookup("lambda1") >> lambda1_;
	transportProperties.lookup("lambda2") >> lambda2_;
	transportProperties.lookup("lambda3") >> lambda3_;
	transportProperties.lookup("lambda4") >> lambda4_;
	transportProperties.lookup("lambda5") >> lambda5_;

	//Update kEff
	//kEff = fvc::interpolate( twoPhaseProperties.lambda() ) + fvc::interpolate( rho*twoPhaseProperties.cp()*turbulence->nut() );
	kEff1 = fvc::interpolate( alpha1*lambda1_ + (1-alpha1)*lambda0_ ) + fvc::interpolate( rho*twoPhaseProperties.cp()*turbulence->nut() );
	kEff2 = fvc::interpolate( alpha1*lambda2_ + (1-alpha1)*lambda0_ ) + fvc::interpolate( rho*twoPhaseProperties.cp()*turbulence->nut() );
	kEff3 = fvc::interpolate( alpha1*lambda3_ + (1-alpha1)*lambda0_ ) + fvc::interpolate( rho*twoPhaseProperties.cp()*turbulence->nut() );
	kEff4 = fvc::interpolate( alpha1*lambda4_ + (1-alpha1)*lambda0_ ) + fvc::interpolate( rho*twoPhaseProperties.cp()*turbulence->nut() );
	kEff5 = fvc::interpolate( alpha1*lambda5_ + (1-alpha1)*lambda0_ ) + fvc::interpolate( rho*twoPhaseProperties.cp()*turbulence->nut() );

    surfaceScalarField alphaEffRho
    (
         "alphaEff",
        fvc::interpolate(rho) * ( twoPhaseProperties.alphaf() + fvc::interpolate( turbulence->nut() ) )
    );

	//Use limited alpha1:
    const volScalarField limAlpha1( min(max(alpha1, scalar(0)), scalar(1)) );

	//Latent heat accounted for in source term
	T1.correctBoundaryConditions();
	H1 = ( (T1 - T_0)*(limAlpha1*rho1*cp1 + (1-limAlpha1)*rho2*cp2) )/rho;
	T2.correctBoundaryConditions();
	H2 = ( (T2 - T_0)*(limAlpha1*rho1*cp1 + (1-limAlpha1)*rho2*cp2) )/rho;
	T3.correctBoundaryConditions();
	H3 = ( (T3 - T_0)*(limAlpha1*rho1*cp1 + (1-limAlpha1)*rho2*cp2) )/rho;
	T4.correctBoundaryConditions();
	H4 = ( (T4 - T_0)*(limAlpha1*rho1*cp1 + (1-limAlpha1)*rho2*cp2) )/rho;
	T5.correctBoundaryConditions();
	H5 = ( (T5 - T_0)*(limAlpha1*rho1*cp1 + (1-limAlpha1)*rho2*cp2) )/rho;

	//Scalar factor for artificial diffusion in the energy equation
	scalar ChillaxFac = 1.0;

    //Energy-temperature loops
    label nEnergyLoops(readLabel(pimple.dict().lookup("nEnergyLoops")));
    for (int EEqnCount=0; EEqnCount < nEnergyLoops; EEqnCount++)
    {
        //Form and solve the energy equation
        fvScalarMatrix EEqn
        (
            fvm::ddt(rho, H1)
            + fvm::div(rhoPhi, H1)
            - fvc::laplacian(kEff1, T1)
			- ChillaxFac*( fvm::laplacian(alphaEffRho, H1) - fvc::laplacian(alphaEffRho, H1) )
			+ phaseChangeModel1->Q_pc()
        );
        EEqn.solve();
        //Now reevaluate T for the updated enthalpy fields
        T1 = T_0 + rho*H1/( limAlpha1*rho1*cp1 + (1-limAlpha1)*rho2*cp2 );
     }

 	for (int EEqnCount=0; EEqnCount < nEnergyLoops; EEqnCount++)
    {
        //Form and solve the energy equation
        fvScalarMatrix EEqn
        (
            fvm::ddt(rho, H2)
            + fvm::div(rhoPhi, H2)
            - fvc::laplacian(kEff2, T2)
			- ChillaxFac*( fvm::laplacian(alphaEffRho, H2) - fvc::laplacian(alphaEffRho, H2) )
			+ phaseChangeModel2->Q_pc()
        );
        EEqn.solve();
        //Now reevaluate T for the updated enthalpy fields
        T2 = T_0 + rho*H2/( limAlpha1*rho1*cp1 + (1-limAlpha1)*rho2*cp2 );
     }
 
	for (int EEqnCount=0; EEqnCount < nEnergyLoops; EEqnCount++)
    {
        //Form and solve the energy equation
        fvScalarMatrix EEqn
        (
            fvm::ddt(rho, H3)
            + fvm::div(rhoPhi, H3)
            - fvc::laplacian(kEff3, T3)
			- ChillaxFac*( fvm::laplacian(alphaEffRho, H3) - fvc::laplacian(alphaEffRho, H3) )
			+ phaseChangeModel3->Q_pc()
        );
        EEqn.solve();
        //Now reevaluate T for the updated enthalpy fields
        T3 = T_0 + rho*H3/( limAlpha1*rho1*cp1 + (1-limAlpha1)*rho2*cp2 );
     }

	 for (int EEqnCount=0; EEqnCount < nEnergyLoops; EEqnCount++)
    {
        //Form and solve the energy equation
        fvScalarMatrix EEqn
        (
            fvm::ddt(rho, H4)
            + fvm::div(rhoPhi, H4)
            - fvc::laplacian(kEff4, T4)
			- ChillaxFac*( fvm::laplacian(alphaEffRho, H4) - fvc::laplacian(alphaEffRho, H4) )
			+ phaseChangeModel4->Q_pc()
        );
        EEqn.solve();
        //Now reevaluate T for the updated enthalpy fields
        T4 = T_0 + rho*H4/( limAlpha1*rho1*cp1 + (1-limAlpha1)*rho2*cp2 );
     }

	 for (int EEqnCount=0; EEqnCount < nEnergyLoops; EEqnCount++)
    {
        //Form and solve the energy equation
        fvScalarMatrix EEqn
        (
            fvm::ddt(rho, H5)
            + fvm::div(rhoPhi, H5)
            - fvc::laplacian(kEff5, T5)
			- ChillaxFac*( fvm::laplacian(alphaEffRho, H5) - fvc::laplacian(alphaEffRho, H5) )
			+ phaseChangeModel5->Q_pc()
        );
        EEqn.solve();
        //Now reevaluate T for the updated enthalpy fields
        T5 = T_0 + rho*H5/( limAlpha1*rho1*cp1 + (1-limAlpha1)*rho2*cp2 );
     }


