
	//Update kEff
	kEff = twoPhaseProperties.lambdaf() + fvc::interpolate( rho*twoPhaseProperties.cp()*turbulence->nut() );

    surfaceScalarField alphaEffRho
    (
         "alphaEff",
        fvc::interpolate(rho) * ( twoPhaseProperties.alphaf() + fvc::interpolate( turbulence->nut() ) )
    );
 
	//Use limited alpha1:
    const volScalarField limAlpha1( min(max(alpha1, scalar(0)), scalar(1)) );

	//Check how many Energy equation sub-cycles to use
	//label nEnergySubCycles( pimple.dict().lookupOrDefault<int>("nEnergySubCycles", 1) );

	//Zero out accumulated fields
	Q_pc_accum       = dimensionedScalar("0", dimensionSet(1, -1, -3, 0, 0), 0);
	PCV_accum        = dimensionedScalar("0", dimensionSet(0,  0, -1, 0, 0), 0);
	alpha1Gen_accum  = dimensionedScalar("0", dimensionSet(0,  0, -1, 0, 0), 0);

	//Keep a reference to the global time step:
	dimensionedScalar totalDeltaT = runTime.deltaT();

//Info<< "T mean 1: " << gAverage( T.internalField() ) << endl;
//Info<< "H mean 1: " << gAverage( H.internalField() ) << endl;


	//Latent heat accounted for in source term
	T.correctBoundaryConditions();
	H = ( (T - T_0)*(limAlpha1*rho1*cp1 + (1-limAlpha1)*rho2*cp2) )/rho;
Info << "number of Energy Sub Cycles: " << numESCycles << endl;
	for ( subCycle<volScalarField> EnergySubCycle(H, numESCycles); !(++EnergySubCycle).end(); )
	{

		//Latent heat accounted for in source term
		T.correctBoundaryConditions();
		H = ( (T - T_0)*(limAlpha1*rho1*cp1 + (1-limAlpha1)*rho2*cp2) )/rho;


		//Update phase change model at new temperature field values
		phaseChangeModel->correct();

//Info<< "T mean 2: " << gAverage( T.internalField() ) << endl;
//Info<< "H mean 2: " << gAverage( H.internalField() ) << endl;

		//Scalar factor for artificial diffusion in the energy equation
		scalar ChillaxFac = 1.0;

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
				- ChillaxFac*( fvm::laplacian(alphaEffRho, H) - fvc::laplacian(alphaEffRho, H) )
				+ phaseChangeModel->Q_pc()
			);
			EEqn.solve();
			//Now reevaluate T for the updated enthalpy fields

//Info<< "T mean 3: " << gAverage( T.internalField() ) << endl;
//Info<< "H mean 3: " << gAverage( H.internalField() ) << endl;


			T = T_0 + rho*H/( limAlpha1*rho1*cp1 + (1-limAlpha1)*rho2*cp2 );

//Info<< "T mean 4: " << gAverage( T.internalField() ) << endl;
//Info<< "H mean 4: " << gAverage( H.internalField() ) << endl;

		}

		//Accumulate Q_PC_accum, PCV_accum, alpha1Gen
		Q_pc_accum       += (runTime.deltaT()/totalDeltaT)*phaseChangeModel->Q_pc();
		PCV_accum        += (runTime.deltaT()/totalDeltaT)*phaseChangeModel->PCV();
		alpha1Gen_accum  += (runTime.deltaT()/totalDeltaT)*phaseChangeModel->alpha1Gen();
	}
