//Implementation based on work by Raeini et al. (2012) and Roberto Castilla Lopez (2013)


{
	//Step 1: smoothing the phase fraction field (2 passes)
	scalar CSK = 0.5;
	//pass 1
	volScalarField alpha1s = CSK * (fvc::average(fvc::interpolate(alpha1))) + (1.0 - CSK) * alpha1;
	//pass 2
	alpha1s = CSK * (fvc::average(fvc::interpolate(alpha1s))) + (1.0 - CSK) * alpha1s;

	//Step 2: initialize interface curvature Kappa
	const volVectorField gradAlpha = fvc::grad(alpha1s);
	const dimensionedScalar deltaN("deltaN", dimensionSet(0,-1,0,0,0,0,0), 1E-16);
	const volVectorField ns(gradAlpha/(mag(gradAlpha) + deltaN));
	volScalarField K = fvc::div(ns);

	//Step 3: smooth curvature field (2 passes)
	volScalarField w = Foam::sqrt( mag( alpha1*(1.0 - alpha1) ) + 1.0E-6);
	volScalarField factor = 2.0*Foam::sqrt( mag( alpha1*(1.0 - alpha1) ) );
	//pass 1
	volScalarField Ks_star = fvc::average(fvc::interpolate(K*w))/fvc::average(fvc::interpolate(w));
	volScalarField Ks = factor * K + (1.0 - factor) * Ks_star;
	//pass 2
	Ks_star = fvc::average(fvc::interpolate(Ks*w))/fvc::average(fvc::interpolate(w));
	Ks = factor * K + (1.0 - factor) * Ks_star;
	//Step 4: compute smoothed curvature on faces
	surfaceScalarField Kf = fvc::interpolate(w*Ks)/fvc::interpolate(w);	

	//Step 5: compute interface delta function from sharpened interface field
	scalar Cpc = 0.5;
	volScalarField alpha1_pc = 1.0/(1.0-Cpc) * (min( max(alpha1,Cpc/2.0), (1.0-Cpc/2.0) ) - Cpc/2.0);
	surfaceScalarField deltasf = fvc::snGrad(alpha1_pc);

	//Step 6: compute surface tension force on faces
	const dimensionedScalar dummyA("DummyA", dimensionSet(0,2,0,0,0,0,0), 1.0);
	//surfaceScalarField fcf = -interface.sigma()*Kf*deltasf;
	fcf = -interface.sigma()*Kf*deltasf;
	//Step 7: filter surface tension forces to only be normal to interfaces - not 100% sure if sure be dotted with face normals or interface normals...
	const scalar filterRelax = 0.9;
	fcf_filter = (deltasf/(mag(deltasf)+deltaN)) * ( filterRelax*fcf_filter + (1.0-filterRelax)*( fvc::interpolate( fvc::grad(pc) - (fvc::grad(pc) & ns)*ns ) & (mesh.Sf()/mesh.magSf()) ) );
	fcf = fcf - fcf_filter;

	//Step 8: produce fc on cell centers
	fc = fvc::average(fcf*mesh.Sf()/mesh.magSf());
	//Step 9: solve for capillary pressure field:
	fvScalarMatrix pcEqn
	(
		fvm::laplacian(pc) == fvc::div(fcf * mesh.magSf() )
	);
	pcEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

	pcEqn.solve();
}
