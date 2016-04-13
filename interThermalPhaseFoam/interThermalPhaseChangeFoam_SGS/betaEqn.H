// Alpha1 on faces
const surfaceScalarField alpha1f = min(max(fvc::interpolate(alpha1), scalar(0)), scalar(1));

// Velocity of bubbles on the faces
const surfaceScalarField Ubetaf = mag(fvc::interpolate(Ubeta));

// Lookup for rate parameters
const dimensionedScalar r = transportProperties.lookup("r");
const dimensionedScalar c1 = transportProperties.lookup("c1");
const dimensionedScalar beta0 = transportProperties.lookup("beta0");

// Cell size definition
const surfaceScalarField Delta =  Foam::sqrt(mesh.magSf());

fvScalarMatrix betaEqn
(
        fvm::ddt(alpha1, beta)
      + fvm::div(Ubetaf, beta)
);

solve
(
    		betaEqn
 		==
			fvc::average(alpha1f)*fvm::laplacian(turbulence->nuEff(), beta)
		  //fvc::average(alpha1f)*fvm::laplacian(nuEff, beta)
		  - fvm::Sp((1.0 - alpha1)*(r/fvc::average(Delta)), beta)		  
		//- beta * (1.0 - alpha1) * (r/fvc::average(Delta))
		  - max((beta - beta0), 0.0) * (c1/fvc::average(Delta))
);


Info<< "SGS beta volume fraction = "
    << beta.weightedAverage(mesh.Vsc()).value()
    << "  Min(" << beta.name() << ") = " << min(beta).value()
    << "  Max(" << beta.name() << ") = " << max(beta).value()
    << endl;
