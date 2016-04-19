
Switch CombinedHSPressure("no");
runTime.controlDict().readIfPresent("CombinedHSPressure", CombinedHSPressure);
const scalar CHSP = CombinedHSPressure ? 1.0 : 0;

    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field alpha1\n" << endl;
    volScalarField alpha1
    (
        IOobject
        (
            "alpha1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field beta\n" << endl;
    volScalarField beta
    (
        IOobject
        (
            "beta",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhi.H"

    
    //Initialize the temperature field
    Info<< "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    //Initialize the enthalpy field
    volScalarField H
    (
        IOobject
        (
            "H",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading transportProperties\n" << endl;
    twoPhaseThermalMixture twoPhaseProperties(U, phi);

    //Store a couple handy constants
    const dimensionedScalar& rho1 = twoPhaseProperties.rho1();
    const dimensionedScalar& rho2 = twoPhaseProperties.rho2();
    const dimensionedScalar& cp1 = twoPhaseProperties.cp1();
    const dimensionedScalar& cp2 = twoPhaseProperties.cp2();
	const dimensionedScalar& nu1 = twoPhaseProperties.nu1();
    //const dimensionedScalar& nu2 = twoPhaseProperties.nu2();


    // Need to store rho for ddt(rho, U)
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        alpha1*rho1 + (scalar(1) - alpha1)*rho2,
        alpha1.boundaryField().types()
    );

    rho.oldTime();

    // Mass flux
    // Initialisation does not matter because rhoPhi is reset after the
    // alpha1 solution before it is used in the U equation.
    surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rho*phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho1*phi
    );

    // Construct incompressible turbulence model
    //RAS moved this above interface properties
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, twoPhaseProperties)
    );

    //The following two fields supports the Kistler contact angle model
    volScalarField muEffKistler
    (
        IOobject
        (
            "muEffKistler",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        twoPhaseProperties.mu() + rho*turbulence->nut()
    );

    volScalarField muEff
    (
        IOobject
        (
            "muEff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        twoPhaseProperties.mu() + rho*turbulence->nut()
     /* mixture.muf()
      + fvc::interpolate(rho*turbulence->nut()) */
    );

    volScalarField nuEff
    (
        IOobject
        (
            "nuEff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        muEff/rho 
    );

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    //Kistler contact angle model needs the sigma vol scalar field to be in the objectRegistry
    // - even though sigma is already in the interface object below
    volScalarField sigma
    (
        IOobject
        (
            "sigmaKistler",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(transportProperties.lookup("sigma"))
    );

    //Reference temperature for calculating enthalpies
    volScalarField T_0
    (
        IOobject
        (
            "T_0",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(transportProperties.lookup("T_0"))
    );

    // Construct interface from alpha1 distribution
    interfaceProperties interface(alpha1, U, twoPhaseProperties);


	//This dictionary contains important props for phase change:
    const dictionary& thermalPhaseChangeProperties
    (
        transportProperties.subDict("thermalPhaseChange")
    );

	//This dictionary contains important props for phase change:
    const dictionary& surfaceTensionForceProperties
    (
        transportProperties.subDict("surfaceTensionForce")
    );


	//Initialize the phase change model:
	autoPtr<thermalPhaseChangeModel> phaseChangeModel
	(
		thermalPhaseChangeModel::New( thermalPhaseChangeProperties.lookup("model"), thermalPhaseChangeProperties, twoPhaseProperties, T, alpha1 )
	);

	//Initialize the surface tension force model:
	autoPtr<surfaceTensionForceModel> stfModel
	(
		surfaceTensionForceModel::New( surfaceTensionForceProperties.lookup("model"), surfaceTensionForceProperties, interface, alpha1 )
	);
	

    #include "readGravitationalAcceleration.H"

    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf());

	//Assigning a value for SGS bubbles diameter
	const dimensionedScalar diam("diam", dimensionSet(0,1,0,0,0,0,0), 1E-7);


	Info<< "Reading field Ubeta\n" << endl;
	//Drift Velocity of SGS Bubbles
	volVectorField Ubeta1
    (
        IOobject
        (
            "Ubeta1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        (1.0/18.0)*g*pow(diam,2.0)*(rho1-rho2)/(nu1*rho1)
    );
	//SGS Bubbles Velovity definition
    volVectorField Ubeta
    (
        IOobject
        (
            "Ubeta",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U + Ubeta1
    );

	surfaceScalarField phibeta
	(
		IOobject
		(
		    "phibeta",
		    runTime.timeName(),
		    mesh,
		    IOobject::NO_READ,
		    IOobject::AUTO_WRITE
		),
		linearInterpolate(Ubeta) & mesh.Sf()
	);

    //Total Pressure Field:
	volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + (1.0-CHSP)*rho*gh + stfModel->pcap()   //Add in surface tension
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        p_rgh,
        mesh.solutionDict().subDict("PIMPLE"),
        pRefCell,
        pRefValue
    );


	if (p_rgh.needReference())
	{
	    p += dimensionedScalar
	    (
	        "p",
	        p.dimensions(),
	        pRefValue - getRefCellValue(p, pRefCell)
	    );
	    p_rgh = p - (1.0-CHSP)*rho*gh - stfModel->pcap();  //Add in surface tension
	}

    fv::IOoptionList fvOptions(mesh);

	//Create kEff here so it can be accessed globally
	surfaceScalarField kEff
	(
        IOobject
        (
            "kEff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
		fvc::interpolate( twoPhaseProperties.lambda() )
		+ fvc::interpolate( rho*twoPhaseProperties.cp()*turbulence->nut() )
	);

