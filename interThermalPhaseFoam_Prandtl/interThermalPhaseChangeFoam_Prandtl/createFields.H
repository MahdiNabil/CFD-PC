
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

    
    //Initialize the temperature fields
    Info<< "Reading field T\n" << endl;
    volScalarField T1
    (
        IOobject
        (
            "T1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

	volScalarField T2
    (
        IOobject
        (
            "T2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

	volScalarField T3
    (
        IOobject
        (
            "T3",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

	volScalarField T4
    (
        IOobject
        (
            "T4",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

	volScalarField T5
    (
        IOobject
        (
            "T5",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    //Initialize the enthalpy fields
    volScalarField H1
    (
        IOobject
        (
            "H1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

	volScalarField H2
    (
        IOobject
        (
            "H2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

	volScalarField H3
    (
        IOobject
        (
            "H3",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

	volScalarField H4
    (
        IOobject
        (
            "H4",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

	volScalarField H5
    (
        IOobject
        (
            "H5",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

	//Initialize the thermal conductivities
	dimensionedScalar lambda0_("lambda0", dimensionSet(1, 1, -3, -1, 0), 0);
	dimensionedScalar lambda1_("lambda1", dimensionSet(1, 1, -3, -1, 0), 0);
	dimensionedScalar lambda2_("lambda2", dimensionSet(1, 1, -3, -1, 0), 0);
	dimensionedScalar lambda3_("lambda3", dimensionSet(1, 1, -3, -1, 0), 0);
	dimensionedScalar lambda4_("lambda4", dimensionSet(1, 1, -3, -1, 0), 0);
	dimensionedScalar lambda5_("lambda5", dimensionSet(1, 1, -3, -1, 0), 0);

    Info<< "Reading transportProperties\n" << endl;
    twoPhaseThermalMixture twoPhaseProperties(U, phi);
    //Store a couple handy constants
    const dimensionedScalar& rho1 = twoPhaseProperties.rho1();
    const dimensionedScalar& rho2 = twoPhaseProperties.rho2();
    const dimensionedScalar& cp1 = twoPhaseProperties.cp1();
    const dimensionedScalar& cp2 = twoPhaseProperties.cp2();

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


	//Initialize the phase change models:
	autoPtr<thermalPhaseChangeModel> phaseChangeModel1
	(
		thermalPhaseChangeModel::New( thermalPhaseChangeProperties.lookup("model"), thermalPhaseChangeProperties, twoPhaseProperties, T1, alpha1 )
	);

	autoPtr<thermalPhaseChangeModel> phaseChangeModel2
	(
		thermalPhaseChangeModel::New( thermalPhaseChangeProperties.lookup("model"), thermalPhaseChangeProperties, twoPhaseProperties, T2, alpha1 )
	);

	autoPtr<thermalPhaseChangeModel> phaseChangeModel3
	(
		thermalPhaseChangeModel::New( thermalPhaseChangeProperties.lookup("model"), thermalPhaseChangeProperties, twoPhaseProperties, T3, alpha1 )
	);

	autoPtr<thermalPhaseChangeModel> phaseChangeModel4
	(
		thermalPhaseChangeModel::New( thermalPhaseChangeProperties.lookup("model"), thermalPhaseChangeProperties, twoPhaseProperties, T4, alpha1 )
	);
	
	autoPtr<thermalPhaseChangeModel> phaseChangeModel5
	(
		thermalPhaseChangeModel::New( thermalPhaseChangeProperties.lookup("model"), thermalPhaseChangeProperties, twoPhaseProperties, T5, alpha1 )
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
	surfaceScalarField kEff1
	(
        IOobject
        (
            "kEff1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
		fvc::interpolate( alpha1*lambda1_ + (1-alpha1)*lambda0_ )  
		+ fvc::interpolate( rho*twoPhaseProperties.cp()*turbulence->nut() )
	);

	surfaceScalarField kEff2
	(
        IOobject
        (
            "kEff2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
		fvc::interpolate( alpha1*lambda2_ + (1-alpha1)*lambda0_ ) 
	    + fvc::interpolate( rho*twoPhaseProperties.cp()*turbulence->nut() )
	);

	surfaceScalarField kEff3
	(
        IOobject
        (
            "kEff3",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
		fvc::interpolate( alpha1*lambda3_ + (1-alpha1)*lambda0_ ) 
	    + fvc::interpolate( rho*twoPhaseProperties.cp()*turbulence->nut() )
	);

	surfaceScalarField kEff4
	(
        IOobject
        (
            "kEff4",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
		fvc::interpolate( alpha1*lambda4_ + (1-alpha1)*lambda0_ ) 
		+ fvc::interpolate( rho*twoPhaseProperties.cp()*turbulence->nut() )
	);

	surfaceScalarField kEff5
	(
        IOobject
        (
            "kEff5",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
		fvc::interpolate( alpha1*lambda5_ + (1-alpha1)*lambda0_ ) 
		+ fvc::interpolate( rho*twoPhaseProperties.cp()*turbulence->nut() )
	);


