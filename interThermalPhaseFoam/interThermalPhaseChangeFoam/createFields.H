
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
        p_rgh + (1.0-CHSP)*rho*gh + stfModel->pcap()  
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
	    p_rgh = p - (1.0-CHSP)*rho*gh - stfModel->pcap(); 
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

	//Modified version of the solver, with multi-energy steps per global time step
	//Need to accumulate PCV, Q_PC, alpha1_gn over time steps

	//First run thermalPhaseChangeModel correct to get initial values
	phaseChangeModel->correct();

	volScalarField Q_pc_accum
	(
		IOobject
		(
			"Q_pc_accum",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		phaseChangeModel->Q_pc()
	);

	volScalarField PCV_accum
	(
		IOobject
		(
			"PCV_accum",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		phaseChangeModel->PCV()
	);

	volScalarField alpha1Gen_accum
	(
		IOobject
		(
			"alpha1Gen_accum",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		phaseChangeModel->alpha1Gen()
	);


