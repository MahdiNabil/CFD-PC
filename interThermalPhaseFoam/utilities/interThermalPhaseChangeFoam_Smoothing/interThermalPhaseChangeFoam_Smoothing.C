/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015, Alex Rattner and Mahdi Nabil
     \\/     M anipulation  | Multiscale Thermal Fluids and Energy (MTFE) Laboratory, PSU 
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    interThermalPhaseChangeFoam

Description
    Solver for 2 incompressible, flow using a VOF (volume of fluid) phase-
    fraction based interface capturing approach.

    The momentum and other fluid properties are of the "mixture" and a single
    momentum equation is solved.

    Thermal transport and thermally driven phase change effects are implemented
    in this code.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

    For a two-fluid approach see twoPhaseEulerFoam.

    Support for the Kistler (1993) dynamic contact angle model based on earlier
    work by Edin Berberovic (2008)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IOobject.H"
#include "MULES.H"
#include "subCycle.H"
#include "interfaceProperties.H"
#include "twoPhaseThermalMixture.H"
#include "turbulenceModel.H"
#include "interpolationTable.H"
#include "pimpleControl.H"
#include "wallFvPatch.H"
#include "fvIOoptionList.H"
#include "MeshGraph.H"
#include "thermalPhaseChangeModel.H"
#include "fixedFluxPressureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pimpleControl pimple(mesh);

    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
    #include "getCellDims.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "FourierNo.H"
        #include "setDeltaT.H"
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

		//Update turbulence and two phase properties
        twoPhaseProperties.correct();

        //Update fields for Kistler model
        muEffKistler = twoPhaseProperties.mu() + rho*turbulence->nut();

	//smoothing operation
	#include "alphaEqn.H"

	Info<< "****Phase change energy: " << gSum( phaseChangeModel->Q_pc()*mesh.V() ) << " W" << endl;

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
