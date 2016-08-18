/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 Alex Rattner and Mahdi Nabil
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "EmpiricalRateParameter.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalPhaseChangeModels
{
    defineTypeNameAndDebug(EmpiricalRateParameter, 0);
    addToRunTimeSelectionTable
    (
        thermalPhaseChangeModel,
        EmpiricalRateParameter,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalPhaseChangeModels::EmpiricalRateParameter::EmpiricalRateParameter
(
        const word& name,
        const dictionary& thermalPhaseChangeProperties,
        const twoPhaseThermalMixture& twoPhaseProperties,
        const volScalarField& T,
        const volScalarField& alpha1
)
:
    thermalPhaseChangeModel
    (
        name,
        thermalPhaseChangeProperties,
        twoPhaseProperties,
        T,
        alpha1
    ),
    Q_pc_
    (
        IOobject
        (
            "PhaseChangeHeat",
            T_.time().timeName(),
            T.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        T.mesh(),
        dimensionedScalar( "dummy", dimensionSet(1,-1,-3,0,0,0,0), 0 )
    )
{

    // reading rl and rv
    thermalPhaseChangeProperties_.lookup("rl") >> rl;   
    thermalPhaseChangeProperties_.lookup("rv") >> rv;   

    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::thermalPhaseChangeModels::EmpiricalRateParameter::calcQ_pc()
{
    const dimensionedScalar& rhol = twoPhaseProperties_.rho1();
    const dimensionedScalar& rhov = twoPhaseProperties_.rho2();

    Q_pc_ = 
          pos(T_ - T_sat_)*h_lv_*rl*alpha1_*rhol*((T_ - T_sat_)/T_sat_)
        + neg(T_ - T_sat_)*h_lv_*rv*(1.0 - alpha1_)*rhov*((T_ - T_sat_)/T_sat_);
}


bool Foam::thermalPhaseChangeModels::EmpiricalRateParameter::
read(const dictionary& thermalPhaseChangeProperties)
{
    thermalPhaseChangeModel::read(thermalPhaseChangeProperties);
    thermalPhaseChangeProperties_.lookup("rl") >> rl;   
    thermalPhaseChangeProperties_.lookup("rv") >> rv;
    return true;
}


// ************************************************************************* //
