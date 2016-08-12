/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 Alex Rattner
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

#include "thermalPhaseChangeModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermalPhaseChangeModel, 0);
    defineRunTimeSelectionTable(thermalPhaseChangeModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalPhaseChangeModel::thermalPhaseChangeModel
(
    const word& name,
    const dictionary& thermalPhaseChangeProperties,
    const twoPhaseThermalMixture& twoPhasePropeties,
    const volScalarField& T,
    const volScalarField& alpha1
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            T.time().constant(),
            T.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    name_(name),
    thermalPhaseChangeProperties_(thermalPhaseChangeProperties),
    twoPhaseProperties_(twoPhasePropeties),
    T_(T),
    alpha1_(alpha1),
    T_sat_(thermalPhaseChangeProperties_.lookup("T_sat")),
    h_lv_(thermalPhaseChangeProperties_.lookup("h_lv")),
    sw_PCV("yes"),
    sw_alpha1Gen("yes")
{
    read(thermalPhaseChangeProperties);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::thermalPhaseChangeModel::PCV() const
{
    if (sw_PCV)
    {
        tmp<volScalarField> Q_pc = this->Q_pc();
        return tmp<volScalarField> 
        ( 
            (Q_pc/h_lv_)*(   (scalar(1.0)/twoPhaseProperties_.rho2()) 
                           - (scalar(1.0)/twoPhaseProperties_.rho1()) ) );
    }
    else
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (           
                IOobject
                (
                    "PCV",
                    T_.time().timeName(),
                    T_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                T_.mesh(),
                dimensionedScalar( "dummy", dimensionSet(0,0,-1,0,0,0,0), 0 )
            )
        );
    }
}

Foam::tmp<Foam::volScalarField> Foam::thermalPhaseChangeModel::alpha1Gen() const
{
    if (sw_alpha1Gen)
    {
        tmp<volScalarField> Q_pc = this->Q_pc();
        
        return tmp<volScalarField>
        ( -Q_pc/(twoPhaseProperties_.rho1()*h_lv_) );
    }
    else
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (           
                IOobject
                (
                    "alpha1Gen",
                    T_.time().timeName(),
                    T_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                T_.mesh(),
                dimensionedScalar( "dummy", dimensionSet(0,0,-1,0,0,0,0), 0 )
            )
        );
    }
}


bool Foam::thermalPhaseChangeModel::
read(const dictionary& thermalPhaseChangeProperties)
{
    thermalPhaseChangeProperties_ = thermalPhaseChangeProperties;
    //Update these properties when the dictionary is re-read
    thermalPhaseChangeProperties_.lookup("h_lv") >> h_lv_;
    thermalPhaseChangeProperties_.lookup("T_sat") >> T_sat_;
    //Update toggles for alpha1Gen, PCV

    //turn on volume change due to phase change
    thermalPhaseChangeProperties.readIfPresent("DilatationSource", sw_PCV);

    //turn on phase volume fraction change
    thermalPhaseChangeProperties.readIfPresent
    (
        "PhaseFractionSource",
        sw_alpha1Gen
    );

    return true;
}


// ************************************************************************* //
