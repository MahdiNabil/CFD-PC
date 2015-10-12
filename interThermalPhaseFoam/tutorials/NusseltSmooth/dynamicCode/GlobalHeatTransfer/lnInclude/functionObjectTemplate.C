/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "functionObjectTemplate.H"
#include "Time.H"
#include "fvCFD.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(GlobalHeatTransferFunctionObject, 0);


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const objectRegistry& GlobalHeatTransferFunctionObject::obr() const
{
    return obr_;
}


const fvMesh& GlobalHeatTransferFunctionObject::mesh() const
{
    return refCast<const fvMesh>(obr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

GlobalHeatTransferFunctionObject::GlobalHeatTransferFunctionObject
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool
)
:
    name_(name),
    obr_(obr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

GlobalHeatTransferFunctionObject::~GlobalHeatTransferFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void GlobalHeatTransferFunctionObject::read(const dictionary& dict)
{
    if (false)
    {
        Info<<"read GlobalHeatTransfer sha1: 079db10f6d63460a7d4a2a67a7f4a094914a0876\n";
    }

//{{{ begin code
    
//}}} end code
}


void GlobalHeatTransferFunctionObject::execute()
{
    if (false)
    {
        Info<<"execute GlobalHeatTransfer sha1: 079db10f6d63460a7d4a2a67a7f4a094914a0876\n";
    }

//{{{ begin code
    
//}}} end code
}


void GlobalHeatTransferFunctionObject::end()
{
    if (false)
    {
        Info<<"end GlobalHeatTransfer sha1: 079db10f6d63460a7d4a2a67a7f4a094914a0876\n";
    }

//{{{ begin code
    
//}}} end code
}


void GlobalHeatTransferFunctionObject::timeSet()
{
    if (false)
    {
        Info<<"timeSet GlobalHeatTransfer sha1: 079db10f6d63460a7d4a2a67a7f4a094914a0876\n";
    }

//{{{ begin codeTime
    
//}}} end code
}


void GlobalHeatTransferFunctionObject::write()
{
    if (false)
    {
        Info<<"write GlobalHeatTransfer sha1: 079db10f6d63460a7d4a2a67a7f4a094914a0876\n";
    }

//{{{ begin code
    #line 105 ".GlobalHeatTransfer"
//Some constants first:                 
                        scalar k = 0.0892;


                        //***********************************************************************
                        //Now get wall heat flux
                        //Get mesh boundary
                        const fvBoundaryMesh& bMesh = mesh().boundary();
                        //Next get wall patch:
                        label WallPatchID  = bMesh.findPatchID("LeftSide");
                        const fvPatch& WallPatch = bMesh[WallPatchID]; 
                        //Get temp gradient on the wall
                        const volScalarField& T = mesh().lookupObject<volScalarField>("T");
                        const surfaceScalarField SnGradT = fvc::snGrad(T);
                        const scalarField GradTWall = SnGradT.boundaryField()[WallPatchID];

                        //Sum up heat flux on the wall:
			const scalar Q_Wall = gSum( WallPatch.magSf() * GradTWall * k ) / gSum( WallPatch.magSf() );
                        //***********************************************************************



                        //***********************************************************************
                        //Finally print out results:

                        //Get t and dt for reference
                        scalar t = mesh().time().value();
                        scalar dt = mesh().time().deltaTValue();

                        //Now write out data:
                        if( Pstream::master() == true )
                        {
                                std::ofstream fs;
                                fs.open ("WallHeatFlux.dat", std::fstream::app);
                                fs.precision(12);
                                fs << t << "\t" << dt << "\t" << Q_Wall << "\n";
                                fs.close();
                        }
//}}} end code
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //

