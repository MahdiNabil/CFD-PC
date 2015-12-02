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

defineTypeNameAndDebug(DataSummaryFunctionObject, 0);


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const objectRegistry& DataSummaryFunctionObject::obr() const
{
    return obr_;
}


const fvMesh& DataSummaryFunctionObject::mesh() const
{
    return refCast<const fvMesh>(obr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DataSummaryFunctionObject::DataSummaryFunctionObject
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

DataSummaryFunctionObject::~DataSummaryFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void DataSummaryFunctionObject::read(const dictionary& dict)
{
    if (false)
    {
        Info<<"read DataSummary sha1: 4e50ff010403fa4f569a4a9f2dfc7b9df5a60f6d\n";
    }

//{{{ begin code
    
//}}} end code
}


void DataSummaryFunctionObject::execute()
{
    if (false)
    {
        Info<<"execute DataSummary sha1: 4e50ff010403fa4f569a4a9f2dfc7b9df5a60f6d\n";
    }

//{{{ begin code
    
//}}} end code
}


void DataSummaryFunctionObject::end()
{
    if (false)
    {
        Info<<"end DataSummary sha1: 4e50ff010403fa4f569a4a9f2dfc7b9df5a60f6d\n";
    }

//{{{ begin code
    
//}}} end code
}


void DataSummaryFunctionObject::timeSet()
{
    if (false)
    {
        Info<<"timeSet DataSummary sha1: 4e50ff010403fa4f569a4a9f2dfc7b9df5a60f6d\n";
    }

//{{{ begin codeTime
    
//}}} end code
}


void DataSummaryFunctionObject::write()
{
    if (false)
    {
        Info<<"write DataSummary sha1: 4e50ff010403fa4f569a4a9f2dfc7b9df5a60f6d\n";
    }

//{{{ begin code
    #line 92 ".DataSummary"
//***********************************************************************
			//First get t and dt
			scalar t = mesh().time().value();
			scalar dt = mesh().time().deltaTValue();
			//***********************************************************************

			//***********************************************************************
			//Next, total phase change heating rate
			const volScalarField& Q_pc = mesh().lookupObject<volScalarField>("PhaseChangeHeat");
			const scalar Q_pcInt = gSum( -mesh().V() * Q_pc.internalField() );  // [W]
			//***********************************************************************


			//***********************************************************************
			//Now get void fraction
			const volScalarField& alpha1 = mesh().lookupObject<volScalarField>("alpha1");

			//Sum up to get hydrostatic DP
			const scalar V_Vapor = gSum( (1.0-alpha1.internalField()) * mesh().V() );
			const scalar V_total  = gSum( mesh().V() );
			const scalar VoidFrac = V_Vapor/V_total;
			//***********************************************************************


			//***********************************************************************
			//Now get void fraction
			const volVectorField& U = mesh().lookupObject<volVectorField>("U");

			//Sum up to get hydrostatic DP
			const vector U_Vapor = gSum( (1.0-alpha1.internalField()) * U.internalField() * mesh().V() );

			const scalar Uy_Vapor = U_Vapor[1]/V_Vapor;
			//***********************************************************************


			//***********************************************************************
			//Finally print out results:
			//Now write out data:
			if( Pstream::master() == true )
			{
				std::ofstream fs;
				fs.open ("Bubble_Condensation.dat", std::fstream::app);
				fs.precision(8);
				fs << t << "\t" << Q_pcInt << "\t" << VoidFrac << "\t" << Uy_Vapor << "\n" ;
				fs.close();
			}
//}}} end code
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //

