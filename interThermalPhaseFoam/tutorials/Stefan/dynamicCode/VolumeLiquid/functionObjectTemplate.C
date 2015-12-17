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

defineTypeNameAndDebug(VolumeLiquidFunctionObject, 0);


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const objectRegistry& VolumeLiquidFunctionObject::obr() const
{
    return obr_;
}


const fvMesh& VolumeLiquidFunctionObject::mesh() const
{
    return refCast<const fvMesh>(obr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

VolumeLiquidFunctionObject::VolumeLiquidFunctionObject
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

VolumeLiquidFunctionObject::~VolumeLiquidFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void VolumeLiquidFunctionObject::read(const dictionary& dict)
{
    if (false)
    {
        Info<<"read VolumeLiquid sha1: b086d6311dfa1c766248101f3ede485f1b974d88\n";
    }

//{{{ begin code
    
//}}} end code
}


void VolumeLiquidFunctionObject::execute()
{
    if (false)
    {
        Info<<"execute VolumeLiquid sha1: b086d6311dfa1c766248101f3ede485f1b974d88\n";
    }

//{{{ begin code
    
//}}} end code
}


void VolumeLiquidFunctionObject::end()
{
    if (false)
    {
        Info<<"end VolumeLiquid sha1: b086d6311dfa1c766248101f3ede485f1b974d88\n";
    }

//{{{ begin code
    
//}}} end code
}


void VolumeLiquidFunctionObject::timeSet()
{
    if (false)
    {
        Info<<"timeSet VolumeLiquid sha1: b086d6311dfa1c766248101f3ede485f1b974d88\n";
    }

//{{{ begin codeTime
    
//}}} end code
}


void VolumeLiquidFunctionObject::write()
{
    if (false)
    {
        Info<<"write VolumeLiquid sha1: b086d6311dfa1c766248101f3ede485f1b974d88\n";
    }

//{{{ begin code
    #line 71 ".VolumeLiquid"
//Sum up liquid volume:
			const volScalarField& alpha1 = mesh().lookupObject<volScalarField>("alpha1");
            		const scalar VolLiquid = gSum( mesh().V() * alpha1.internalField() );
                        //***********************************************************************

                        //Finally print out results:

                        //Get t
                        scalar t = mesh().time().value();           

                        //Now write out data:
                        if( Pstream::master() == true )
                        {
                                std::ofstream fs;
                                fs.open ("LiquidAccumulation.dat", std::fstream::app);
                                fs.precision(12);
                                fs << t << "\t" << VolLiquid << "\n";
                                fs.close();
                        }
//}}} end code
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //

