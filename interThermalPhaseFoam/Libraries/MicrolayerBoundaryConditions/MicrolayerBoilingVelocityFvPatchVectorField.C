/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "MicrolayerBoilingVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include <stdio.h>
#include "fvcGrad.H"



//#include "setRootCase.H"
//#include "createTime.H"
//#include "createMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MicrolayerBoilingVelocityFvPatchVectorField::
MicrolayerBoilingVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    transportProperties(
        IOobject
        (
            "transportProperties",
            db().time().constant(),
            db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    OldMicrolayerThickness(p.size(), pTraits<scalar>::zero)
{
    MicrolayerThicknessInitialized = false;
}


Foam::MicrolayerBoilingVelocityFvPatchVectorField::
MicrolayerBoilingVelocityFvPatchVectorField
(
    const MicrolayerBoilingVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    transportProperties(
        IOobject
        (
            "transportProperties",
            db().time().constant(),
            db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    OldMicrolayerThickness(ptf.OldMicrolayerThickness, mapper)
{

    MicrolayerThicknessInitialized = true;
}


Foam::MicrolayerBoilingVelocityFvPatchVectorField::
MicrolayerBoilingVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF),
    transportProperties(
        IOobject
        (
            "transportProperties",
            db().time().constant(),
            db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    OldMicrolayerThickness(p.size(), pTraits<scalar>::zero)
{
    //Set initial microlayer thickness from dictionary if it exists

    if (dict.found("MicrolayerThickness"))
    {
	//const scalarField tmpMicrolayerThickness("gradient", dict, p.size());
        //OldMicrolayerThickness = tmpMicrolayerThickness;
        OldMicrolayerThickness = scalarField("MicrolayerThickness", dict, p.size());
	MicrolayerThicknessInitialized = true;
    }
    else
    {
	MicrolayerThicknessInitialized = false;
    }


	//Initialize the patch field for the first run:	
	updateCoeffs();
}


Foam::MicrolayerBoilingVelocityFvPatchVectorField::
MicrolayerBoilingVelocityFvPatchVectorField
(
    const MicrolayerBoilingVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    transportProperties(
        IOobject
        (
            "transportProperties",
            db().time().constant(),
            db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    OldMicrolayerThickness(ptf.OldMicrolayerThickness)
{
    MicrolayerThicknessInitialized = true;
}


Foam::MicrolayerBoilingVelocityFvPatchVectorField::
MicrolayerBoilingVelocityFvPatchVectorField
(
    const MicrolayerBoilingVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    transportProperties(
        IOobject
        (
            "transportProperties",
            db().time().constant(),
            db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    OldMicrolayerThickness(ptf.OldMicrolayerThickness)
{
    MicrolayerThicknessInitialized = true;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



void Foam::MicrolayerBoilingVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (db().time().timeOutputValue() != 0)
    {
	if (oldTimeStep == db().time().timeOutputValue())
	{
	return;
	}
    }
    //const scalar t = db().time().timeOutputValue();

    // a simpler way of doing this would be nice
    //const scalar avgU = -flowRate_->value(t)/gSum(patch().magSf());

	//operator==(n*avgU);

	//Cell face centers on the patch
	tmp<vectorField> Cfcs = -patch().Cf();

	//Inward pointing normals into the domain
	tmp<vectorField> n = -patch().nf();

	// Extracting transport value
	const dimensionedScalar& T_sat = transportProperties.subDict("thermalPhaseChange").lookup("T_sat");
	const dimensionedScalar& h_lv = transportProperties.subDict("thermalPhaseChange").lookup("h_lv");
	const dimensionedScalar& R_g = transportProperties.subDict("thermalPhaseChange").lookup("R_g");
	const dimensionedScalar& ThermConductivityL = transportProperties.subDict("phase1").lookup("lambda");
	const dimensionedScalar& rho1 = transportProperties.subDict("phase1").lookup("rho");
	const dimensionedScalar& rho2 = transportProperties.subDict("phase2").lookup("rho");

	// Extracting scalar fields
	const fvPatchField<scalar>& alpha1p = patch().lookupPatchField<volScalarField, scalar>("alpha1");
	const fvPatchField<scalar>& Tp = patch().lookupPatchField<volScalarField, scalar>("T");

	// Declare Microlayer Thickness. Using 'Tp' to set size
	scalarField MicrolayerThickness = Tp*0;// = (Tp*0);

	// Update Microlayer Thickness from previous timestep (if t != 0)
	if (MicrolayerThicknessInitialized)
	{
		MicrolayerThickness = OldMicrolayerThickness;
	}
	else
	{
		//scalarField MicrolayerThickness = Tp*0;
		MicrolayerThicknessInitialized = true;
	}

	// Calculate Microlayer Thickness
	forAll( alpha1p, celli)
	{	

		if ((MicrolayerThickness[celli] == 0) && (alpha1p[celli] < 0.01))
		{		
			//ASR - radial distance from bubble center [Utaka 2013]		
			MicrolayerThickness[celli] = mag(Cfcs()[celli])*4.46E-3;

			//multibubble approx radially-independent initial thickness
			//MicrolayerThickness[celli] = 4.46E-6;	
		}
		else if (alpha1p[celli] > 0.9)	//reset neg values if new liquid over patch
		{
			MicrolayerThickness[celli] = 0;
		}
		//else if ((MicrolayerThickness[celli] < 0) && (alpha1p[celli] < 0.01))  //neg value signifying microlayer already depleted. Do nothing
		//{}
		
	}

	//Calculate heat flux through microlayer
	scalarField Q_Microlayer(Tp*0);
	forAll( patch().Cf(), celli)
	{
		if (MicrolayerThickness[celli] > 0)
		{
			//Guo & El-Genk [1993]
			Q_Microlayer[celli] = ((Tp[celli] - T_sat.value())/((MicrolayerThickness[celli]/ThermConductivityL.value()) + ((1/(1.0*rho2.value()*h_lv.value()))  * pow(2*3.141593*R_g.value()*pow(T_sat.value(),3),0.5))));
		}

	}
	//Microlayer Evaporated per timestep
	scalarField MicrolayerEvap
	(
		//Guo & El-Genk [1993]
		((Q_Microlayer/(rho1.value()*h_lv.value()))*db().time().deltaTValue())
	);

	//limit Q max to total Microlayer evap
	forAll( patch().Cf(), celli)
	{
		if ((MicrolayerThickness[celli] > 0) && (MicrolayerEvap[celli] > MicrolayerThickness[celli]))
		{
			Q_Microlayer[celli] = (MicrolayerThickness[celli] * (rho1.value()*h_lv.value())) / db().time().deltaTValue();
		}
	}
//if (first_iteration)
//{

	//Update MicrolayerThickness
	forAll( patch().Cf(), celli)
	{
		if (MicrolayerThickness[celli] > 0)
		{
			MicrolayerThickness[celli] = (MicrolayerThickness[celli] - MicrolayerEvap[celli]);
			if (MicrolayerThickness[celli] <= 0)
			{
				//force negative value once MC depleted
				MicrolayerThickness[celli] = -1;
			}		
		}
	}
//}



//Calculate vectorfield of gas
    	const vectorField VaporInletVelocityp
    	(
		//n*MicrolayerThickness
		n*(Q_Microlayer/(rho2.value()*h_lv.value()))
		//n*(VaporInletVelocity*pos(0.01 - alpha1p))
    	);

//scalarField gradAlpha1p = fvc::grad(alpha1p);
/*
if (db().time().value() != 0)
{
	forAll(VaporInletVelocityp, celli)
	{
Info << "alpha1 - " << alpha1p[celli] << " | OldMicroThickness - " << OldMicrolayerThickness[celli]  << " | MicroThickness - " << MicrolayerThickness[celli] << " | LayerEvap - " << MicrolayerEvap[celli] << " | Q - " << Q_Microlayer[celli] << endl;
Info << "Vapor Velocity - " << mag(VaporInletVelocityp[celli]) << " | Time - " << db().time().timeOutputValue() << " | DTime - " << db().time().deltaTValue()<< endl;

	}
}
*/

	//Update OldMicrolayerThickness
	OldMicrolayerThickness = MicrolayerThickness;

	
	operator==(VaporInletVelocityp);

    fixedValueFvPatchVectorField::updateCoeffs();

	oldTimeStep = db().time().timeOutputValue();

Info << "************************************" << endl;
Info << "Tot ML_Thickness" << gSum(MicrolayerThickness) << endl;
Info << "************************************" << endl;

}


void Foam::MicrolayerBoilingVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
//    flowRate_->writeData(os);
//    if (!volumetric_)
//    {
//        writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
//        writeEntryIfDifferent<scalar>(os, "rhoInlet", -VGREAT, rhoInlet_);
//    }
    writeEntry("value", os);
    OldMicrolayerThickness.writeEntry("MicrolayerThickness", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       MicrolayerBoilingVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
