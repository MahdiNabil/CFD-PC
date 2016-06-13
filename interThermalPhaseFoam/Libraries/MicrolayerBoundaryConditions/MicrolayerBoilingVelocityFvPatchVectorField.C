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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MicrolayerBoilingVelocityFvPatchVectorField::
MicrolayerBoilingVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF)//,
	//rhoInlet_(0.0)
{}


Foam::MicrolayerBoilingVelocityFvPatchVectorField::
MicrolayerBoilingVelocityFvPatchVectorField
(
    const MicrolayerBoilingVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper)//,
	//rhoInlet_(ptf.rhoInlet_)
{}


Foam::MicrolayerBoilingVelocityFvPatchVectorField::
MicrolayerBoilingVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF)//,
//    rhoInlet_(dict.lookupOrDefault<scalar>("rhoInlet", -VGREAT))
{
	//Initialize the patch field for the first run:	
	updateCoeffs();

/*
    if (dict.found("volumetricFlowRate"))
    {
        volumetric_ = true;
        flowRate_ = DataEntry<scalar>::New("volumetricFlowRate", dict);
        rhoName_ = "rho";
    }
    else if (dict.found("massFlowRate"))
    {
        volumetric_ = false;
        flowRate_ = DataEntry<scalar>::New("massFlowRate", dict);
        rhoName_ = word(dict.lookupOrDefault<word>("rho", "rho"));
    }
    else
    {
        FatalIOErrorIn
        (
            "MicrolayerBoilingVelocityFvPatchVectorField::"
            "MicrolayerBoilingVelocityFvPatchVectorField"
            "(const fvPatch&, const DimensionedField<vector, volMesh>&,"
            " const dictionary&)",
            dict
        )   << "Please supply either 'volumetricFlowRate' or"
            << " 'massFlowRate' and 'rho'" << exit(FatalIOError);
    }


    // Value field require if mass based
    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        evaluate(Pstream::blocking);
    }
*/

}


Foam::MicrolayerBoilingVelocityFvPatchVectorField::
MicrolayerBoilingVelocityFvPatchVectorField
(
    const MicrolayerBoilingVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf)//,
//    rhoInlet_(ptf.rhoInlet_)
{}


Foam::MicrolayerBoilingVelocityFvPatchVectorField::
MicrolayerBoilingVelocityFvPatchVectorField
(
    const MicrolayerBoilingVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF)//,
    //rhoInlet_(ptf.rhoInlet_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::MicrolayerBoilingVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

//    const scalar t = db().time().timeOutputValue();

    // a simpler way of doing this would be nice
    //const scalar avgU = -flowRate_->value(t)/gSum(patch().magSf());

    //tmp<vectorField> n = patch().nf();

	//operator==(n*avgU);

	//Inward pointing normals into the domain
	tmp<vectorField> n = -patch().nf();


	//Find out which cell faces are liquid vs. vapor:
    const fvPatchScalarField& alpha1p = patch().lookupPatchField<volScalarField, scalar>("alpha1");

	const scalar VaporInletVelocity = 0.01;

    const vectorField VaporInletVelocityp
    (
		n*VaporInletVelocity*pos(0.01 - alpha1p)
    );


	operator==(VaporInletVelocityp);

    fixedValueFvPatchVectorField::updateCoeffs();
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
