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

Class
    Foam::MicrolayerBoilingVelocityFvPatchVectorField

Group
    grpInletBoundaryConditions

Description


SourceFiles
    MicrolayerBoilingVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef MicrolayerBoilingVelocityFvPatchVectorField_H
#define MicrolayerBoilingVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "DataEntry.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
           Class MicrolayerBoilingVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class MicrolayerBoilingVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data
    IOdictionary transportProperties;


public:

   //- Runtime type information
   TypeName("MicrolayerBoilingVelocity");

scalarField OldMicrolayerThickness;
scalar oldTimeStep;
scalarField Q_Microlayer;

   // Constructors

        //- Construct from patch and internal field
        MicrolayerBoilingVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        MicrolayerBoilingVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  MicrolayerBoilingVelocityFvPatchVectorField
        //  onto a new patch
        MicrolayerBoilingVelocityFvPatchVectorField
        (
            const MicrolayerBoilingVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        MicrolayerBoilingVelocityFvPatchVectorField
        (
            const MicrolayerBoilingVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new MicrolayerBoilingVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        MicrolayerBoilingVelocityFvPatchVectorField
        (
            const MicrolayerBoilingVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new MicrolayerBoilingVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
