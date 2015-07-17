/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 Alex Rattner
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
    Foam::fixedWallHeatFluxFvPatchScalarField

Group
    PhaseChangeBCs

Description
	This BC is used to specify a fixed wall heat flux (in W/m^2). It is based
	on the fixedGradient BC described below:

    This boundary condition supplies a fixed gradient condition, such that
    the patch values are calculated using:

        \f[
            x_p = x_c + \frac{\nabla(x)}{\Delta}
        \f]

    where
    \vartable
        x_p      | patch values
        x_c      | internal field values
        \nabla(x)| gradient (user-specified)
        \Delta   | inverse distance from patch face centre to cell centre
    \endvartable

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        gradient     | gradient                | yes         | 
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            fixedGradient;
        gradient        uniform 0;
    }
    \endverbatim

SourceFiles
    fixedWallHeatFluxFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedWallHeatFluxFvPatchScalarField_H
#define fixedWallHeatFluxFvPatchScalarField_H

//#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class fixedWallHeatFluxFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fixedWallHeatFluxFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private data

        //scalarField gradient_;


public:

    //- Runtime type information
    TypeName("fixedWallHeatFlux");


    // Constructors

        //- Construct from patch and internal field
        fixedWallHeatFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedWallHeatFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given fixedWallHeatFluxFvPatchScalarField
        //  onto a new patch
        fixedWallHeatFluxFvPatchScalarField
        (
            const fixedWallHeatFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedWallHeatFluxFvPatchScalarField
        (
            const fixedWallHeatFluxFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedWallHeatFluxFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedWallHeatFluxFvPatchScalarField
        (
            const fixedWallHeatFluxFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedWallHeatFluxFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Return defining fields
/*
            //- Return gradient at boundary
            virtual scalarField& gradient()
            {
                return gradient_;
            }

            virtual const scalarField& gradient() const
            {
                return gradient_;
            }

*/


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );

/*
        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<scalarField> snGrad() const
            {
                return gradient_;
            }

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::blocking
            );

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<scalarField> valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<scalarField > valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<scalarField> gradientInternalCoeffs() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<scalarField> gradientBoundaryCoeffs() const;
*/

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//#   include "fixedWallHeatFluxFvPatchScalarField.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
