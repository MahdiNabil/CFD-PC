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

Class
    Foam::surfaceTensionForceModel::SST

Description
    Sharp surface tension force model, from Raeini, A. Q., Blunt, M. J., & 
    Bijeljic, B. (2012). Modelling two-phase flow in porous media at the pore 
    scale using the volume-of-fluid method. Journal of Computational Physics, 
    231(17), 5653-5668.


SourceFiles
    SST.C

\*--------------------------------------------------------------------*/

#ifndef SST_H
#define SST_H

#include "surfaceTensionForceModel.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceTensionForceModels
{

/*--------------------------------------------------------------------*\
                             Class SST
\*--------------------------------------------------------------------*/

class SST
:
    public surfaceTensionForceModel
{

// Private Member Functions

    // This phase change model needs access to the mesh:
    const fvMesh& mesh_;

    //Some parameters
    const scalar Cpc; //Sharpening factor
    //Cut off forces below this factor
    const scalar phi_c_thresholdFactor;
    //Component of non-interface-normal alpha1 to remove
    const scalar filterFactor; 

    //- Surface tension forces on cell faces
    volScalarField pc;
    surfaceScalarField Fstffv;
    volVectorField fc;
    surfaceScalarField fcf;
    surfaceScalarField fcf_filter;
    label pcRefCell;
    scalar pcRefValue;




public:

    //- Runtime type information
    TypeName("SST");

    // Constructors

    //- construct from components
    SST
    (
        const word& name,
        const dictionary& surfaceTensionForceProperties,
        const interfaceProperties& interface,
        const volScalarField& alpha1
    );


    //- Destructor
    ~SST()
    {}


    // Member Functions

        //- Correct the surface tension model
        void correct();
        
        //- Returns surface tension force
        tmp<surfaceScalarField> Fstff() const
        {
            return Fstffv;
        }

        //- Return actual capillary pressure
        tmp<volScalarField> pcap() const
        {
            return pc*0.0;
        }

        //- Calculates and Return phi_c
        tmp<surfaceScalarField> phi_c(const surfaceScalarField& rAUf_) const;

        //- Read the transportProperties dictionary and update
        bool read(const dictionary& surfaceTensionForceProperties);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceTensionForceModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
