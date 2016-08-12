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
    Foam::surfaceTensionForceModel

Description
    Parent class for general surface tension force models in VOF solver.
    
SourceFiles
    surfaceTensionForceModel.C
    newSurfaceTensionForceModel.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceTensionForceModel_H
#define surfaceTensionForceModel_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "tmp.H"
#include "autoPtr.H"
#include "interfaceProperties.H"
#include "geometricZeroField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class surfaceTensionForceModel Declaration
\*---------------------------------------------------------------------------*/

class surfaceTensionForceModel
:
//    public transportModel,
    public IOdictionary
{

protected:

    // Protected data
        word name_;
        dictionary surfaceTensionForceProperties_;
        //Thermal mixture props
        const interfaceProperties& interface_;
        const volScalarField& alpha1_;

    // Private Member Functions

        //- Disallow copy construct
        surfaceTensionForceModel(const surfaceTensionForceModel&);

        //- Disallow default bitwise assignment
        void operator=(const surfaceTensionForceModel&);


public:

    //- Runtime type information
    TypeName("surfaceTensionForceModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            surfaceTensionForceModel,
            dictionary,
            (
                const word& name,
                const dictionary& surfaceTensionForceProperties,
                const interfaceProperties& interface,
                const volScalarField& alpha1
            ),
            (name, surfaceTensionForceProperties, interface, alpha1)
        );


    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<surfaceTensionForceModel> New
        (
            const word& name,
            const dictionary& surfaceTensionForceProperties,
            const interfaceProperties& interface,
            const volScalarField& alpha1
        );


    // Constructors

        //- Construct from components
        surfaceTensionForceModel
        (
            const word& name,
            const dictionary& surfaceTensionForceProperties,
            const interfaceProperties& interface,
            const volScalarField& alpha1
        );


    //- Destructor
    virtual ~surfaceTensionForceModel()
    {}


    // Member Functions
        
        //- Returns the surface tension force on cell faces
        virtual tmp<surfaceScalarField> Fstff() const = 0;

        //- Correct the surface tension force model
        virtual void correct() = 0;
        
        //  Returns the capillary pressure
        virtual tmp<volScalarField> pcap() const;

        //  Returns the capillary flux
        virtual tmp<surfaceScalarField> 
        phi_c(const surfaceScalarField& rAUf_) const;

        //- Read the thermalPhaseChange dictionary and update
        virtual bool read(const dictionary& surfaceTensionForceProperties) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
