/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation and
     \\/     M anipulation  | Copyright (C) 2016 Alex Rattner and Mahdi Nabil
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
    Foam::twoPhaseThermalMixture

Description
    A two-phase incompressible transportModel with thermal properties

SourceFiles
    twoPhaseThermalMixture.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhaseThermalMixture_H
#define twoPhaseThermalMixture_H

#include "transportModel.H"
#include "viscosityModel.H"
#include "conductivityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "IOdictionary.H"
#include "fvMesh.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class twoPhaseThermalMixture Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseThermalMixture
:
    public transportModel,
    public IOdictionary
{
protected:

    // Protected data
    word phase1Name_;
    word phase2Name_;

    autoPtr<viscosityModel> nuModel1_;
    autoPtr<viscosityModel> nuModel2_;

    // Added by ASR
    autoPtr<conductivityModel> lambdaModel1_;
    autoPtr<conductivityModel> lambdaModel2_;

    dimensionedScalar rho1_;
    dimensionedScalar rho2_;

    // Added by ASR
    dimensionedScalar cp1_;
    dimensionedScalar cp2_;

    const volVectorField& U_;
    const surfaceScalarField& phi_;

    const volScalarField& alpha1_;
    volScalarField nu_;

    // Added by ASR
    volScalarField lambda_;
    volScalarField cp_;
    volScalarField rho_;
    // Added by ASR
    const fvMesh& mesh_;

    // Whether or not to use improved transport property blending models
    Switch ImprovedTransportBlending;

    // Private Member Functions

    //- Calculate and return the laminar viscosity
    void calcNu();

    //- ASR Calculate and return the laminar conductivity
    void calcLambda();

    //- ASR Calculate specific heat
    void calcCp();

    //- ASR Calculate local density
    void calcRho();

public:


    // Constructors

    //- Construct from components
    twoPhaseThermalMixture
    (
        const volVectorField& U,
        const surfaceScalarField& phi,
        const word& alpha1Name = "alpha1"
    );


    //- Destructor
    ~twoPhaseThermalMixture()
    {}


    // Member Functions

        const word phase1Name() const
        {
            return phase1Name_;
        }

        const word phase2Name() const
        {
            return phase2Name_;
        }

        //- Return const-access to phase1 viscosityModel
        const viscosityModel& nuModel1() const
        {
            return nuModel1_();
        }

        //- Return const-access to phase2 viscosityModel
        const viscosityModel& nuModel2() const
        {
            return nuModel2_();
        }

        //- RAS Return const-access to phase1 conductivityModel
        const conductivityModel& lambdaModel1() const
        {
            return lambdaModel1_();
        }

        //- RAS Return const-access to phase2 conductivityModel
        const conductivityModel& lambdaModel2() const
        {
            return lambdaModel2_();
        }

        //- Return const-access to phase1 density
        const dimensionedScalar& rho1() const
        {
            return rho1_;
        }

        //- Return const-access to phase2 density
        const dimensionedScalar& rho2() const
        {
            return rho2_;
        };

        //- Return const-access to phase1 specific heat
        const dimensionedScalar& cp1() const
        {
            return cp1_;
        }

        //- Return const-access to phase2 specific heat
        const dimensionedScalar& cp2() const
        {
            return cp2_;
        };

        //- Return the dynamic laminar viscosity
        tmp<volScalarField> mu() const;

        //- RAS Return the conductivity
        tmp<volScalarField> lambda() const;

        //- Return the face-interpolated dynamic laminar viscosity
        tmp<surfaceScalarField> muf() const;

        //- RAS Return the face-interpolated dynamic laminar conductivity
        tmp<surfaceScalarField> lambdaf() const;

        //- Return the kinematic laminar viscosity
        tmp<volScalarField> nu() const
        {
            return nu_;
        }

        tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        //- Return the face-interpolated kinematic laminar viscosity
        tmp<surfaceScalarField> nuf() const;

        //- ASR Return the specific heats
        tmp<volScalarField> cp() const;

        //- ASR Return the specific heats at faces
        tmp<surfaceScalarField> cpf() const;

        //- ASR for calculating thermal diffusivity
        tmp<volScalarField> alpha() const;

        //- ASR for calculating alpha at faces
        tmp<surfaceScalarField> alphaf() const;

        //- ASR return the density
        tmp<volScalarField> rho() const;

        //- ASR Correct the laminar viscosity and conductivity
        virtual void correct()
        {
            calcNu();
            calcLambda();
            calcCp();
            calcRho();
        }

        //- Read base transportProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
