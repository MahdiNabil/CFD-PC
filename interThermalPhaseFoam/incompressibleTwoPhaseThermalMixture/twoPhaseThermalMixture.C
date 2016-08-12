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

    Extended by Alex Rattner to support thermal transport


\*---------------------------------------------------------------------------*/

#include "twoPhaseThermalMixture.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"
#include "wallFvPatch.H"
#include "fvMesh.H"
// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

//- Calculate and return the laminar viscosity
void Foam::twoPhaseThermalMixture::calcNu()
{
    nuModel1_->correct();
    nuModel2_->correct();

    const volScalarField limitedAlpha1
    (
        "limitedAlpha1",
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    // Average kinematic viscosity calculated from dynamic viscosity
    nu_ = mu()/(limitedAlpha1*rho1_ + (scalar(1) - limitedAlpha1)*rho2_);
}

// Calculate and return the thermal conductivity
void Foam::twoPhaseThermalMixture::calcLambda()
{
    // Apply thermal conductivity corrections    
    lambdaModel1_->correct();
    lambdaModel2_->correct();
    
    // We may need to calculate lambda here, somehow
    lambda_ = lambda();
}

//- Calculate specific heat
void Foam::twoPhaseThermalMixture::calcCp()
{
    cp_ = cp();
}

//- Calculate density
void Foam::twoPhaseThermalMixture::calcRho()
{
    rho_ = rho();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseThermalMixture::twoPhaseThermalMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const word& alpha1Name
)
:
    // transportModel(U, phi), -ASR: the old version of OpenFOAM had transport
    // model inherit from IODictionary, and had a corresponding constructor
    // with this signature. In the new version they are decoupled, and have
    // to be inherited separately
    transportModel(),
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),


    phase1Name_("phase1"),
    phase2Name_("phase2"),

    nuModel1_
    (
        viscosityModel::New
        (
            "nu1",
            subDict(phase1Name_),
            U,
            phi
        )
    ),
    nuModel2_
    (
        viscosityModel::New
        (
            "nu2",
            subDict(phase2Name_),
            U,
            phi
        )
    ),
    lambdaModel1_
    (
        conductivityModel::New
        (
            "lambda1",
            subDict(phase1Name_),
            U,
            phi
        )
    ),
    lambdaModel2_
    (
        conductivityModel::New
        (
            "lambda2",
            subDict(phase2Name_),
            U,
            phi
        )
    ),
    rho1_(nuModel1_->viscosityProperties().lookup("rho")),
    rho2_(nuModel2_->viscosityProperties().lookup("rho")),
    //First get the specific heats from the dictionary
    cp1_( subDict(phase1Name_).lookup("cp") ),
    cp2_( subDict(phase2Name_).lookup("cp") ),
    U_(U),
    phi_(phi),
    alpha1_(U_.db().lookupObject<const volScalarField> (alpha1Name)),
    nu_
    (
        IOobject
        (
            "nu",
            U_.time().timeName(),
            U_.db()
        ),
        U_.mesh(),
        dimensionedScalar("nu", dimensionSet(0, 2, -1, 0, 0), 0),
        calculatedFvPatchScalarField::typeName
    ),
    //for the thermal conductivity
    lambda_
    (
        IOobject
        (
            "lambda",
            U_.time().timeName(),
            U_.db()
        ),
        U_.mesh(),
        dimensionedScalar("lambda", dimensionSet(1, 1, -3, -1, 0), 0),
        calculatedFvPatchScalarField::typeName
    ),
    //Now define the cp field
    cp_
    (
        IOobject
        (
            "cp",
            U_.time().timeName(),
            U_.db()
        ),
        U_.mesh(),
        dimensionedScalar("cp", dimensionSet(0, 2, -2, -1, 0), 0),
        calculatedFvPatchScalarField::typeName
    ),
    //Now define the rho field
    rho_
    (
        IOobject
        (
            "rho",
            U_.time().timeName(),
            U_.db()
        ),
        U_.mesh(),
        dimensionedScalar("rho", dimensionSet(1, -3, 0, 0, 0), 0),
        calculatedFvPatchScalarField::typeName
    ),
    mesh_(alpha1_.mesh() ),
    ImprovedTransportBlending("no")

{
    //Read fluid properties, possibly switch on ImprovedTransportBlending
    read();

    calcNu();
    //Initial calculation
    calcLambda();
    calcCp();
    calcRho();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::twoPhaseThermalMixture::mu() const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    //Improved transport blending Marschall et al. (2012) model
    if (ImprovedTransportBlending)
    {
        //Proper interface blending for viscosity model:
        //Cell face normals
        const surfaceVectorField n_face = mesh_.Sf()/mesh_.magSf();
        //Small unit of 1/L to avoid divide by zero issues
        const dimensionedScalar eps_1_L("eps_1_L", dimless/dimLength, SMALL );
        //Get the phase fraction cell face gradients
        const surfaceScalarField alpha1f = fvc::interpolate( limitedAlpha1 );
        const surfaceVectorField grad_alpha1f =
            fvc::interpolate( fvc::grad(limitedAlpha1) );
        //Get the cell-face interface normal
        const surfaceVectorField n_alpha1 =
            grad_alpha1f / ( mag( grad_alpha1f ) + eps_1_L );
        //Weighting between serial and parallel blending of viscosity:
        const surfaceScalarField eta_mu = mag( n_alpha1 & n_face );
        const surfaceScalarField mu_ser =
              alpha1f*(rho1_*fvc::interpolate( nuModel1().nu() ) ) 
            + (1.0 - alpha1f)*(rho2_*fvc::interpolate( nuModel2().nu() ) );
        const surfaceScalarField mu_par =
              1.0 /( alpha1f/(rho1_*fvc::interpolate( nuModel1().nu() ) ) 
            + (1.0 - alpha1f)/(rho2_*fvc::interpolate( nuModel2().nu() ) ) );
        const surfaceScalarField mu_impf =
            (1.0 - eta_mu)*mu_ser + eta_mu*mu_par;
        
        return tmp<volScalarField>
        (
            new volScalarField
                (
                "mu",
                 fvc::average(mu_impf)
            )
        );

    }
    else
    {
        return tmp<volScalarField>
        (
            new volScalarField
                (
                "mu",
                   limitedAlpha1*rho1_*nuModel1_->nu() 
                 + (scalar(1) - limitedAlpha1)*rho2_*nuModel2_->nu()
            )
        );
    }
}

Foam::tmp<Foam::volScalarField> Foam::twoPhaseThermalMixture::lambda() const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    //Improved thermal transport blending
    if (ImprovedTransportBlending)
    {
        //Proper interface blending for viscosity model:
        //Cell face normals
        const surfaceVectorField n_face = mesh_.Sf()/mesh_.magSf();
        //Small unit of 1/L to avoid divide by zero issues
        const dimensionedScalar eps_1_L("eps_1_L", dimless/dimLength, SMALL );
        //Get the phase fraction cell face gradients
        const surfaceScalarField alpha1f = fvc::interpolate( limitedAlpha1 );
        const surfaceVectorField grad_alpha1f =
            fvc::interpolate( fvc::grad(limitedAlpha1) );
        //Get the cell-face interface normal
        const surfaceVectorField n_alpha1 =
            grad_alpha1f/( mag( grad_alpha1f ) + eps_1_L );
        //Weighting between serial and parallel blending of conductivity:
        const surfaceScalarField eta_lambda = mag( n_alpha1 & n_face );
        const surfaceScalarField lambda_ser =
              alpha1f * fvc::interpolate( lambdaModel1_->lambda() )
            + (1.0 - alpha1f)*fvc::interpolate( lambdaModel2_->lambda() );
        const surfaceScalarField lambda_par =
            1.0 /( alpha1f / fvc::interpolate( lambdaModel1_->lambda() )
                   + (1.0 - alpha1f)/fvc::interpolate(lambdaModel2_->lambda()));

        const surfaceScalarField lambda_impf = 
            eta_lambda*lambda_par + (1.0 - eta_lambda)*lambda_ser;

        return tmp<volScalarField>
        (
            new volScalarField
                (
                "lambda",
                 fvc::average(lambda_impf)
            )
        );

    }
    else
    {
        return tmp<volScalarField>
        (
            //This is a pretty naive model for mixture conductivity
            // (linear interpolation), but easy to evaluate
            new volScalarField
            (
                "lambda",
                limitedAlpha1*lambdaModel1_->lambda()
                + (scalar(1) - limitedAlpha1)*lambdaModel2_->lambda()
            )
        );
    }

}

Foam::tmp<Foam::volScalarField> Foam::twoPhaseThermalMixture::cp() const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    //Calculate average cp
    return tmp<volScalarField>
    (
        new volScalarField
        (
            "cp",
             ( cp1_*rho1_*limitedAlpha1 + cp2_*rho2_*(scalar(1) - limitedAlpha1) )
            /( rho1_*limitedAlpha1 + rho2_*(scalar(1) - limitedAlpha1) )
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseThermalMixture::rho() const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    //Calculate average density
    return tmp<volScalarField>
    (
        new volScalarField
        (
            "rho",
            limitedAlpha1*rho1_ + (scalar(1) - limitedAlpha1)*rho2_
        )
    );
}

Foam::tmp<Foam::volScalarField> Foam::twoPhaseThermalMixture::alpha() const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    //Calculate average alpha
    return tmp<volScalarField>
    (
        new volScalarField
        (
            "alpha",
            lambda()/(cp()*(   rho1_*limitedAlpha1 
                             + rho2_*(scalar(1) - limitedAlpha1) ) )
        )
    );
}

Foam::tmp<Foam::surfaceScalarField> Foam::twoPhaseThermalMixture::muf() const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    //Improved transport blending Marschall et al. (2012) model
    if (ImprovedTransportBlending)
    {
        //Proper interface blending for viscosity model:
        //Cell face normals
        const surfaceVectorField n_face = mesh_.Sf()/mesh_.magSf();
        //Small unit of 1/L to avoid divide by zero issues
        const dimensionedScalar eps_1_L("eps_1_L", dimless/dimLength, SMALL );
        //Get the phase fraction cell face gradients
        const surfaceScalarField alpha1f = fvc::interpolate( limitedAlpha1 );
        const surfaceVectorField grad_alpha1f =
             fvc::interpolate( fvc::grad(limitedAlpha1) );
        //Get the cell-face interface normal
        const surfaceVectorField n_alpha1 =
            grad_alpha1f / ( mag( grad_alpha1f ) + eps_1_L );
        //Weighting between serial and parallel blending of viscosity:
        const surfaceScalarField eta_mu = mag( n_alpha1 & n_face );
        const surfaceScalarField mu_ser =
              alpha1f*(rho1_*fvc::interpolate( nuModel1().nu() ) ) 
            + (1.0 - alpha1f)*(rho2_*fvc::interpolate( nuModel2().nu() ) );
        const surfaceScalarField mu_par =
              1.0 /( alpha1f/(rho1_*fvc::interpolate( nuModel1().nu() ) ) 
            + (1.0 - alpha1f)/(rho2_*fvc::interpolate( nuModel2().nu() ) ) );
        const surfaceScalarField mu_impf =
            (1.0 - eta_mu)*mu_ser + eta_mu*mu_par;
                
        return tmp<surfaceScalarField>
        (
            new surfaceScalarField
                (
                "mu",
                 mu_impf
            )
        );

    }
    else
    {
        const surfaceScalarField alpha1f = fvc::interpolate( limitedAlpha1 );
        return tmp<surfaceScalarField>
        (
            new surfaceScalarField
            (
                "mu",
                alpha1f*rho1_*fvc::interpolate(nuModel1_->nu())
              + (scalar(1) - alpha1f)*rho2_*fvc::interpolate(nuModel2_->nu())
            )
        );
    }
}

Foam::tmp<Foam::surfaceScalarField> Foam::twoPhaseThermalMixture::
lambdaf() const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );


    const surfaceScalarField alpha1f
    (
        fvc::interpolate(limitedAlpha1)
    );



    //Improved thermal transport blending
    if (ImprovedTransportBlending)
    {
        //Proper interface blending for viscosity model:
        //Cell face normals
        const surfaceVectorField n_face = mesh_.Sf()/mesh_.magSf();
        //Small unit of 1/L to avoid divide by zero issues
        const dimensionedScalar eps_1_L("eps_1_L", dimless/dimLength, SMALL );
        //Get the phase fraction cell face gradients
        const surfaceVectorField grad_alpha1f =
            fvc::interpolate( fvc::grad(limitedAlpha1) );
        //Get the cell-face interface normal
        const surfaceVectorField n_alpha1 =
            grad_alpha1f / ( mag( grad_alpha1f ) + eps_1_L );
        //Weighting between serial and parallel blending of conductivity:
        const surfaceScalarField eta_lambda = mag( n_alpha1 & n_face );
        const surfaceScalarField lambda_ser =
              alpha1f * fvc::interpolate( lambdaModel1_->lambda() )
            + (1.0 - alpha1f) * fvc::interpolate( lambdaModel2_->lambda() );
        const surfaceScalarField lambda_par =
              1.0/( alpha1f / fvc::interpolate( lambdaModel1_->lambda() )
            + (1.0 - alpha1f) / fvc::interpolate( lambdaModel2_->lambda() ) );
        const surfaceScalarField lambda_impf =
            eta_lambda*lambda_par + (1.0 - eta_lambda)*lambda_ser;


        return tmp<surfaceScalarField>
        (
            new surfaceScalarField
            (
                "lambda",
                lambda_impf
            )
        );

    }
    else
    {
        return tmp<surfaceScalarField>
        (
            new surfaceScalarField
            (
                "lambdaf",
                alpha1f*fvc::interpolate(lambdaModel1_->lambda())
                +   (scalar(1) - alpha1f)
                   *fvc::interpolate(lambdaModel2_->lambda())
            )
        );
    }

}

Foam::tmp<Foam::surfaceScalarField> Foam::twoPhaseThermalMixture::alphaf() const
{
    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "alphaf",
            fvc::interpolate(alpha())
        )
    );
}

//Get cp at faces
Foam::tmp<Foam::surfaceScalarField> Foam::twoPhaseThermalMixture::cpf() const
{
    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "cpf",
            fvc::interpolate(cp())
        )
    );
}


Foam::tmp<Foam::surfaceScalarField> Foam::twoPhaseThermalMixture::nuf() const
{
    const surfaceScalarField alpha1f
    (
        min(max(fvc::interpolate(alpha1_), scalar(0)), scalar(1))
    );

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "nuf",
            (
                alpha1f*rho1_*fvc::interpolate(nuModel1_->nu())
              + (scalar(1) - alpha1f)*rho2_*fvc::interpolate(nuModel2_->nu())
            )/(alpha1f*rho1_ + (scalar(1) - alpha1f)*rho2_)
        )
    );
}



bool Foam::twoPhaseThermalMixture::read()
{
    if (transportModel::read())
    {
        if  //Check if the proper data is in the dictionaries
        (
            nuModel1_().read(subDict(phase1Name_))
            && nuModel2_().read(subDict(phase2Name_))
            && lambdaModel1_().read(subDict(phase1Name_))
            && lambdaModel2_().read(subDict(phase2Name_))
        )
        {
            nuModel1_->viscosityProperties().lookup("rho") >> rho1_;
            nuModel2_->viscosityProperties().lookup("rho") >> rho2_;
            subDict(phase1Name_).lookup("cp") >> cp1_;
            subDict(phase2Name_).lookup("cp") >> cp2_;
            readIfPresent
            (
                "ImprovedTransportBlending",
                ImprovedTransportBlending
            );

            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
