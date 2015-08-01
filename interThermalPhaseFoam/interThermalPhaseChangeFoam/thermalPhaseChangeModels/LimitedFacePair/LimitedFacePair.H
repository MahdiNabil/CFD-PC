/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 Alex Rattner
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
	Foam::thermalPhaseChangeModel::LimitedFacePair

Description
	This model forces internal interface cells to the equilibrium state,
	causing evaporation or condensation as appropriate. The model finds
    interface cell faces, and applies phase change to the cell on either
    side (i.e. evaporation in cell with lower alpha, condensation in the
    cell with higher alpha)

    This model works well for condensation, but seems to blow up for
    evaporation.


SourceFiles
    LimitedFacePair.C

\*--------------------------------------------------------------------*/

#ifndef LimitedFacePair_H
#define LimitedFacePair_H

#include "thermalPhaseChangeModel.H"
#include "MeshGraph.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalPhaseChangeModels
{

/*--------------------------------------------------------------------*\
                              Class LimitedInterfaceEquilibrium
\*--------------------------------------------------------------------*/

class LimitedFacePair
:
    public thermalPhaseChangeModel
{

// Private Member Functions

	//- Calculate and return the phase change heat rate
	void calcQ_pc();

	//- This phase change model needs access to the mesh:
	const fvMesh& mesh_;

	//- Phase change heat
	volScalarField Q_pc_;

	//- Graph of mesh grid, used for finding interface cells
	MeshGraph InterfaceMeshGraph;

	//- Field for tracking Interface cells
    volScalarField InterfaceField;

	//- Field for blocking evaporation on wall cells
	volScalarField WallField;

public:

    //- Runtime type information
    TypeName("LimitedFacePair");


    // Constructors

	//- construct from components
	LimitedFacePair
	(
		const word& name,
		const dictionary& thermalPhaseChangeProperties,
		const twoPhaseThermalMixture& twoPhaseProperties,
		const volScalarField& T,
		const volScalarField& alpha1
	);


    //- Destructor
    ~LimitedFacePair()
    {}


    // Member Functions

        //- Return phase change heat
		tmp<volScalarField> Q_pc() const
		{
			return Q_pc_;
		}

        //- Correct the thermalPhaseChange model
        void correct()
		{
			calcQ_pc();
		}
		
		//Overwrite the default PCV/alpha1Gen functions
		//- Return Volumetric generation rate of volume due to phase change
		//  i.e. how much fluid volume is generated due to evaporation
		//  or destroyed due to condensation, units (1/s)
        virtual tmp<volScalarField> PCV() const;

        //- Return volumetric production rate of the alpha1 field due to phase change
		//  This is positive in the case of condensation, or negative
		//  for evaporation, units (1/s)
        virtual tmp<volScalarField> alpha1Gen() const;


        //- Read the transportProperties dictionary and update
        bool read(const dictionary& thermalPhaseChangeProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermalPhaseChangeModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
