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
	Foam::thermalPhaseChangeModel::LimitedInterfaceEquilibrium

Description
	This model forces wall and interface cells to the equilibrium state,
	causing evaporation or condensation as appropriate. The phase change
	heat term is limited to the amount of evaporating/condensing material
	in the cell.

    This model works well for condensation, but seems to blow up for
    evaporation.


SourceFiles
    LimitedInterfaceEquilibrium.C

\*--------------------------------------------------------------------*/

#ifndef LimitedInterfaceEquilibrium_H
#define LimitedInterfaceEquilibrium_H

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

class LimitedInterfaceEquilibrium
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
    TypeName("LimitedInterfaceEquilibrium");


    // Constructors

	//- construct from components
	LimitedInterfaceEquilibrium
	(
		const word& name,
		const dictionary& thermalPhaseChangeProperties,
		const twoPhaseThermalMixture& twoPhaseProperties,
		const volScalarField& T,
		const volScalarField& alpha1
	);


    //- Destructor
    ~LimitedInterfaceEquilibrium()
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
		

        //- Read the transportProperties dictionary and update
        bool read(const dictionary& thermalPhaseChangeProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermalPhaseChangeModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
