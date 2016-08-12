/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015, Alex Rattner and Mahdi Nabil
     \\/     M anipulation  | Multiscale Thermal Fluids and Energy (MTFE)
                            | Laboratory, PSU 
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

Global
    getCellDims

Description
    Finds the characteristic size (minimum edge length) for mesh cells.

\*---------------------------------------------------------------------------*/

volScalarField CellDims
(
    IOobject
    (
        "CellDims",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
     ),
     mesh,
     dimensionedScalar("dummy", dimensionSet(0,1,0,0,0,0,0), 0)
);


// Like above, but longest cell dimensions - useful for checking interface
// location
volScalarField MaxCellDims
(
    IOobject
    (
        "MaxCellDims",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
     ),
     mesh,
     dimensionedScalar("dummy", dimensionSet(0,1,0,0,0,0,0), 0)
);

const cellList& cells = mesh.cells();

forAll(cells, c) // over all cells
{
    scalar Dim = GREAT;
    scalar MaxDim = SMALL;

    const labelList& curEdges = mesh.cellEdges()[c];

    // Now go through cell edges and get min
    forAll (curEdges, e)
    {
        scalar len = mesh.edges()[ curEdges[e] ].mag( mesh.points() );
        if (len < Dim)
        {  
            Dim = len;
        }
        if (len > MaxDim)
        {
            MaxDim = Dim;
        }
    }
    CellDims[c] = Dim;
    MaxCellDims[c] = MaxDim;
}

Info<< "Minimum cell dimension: " <<  gMin(CellDims) << endl;
