/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    RiddersRoot

Description
    Ridder's method of root finiding given a function, bracketed root
    and accuracy.  Based on Numerical Recipes in C++, Section 9.2,
    page 362.

    Function is provided as a template parameter function object, evaluated
    using operator()(const scalar x)

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    RiddersRoot.C

\*---------------------------------------------------------------------------*/

#ifndef RiddersRoot_H
#define RiddersRoot_H

#include "scalar.H"
#include "error.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//Dummy base class
class FuncBase
{
public:
    // operator() on a scalar
    virtual scalar operator()(scalar x) const
    {   return 0;    }
    //virtual destructor for reasons I don't understand
    virtual ~FuncBase(){}
};

/*---------------------------------------------------------------------------*\
                           Class RiddersRoot Declaration
\*---------------------------------------------------------------------------*/

class RiddersRoot
{
private:

    // Private data

        //- Reference to a function
        const FuncBase& f_;

        //- Required accuracy
        const scalar eps_;


    // Private member functions

        //- Return Numerical Recipes style sign check
        inline static scalar checkSign(const scalar& a, const scalar& b)
        {
            if (b >= 0)
            {
                if (a >= 0)
                {
                    return false;
                }
                else
                {
                    return true;
                }
            }
            else
            {
                if (a >= 0)
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
        }


public:

    // Static data members

        //- Maximum number of iterations
        static const label maxIter;


    // Constructors

        //- Construct given a function
        RiddersRoot(const FuncBase& f, const scalar eps);


    // Destructor - default


    // Member Functions

        //- Return root
        scalar root(const scalar x0, const scalar x1) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
