/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::incompressible::LESModels::SIGMA

Group
    grpIcoLESTurbulence

Description
    The  SIGMA Model for incompressible flows.
    Reference:
    \verbatim
        Nicoud, Franck, et al. (2011)
        Using singular values to build a subgrid-scale model for large eddy
        simulations.
        Physics of Fluids 23.8, 085106.
    \endverbatim

    The default model coefficients are
    \verbatim
        SIGMACoeffs
        {
            Ck                  0.094;
            Ce                  1.048;
            Csigma              1.35;
        }
    \endverbatim

SourceFiles
    SIGMA.C

\*---------------------------------------------------------------------------*/

#ifndef SIGMA_H
#define SIGMA_H

#include "GenEddyVisc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class SIGMA Declaration
\*---------------------------------------------------------------------------*/

class SIGMA
:
    public GenEddyVisc
{
    // Private data

        dimensionedScalar Ck_;
        dimensionedScalar Csigma_;


    // Private Member Functions

        //- Update sub-grid scale fields
        void updateSubGridScaleFields(const volTensorField& gradU);

        // Disallow default bitwise copy construct and assignment
        SIGMA(const SIGMA&);
        SIGMA& operator=(const SIGMA&);


public:

    //- Runtime type information
    TypeName("SIGMA");

    // Constructors

        //- Construct from components
        SIGMA
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~SIGMA()
    {}


    // Member Functions
    
        //- Return SGS kinetic energy
        //  calculated from the given velocity gradient
        //volScalarField& k(const volTensorField& gradU) const;
        tmp<volScalarField> k(const volTensorField& gradU) const;

        //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k(fvc::grad(U()));
        }


        //- Correct Eddy-Viscosity and related properties
        virtual void correct(const tmp<volTensorField>& gradU);

        //- Read LESProperties dictionary
        virtual bool read();
        volScalarField& boundMinMax(volScalarField& , const dimensionedScalar& lowerBound, const dimensionedScalar& upperBound) const; 

        vector cubitEqnRoot(const scalar& a,const scalar& b, const scalar& c, const scalar& d) const;

        vector quadraticEqnRoot( const scalar& a, const scalar& b,const scalar& c) const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
