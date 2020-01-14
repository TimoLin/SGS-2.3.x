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
    Foam::compressible::LESModels::WALE

Group
    grpCmpLESTurbulence

Description
    The choric WALE Model for compressible flows.

    Algebraic eddy viscosity SGS model founded on the assumption that
    local equilibrium prevails.
    Thus,
    \verbatim

        B = 2/3*k*I - 2*nuSgs*dev(D)

    where

        D = symm(grad(U));
        k from rho*D:B + ce*rho*k^3/2/delta = 0
        muSgs = ck*rho*sqrt(k)*delta
    \endverbatim

SourceFiles
    WALE.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleWALE_H
#define compressibleWALE_H

#include "GenEddyVisc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class WALE Declaration
\*---------------------------------------------------------------------------*/

class WALE
:
    public GenEddyVisc
{
    // Private data

        dimensionedScalar Ck_;
        dimensionedScalar Cw_;


    // Private Member Functions

        //- Update sub-grid scale fields
        void updateSubGridScaleFields(const volTensorField& gradU);

        // Disallow default bitwise copy construct and assignment
        WALE(const WALE&);
        WALE& operator=(const WALE&);


public:

    //- Runtime type information
    TypeName("WALE");

    // Constructors

        //- Construct from components
        WALE
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const fluidThermo& thermoPhysicalModel,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~WALE()
    {}


    // Member Functions
        
        // Sij^d
        tmp<volSymmTensorField> Sd(const tmp<volTensorField>& gradU) const;

        //- Return SGS kinetic energy
        //  calculated from the given velocity gradient
        tmp<volScalarField> k(const tmp<volTensorField>& gradU) const;
        /*
        const
        {
            volSymmTensorField D(symm(gradU));

            volScalarField a(ce_/delta());
            volScalarField b((2.0/3.0)*tr(D));
            volScalarField c(2*ck_*delta()*(dev(D) && D));

            return sqr((-b + sqrt(sqr(b) + 4*a*c))/(2*a));
        }
        */

        //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k(fvc::grad(U()));
        }


        //- Correct Eddy-Viscosity and related properties
        virtual void correct(const tmp<volTensorField>& gradU);

        //- Read LESProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //