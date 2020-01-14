/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "WALE.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(WALE, 0);
addToRunTimeSelectionTable(LESModel, WALE, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

tmp<volSymmTensorField> WALE::Sd(const tmp<volTensorField>& gradU) const
{
    // Calculate S_ij^d
    return dev(symm(gradU & gradU));
}

tmp<volScalarField> WALE::k(const tmp<volTensorField>& gradU) const
{
    volScalarField magSqrSd(magSqr(Sd(gradU)));
    
    return
    (
        sqr(sqr(Cw_)*delta()/Ck_)*
        (
            pow3(magSqrSd) //|Sij^d|^3
            /(
                sqr
                (
                    pow(magSqr(symm(gradU)), 5.0/2.0) //|Sij|^(5/2)
                  + pow(magSqrSd,5.0/4.0) //|Sij^d|^(5/4)
                )
              + dimensionedScalar
                (
                    "small",
                    dimensionSet(0,0,-10,0,0),
                    SMALL
                )
             )
        )
     );
}

void WALE::updateSubGridScaleFields(const volTensorField& gradU)
{
    nuSgs_ = Ck_*delta()*sqrt(k(gradU));
    nuSgs_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

WALE::WALE
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    LESModel(modelName, U, phi, transport, turbulenceModelName),
    GenEddyVisc(U, phi, transport),

    Ck_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ck",
            coeffDict_,
            0.094
        )
    ),
    Cw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw",
            coeffDict_,
            0.325
        )
    )
{
    updateSubGridScaleFields(fvc::grad(U));

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void WALE::correct(const tmp<volTensorField>& gradU)
{
    GenEddyVisc::correct(gradU);
    updateSubGridScaleFields(gradU());
}


bool WALE::read()
{
    if (GenEddyVisc::read())
    {
        Ck_.readIfPresent(coeffDict());
        Cw_.readIfPresent(coeffDict());
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
