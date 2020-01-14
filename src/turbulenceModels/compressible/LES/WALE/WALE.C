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

\*---------------------------------------------------------------------------*/

#include "WALE.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
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

    muSgs_ = Ck_*rho()*delta()*sqrt(k(gradU));
    muSgs_.correctBoundaryConditions();

    alphaSgs_ = muSgs_/Prt_;
    alphaSgs_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

WALE::WALE
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const fluidThermo& thermoPhysicalModel,
    const word& turbulenceModelName,
    const word& modelName
)
:
    LESModel(modelName, rho, U, phi, thermoPhysicalModel, turbulenceModelName),
    GenEddyVisc(rho, U, phi, thermoPhysicalModel),

    Ck_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ck",
            coeffDict_,
            0.02
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
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
