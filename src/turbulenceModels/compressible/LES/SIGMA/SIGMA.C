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

#include "SIGMA.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(SIGMA, 0);
addToRunTimeSelectionTable(LESModel, SIGMA, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

tmp<volScalarField> SIGMA::k(const tmp<volTensorField>& gradU) const
{   
    // Weired compile error. Maybe because gradU is tmp.
    
    const volTensorField& gradUClone = gradU();
    const volTensorField& gradUT = gradUClone.T();

    const volTensorField& G(gradUT & gradU);

    // Sigma vector
    volScalarField sigma1 = tr(G);
    volScalarField sigma2 = 0.5*(sqr(tr(G)))-(tr(G & G));
    volScalarField sigma3 = det(G);
    
    return
    (
        sqr(sqr(Csigma_)*delta()/Ck_)*
        (
            (
                sigma3*(sigma1-sigma2)*(sigma2-sigma3)
            )/(
                 sqr(sigma1)
               + dimensionedScalar
                 (
                    "small",
                    dimensionSet(0,0,-2,0,0),
                    SMALL
                 )
              )
        )
     );
}

 
void SIGMA::updateSubGridScaleFields(const volTensorField& gradU)
{

    muSgs_ = Ck_*rho()*delta()*sqrt(k(gradU));
    muSgs_.correctBoundaryConditions();

    alphaSgs_ = muSgs_/Prt_;
    alphaSgs_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

SIGMA::SIGMA
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
    Csigma_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Csigma",
            coeffDict_,
            1.35
        )
    )    
{
    updateSubGridScaleFields(fvc::grad(U));

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void SIGMA::correct(const tmp<volTensorField>& gradU)
{
    GenEddyVisc::correct(gradU);
    updateSubGridScaleFields(gradU());
}


bool SIGMA::read()
{
    if (GenEddyVisc::read())
    {
        Ck_.readIfPresent(coeffDict());
        Csigma_.readIfPresent(coeffDict());

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
