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

#include "SIGMA.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(SIGMA, 0);
addToRunTimeSelectionTable(LESModel, SIGMA, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

tmp<volScalarField> SIGMA::k(const volTensorField& gradU) const
{

    // Step 1: Build matrix G = g^t g :
    // (G is symmetric semi-definite positive with positive eigen value)
    volTensorField G(dev(gradU.T()) & dev(gradU));
    
    
    // Step 2: Compute matrix G's invariants:
     volScalarField I1 = tr(G);
     volScalarField I2 = 0.5*(sqr(I1)-tr(G & G));
     volScalarField I3 = det(G);

    // Step 3: Compute angles from invariants:
    volScalarField A1 = sqr(I1)/9.0 - I2/3.0;
    volScalarField A2 = pow3(I1)/27.0 - I1*I2/6.0 + I3/2.0;
    volScalarField a21 = A2/(pow(A1,1.5));
    //boundMinMax(a21, -1.0, 1.0);
    volScalarField A3 = 1.0/3.0 * acos(a21);

    //bound(A1, Amin);
    Info<<"This is 2:"<<min(A1).value()<<","<<max(A1).value()<<","<<min(I1).value()<<endl;

    Info<<"This is 3:"<<min(A2).value()<<endl;

    Info<<"This is 5:"<<min(a21).value()<<"\t"<<max(a21).value()<<endl;

    //volScalarField A3 = 1.0/3.0 * acos(A2/(pow(A1,1.5)));
    Info<<"This is 4"<<min(A3).value()<< "\t"<<max(A3).value()<<endl;
    
    // Step 4: Compute Sigular values:
    const scalar& PI = constant::mathematical::pi;
    boundMinMax(A3, 0, 1.0/3.0*PI);
    //volScalarField temp = I1/3.0+2*sqrt(A1)*cos(A3);
    //Info<<"Eig1:"<<min(temp).value()<<"\t"<<max(temp).value()<<endl;
    //temp = I1/3.0-2*sqrt(A1)*cos(PI/3.0+A3);
    //Info<<"Eig1:"<<min(temp).value()<<"\t"<<max(temp).value()<<endl;
    //temp = I1/3.0-2*sqrt(A1)*cos(PI/3.0-A3);
    //Info<<"Eig1:"<<min(temp).value()<<"\t"<<max(temp).value()<<endl;
    
    //volVectorField E(this->U_);

    //Info<<"This is 51:"<<min(E.component(1)).value()<< "\t"<<max(A3).value()<<endl;
    //forAll(E, cellI)
    //{
        ////Info<<"Hah:"<<cellI<<"\t"<<E[cellI]<<"\t"<<gradU[cellI]<<"\t"<<endl;

        //E[cellI] = cubitEqnRoot(1.0, -I1[cellI], I2[cellI], -I3[cellI]);
        //vector temp = cubitEqnRoot(1.0, -I1[cellI], I2[cellI], -I3[cellI]);
        ////Info<<"Hah2:"<<cellI<<"\t"<<E[cellI]<<"\t"<<gradU[cellI]<<"\t"<<endl;
        //if (temp[0] < 0 || temp[1]<0 || temp[2]<0 ) 
        //{
            //Info<<"Hah:"<<cellI<<"\t"<<E[cellI]<<"\t"<<gradU[cellI]<<"\t"<<endl;
        //}
    //}

    //Info<<"This is 52:"<<min(E.component(1)).value()<< "\t"<<max(A3).value()<<endl;
    volScalarField lambda1 = I1/3.0 + 2*sqrt(A1)*cos(A3);
    volScalarField lambda2 = I1/3.0 - 2*sqrt(A1)*cos(PI/3.0+A3);
    volScalarField lambda3 = I1/3.0 - 2*sqrt(A1)*cos(PI/3.0-A3);
    Info<<"This is 22:"<<min(lambda1).value()<<","<<max(lambda1).value()<<endl;

    Info<<"This is 33:"<<min(lambda2).value()<<endl;

    Info<<"This is 55:"<<min(lambda3).value()<<"\t"<<max(lambda3).value()<<endl;
    volScalarField sigma1 = sqrt(lambda1);
    volScalarField sigma2 = sqrt(lambda2);
    volScalarField sigma3 = sqrt(lambda3);
    //volScalarField sigma1 = sqrt(E.component(0));
    //Info<<"This is 6:"<<min(A3).value()<< "\t"<<max(A3).value()<<endl;
    //volScalarField sigma2 = sqrt(E.component(1));
    //Info<<"This is 7:"<<min(A3).value()<< "\t"<<max(A3).value()<<endl;
    //volScalarField sigma3 = sqrt(E.component(2));
    //Info<<"This is 8:"<<min(A3).value()<< "\t"<<max(A3).value()<<endl;

    /*
    // Yet Another Implementation
    // OF-2.3.x may give wrong eigenValues sometimes. The tensor libary has been rewritten in the latter
    // version and eigenValues calculations are improved and known bugs are fixed in commit 82b3c0c(of-6)
    // Therefore in our version eigenValues function is not recommended.
    volVectorField E = eigenValues(G);
    volScalarField sigma1 = sqrt(E.component(0));
    volScalarField sigma2 = sqrt(E.component(1));
    volScalarField sigma3 = sqrt(E.component(2));
    */

    return
        (
         sqr(sqr(Csigma_)*delta()/Ck_)*
         (
          sqr(
              sigma3*(sigma1-sigma2)*(sigma2-sigma3)
             )/(
                 sqr(sqr(sigma1))
                 + dimensionedScalar
                 (
                  "small",
                  dimensionSet(0,0,-8,0,0),
                  SMALL
                 )
               )
         )
        );

}

void SIGMA::updateSubGridScaleFields(const volTensorField& gradU)
{   

    Info<<gradU.T()()[100]<<gradU[100]<<(gradU.T()()[100] & gradU[100])<<endl;

    nuSgs_ = Ck_*delta()*sqrt(k(gradU));
    nuSgs_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

SIGMA::SIGMA
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

volScalarField& SIGMA::boundMinMax(volScalarField& vsf, const dimensionedScalar& lowerBound, const dimensionedScalar& upperBound) const
{
    //bound field in the range of [lowerBound, upperBound]
    const scalar minVsf = min(vsf).value();
    const scalar maxVsf = max(vsf).value();
    if (minVsf < lowerBound.value())
    { 
        vsf.internalField() = max
        (
            max
            (
                vsf.internalField(),
                fvc::average(max(vsf, lowerBound))().internalField()
              * pos(-vsf.internalField())
            ),
            lowerBound.value()
        );

        vsf.boundaryField() = max(vsf.boundaryField(), lowerBound.value());
    }
    if (maxVsf > lowerBound.value())
    { 
        vsf.internalField() = min
        (
            min
            (
                vsf.internalField(),
                fvc::average(min(vsf, lowerBound))().internalField()
              * pos(-vsf.internalField())
            ),
            upperBound.value()
        );
        vsf.boundaryField() = min(vsf.boundaryField(), upperBound.value());
    }
    return vsf;
}

vector SIGMA::cubitEqnRoot(const scalar& a, const scalar& b,const scalar& c,const scalar& d) const
{
    /*
     Taken from OF-6
        This function solves a cubic equation of the following form:
        a*x^3 + b*x^2 + c*x + d = 0
          x^3 + B*x^2 + C*x + D = 0
    */
    const scalar p = c*a - b*b/3;
    const scalar q = b*b*b*(2.0/27.0) - b*c*a/3 + d*a*a;
    const scalar disc = p*p*p/27 + q*q/4;

    // How many roots of what types are available?
    const bool oneReal = disc == 0 && p == 0;
    const bool twoReal = disc == 0 && p != 0;
    const bool threeReal = disc < 0;
    // const bool oneRealTwoComplex = disc > 0;

    static const scalar sqrt3 = sqrt(3.0);

    scalar x;

    if (oneReal)
    {
        //const Roots<1> r = linearEqn(a, b/3).roots();
        const scalar r = -b/3 / a;
        return vector(r,r,r);
    }
    else if (twoReal)
    {
        if (q*b > 0)
        {
            x = - 2*cbrt(q/2) - b/3;
        }
        else
        {
            x = cbrt(q/2) - b/3;
            //const Roots<1> r = linearEqn(- a, x).roots();
            //return Roots<3>(Roots<2>(r, r), linearEqn(x*x, a*d).roots());
            const scalar r1 = x/a;
            const scalar r2 = -a*d/(sqr(x));
            return vector(r1,r1,r2);
        }
    }
    else if (threeReal)
    {
        const scalar wCbRe = - q/2, wCbIm = sqrt(- disc);
        const scalar wAbs = cbrt(hypot(wCbRe, wCbIm));
        const scalar wArg = atan2(wCbIm, wCbRe)/3;
        const scalar wRe = wAbs*cos(wArg), wIm = wAbs*sin(wArg);
        if (b > 0)
        {
            x = - wRe - mag(wIm)*sqrt3 - b/3;
        }
        else
        {
            x = 2*wRe - b/3;
        }
    }
    //return
    //Roots<3>
    //(
    //linearEqn(- a, x).roots(),
    //quadraticEqn(- x*x, c*x + a*d, d*x).roots()
    //);
    const scalar r1 = x/a;
    const vector r2 = quadraticEqnRoot(- x*x, c*x + a*d, d*x);

    return vector(r1, r2[0], r2[1]);
}
vector SIGMA::quadraticEqnRoot(const scalar& a, const scalar& b, const scalar& c) const
{
    /*
        This function solves a quadraticEqn equation of the following form:
        a*x^2 + b*x + c = 0
          x^2 + B*x + C = 0
    */
    
    if(a==0)
    {
        const scalar r = -c/b;
        return vector(r, r, 0);
    }

    const scalar disc = b*b/4 - a*c;

    const bool oneReal = disc == 0;
    const bool twoReal = disc > 0;   

    if (oneReal)
    {
        const scalar r = -b/2 /a;
        return vector (r, r, 0);
    }
    else if (twoReal)
    {
        const scalar x = - b/2 - sign(b)*sqrt(disc);
        //return Roots<2>(linearEqn(- a, x).roots(), linearEqn(- x, c).roots());
        return vector (x/a,  c/x, 0.0);
    }
    else // if (twoComplex)
    {
        //return Roots<2>(roots::complex, 0);
        FatalErrorIn
        (
            "SIGMA::quadraticEqnRoot"
            "("
                "scalar& a,"
                "scalar& b,"
                "scalar& c"
            ")"
        )   << " Complex roots for coefficient: "
            << a<<" x^2 + "<<b<<" x + "<<c<<"=0"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
