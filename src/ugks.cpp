// ----------------------------------------------------------------------------- //
//  __    __      _____     __   ___    _____                ____     ______     //
//  ) )  ( (     / ___ \   () ) / __)  / ____\              /   /    (_  __ \    // 
// ( (    ) )   / /   \_)  ( (_/ /    ( (___    ________   / /) )      ) ) \ \   //
//  ) )  ( (   ( (  ____   ()   (      \___ \  (________) /_/( (      ( (   ) )  // 
// ( (    ) )  ( ( (__  )  () /\ \         ) )                ) )      ) )  ) )  //
//  ) \__/ (    \ \__/ /   ( (  \ \    ___/ /                ( (      / /__/ /   // 
//  \______/     \____/    ()_)  \_\  /____/                 /__\    (______/    //
//                                                                               //
// ----------------------------------------------------------------------------- //

#include <iostream>
#include <cmath>

// Output part
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

#include "ugks.hpp"

// ----------------------------------------------------------------------------- //
// Description:                                                                  //
// Tools for basic kinetic theory                                                //
// ----------------------------------------------------------------------------- //

namespace Tools
{

scalar  abs(const scalar num)
{
    return num > 0 ? num : -num;
}

tensor1 sign(const tensor1 &array)
{
    tensor1 sign = array;
    sign[ sign >= 0.0 ] =  1.0;
    sign[ sign <  0.0 ] = -1.0;

    return sign;
}

scalar max(const scalar a, const scalar b)
{
    return a > b ? a : b;
}

scalar min(const scalar a, const scalar b)
{
    return a < b ? a : b;
}

scalar getGamma(const scalar Kr, const scalar Kv)
{
    return ( 5+Kr+Kv) / ( 3+Kr+Kv );
}

scalar getSoundSpeed(const scalar gamma, const scalar T)
{
    return sqrt(gamma * Const::R * T);
}

scalar getMostProbableSpeed(const scalar T)
{
    return sqrt(2 * Const::R * T);
}

scalar getTau(const scalar rho, const scalar lambda)
{
    using namespace Const;
    return muRef * pow( 2*lambda, 1-omega) / ( rho* pow( R*Tref, omega ) );
}

scalar TtoLambda(const scalar T)
{
    using Const::R;
    return 1 / ( 2.0*R*T );
}

scalar lambdaToT(const scalar lambda)
{
    using Const::R;
    return 1 / ( 2.0*R*lambda );
}

tensor1 conToPrim(const tensor1 &con)
{
    const scalar rho    = con[0];
    const scalar U      = con[1] / con[0];
    const scalar lambda = ( 3+Const::Kr ) * rho
                        / ( 4*( con[2] - 0.5*rho*U*U ) );

    return {rho, U, lambda};
}

tensor1 primToCon(const tensor1 &prim) {
    const scalar rho  = prim[0];
    const scalar rhoU = prim[0] * prim[1];
    const scalar rhoE = 0.5*prim[0] * ( prim[1]*prim[1] + ( 3+Const::Kr )/( 2*prim[2] ) ); 

    return { rho, rhoU, rhoE };
}

void  setMaxWell(const tensor1 &prim, tensor1 &Heq, tensor1 &Beq)
{
    using Const::pi;
    using Const::Kr;

    const scalar rho = prim[0], U = prim[1], lambda = prim[2];

    tensor1 c = Const::uSpace - U;

    Heq = rho*sqrt( lambda/pi )*exp( -lambda*c*c );
    Beq = Heq*(Kr + 2) / ( 2.0*lambda );
}

tensor1 getUspace(const label uSize, const scalar T)
{
    tensor1 u(uSize);

    const scalar C    = getMostProbableSpeed(T);

    const scalar uMin = -5.0 * C;
    const scalar uMax =  5.0 * C;
    const scalar du   = ( uMax - uMin ) / uSize;

    for(label i = 0; i < uSize; ++i)
    {
        u[i] =  uMin + du*(i+0.5);
    }

    return u;
}

tensor1 getUweight(const label uSize, const scalar T)
{
    const scalar C    = getMostProbableSpeed(T);

    const scalar uMin = -5.0 * C;
    const scalar uMax =  5.0 * C;
    const scalar du   = ( uMax - uMin ) / uSize;

    tensor1 uW(uSize);

    for(label i = 0; i < uSize; ++i)
    {
        uW[i] = du;
    }

    return uW;
}
}

// ----------------------------------------------------------------------------- //
// Description:                                                                  //
// Calculate moment to connect micro and macro variables                         //
// ----------------------------------------------------------------------------- //
namespace Moment
{

void calcMomentPositive(const tensor1 &prim, tensor1 &MuL, tensor1 &Mxi)
{
    using Const::pi;

    const scalar U = prim[1], lambda = prim[2];

    // Moments of normal velocity
    MuL[0] = 0.5 * erfc( -sqrt(lambda)* U ) ;
    MuL[1] = U * MuL[0] + 0.5*exp( -lambda*U*U ) / sqrt( pi*lambda );

    for (std::size_t i = 2; i < (MuL).size(); ++i)
    {
        MuL[i] = U * MuL[i - 1] + 0.5 * (i - 1) * MuL[i - 2] / lambda;
    }

    // Moments of \xi
    calcMomentXi(lambda, Mxi);
}

void calcMomentNegative(const tensor1 &prim, tensor1 &MuR, tensor1 &Mxi)
{
    using Const::pi;

    const scalar U = prim[1], lambda = prim[2];

    // Moments of normal velocity
    MuR[0] = 0.5 * erfc( sqrt(lambda)* U ) ;
    MuR[1] = U * MuR[0] - 0.5*exp( -lambda*U*U ) / sqrt( pi*lambda );

    for (std::size_t i = 2; i < (MuR).size(); ++i)
    {
        MuR[i] = U * MuR[i - 1] + 0.5 * (i - 1) * MuR[i - 2] / lambda;
    }

    calcMomentXi(lambda, Mxi);
}

void calcMomentComplete(const tensor1 &prim, tensor1 &Mu, tensor1 &Mxi)
{
    const scalar lambda = prim[2];

    // Moments of normal velocity
    calcMomentComplete_u(prim, Mu);

    // Moments of \xi
    calcMomentXi(lambda, Mxi);
}

void calcMomentComplete_u(const tensor1 &prim, tensor1 &Mu)
{
    const scalar U = prim[1], lambda = prim[2];

    // Moments of normal velocity
    Mu[0] = 1.0;
    Mu[1] = U;
    for (std::size_t i = 2; i < (Mu).size(); ++i)
    {
        Mu[i] = U * Mu[i - 1] + 0.5 * (i - 1) * Mu[i - 2] / lambda;
    }
}

void calcMomentXi(const scalar &lambda, tensor1 &Mxi)
{
    using Const::Kr;

    const scalar K = Kr + 2.0; 
    Mxi[0] = 1.0;                                           // <\xi^0>
    Mxi[1] = 0.5 * K / lambda;                              // <\xi^2>
    Mxi[2] = K * ( K+2.0 )/( 4.0*lambda*lambda );           // <\xi^4>
}

tensor1 Moment_psi(const tensor1 &Mu, const tensor1 &Mxi, label alpha, label delta)
{
    tensor1 moment_psi(3);

    moment_psi[0] = Mu[alpha]   * Mxi[delta/2];
    moment_psi[1] = Mu[alpha+1] * Mxi[delta/2];
    moment_psi[2] 
                  = 0.5 * ( Mu[alpha+2] * Mxi[delta/2] 
                          + Mu[alpha]   * Mxi[(delta+2)/2] );

    return moment_psi;
}

tensor1 Moment_slope_psi_psi
(
    const tensor1 &slope, const tensor1 &Mu, const tensor1 &Mxi, label alpha
)
{
    tensor1 a_psi_psi(3);
    a_psi_psi = slope[0] *  Moment_psi(Mu, Mxi, alpha + 0, 0) +
                slope[1] *  Moment_psi(Mu, Mxi, alpha + 1, 0) +
                slope[2] * (Moment_psi(Mu, Mxi, alpha + 2, 0) +
                            Moment_psi(Mu, Mxi, alpha + 0, 2));
    return a_psi_psi;
}

tensor1 microSlope(const tensor1 &slope, const tensor1 &prim)
{
    scalar srho  = slope[0];
    scalar srhoU = slope[1];
    scalar srhoE = slope[2];

    scalar rho = prim[0], U = prim[1], lambda = prim[2];

    scalar a1, a2, a3;

    const scalar K = Const::Kr;

    const scalar A = srhoU - U*srho;
    scalar B = 2.0*srhoE - ( U*U + ( 3.0+K )/( 2.0*lambda ) )*srho;

    a3 = 2.0*lambda*lambda / (3.0+K) * ( B - 2.0*( U*A ) );
    a2 = 2.0*( lambda*A - U*a3 );
    a1 = srho - (U*a2) - a3* ( U*U + ( 3.0+K )/( 2.0*lambda ) );

    return { a1/rho, a2/rho, a3/rho };
}

}


// ----------------------------------------------------------------------------- //
// Description:                                                                  //
// Constructor and deconstructor for cell class                                  //
// ----------------------------------------------------------------------------- //

FVM::Cell::Cell(const label i)
:
    cellID(i), W(3), prim(3), sW(3),
    H(Const::uSize), B(Const::uSize), sH(Const::uSize), sB(Const::uSize)
{}

FVM::Cell::~Cell()
{}

// ----------------------------------------------------------------------------- //
// Description:                                                                  //
// Constructor and deconstructor for face class                                  //
// ----------------------------------------------------------------------------- //

FVM::Face::Face(const Cell& Lcell, const Cell& Rcell)
:
    cellL(Lcell), cellR(Rcell), faceID(Lcell.cellID),
    Wl(3), Wr(3), primL(3), primR(3),
    W0(3), prim0(3),
    Hl(Const::uSize), Hr(Const::uSize), Bl (Const::uSize),  Br(Const::uSize),
    H0(Const::uSize), B0(Const::uSize), sH0(Const::uSize), sB0(Const::uSize),
    flux(3), fluxH(Const::uSize), fluxB(Const::uSize)
{}

FVM::Face::~Face()
{}

tensor1 FVM::Face::getTimeCoef(const scalar dt, const scalar tau)
{
    const scalar rr = dt/tau;
    const scalar eta = exp( -dt/tau );
    
    // Time coefficients for free transport, i.e. q4 and q5
    scalar tFr, tFrX;

    // Time coefficients for 2nd GKS Ffr term in time
    scalar tFrT;  

    if( rr < 1.0E-3)
    {
        scalar tmp_a = rr * ( rr * ( rr * ( rr*( rr*( 7.0-rr ) - 42.0 ) +  210.0 ) - 840.0 ) + 2520.0 ) 
                     / 5040.0;
        scalar tmp_b = rr * ( rr * ( rr * ( rr*( rr*( 48.0-7.0*rr ) - 280.0 ) + 1344.0) -5040.0 ) + 13440.0) / 40320.0;

        tFr  = dt * ( 1.0-tmp_a );
        tFrX = dt * dt * ( tmp_b-0.5 );
        tFrT = -tau*dt*( 1.0-tmp_a );
    }
    else
    {
        tFr  = tau * ( 1-eta );
        tFrX = tau*dt*eta - tau*tau*( 1-eta );
        tFrT = - tau*tau*( 1-eta );
    }

    // Time coefficients for equilibrium transport, i.e. q1, q2, q3
    const scalar tEq  = dt - tau*( 1-eta );
    const scalar tEqX = 2.0*tau*tau*( 1-eta ) - tau*dt - tau*dt*eta;
    const scalar tEqT = 0.5*dt*dt - tau*dt + tau*tau*( 1-eta );

    if(Const::isFirstOrder)
    {
        return {tEq, 0, 0, tFr, 0, 0};
    }
    else
    {
        return { tEq, tEqX, tEqT, tFr, tFrX, tFrT};
    }
    
}

void FVM::Face::interpolate()
{   
    using namespace Const;
    using namespace Tools;

    const tensor1 delta = ( sign(uSpace) + 1.0 ) * 0.5;

    // A macro by cpp code to interpolate H and B
    #define interpolateDisc(F)                                  \
    F##l = cellL.F + 0.5*dx*cellL.s##F;                         \
    F##r = cellR.F - 0.5*dx*cellR.s##F;                         \
    F##0 = F##l * delta + F##r * ( 1-delta );                   \
    s##F##0 = cellL.s##F * delta + cellR.s##F * ( 1-delta );

    // Get Hl, Hr, H0 and sH0
    interpolateDisc(H)

    // Get Bl, Br, B0 and sB0
    interpolateDisc(B)

    // Get Wl, Wr, primL, primR
    Wl = cellL.W + 0.5*dx*cellL.sW;
    Wr = cellR.W - 0.5*dx*cellR.sW;
    primL = conToPrim(Wl);
    primR = conToPrim(Wr);
}

tensor1 FVM::Face::get2ndOrderMicroFlux(const tensor1 a, const tensor1 Xeq, const tensor1 Yeq)
{
    const tensor1 u = Const::uSpace;
    // Need to take a look
    return a[0]*Xeq + a[1]*u*Xeq + a[2]*( u*u*Xeq + Yeq );
}

void FVM::Face::getFlux(const scalar dt)
{
    using namespace Const;
    using namespace Tools;
    using namespace Moment;

    interpolate();

    // Get W0
    const scalar rho0  = ( uWeight*H0 ).sum();
    const scalar rhoU0 = ( uWeight*H0*uSpace ).sum();
    const scalar rhoE0 = 0.5 * ( uWeight*H0*uSpace*uSpace ).sum() + 0.5 * ( uWeight*B0 ).sum();

    const tensor1 W0 = { rho0, rhoU0, rhoE0 };

    // Get Prim0
    const tensor1 prim0   = conToPrim(W0); 
    const scalar  lambda0 = prim0[2]; 

    // Get tau on face by W0
    scalar tau = getTau(rho0, lambda0);

    // Calculate numerical viscosity
    const scalar rhoL = primL[0], lambdaL = primL[2];
    const scalar rhoR = primR[0], lambdaR = primR[2];

    scalar muNum  = Tools::abs( rhoL/lambdaL - rhoR/lambdaR ) 
                  / Tools::abs( rhoL/lambdaL + rhoR/lambdaR );
    
    scalar tauNum = muNum * dt/CFL;

    tau += tauNum;

    // Get time coefficient
    const tensor1 tCoef = getTimeCoef(dt, tau);

    // Get micro slope aL, aR by the slope of H and B
    const tensor1 sWl = ( W0-cellL.W ) / ( 0.5*dx );
    const tensor1 sWr = ( cellR.W-W0 ) / ( 0.5*dx );

    const tensor1 aL = microSlope(sWl, prim0);
    const tensor1 aR = microSlope(sWr, prim0);

    // Get ∂W0/∂t and A0
    tensor1 MuL(7), MuR(7), MxiL(3), MxiR(3);

    calcMomentPositive(prim0, MuL, MxiL);
    calcMomentNegative(prim0, MuR, MxiR);

    const tensor1 MauL = Moment_slope_psi_psi(aL, MuL, MxiL, 1);
    const tensor1 MauR = Moment_slope_psi_psi(aR, MuR, MxiR, 1);
    const tensor1 sWt  = -rho0 * ( MauL + MauR );

    const tensor1 A0   = microSlope(sWt, prim0);

     /* * * * * * * * * * * * Macro Flux * * * * * * * * * * * * */

    tensor1 Feq(3), Ffr(3);

    tensor1 Mu(6), Mxi(3);
    calcMomentComplete(prim0, Mu, Mxi);

    //- First order Feq
    Feq = tCoef[0] * rho0 * Moment_psi(Mu, Mxi, 1, 0);

    //- Second order Feq in space
    Feq += tCoef[1] * rho0 
         * ( Moment_slope_psi_psi(aL, MuL, Mxi, 2) 
           + Moment_slope_psi_psi(aR, MuR, Mxi, 2) );
    
    //- Second order Feq in time
    Feq += tCoef[2] * rho0 * Moment_slope_psi_psi(A0, Mu, Mxi, 1);

    //- First order Ffr for mass, momentum and energy
    const tensor1 uW = uWeight, u = uSpace;
    Ffr[0] = tCoef[3] * ( uW*u*H0 ).sum();
    Ffr[1] = tCoef[3] * ( uW*u*u*H0).sum();
    Ffr[2] = tCoef[3] * 0.5 * ( ( uW*u*u*u*H0 ).sum() + ( uW*u*B0 ).sum() );

    //- Second order Ffr for mass, momentum and energy
    Ffr[0] += tCoef[4] * ( uW*u*u*sH0 ).sum();
    Ffr[1] += tCoef[4] * ( uW*u*u*u*sH0 ).sum();
    Ffr[2] += tCoef[4] * 0.5 * ( (uW*u*u*u*u*sH0).sum() + (uW*u*u*sB0).sum() );
    
    flux = Feq + Ffr;

    /* * * * * * * * * * * * Micro Flux * * * * * * * * * * * * */

    const tensor1 delta = ( sign(uSpace) + 1.0 ) * 0.5;
    Tools::setMaxWell(prim0, Heq, Beq);

    fluxH = 0.0; fluxB = 0.0;

    fluxH = tCoef[0] * u * Heq
          + tCoef[1] * u * u * get2ndOrderMicroFlux(aL, Heq, Beq) * delta
          + tCoef[1] * u * u * get2ndOrderMicroFlux(aR, Heq, Beq) * ( 1-delta )
          + tCoef[2] * u     * get2ndOrderMicroFlux(A0, Heq, Beq)
          + tCoef[3] * u * H0
          + tCoef[4] * u * u * sH0;
    
    const scalar K = Const::Kr + 2.0;
    const scalar Mxi4 = K * ( K+2.0 )/( 4.0*lambda0*lambda0 );

    fluxB = tCoef[0] * u * Beq
          + tCoef[1] * u * u * get2ndOrderMicroFlux(aL, Beq, Mxi4*Heq) * delta
          + tCoef[1] * u * u * get2ndOrderMicroFlux(aR, Beq, Mxi4*Heq) * ( 1-delta )
          + tCoef[2] * u     * get2ndOrderMicroFlux(A0, Beq, Mxi4*Heq)
          + tCoef[3] * u * B0
          + tCoef[4] * u * u * sB0;

}

void FVM::Face::getGKSFlux(const scalar dt)
{
    using namespace Const;
    using namespace Tools;
    using namespace Moment;

    // Get Wl, Wr and their moments
    interpolate();

    const scalar rhoL = primL[0], lambdaL = primL[2];
    const scalar rhoR = primR[0], lambdaR = primR[2];

    tensor1 MuL(7), MxiL(3);
    tensor1 MuR(7), MxiR(3);

    calcMomentPositive(primL, MuL, MxiL);
    calcMomentNegative(primR, MuR, MxiR);

    // Get W0 and its moments
    const tensor1 W0    = rhoL * Moment_psi(MuL, MxiL, 0, 0) 
                        + rhoR * Moment_psi(MuR, MxiR, 0, 0) ;
    const tensor1 prim0 = conToPrim(W0);
    const scalar  rho0  = prim0[0];

    tensor1 Mu(6), Mxi(3);
    calcMomentComplete(prim0, Mu, Mxi);

    // Get microslope aL, aR, Al, Ar, a0, A0
    const tensor1 aL = microSlope(cellL.sW, primL);
    const tensor1 aR = microSlope(cellR.sW, primR);

    const tensor1 sW0X = rhoL*Moment_slope_psi_psi(aL, MuL, MxiL, 0)
                       + rhoR*Moment_slope_psi_psi(aR, MuR, MxiR, 0);
    const tensor1 a0 = microSlope(sW0X, prim0);

    const tensor1 sWlT = - rhoL*Moment_slope_psi_psi(aL, MuL, MxiL, 1);
    const tensor1 sWrT = - rhoR*Moment_slope_psi_psi(aR, MuR, MxiR, 1);
    const tensor1 sW0T = - rho0*Moment_slope_psi_psi(a0, Mu , Mxi , 1);

    const tensor1 Al = microSlope(sWlT, primL);
    const tensor1 Ar = microSlope(sWrT, primR);
    const tensor1 A0 = microSlope(sW0T, prim0);


    // Get time coefficient

    //- Get physical tau
    scalar tau = getTau(prim0[0], prim0[2]);

    //- Get artifical viscosity
    scalar muNum  = Tools::abs( rhoL/lambdaL - rhoR/lambdaR ) 
                  / Tools::abs( rhoL/lambdaL + rhoR/lambdaR );
    scalar tauNum = muNum * dt / Const::CFL;

    tau += tauNum;
    
    //- Get the coefficient
    const tensor1 tCoef = getTimeCoef(dt, tau);

    // Get flux

    // Free transport part
    tensor1 Ffr(3);

    //- First order part for free transport part of GKS flux 
    Ffr = tCoef[3] * ( rhoL * Moment_psi(MuL, MxiL, 1, 0) + rhoR * Moment_psi(MuR, MxiR, 1, 0) );
    
    //- Second order in space for free transport part of GKS flux
    Ffr += tCoef[4] * ( rhoL*Moment_slope_psi_psi(aL, MuL, MxiL, 2) 
                      + rhoR*Moment_slope_psi_psi(aR, MuR, MxiR, 2) );

    //- Second order in time for free transport part of GKS flux
    Ffr += tCoef[5] * ( rhoL*Moment_slope_psi_psi(Al, MuL, MxiL, 1) 
                      + rhoR*Moment_slope_psi_psi(Ar, MuR, MxiR, 1) );

    // Collision part
    tensor1 Feq(3);

    //- First order part for collision part of GKS flux 
    Feq = tCoef[0] * rho0 * Moment_psi(Mu, Mxi, 1, 0);
    
    //- Second order in space for collision part of GKS flux
    Feq += tCoef[1] * rho0 * Moment_slope_psi_psi(a0, Mu, Mxi, 2);

    //- Second order in time for collision part of GKS flux
    Feq += tCoef[2] * rho0 * Moment_slope_psi_psi(A0, Mu, Mxi, 1);

    flux = Feq + Ffr;
}

FVM::Solver::Solver()
:
    step(0), runTime(0.0), dt(0.0)
{   
    using Const::numMesh;

    // Initialize cells: 
    // cells[numMesh] is Right ghost cell, [numMesh+1] is left ghost cell
    for (label i = 0; i < numMesh + 2; i++)
    {
        cells.emplace_back(i);
    }

    // Initialize faces
    for (label i = 0; i <= numMesh; i++)
    {
        if( i == 0 )  // cell at lift side
        {
            faces.emplace_back( cells.at( numMesh+1 ), cells.at( i ) );
        }
        else
        {
            faces.emplace_back( cells.at( i-1 ), cells.at(i) );
        }
    }
}

void FVM::Solver::initlize()
{
    using namespace Const;
    using Tools::TtoLambda;
    using Tools::primToCon;
    using Tools::setMaxWell;

    // Set left side
    for ( label i = 0; i < label( numMesh/2 ); ++i )
    {
        cells[i].prim = { rhoL, Ul, TtoLambda(Tl) };
        setMaxWell(cells[i].prim, cells[i].H, cells[i].B);
    }

    // Set right side
    for ( label i = label( numMesh/2 ); i < numMesh; ++i )
    {
        cells.at(i).prim = { rhoR, Ur, TtoLambda(Tr) };
        setMaxWell(cells[i].prim, cells[i].H, cells[i].B);
    }

    // Set Ghost Cells
    Cell &cellR = cells[numMesh];
    Cell &cellL = cells[ numMesh+1 ];
    cellR.prim = { rhoR, Ur, TtoLambda(Tr) };  // right
    setMaxWell(cellR.prim, cellR.H, cellR.B);

    cellL.prim = { rhoL, Ul, TtoLambda(Tl) };  // left
    setMaxWell(cellL.prim, cellL.H, cellL.B);

    // Get Conservative variables at each cell
    for ( auto& cell : cells )
    {
        cell.W = primToCon(cell.prim);
        cell.sH = 0.0; cell.sB = 0.0;
    }
}

void FVM::Solver::getDt()
{
    using Const::CFL;
    using Const::gamma;
    using Tools::min;

    dt = 1e150;

    for( auto& cell: cells )
    {
        const scalar lambda = cell.prim[2];
        const scalar U      = cell.prim[1];

        // Get sound speed
        const scalar sos = sqrt( 0.5*gamma/lambda );

        // Get uMax 
        const scalar Umax = Tools::max( Const::uMax, U ) + sos;

        // Get dtLocal
        const scalar dtLocal = Const::dx / Umax;

        dt = min( dt, dtLocal );
    }

    dt *= CFL;
}

void FVM::Solver::reconstruct()
{
    using Const::numMesh;
    using Const::dx;
    using Const::uSize;

    using Tools::sign;
    using Tools::abs;

    
    for ( label i = 0; i < numMesh; ++i )
    {   
        if (Const::isFirstOrder) break;

        Cell *cellL = &cells[i], *cellR = &cells[ i+1 ], *cellC = &cells[i];

        if ( i == 0)
        {
            cellL = &cells[ numMesh+1 ];
        }
        else
        {
            cellL = &cells[ i-1 ];
        }
        
        tensor1 sL(uSize), sR(uSize);

        #define reconst(F)                                          \
        sL = ( cellC->F - cellL->F ) / dx;                          \
        sR = ( cellR->F - cellC->F ) / dx;                          \
        cellC->s##F = ( sign(sL) + sign(sR) ) * abs(sR) * abs(sL)   \
                  / ( abs(sR) + abs(sL) + 1e-50 );  

        // Reconstruct micro variables
        reconst(H)
        reconst(B)

        // Reconstruct macro variables
        sL.resize(3); sR.resize(3);
        reconst(W)

                          
    }
}

void FVM::Solver::evolve()
{
    for(auto &face : faces)
    {   
        using Const::fluxType;

        if(Const::fluxType == "UGKS" )
        {
            face.getFlux(dt);
        }
        else if (fluxType == "GKS" )
        {
            face.getGKSFlux(dt);
        }
    }
}

void FVM::Solver::update()
{
    using Const::numMesh;
    using Const::dx;
    using Const::uSize;
    using Tools::conToPrim;
    using Tools::getTau;
    using Tools::setMaxWell;

    for ( label i = 0; i < numMesh; ++i )
    {   

        tensor1 HeqOld(uSize), BeqOld(uSize);
        setMaxWell(cells[i].prim, HeqOld, BeqOld);
        const scalar tauOld = getTau(cells[i].prim[0], cells[i].prim[2]);
        
        // Update macro variables
        cells[i].W += (  faces[i].flux - faces[i+1].flux  ) / dx;
        
        cells[i].prim = conToPrim( cells[i].W );

        const scalar tau = getTau(cells[i].prim[0], cells[i].prim[2]);
        tensor1 Heq(Const::uSize), Beq(Const::uSize);
        Tools::setMaxWell(cells[i].prim, Heq, Beq);

        cells[i].H = ( cells[i].H + ( faces[i].fluxH - faces[i+1].fluxH ) / dx 
                       + 0.5 * dt * ( Heq/tau + ( HeqOld-cells[i].H )/tauOld ) ) 
                   / (1 + 0.5*dt/tau);

        cells[i].B = ( cells[i].B + ( faces[i].fluxB - faces[i+1].fluxB ) / dx 
                       + 0.5*dt*( Beq/tau + ( BeqOld-cells[i].B )/tauOld ) ) 
                   / (1 + 0.5*dt/tau);
    }

    runTime += dt;
    ++step;

}

void FVM::Solver::write(std::string str)
{
    std::stringstream fileName;

    // File name like "prim1000.dat"
    fileName << "./result/prim" << str << ".dat" << std::endl;

    std::string s;

    fileName >> s;

    std::ofstream out(s.c_str());

    // Write as Tecplot format
    out << "VARIABLES = X, Density, U, Temperature, P, rhoU, rhoE" << std::endl;

    out << "ZONE I = " << Const::numMesh << ", DATAPACKING=POINT"
        << std::endl;

    for ( label i = 0; i < Const::numMesh; i++ )
    {
        const scalar cellCenter = -0.5*Const::L + ( cells[i].cellID + 0.5 ) * Const::dx;

        const scalar rho = cells[i].prim[0];

        const scalar T = Tools::lambdaToT(cells[i].prim[2]);
        const scalar P = rho * Const::R * T;

        out << std::setprecision(10)
            << cellCenter           << " "
            << cells[i].prim[0]     << " "
            << cells[i].prim[1]     << " "
            << T                    << " "
            << P                    << " "
            << cells[i].W[1]        << " "
            << cells[i].W[2]        << " "
            << std::endl;
    }
    out.close();
}

void FVM::Solver::info()
{
    using std::cout;
    using std::endl;

    cout << "step = "    << step     << "\t" 
         << "runTime = " << runTime  << "\t" 
         << "dt = "      << dt      
         << endl;
}

int main(int argc, char* argv[])
{
    using namespace Const;

    FVM::Solver ugks;
    
    ugks.initlize();

    while (ugks.runTime < Const::stopTime && ugks.step < Const::stopStep)
    {
        ugks.getDt();
        ugks.reconstruct();
        ugks.evolve();
        ugks.update();

        if (ugks.step % Const::writeInterval == 0)
        {
            ugks.write( std::to_string(ugks.step) );
        }

        ugks.info();
    }
    
    ugks.write("End");
    
    return 0;
}

