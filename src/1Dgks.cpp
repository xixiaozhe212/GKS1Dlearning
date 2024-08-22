// 头文件
#include <iostream> // 标准库，定义输入输出
#include <cmath>

#include <fstream>// 输出
#include <iomanip>
#include <sstream>
#include <string>

// 引入自定义头文件：声明函数和类。
#include "kfvs.hpp"



////////////////////////////////////////////////////////////////////////////////////////////////////
// 第一步，定义工具：定义之后求解中用的到的函数，变量，中间值等。namespace定义命名空间，tools是名称，后面是内容。
namespace Tools
{

// 绝对值函数。返回类型和参数类型为标量；函数体包含一个三元操作符，格式是 condition ? value_if_true : value_if_false；。
scalar  abs(const scalar num)
{
    return num > 0 ? num : -num;
}

// // sign函数。返回类型和参数类型为一维张量。用于2阶格式 van Leer limiter.
tensor1 sign(const tensor1 &array)
{
    tensor1 sign = array;
    sign[ sign >= 0.0 ] =  1.0;
    sign[ sign <  0.0 ] = -1.0;

    return sign;
}

// 最值函数。返回类型和参数类型为标量。
scalar max(const scalar a, const scalar b)
{
    return a > b ? a : b;
}

scalar min(const scalar a, const scalar b)
{
    return a < b ? a : b;
}

//最大慨然速度
scalar getMostProbableSpeed(const scalar T)
{
    return sqrt(2 * Const::R * T);
}

//λ，T to λ，用温度来计算λ
scalar T2lambda(const scalar T)
{
    using Const::R;
    return 1 / ( 2.0*R*T );
}

//温度T，λ to T，用λ来计算温度
scalar lambda2T(const scalar lambda)
{
    using Const::R;
    return 1 / ( 2.0*R*lambda );
}

//守恒变量to原始变量
//λ，g中参数，与温度、摩尔质量、玻尔兹曼常数有关
tensor1 con2prim(const tensor1 &con)
{
    const scalar N = Const::Kr + Const::Kv;
    const scalar K = N + 2.0;

    const scalar rho    = con[0];
    const scalar U      = con[1] / con[0];
    const scalar lambda = ( K+1 ) * rho
                        / ( 4*( con[2] - 0.5*rho*U*U ) );

    return {rho, U, lambda};
}

//原始变量 to 守恒变量W
tensor1 prim2con(const tensor1 &prim) {
    const scalar rho  = prim[0];
    const scalar rhoU = prim[0] * prim[1];
    const scalar rhoE = 0.5*prim[0] * ( prim[1]*prim[1] + ( 3+Const::Kr )/( 2*prim[2] ) ); 

    return { rho, rhoU, rhoE };
}


// //无返回值函数，如果硬要加就用，int，return 0；没看懂1
// void  setMaxWell(const tensor1 &prim, tensor1 &Heq, tensor1 &Beq)
// {
//     using Const::pi;
//     using Const::Kr;

//     const scalar rho = prim[0], U = prim[1], lambda = prim[2];

//     tensor1 c = Const::uSpace - U;

//     Heq = rho*sqrt( lambda/pi )*exp( -lambda*c*c );//为什么是1/2次幂
//     Beq = Heq*(Kr + 2) / ( 2.0*lambda );//不知道是什么
// }

//求粒子的速度
tensor1 getUspace(const label uSize, const scalar T)
{
    tensor1 u(uSize);//一维数组

    const scalar uMin = -5.0 * getMostProbableSpeed(T);//5*sigma原则
    const scalar uMax =  5.0 * getMostProbableSpeed(T);
    const scalar du   = ( uMax - uMin ) / uSize;

    for(label i = 0; i < uSize; ++i)
    {
        u[i] =  uMin + du*(i+0.5);
    }

    return u;
}

//不知道啥玩意2
// tensor1 getUweight(const label uSize, const scalar T)
// {
//     tensor1 uW(uSize);//一维数组

//     const scalar uMin = -5.0 * getMostProbableSpeed(T);//5*sigma原则
//     const scalar uMax =  5.0 * getMostProbableSpeed(T);
//     const scalar du   = ( uMax - uMin ) / uSize;

//     for(label i = 0; i < uSize; ++i)
//     {
//         uW[i] =  du;
//     }

//     return uW;
// }

}

////////////////////////////////////////////////////////////////////////////////////////////////////
// 第二步，计算通量中会用到的矩阵

namespace Moment
{

void cal_Moment_xi(const scalar &lambda, tensor1 &Mxi)
{
    const scalar N = Const::Kr + Const::Kv;
    const scalar K = N + 2.0;

    Mxi[0] = 1.0;                                     // <\xi^0>
    Mxi[1] = K/(2.0*lambda);                          // <\xi^2>
    Mxi[2] = (3.0*K + K*(K-1.0))/(4*lambda*lambda);   // <\xi^4>
}

void cal_Moment_u(const tensor1 &prim, tensor1 &Mu, tensor1 &Mxi)
{
    
    const scalar U = prim[1], lambda = prim[2];
    
    Mu[0] = 1.0;                                     // <u^0>
    Mu[1] = U;                                       // <u^1>

    for (std::size_t i = 2; i < (Mu).size(); ++i)//循环变量 i 从 2 开始，循环继续执行直到 i 达到或超过 MuL 数组的大小 MuL.size()。
    {
        Mu[i] = U * Mu[i - 1] + 0.5 * (i - 1) * Mu[i - 2] / lambda;
    }

    cal_Moment_xi(lambda, Mxi);
}

void cal_Moment_uPositive(const tensor1 &prim, tensor1 &MuL, tensor1 &Mxi)// x<0, left
{
    using Const::pi;

    const scalar U = prim[1], lambda = prim[2];

    MuL[0] = 0.5 * erfc( -sqrt(lambda)* U ) ;
    MuL[1] = U * MuL[0] + 0.5*exp( -lambda*U*U ) / sqrt( pi*lambda );

    for (std::size_t i = 2; i < (MuL).size(); ++i)
    {
        MuL[i] = U * MuL[i - 1] + 0.5 * (i - 1) * MuL[i - 2] / lambda;
    }

    cal_Moment_xi(lambda, Mxi);

}

void cal_Moment_uNegative(const tensor1 &prim, tensor1 &MuR, tensor1 &Mxi)// x>0, right
{
    using Const::pi;

    const scalar U = prim[1], lambda = prim[2];

    MuR[0] = 0.5 * erfc( sqrt(lambda)* U ) ;
    MuR[1] = U * MuR[0] - 0.5*exp( -lambda*U*U ) / sqrt( pi*lambda );

    for (std::size_t i = 2; i < (MuR).size(); ++i)
    {
        MuR[i] = U * MuR[i - 1] + 0.5 * (i - 1) * MuR[i - 2] / lambda;
    }
   
    cal_Moment_xi(lambda, Mxi);

}

//两种求斜率的方法
tensor1 micro_slope(const tensor1 &slope, const tensor1 &prim)
{
    const scalar N = Const::Kr + Const::Kv;
    const scalar K = N + 2.0;
    
    scalar a1, a2, a3;
    
    scalar srho  = slope[0];// partia rho/partial x
    scalar srhoU = slope[1];
    scalar srhoE = slope[2];

    scalar rho = prim[0], U = prim[1], lambda = prim[2];

    scalar rho_d = 1/rho;
    scalar A = rho_d * (srhoU - U*srho);
    scalar B = rho_d * (2.0*srhoE - ( U*U + ( 1.0+K )/( 2.0*lambda ) )*srho);

    a3 = 4.0*lambda*lambda / (1.0+K) * ( B - 2.0*( U*A ) );
    a2 = 2.0*lambda*A - U*a3;
    a1 = rho_d *srho - (U*a2) - a3* ( 0.5 *U*U + ( 1.0+K )/( 4.0*lambda ) );

    return { a1, a2, a3 };
}

// tensor1 micro_slope(const tensor1 &slope, const tensor1 &prim, const tensor1 &con)
// {
//     const scalar N = Const::Kr + Const::Kv;
//     const scalar K = N + 2.0;
    
//     scalar a1, a2, a3;
    
//     scalar srho  = slope[0];// partia rho/partial x
//     scalar srhoU = slope[1];
//     scalar srhoE = slope[2];

//     scalar rho = prim[0], U = prim[1], lambda = prim[2];
//     scalar rho_d = 1/rho;
    
//     scalar sU = rho_d * (srhoU - U * srho);
//     scalar sLambda = (K+1.0)/4.0*(1.0/pow(con[2]/rho-0.5*U*U,2))*(-rho_d*srhoE+con[2]/(rho*rho)*srho+U*sU);

//     a1 = rho_d * srho - 2.0 * lambda * U * sU + (0.5*(K+1)/lambda-U*U)*sLambda;
//     a2 = 2.0*(lambda*sU+U*sLambda);
//     a3 = -2*sLambda;

//     return { a1, a2, a3 };
// }

tensor1 Moment_half(const tensor1 &slope, label m, label n, const tensor1 &Mu, const tensor1 &Mxi)
{
    tensor1 moment(3);

    moment[0] = slope[0] * Mu[m] * Mxi[n] \
                + slope[1] * Mu[m+1] * Mxi[n] \
                + slope[2] * 0.5 * (Mu[m+2] * Mxi[n] + Mu[m] * Mxi[n+1]);

    moment[1] = slope[0] * Mu[m+1] * Mxi[n] \
                + slope[1] * Mu[m+2] * Mxi[n] \
                + slope[2] * 0.5 * (Mu[m+3] * Mxi[n] + Mu[m+1] * Mxi[n+1]);

    moment[2] = 0.5 * (slope[0] * (Mu[m+2] * Mxi[n] + Mu[m] * Mxi[n+1]) \
                + slope[1] * (Mu[m+3] * Mxi[n] + Mu[m+1] * Mxi[n+1]) \
                + slope[2] * 0.5 * (Mu[m+4] * Mxi[n] + Mu[m] * Mxi[n+2] + 2 * Mu[m+2] * Mxi[n+1]));

    return moment;
}

}


////////////////////////////////////////////////////////////////////////////////////////////////////
// 第三步，构造函数，析构函数（可省略），这块没有理解
FVM::Cell::Cell(const label i)
:
    cellID(i), W(3), prim(3), sW(3)
    // H(Const::uSize), B(Const::uSize), sH(Const::uSize), sB(Const::uSize)
{
    sW[0] = 0.0;
    sW[1] = 0.0;
    sW[2] = 0.0;
}

FVM::Cell::~Cell()
{}



FVM::Face::Face(const Cell& Lcell, const Cell& Rcell)
:
    cellL(Lcell), cellR(Rcell), faceID(Lcell.cellID),
    Wl(3), Wr(3), primL(3), primR(3),
    W_c(3), prim_c(3),
    flux(3)
{}

FVM::Face::~Face()
{}


// 插值函数,Get Wl, Wr, primL, primR
void FVM::Face::interpolate()
{   
    using namespace Const;
    using namespace Tools;

    if (isFirstOrder)
    {
        Wl = cellL.W;
        Wr = cellR.W;
    }
    else
    {
        Wl = cellL.W + 0.5*dx*cellL.sW;
        Wr = cellR.W - 0.5*dx*cellR.sW;
    }
    
    primL = con2prim(Wl);
    primR = con2prim(Wr);
}




////////////////////////////////////////////////////////////////////////////////////////////////////
//第三步，计算通量。

void FVM::Face::get_KFVS_flux(const scalar dt)
{
    using namespace Const;
    using namespace Tools;
    using namespace Moment;

    // Get Wl, Wr and their moments
    interpolate();

    const scalar rhoL = primL[0], lambdaL = primL[2];
    const scalar rhoR = primR[0], lambdaR = primR[2];

    tensor1 MuL(7), MuR(7), MxiL(3), MxiR(3);
    
    // first order part
    cal_Moment_uPositive(primL, MuL, MxiL);
    cal_Moment_uNegative(primR, MuR, MxiR);
    
    // second order part
    tensor1 a_L = micro_slope(cellL.sW, primL);
    tensor1 a_R = micro_slope(cellR.sW, primR);

    tensor1 t(2);
    tensor1 unit = {1, 0, 0};
    
    t[0] = dt;
    t[1] = - 0.5 * dt * dt;
    
    flux = t[0] * ( rhoL * Moment_half(unit, 1, 0, MuL, MxiL) + rhoR * Moment_half(unit, 1, 0, MuR, MxiR));
    flux += t[1] * ( rhoL * Moment_half(a_L, 2, 0, MuL, MxiL) + rhoR * Moment_half(a_R, 2, 0, MuR, MxiR));
    
}


void FVM::Face::get_GKS_flux(const scalar dt)
{
    using namespace Const;
    using namespace Tools;
    using namespace Moment;
    tensor1 unit = {1, 0, 0};

    // Get Wl, Wr and their moments
    interpolate();

    const scalar rhoL = primL[0], lambdaL = primL[2];
    const scalar rhoR = primR[0], lambdaR = primR[2];

    tensor1 MuL(7), MuR(7), MxiL(3), MxiR(3);
   
    cal_Moment_uPositive(primL, MuL, MxiL);
    cal_Moment_uNegative(primR, MuR, MxiR);
    
    //macroscopic variables at the interface, j+1/2, collapsed left and right
    const tensor1 W_c = rhoL * Moment_half(unit, 0, 0, MuL, MxiL) + rhoR * Moment_half(unit, 0, 0, MuR, MxiR);
    
    const tensor1 prim_c = con2prim(W_c);
    const scalar rhoC = prim_c[0];

    tensor1 MuC(6), MxiC(3);
    cal_Moment_u(prim_c, MuC, MxiC);

    // second order part
    tensor1 a_L = micro_slope(cellL.sW, primL);
    tensor1 a_R = micro_slope(cellR.sW, primR);
    
    tensor1 slope_t_L = -rhoL * Moment_half(a_L, 1, 0, MuL, MxiL);
    tensor1 slope_t_R = -rhoR * Moment_half(a_R, 1, 0, MuR, MxiR);

    tensor1 a_t_L = micro_slope(slope_t_L, primL);
    tensor1 a_t_R = micro_slope(slope_t_R, primR);

    tensor1 slope_C = rhoL * Moment_half(a_L, 0, 0, MuL, MxiL) + rhoR * Moment_half(a_R, 0, 0, MuR, MxiR);
    tensor1 a_C = micro_slope(slope_C, prim_c);

    tensor1 slope_t_C = -rhoC * Moment_half(a_C, 1, 0, MuC, MxiC);
    tensor1 a_t_C = micro_slope(slope_t_C, prim_c);

    //inviscid flow, - Get artifical viscosity
    double C = fabs(rhoL / lambdaL - rhoR / lambdaR) / fabs(rhoL / lambdaL + rhoR / lambdaR);//不加f可能求出来的是整型
    const scalar C1 = 0.05, C2 = 1;
    const scalar tau = C1 * dt;
    const scalar taunum = C1 * dt + C2 * C * dt;

    scalar eta = exp(-dt/taunum);

    tensor1 t(6);
    tensor1 Ffr(3), Feq(3);
    
    t[0] = taunum * (1 - eta);
    t[1] = taunum * taunum * (eta - 1) + taunum * tau * (eta - 1) + taunum * eta * dt;
    t[2] = taunum * tau * (eta - 1);
    t[3] = dt + taunum * eta - taunum;
    t[4] = taunum * taunum * (1 - eta) + taunum * tau * (1 - eta) - taunum * eta * dt - tau * dt;
    t[5] = 0.5 * dt * dt + taunum * tau * (1 - eta) - tau * dt;

    // Free transport part, F0
    Ffr = t[0] * (rhoL * Moment_half(unit, 1, 0, MuL, MxiL) + rhoR * Moment_half(unit, 1, 0, MuR, MxiR));
    Ffr += t[1] * (rhoL * Moment_half(a_L, 2, 0, MuL,  MxiL) + rhoR * Moment_half(a_R, 2, 0, MuR,  MxiR));
    Ffr += t[2] * (rhoL * Moment_half(a_t_L, 1, 0, MuL,  MxiL) + rhoR * Moment_half(a_t_R, 1, 0, MuR,  MxiR));
    
    //sencond part of the flux, F1
    Feq = t[3] * rhoC * Moment_half(unit, 1, 0, MuC,  MxiC);
    Feq += t[4] * rhoC * Moment_half(a_C, 2, 0, MuC,  MxiC);
    Feq += t[5] * rhoC * Moment_half(a_t_C, 1, 0, MuC,  MxiC);

    flux = Feq + Ffr;
}

// 构造数据结构，cell和face
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
        M.push_back(0.0);
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

//数据初始化函数
void FVM::Solver::initlize()
{
    using namespace Const;
    using Tools::T2lambda;
    using Tools::prim2con;
    

    // Set left side
    for ( label i = 0; i < label( numMesh/2 ); ++i )
    {
        cells[i].prim = { rhoL, Ul, T2lambda(Tl) };
        // setMaxWell(cells[i].prim, cells[i].H, cells[i].B);
    }

    // Set right side
    for ( label i = label( numMesh/2 ); i < numMesh; ++i )
    {
        cells.at(i).prim = { rhoR, Ur, T2lambda(Tr) };
        // setMaxWell(cells[i].prim, cells[i].H, cells[i].B);
    }

    // Set Ghost Cells
    Cell &cellR = cells[numMesh];
    Cell &cellL = cells[ numMesh+1 ];
    cellR.prim = { rhoR, Ur, T2lambda(Tr) };  // right
    // setMaxWell(cellR.prim, cellR.H, cellR.B);

    cellL.prim = { rhoL, Ul, T2lambda(Tl) };  // left
    // setMaxWell(cellL.prim, cellL.H, cellL.B);

    // Get Conservative variables at each cell
    for ( auto& cell : cells )
    {
        cell.W = prim2con(cell.prim);
        // cell.sH = 0.0; cell.sB = 0.0;
    }
}

void FVM::Solver::getDt()// 得到时间步长 dt
{
    using Const::CFL;
    using Const::gamma;
    using Tools::min;

    dt = 1e150;// 发散的时候

    for( auto& cell: cells )//使用 auto 自动推导 cell 的类型，并通过引用（&）访问 cells 容器中的每个元素。这样可以避免复制每个元素，提高性能。
    {
        const scalar lambda = cell.prim[2];
        const scalar U      = cell.prim[1];

        const scalar sos = sqrt( 0.5*gamma/lambda );// Get sound speed
        const scalar Umax = Tools::max( Const::uMax, U ) + sos;// Get uMax 
        const scalar dtLocal = Const::dx / Umax;// Get dtLocal

        dt = min( dt, dtLocal );
    }

    dt *= CFL;
}

void FVM::Solver::reconstruct()// 重构函数
{
    using Const::numMesh;
    using Const::dx;
    using Const::uSize;

    using Tools::sign;
    using Tools::abs;

    
    for ( label i = 0; i < numMesh+2; ++i )
    {   
        if (Const::isFirstOrder) break;//语句将退出当前的循环

        if (i == 0) {
            tensor1 s = (cells[1].W - cells[0].W)/dx;
            tensor1 r = {0,0,0};
            cells[i].sW = (sign(s)+sign(r))*(s*r)/(abs(s)+abs(r)+1e-40);
        } else if (i == numMesh+1) {
            tensor1 s = {0,0,0};
            tensor1 r = (cells[numMesh+1].W - cells[numMesh].W)/dx;
            cells[i].sW = (sign(s)+sign(r))*(s*r)/(abs(s)+abs(r)+1e-40);
        } else {
            tensor1 s = (cells[i+1].W - cells[i].W)/dx;
            tensor1 r = (cells[i].W - cells[i-1].W)/dx;
            cells[i].sW = (sign(s)+sign(r))*abs(s*r)/(abs(s)+abs(r)+1e-40);
            // cells[i].sW = min(0.5*abs(s+r), 2*abs(s), 2*abs(r));
        }
    }
}

//演化，求每个face的flux
void FVM::Solver::evolve()
{
    for(auto &face : faces)
    {   
        using Const::fluxType;

        if(Const::fluxType == "KFVS" )
        {
            face.get_KFVS_flux(dt);
        }
        else if (fluxType == "GKS" )
        {
            face.get_GKS_flux(dt);
        }
    }
}





////////////////////////////////////////////////////////////////////////////////////////////////////
//第六步，更新。
void FVM::Solver::update()
{
    using Const::numMesh;
    using Const::stopTime;
    using Const::dx;
    using Const::uSize;
    using Tools::con2prim;
    // using Tools::setMaxWell;

    if (runTime + dt > stopTime)
        dt = stopTime - runTime;

    for ( label i = 0; i < numMesh; ++i )
    {   

        
        // Update macro variables
        cells[i].W += (  faces[i].flux - faces[i+1].flux  ) / dx;
        
        cells[i].prim = con2prim( cells[i].W );

        
    }

    runTime += dt;
    ++step;

}

////////////////////////////////////////////////////////////////////////////////////////////////////
//第七步，写输出文件
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

        const scalar T = Tools::lambda2T(cells[i].prim[2]);
        const scalar P = rho * Const::R * T;

        out << std::setprecision(10)
            << cellCenter           << " "
            << cells[i].prim[0]     << " "
            << cells[i].prim[1]     << " "
            << T                    << " "
            << P                    << " "
            << cells[i].W[1]        << " "
            << cells[i].W[2]        << std::endl;
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


////////////////////////////////////////////////////////////////////////////////////////////////////
//主函数
int main(int argc, char* argv[])
{
    using namespace Const;

    FVM::Solver kfvs;//定义求解器
    
    kfvs.initlize();//初始化

    while (kfvs.runTime < Const::stopTime && kfvs.step < Const::stopStep)//为啥要搞两个判断
    {
        kfvs.getDt();// dt
        kfvs.reconstruct();// 1阶不需要重构
        kfvs.evolve();
        kfvs.update();

        if (kfvs.step % Const::writeInterval == 0)
        {
            kfvs.write( std::to_string(kfvs.step) );
        }

        kfvs.info();
    }
    
    kfvs.write("End");
    
    return 0;
}

