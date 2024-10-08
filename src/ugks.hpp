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
//                                                                               //
// Description:                                                                  //
// Contains Const, Solver, Tools                                                 //
//                                                                               //
// ----------------------------------------------------------------------------- //

#include <valarray>
#include <vector>
#include <string>

// Typedefs

using label   = long long int;
using scalar  = double;
using tensor1 = std::valarray<scalar>;
using tensor2 = std::valarray<std::valarray<scalar>>;
using std::vector;

namespace Tools
{

// Mathematic

    scalar  abs(const scalar num);

    tensor1 sign(const tensor1 &array);

    scalar max(const scalar a, const scalar b);

    scalar min(const scalar a, const scalar b);
    
// Basic Kinetic Formula

    scalar  getGamma(const scalar Kr, const scalar Kv);

    scalar  getSoundSpeed(const scalar gamma, const scalar T);

    scalar  getMostProbableSpeed(const scalar T);

    scalar  getTau(const scalar rho, const scalar lambda);

    scalar  TtoLambda(const scalar T);

    scalar  lambdaToT(const scalar lambda);

    tensor1 conToPrim(const tensor1 &con);  

    tensor1 primToCon(const tensor1 &prim);

    void    setMaxWell(const tensor1 &prim, tensor1 &Heq, tensor1 &Beq);

// Generate velocity space

    //- Generate Discrete velocity point by mid point rule
    tensor1 getUspace(const label uSize, const scalar T);

    //- Generate Discrete velocity weight by mid point rule
    tensor1 getUweight(const label uSize, const scalar T);


}

namespace Moment
{
// Related to Moment Calculation

    void calcMomentPositive(const tensor1 &prim, tensor1 &MuL, tensor1 &Mxi);

    void calcMomentNegative(const tensor1 &prim, tensor1 &MuR, tensor1 &Mxi);

    void calcMomentComplete(const tensor1 &prim, tensor1 &Mu,  tensor1 &Mxi);

    void calcMomentComplete_u(const tensor1 &prim, tensor1 &Mu);

    void calcMomentXi(const scalar &lambda, tensor1 &Mxi);

    tensor1 Moment_psi(const tensor1 &Mu, const tensor1 &Mxi, label alpha, label delta);

    tensor1 Moment_slope_psi_psi(const tensor1 &slope, const tensor1 &Mu, const tensor1 &Mxi, label alpha);

    tensor1 microSlope(const tensor1 &slope, const tensor1 &prim);

}

namespace Const
{

// ----------------------------------------------------------------------------- //
// Description:                                                                  //
// Save initial condition, boundary condition, macro/micro propertites,          //
// and tStop, CFL                                                                // 
// ----------------------------------------------------------------------------- //

// Mesh Infomation 

    //- Mesh number
    inline constexpr label    numMesh = 100;

    //- Length, m
    inline constexpr scalar   L       = 10;

    //- Cell size
    inline constexpr scalar   dx      = L/numMesh;

// Initial Condition

    inline constexpr scalar rhoL = 6.01888398652e-04;
    inline constexpr scalar Ul   = 0.0;
    inline constexpr scalar Tl   = 696.864;   

    inline constexpr scalar rhoR = 7.52360498315e-05;
    inline constexpr scalar Ur   = 0.0;
    inline constexpr scalar Tr   = 557.4912;   
    
// Macro Properity

    //- PI
    inline constexpr scalar pi(M_PI);

    //- Boltzmann constant, m^2 kg s^-2 K^-1
    inline constexpr scalar k         = 1.380649e-23;

    //- Molecule mass, kg (N2)
    inline constexpr scalar m         = 4.65e-26;

    //- Molecule hara sphere diameter, m (N2)
    inline constexpr scalar d         = 4.17e-10;

    //- Nitrogen gas specific gas constant, m^2 s^-2 K^-1 (N2)
    inline constexpr scalar R         = 287;

    //- Mean free path length, m
    inline constexpr scalar l         = 1.0;

    //- Parameters for VHS
    inline constexpr scalar alpha     = 1.0;

    //- Parameters for VHS model
    inline constexpr scalar omega     = 0.81;

    //- Reference Temperature
    static const scalar Tref          = 348.432;
    
    //- Reference Viscosity by VSS
    static const scalar muRef       = 5.0*( alpha+1.0 )*( alpha+2.0 )*sqrt( m*k*Tref/pi ) 
                                    / ( 4.0*alpha*( 5.0-2.0*omega )*( 7.0-2.0*omega )*d*d);

    //- Rotational degree of freedom
    static const scalar Kr            = 2;

    //- Gamma 
    static const scalar gamma         = ( 5.0+Kr ) / ( 3.0+Kr );

    //- Vibrational degree of freedom
    static const bool   isVib         = false;

    //- Is first order or not
    static const bool   isFirstOrder  = false;

    //- Flux type
    static const std::string fluxType = "GKS";


// Micro Properity

    //- Discrete velocity size
    inline constexpr label uSize      = 101;

    //- Discrete velocity point 
    static const tensor1 uSpace       = Tools::getUspace(uSize, Tref);

    //- Max micro velocity speed
    static const scalar uMax          = abs( uSpace ).max();

    //- Discrete velocity size
    static const tensor1 uWeight      = Tools::getUweight(uSize, Tref);  

// Solver

    //- Stop Time
    inline constexpr scalar stopTime     = 0.00263810493558;

    //- Stop Step
    inline constexpr label  stopStep     = 1e5;

    //- CFL number
    inline constexpr scalar CFL          = 0.5;

    //- Write Interval
    inline constexpr label writeInterval = 100;

}

namespace FVM
{

// Basic Element for Finite Volume Element

class Cell
{

public:

// Member Variables

    // Geometry 

    //- Cell ID
    label cellID;

    // Macro Variables

    //- Conservative Variables
    tensor1 W;
    
    //- Primitive Variabls
    tensor1 prim;

    //- Slope of macro veriables
    tensor1 sW;

    //- Temperature and Pressure
    scalar P, T;

    // Micro Variables

    //- Reduced distribution function for conservative part
    tensor1 H;

    //- Reduced distribution function for internal degree of freedom
    tensor1 B;

    //- Slope of micro variables
    tensor1 sH, sB;

// Constructor

    //- Construct by id
    explicit Cell(const label i);

// Destructor
    ~Cell();
};

class Face
{

public: 
// Member Variables
    
    //- Cells at left and right side
    const Cell& cellL, &cellR;

    //- Face ID
    label faceID;

    //- Macro variables at left and right sides
    tensor1 Wl, Wr, primL, primR;

    //- Macro variables after collision
    tensor1 W0, prim0;

    //- Reduced distribution function at left and right side
    tensor1 Hl, Hr, Bl, Br;

    //- Reduced distribution function at cell interface
    tensor1 H0, B0;

    //- Slope of reduced distribution function at cell interface
    tensor1 sH0, sB0;

    //- Equilibrium Maxwellian distribution function at cell interface
    tensor1 Heq, Beq;

    //- Macro flux
    tensor1 flux;

    //- Reduced micro flux
    tensor1 fluxH, fluxB; 

// Constructor

    //- Construct by cells
    explicit Face(const Cell& Lcell, const Cell& Rcell);

// Destructor
    
    ~Face();

// Constructor 

// Member function

    //- Interpolate to face
    void interpolate();

    //- Set time coefficient
    tensor1 getTimeCoef(const scalar dt, const scalar tau);

    //- Calculate 2nd order micro flux
    tensor1 get2ndOrderMicroFlux(const tensor1 a, const tensor1 Xeq, const tensor1 Yeq);
    
    //- Calculate Flux
    void getFlux(const scalar dt);
    void getGKSFlux(const scalar dt);
};

class Solver
{
public:

// Member Variables
    
    //- Basic FVM for cell
    vector<Cell> cells;

    //- Basic FVM for face
    vector<Face> faces;

    //- nth step for iteration 
    label step;

    //- Run time
    scalar runTime;

    //- Time step
    scalar dt;

// Constructor

    Solver();

// Member function 

    //- Initialize the flow region
    void initlize();

    //- Iteration to solve
    void iterate();

    //- Calculate time step
    void getDt();

    //- Reconstruction
    void reconstruct();

    //- Evolution 
    void evolve();

    //- Update
    void update();

    //- Write
    void write(std::string str);

    //- Print information
    void info();

};

}



