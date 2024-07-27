//head documents
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>//专门用于文件的输入和输出

using namespace std;

// constant expression
constexpr double u1 = 0, rho1 = 1, p1 = 1;
constexpr double u2 = 0, rho2 = 0.125, p2 = 0.1;
constexpr double gama = 1.4;
constexpr double delta = 1e-6;//小量
constexpr double L = 1.0;
constexpr double t_end = 0.14;
constexpr int N = 200;

constexpr double Coef1 = (gama+1)/(2*gama);
constexpr double Coef2 = (gama-1)/(2*gama);

vector<vector<double>> x(N, vector<double>(4));//二维数组
vector<double> t(10);



double soundspeed(double p, double rho)
{
    return sqrt(gama*p/rho);
}

//mass flux
double Ai(double rhoi, double pi, double pstar)
{
    return rhoi*soundspeed(pi, rhoi)*sqrt(Coef1*(pstar/pi)+Coef2);
}

double fp(double pstar, double pi, double rhoi)
{
    if (pstar>pi)
    {
        return (pstar-pi)/Ai(rhoi,pi,pstar);
    }
    else
    {
        return 2*soundspeed(pi, rhoi)/(gama-1)*(pow((pstar/pi),Coef2)-1);
    }
}

double F(double pstar)
{
    return fp(pstar, p1, rho1)+fp(pstar,p2,rho2);
}

double diff_F(double pstar)
{
    return (F(pstar+delta)-F(pstar-delta))/(2*delta);
}

double solve(double pstar)
{
    while (abs(F(pstar))>delta)
    {
        pstar = pstar - F(pstar)/diff_F(pstar);
    }
    return pstar;
}

void initial()
{
    for (int i=0; i < t.size(); i++) {
        t[i] = i*t_end/(t.size()-1);
    }

    for (int i = 0; i<N; i++)
    {
        x[i][0] = i*L/(N-1)-0.5;
        if (x[i][0]<0)
        {
            x[i][1] = rho1;
            x[i][2] = u1;
            x[i][3] = p1;
        }
        else
        {
            x[i][1] = rho2;
            x[i][2] = u2;
            x[i][3] = p2;
        }
    }
}

int main()
{
    initial();
    double p34 = solve(0.5);
    double u34 = u1- fp(p34, p1, rho1);
    double Z_r = u2-(p34-p2)/(rho2*(u2-u34));
    double rho4 = rho2*(u2-Z_r)/(u34-Z_r);

    double c_L = (u1+2*soundspeed(p1,rho1)/(gama-1)-u34)*(gama-1)/2;
    double rho3 = gama*p34/pow(c_L,2);

    double Z_h = u1-soundspeed(p1,rho1);
    double Z_t = u34-c_L;

    for (int k = 0; k<t.size(); k++)
    {
        for (int i = 0; i<N; i++)
        {
            if (x[i][0] >= Z_h*t[k] && x[i][0] < Z_t*t[k])
            {
                double c_expand = (gama-1)/(gama+1)*(u1-x[i][0]/t[k])+2/(gama+1)*soundspeed(p1,rho1);
                x[i][2] = c_expand+x[i][0]/t[k];
                x[i][3] = p1*pow((c_expand/soundspeed(p1,rho1)),1/Coef2);
                x[i][1] = gama*x[i][3]/pow(c_expand,2);
            }
            else if(x[i][0] >= Z_t*t[k] && x[i][0] < u34*t[k])
            {
                x[i][1] = rho3;
                x[i][2] = u34;
                x[i][3] = p34;
            }
            else if(x[i][0] >= u34*t[k] && x[i][0] < Z_r*t[k])
            {
                x[i][1] = rho4;
                x[i][2] = u34;
                x[i][3] = p34;
            }
        }
        

        string filename = "sod"+ to_string(t[k]) +".dat";
        ofstream file(filename, ios::out);
        file << "VARIABLES=x,rho,u,p\n";
        file << "ZONE I=" << N << "\n";

        for (int i = 0; i<N; i++)
        {
            for (int j=0; j<4; j++)
            {
                file << x[i][j] <<" ";
            }
            file << "\n";
        }

        file.close();
        
    }
    
    return 0;
}