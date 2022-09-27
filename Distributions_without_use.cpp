#include<cstdlib>
#include<fstream> 
#include<iostream>
#include<random>
#include<cmath>
#include <math.h> 
using namespace std;


//ERF Inverse Function
float myErfInv2(float x){

   float sgn = (x < 0) ? -1.0f : 1.0f;

   x = (1 - x)*(1 + x);        // x = 1 - x*x;
   float lnx = logf(x);

   float aux1 = 2/(M_PI*0.147) + 0.5f * lnx;
   float aux2 = 1/(0.147) * lnx;

   return(sgn*sqrtf(-aux1 + sqrtf(aux1*aux1 - aux2)));
}

//double ICDF_normal_distribution(double beta, double gamma)
double ICDF_normal_distribution(double p, double mu,double sigma)
{   
    float aux = 2*p-1;
    return sqrt(2)*sigma*myErfInv2(aux)+mu;

}


double ICDF_Birnbaum_Sanders_2(double p, double mu, double k, double sigma)
{
    double aux1 = pow(k,2);
    double mu_normal = sigma*(1+(aux1*0.5));
    double aux2 = pow(sigma*k,2);
    double sigma_normal =  aux2*(1+(1.25*aux2));

    double normal = ICDF_normal_distribution(p,mu_normal,sigma_normal);

    double aux3 = k*normal;
    double aux4 = pow(aux3,2)+4;
    double aux5 = sqrt(aux4) + aux3;
        return (0.25*sigma)*pow(aux5,2);
}