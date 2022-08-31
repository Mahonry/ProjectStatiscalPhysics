#include<cstdlib>
#include<fstream> 
#include<iostream>
#include<random>
#include<cmath>
using namespace std;

double p; //Valor aleatorio con distribuion uniforme (0,1)
double mu = 2.76958;//Location parameter
double k = 0.476206;//Shape parameter
double sigma = 2.1753;// Scale Parameter

//Definimos la funcion inversa de la distribucion generalizado de valores extremos acumulada
double ICDF_ExtremeValues(double p,double mu,double k,double sigma)
{
    if(k == 0){
        return mu-sigma*log(-log(p));
    }
    else{
        return ((pow(-1*log(p),-k)-1)*(sigma/k))+mu;
    }
}


int main(){
  
  //Definimos el generador de numeros aleatorios para una distribuicon uniforme
  default_random_engine generator(time(NULL));
  uniform_real_distribution<double> distribution(0.0,1.0);
  
  //Creamos el fichero para guardar los archivos 
  ofstream Data ("ICDF_Distribution.dat");


  //Ciclo For para Validar la distribucion
  for (int i = 0; i<100000; i++){
      double p = distribution(generator);
      double time = ICDF_ExtremeValues(p,mu,k,sigma);
      Data<<time<<endl;
  }



    return 0;
}