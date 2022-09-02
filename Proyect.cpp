#include<cstdlib>
#include<fstream> 
#include<iostream>
#include<random>
#include<cmath>
using namespace std;



double p; //Valor aleatorio con distribuion uniforme (0,1)

//Parametros para la distribucion Extreme Values (Valido para las regiones North, South y West)
double mu = 0;//Location parameter
double k = 0.5;//Shape parameter
double sigma = 1;// Scale Parameter

//Parametros para la distribuicon Birnabum Sanders(Valido para la region East)
double beta = 0; //Scale parameter
double gamma = 0; //Shape parameter



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

//double ICDF_Birnbaum_Sanders(double beta, double gamma)




//GENERADOR DE LA COLA

void queue_generaration(double time_max, 
                            ofstream data, 
                            double (*distribution_taken)(double,double,double,double), 
                            double location_parameter, 
                            double shape_parameter, 
                            double scale_parameter, 
                            )
{
    //Definimos el generador de numeros aleatorios para una distribuicon uniforme
    default_random_engine generator(time(NULL));
    uniform_real_distribution<double> distribution(0.0,1.0);

    //Seteamos los parametros iniciales

    double t = 0; //contador de tiempo
    int N = 1; //contador de numero de vehiculos en la cola

    while(t <= time_max)
    {

        double p = distribution(generator);
        double delta_t = distribution_taken(p, location_parameter, shape_parameter, scale_parameter);
        t += delta_t;
        N++;

    }
}


int main(){
  

  
  //Creamos el fichero para guardar los archivos 
  ofstream Data ("ICDF_Distribution.dat");

  Data<<"P T"<<endl;
      Data<<p<<" "<<time<<endl;
  



    return 0;
}