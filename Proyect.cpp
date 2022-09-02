#include<cstdlib>
#include<fstream> 
#include<iostream>
#include<random>
#include<cmath>
using namespace std;


//Valor aleatorio con distribuion uniforme (0,1)
double p; 

//Tiempo maximo de la simulacion (1 Hora)
double time_max = 3600;

//Parametros para la distribucion Extreme Values (Valido para las regiones North, South y West)
double mu = 2.76958;//Location parameter
double k = 0.476206;//Shape parameter
double sigma = 2.1753;// Scale Parameter

//Parametros para la distribuicon Birnabum Sanders(Valido para la region East)
double beta = 0; //Scale parameter
double gamma_ = 0; //Shape parameter



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
                            ofstream& data, 
                            double (*distribution_taken)(double,double,double,double), 
                            double location_parameter, 
                            double shape_parameter, 
                            double scale_parameter
                            )
{
    //Definimos el generador de numeros aleatorios para una distribuicon uniforme
    default_random_engine generator(time(NULL));
    uniform_real_distribution<double> distribution(0.0,1.0);

    //Seteamos los parametros iniciales
    double t = 0; //contador de tiempo
    int N = 1; //contador de numero de vehiculos en la cola

    //Creamos el encabezado del archivo de salida
    data<<"Delta_t"<<" "<<"Tiempo_acumulado"<<" "<<"#Autos"<<endl;

    //Definimos la simulacion
    while(t <= time_max)
    {
        double p = distribution(generator);
        double delta_t = distribution_taken(p, location_parameter, shape_parameter, scale_parameter);
        t += delta_t;
        N++;

        //Guardamos los datos
        data<<delta_t<<" "<<t<<" "<<N<<endl;
    }
}


int main(){

  //Creamos el fichero para guardar los archivos 
  ofstream Data ("ICDF_Distribution.dat");

  //Prueba simulacion
  queue_generaration(time_max,Data,ICDF_ExtremeValues, mu, k, sigma);

  //Cerramos el archivo
  Data.close();

  



    return 0;
}