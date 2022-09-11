#include<cstdlib>
#include<fstream> 
#include<iostream>
#include<random>
#include<cmath>
#include <math.h> 
using namespace std;

//Numero de veces a efectuar la simulacion
int maximum_repetitions = 50000;

//Valor aleatorio con distribuion uniforme (0,1)
double p; 

//Declaramos la variable de nombre
string namef;

//Tiempo maximo de la simulacion (1 Hora)
double time_max = 3600;

//Parametros para la distribucion Extreme Values (Valido para las regiones North, South y West)
double mu_north = 2.76;//Location parameter
double k_north = 0.47;//Shape parameter
double sigma_north = 2.17;// Scale Parameter

//Parametros para la distribuicon Birnabum Sanders(Valido para la region East)
double sigma_east = 11.2647; //Scale parameter
double mu_east = 0; //Location parameter
double k_east = 1.0928;//Shape parameter 



//Definimos la funcion inversa de la distribucion generalizado de valores extremos acumulada
double ICDF_ExtremeValues(double p,double mu,double k,double sigma)
{
    if(k == 0){
        return mu-sigma*log(-log(p));
    }
    else{
        double a = pow(-1*log(p),-1*k);
        double b = -1*k*mu*pow(-log(p),k);
        double c = sigma * pow(-log(p),k) - sigma;
        return (-1)*(a*(b+c) )/k;
    }
}

//ERF Inverse Function
float myErfInv2(float x){
   float tt1, tt2, lnx, sgn;
   sgn = (x < 0) ? -1.0f : 1.0f;

   x = (1 - x)*(1 + x);        // x = 1 - x*x;
   lnx = logf(x);

   tt1 = 2/(M_PI*0.147) + 0.5f * lnx;
   tt2 = 1/(0.147) * lnx;

   return(sgn*sqrtf(-tt1 + sqrtf(tt1*tt1 - tt2)));
}



//double ICDF_normal_distribution(double beta, double gamma)
double ICDF_normal_distribution(double p, double mu, double k, double sigma)
{   
    float aux = 2*p-1;
    return sqrt(2)*myErfInv2(aux)+mu;

}

double ICDF_Birnbaum_Sanders(double p, double mu, double k, double sigma)
{
    double z = 1.0*rand()/RAND_MAX;
    double aux1 = sigma*0.25;
    double aux2 = k*z;
    double aux3 = pow(aux2,2) + 4;
    double aux4 = sqrt(aux3) + aux2;
    return (0.25)*pow(aux4,2);
}





//GENERADOR DE LA COLA
float queue_generation(double time_max, 
                            ofstream& data, 
                            double (*distribution_selected)(double,double,double,double), 
                            double location_parameter, 
                            double shape_parameter, 
                            double scale_parameter,
                            int iteration
                            )
{
    //Definimos el generador de numeros aleatorios para una distribuicon uniforme
    default_random_engine generator(time(NULL));
    uniform_real_distribution<double> distribution(0.0,1.0);

    //Seteamos los parametros iniciales
    double t = 0; //contador de tiempo
    int N = 1; //contador de numero de vehiculos en la cola

    //Definimos la simulacion
    while(t < time_max)
    {
        double p = distribution(generator);
        double delta_t = distribution_selected(p, 
                                            location_parameter, 
                                            shape_parameter, 
                                            scale_parameter);
        t += delta_t;
        N++;

        //Guardamos los datos
        data<<delta_t<<" "<<t<<" "<<N<<" "<<iteration<<endl;
    }

    return N;
}

//Analisis de tendencia 
void queue_trend_analysis(int maximum_repetitions,
                        double (*distribution_selected)(double,double,double,double),
                        double mu,
                        double k,
                        double sigma)
{
    //Nombre del archivo
    namef = "Results/Easth_Complete_trend_analysis.dat";
    
    //Creamos el fichero para guardar los archivos 
    ofstream Data (namef);
    
    //Creamos el encabezado del archivo de salida
    Data<<"Delta_t"<<" "<<"Tiempo_acumulado"<<" "<<"#Autos"<<" "<<"Iteration"<<endl;

    for(int i = 0; i < maximum_repetitions; i++)
    {

        //Generamos la simulacion
        queue_generation(time_max,Data,distribution_selected, mu, k, sigma, i);


    }

    //Cerramos el archivo
    Data.close();
}

//Formador de la cola 4 links 
float queue_formation(  float Q_0,
                        float V_first,
                        float t_g, //En este caso tomaremos time green (t_g) and time red (t_r) iguales
                        float t_r,
                        float S_flow,
                        float V_g, 
                        float time_max
                        )
{
    //Declaramos el valor final de la cola
    float Q_total;

    //Declaramos Suturated flow rate at green
    float S_g = (t_g*S_flow)/3600;

    //Declaramos el contador de tiempo
    float t = 0;

    //Probamos con un cycle lenght de 60as
    float  cycle_lenght = 60;

    //Creamos el archivo donde se van a guardar los datos 
    //Nombre del archivo
    namef = "Results/queue_formation.dat";
    
    //Creamos el fichero para guardar los archivos 
    ofstream Data (namef);
    
    //Creamos el encabezado del archivo de salida
    Data<<"Delta_t"<<" "<<"Tiempo_acumulado"<<" "<<"#Autos"<<" "<<"Iteration"<<endl;

    float i = 0; //Variable para la funcion queue_generation 



    while(t <= time_max)
    {


        //Calculamos el arrivo en rojo
        V_first = queue_generation(t_r,Data,ICDF_ExtremeValues, mu_north, k_north, sigma_north, i);
        t += t_r;
        
        //Calculamos el arribo en verde
        V_g = queue_generation(t_g,Data,ICDF_ExtremeValues, mu_north, k_north, sigma_north, i);
        t += t_g;
        
        //Calculamos la cola final
        float control = Q_0 + V_first + V_g - S_g;

        if (control <= 0)
        {
            Q_0 = 0;
        }
        else
        {
            Q_0 = Q_0 + V_first + V_g - S_g;
        }

        Q_total += Q_0;
            
        i += 1;

        //Probamos imprimir
        cout<<"Arribo rojo "<<V_first<<endl;
        cout<<"Arribo verde "<<V_g<<endl;
        cout<<"Control "<<control<<endl;
        cout<<"Tiempo "<<t<<endl;
        cout<<"_________"<<endl;

    }

    return Q_total;
}


int main(){

    //Definimos el generador de numeros gaussianos
    srand(time(NULL));

//Generamos los archivos para hacer el trend analysis
    queue_trend_analysis(maximum_repetitions,ICDF_Birnbaum_Sanders,mu_east,k_east,sigma_east);
    //cout<<queue_formation(0,0,60,180,3600,0,time_max)<<endl;

  



    return 0;
}