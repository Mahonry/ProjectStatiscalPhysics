#include<cstdlib>
#include<fstream> 
#include<iostream>
#include<random>
#include<cmath>
#include <math.h> 
using namespace std;

//Numero de veces a efectuar la simulacion
int maximum_repetitions = 50000;

//Declaramos la variable de nombre
string namef;

//Tiempo maximo de la simulacion (1 Hora)
double time_max = 3600;

//Parametros NORTH para la distribucion Extreme Values 
double mu_north = 2.76;//Location parameter
double k_north = 0.47;//Shape parameter
double sigma_north = 2.17;// Scale Parameter

//Parametros EAST para la distribuicon Birnabum Sanders
double sigma_east = 11.26; //Scale parameter
double mu_east = 0; //Location parameter
double k_east = 1.1;//Shape parameter 

//Parametros SOUTH para la distribucion Extreme Values
double mu_south = 2.55;//731;//Location parameter
double k_south = 0.4128;//Shape parameter
double sigma_south = 2;//1.99457;// Scale Parameter

//Parametros WEST para la distribucion Extreme Values
double mu_west = 3.9742;//Location parameter
double k_west = 0.8642;//Shape parameter
double sigma_west = 3.019;// Scale Parameter




//Distribucion Generalized Extreme Values
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

//Distribucion Birnbaum Sanders
double ICDF_Birnbaum_Sanders(double p, double mu, double k, double sigma)
{
    double z = 1.0*rand()/RAND_MAX;
    double aux1 = 0.25*sigma;
    double aux2 = k*z;
    double aux3 = pow(aux2,2) + 4;
    double aux4 = sqrt(aux3) + aux2;
    return (aux1)*pow(aux4,2);
}


//GENERADOR DE LA COLA
double queue_generation(double time_max, 
                            ofstream& data, 
                            double (*distribution_selected)(double,double,double,double), 
                            double location_parameter, 
                            double shape_parameter, 
                            double scale_parameter,
                            int iteration,
                            double max_distribution
                            )
{
    //Definimos el generador de numeros aleatorios para una distribuicon uniforme
    default_random_engine generator(time(NULL));
    uniform_real_distribution<double> distribution(0.0,max_distribution);

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
                        double sigma,
                        string region,
                        double max_distribution)
{
    //Nombre del archivo
    namef = "Results/_" + region + "_Complete_trend_analysis.dat";
    
    //Creamos el fichero para guardar los archivos 
    ofstream Data (namef);
    
    //Creamos el encabezado del archivo de salida
    Data<<"Delta_t"<<" "<<"Tiempo_acumulado"<<" "<<"#Autos"<<" "<<"Iteration"<<endl;

    for(int i = 0; i < maximum_repetitions; i++)
    {

        //Generamos la simulacion
        queue_generation(time_max,Data,distribution_selected, mu, k, sigma, i,max_distribution);


    }

    //Cerramos el archivo
    Data.close();
}





//Formador de la cola 4 links 
double queue_formation(  double Q_0,
                        double V_first,
                        double t_g, 
                        double t_r,
                        double S_flow,
                        double V_g, 
                        double time_max,
                        double max_distribution,
                        double (*distribution_selected)(double,double,double,double),
                        double mu,
                        double k, 
                        double sigma
                        )
{
    //Declaramos el valor final de la cola
    float Q_total;

    //Declaramos Suturated flow rate at green
    float S_g = (t_g*S_flow)/3600;

    //Declaramos el contador de tiempo
    float t = 0;

    //Creamos el archivo donde se van a guardar los datos 
    //Nombre del archivo
    namef = "Results/queue_formation.dat";
    
    //Creamos el fichero para guardar los archivos 
    ofstream Data (namef);
    
    //Creamos el encabezado del archivo de salida
    Data<<"Delta_t"<<" "<<"Tiempo_acumulado"<<" "<<"#Autos"<<" "<<"Iteration"<<endl;

    float i = 0; //Variable para la funcion queue_generation 



    while(t < time_max)
    {


        //Calculamos el arribo en rojo
        V_first = queue_generation(t_r,Data,distribution_selected, mu, k, sigma, i,max_distribution);
        t += t_r;
        
        //Calculamos el arribo en verde
        V_g = queue_generation(t_g,Data,distribution_selected, mu, k, sigma, i,max_distribution);
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


    }

    return Q_total;
}


void for_cola_four_links(int maximum_repetitions)
{
    ofstream Data ("Cola_final.dat");
    Data<<"N"<<" "<<"S"<<" "<<"W"<<" "<<"E"<<" "<<"Tg"<<endl;
    for(int cycle = 60; cycle <= 150; cycle += 10)
    {
        for(int j = 0; j < maximum_repetitions; j++)
        {
            if(j%10000 == 0)//Ciclo para checar que esta funcionando 
            {
                cout<<cycle<<" "<<j<<endl;
            }

            
            Data<<queue_formation(0,0,cycle,cycle*3,1800,0,time_max,1,ICDF_ExtremeValues,mu_north,k_north,sigma_north)<<" "
                <<queue_formation(0,0,cycle,cycle*3,1800,0,time_max,1,ICDF_ExtremeValues,mu_south,k_south,sigma_south)<<" "
                <<queue_formation(0,0,cycle,cycle*3,900,0,time_max,0.999,ICDF_ExtremeValues,mu_west,k_west,sigma_west)<<" "
                <<queue_formation(0,0,cycle,cycle*3,900,0,time_max,1,ICDF_Birnbaum_Sanders,mu_east,k_east,sigma_east)<<" "
                <<cycle<<endl;
            

        }
    }
    Data.close();
}


double vehicle_delay()
{

}


int main(){

//Definimos el generador de numeros gaussianos
    srand(time(NULL));

//Generamos los archivos para hacer el trend analysis
    //string region = "South_180";
    //queue_trend_analysis(maximum_repetitions,ICDF_ExtremeValues,mu_south,k_south,sigma_south,region,1); //.99 para WEST .99999 North//ICDF_Birnbaum_Sanders EAST

//Generamos la cola final para los 4 links 
    for_cola_four_links(maximum_repetitions);



    return 0;
}