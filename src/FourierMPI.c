#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <mpi.h>
/*
    IMPORTANTE
Las muestras que se les da a las funciones son los valores x[n] de la DFT y x(t) de la CFT, si no por ejemplo
en la CFT necesitariamos definir una funcion f que tome valores en el tiempo para poder hacer el CFT, para más
simplicidad, he supuesto que estos valores son los que tienen mis muestras, ademas de acotarlas entre [-1,1].

*/

#define PI 3.14159265358979323846
//Las muestras están entre -1 y 1
#define T_MAX 1
#define T_MIN -1



//Fourier discreto 😇😇
/*void DFT(double complex *Fourier,const double *muestras, const int n){
    for (int i=0;i<n;i++){
        Fourier[i] = 0.0 + 0.0*I;
        for (int j=0;j<n;j++){
            double angle = -2.0*PI*i*j/n;
            Fourier[i] += muestras[j]*cexp(I*angle);
        }
    }    
}*/

//Fourier discreto pero la version MPI

void DFTMPI(double complex *Fourier, const double *muestras, const int n, const int rank, const int size){
    
    int inicio = rank*(n/size);
    int fin = (rank == size-1)? n : inicio + (n/size);
    
    for (int i=inicio;i<fin;i++){
        Fourier[i] = 0.0 + 0.0*I;
        for (int j=0;j<n;j++){
            double angle = -2.0*PI*i*j/n;
            Fourier[i] += muestras[j]*cexp(I*angle);
        }
    }

    //Ahora una vez ejecutado todo esto, hacemos un gather para que le den los resultados de todos los procesos a el 0
    MPI_Gather(Fourier + inicio,(fin-inicio),MPI_DOUBLE_COMPLEX,Fourier,(fin-inicio),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD);

}


//Fourier continuo 🤯🤯
/*void CFT(double complex *Fourier, const double *muestras, const int n, const double paso_temporal){
    for (int i=0;i<n;i++){
        Fourier[i] = 0.0 + 0.0*I;
        double omega = 2.0*PI*i/(T_MAX-T_MIN);//La w de la formula que es el omega
        //Aqui ya entra en juego tanto el intervalo del tiempo como los valores que dan la funcion, es decir el tiempo esta entre -1 y 1
        //y los valores de la funcion están en mis muestras
        //Vamos desde el minimo hasta el maximo pero con nuestro paso temporal para tomar fourier lo más preciso posible
        for (double j = T_MIN;j<T_MAX; j= j+paso_temporal){
            //printf("Estoy en el segundo FOR");
            //Dudas a como acceder con el paso temporal, preguntar a Pablo
            int indice = (int)((j - T_MIN)/paso_temporal);
            if (indice <= 0 ){
                //Si es menor o igual a 0 suponemos que coge el primer elemento
                Fourier[i] += muestras[0]*cexp(-I*omega*j) * paso_temporal;
            }else if (indice >= n-1){ 
                //Si es mayor o igual al numero de elemnto que hay, se coge el ultimo
                Fourier[i] += muestras[n-1]*cexp(-I*omega*j) * paso_temporal;
            }else{
                //Si no es valido y cogemos el valor calculado
                Fourier[i] += muestras[indice]*cexp(-I*omega*j) * paso_temporal;
            }
            

        }
    }    
}*/


//Fourier continuo 🤯🤯
void CFTMPI(double complex *Fourier, const double *muestras, const int n, const double paso_temporal, const int rank, const int size){
    int inicio = rank*(n/size);
    int fin = (rank == size-1)? n : inicio + (n/size);
    for (int i=inicio;i<fin;i++){
        Fourier[i] = 0.0 + 0.0*I;
        double omega = 2.0*PI*i/(T_MAX-T_MIN);//La w de la formula que es el omega
        //Aqui ya entra en juego tanto el intervalo del tiempo como los valores que dan la funcion, es decir el tiempo esta entre -1 y 1
        //y los valores de la funcion están en mis muestras
        //Vamos desde el minimo hasta el maximo pero con nuestro paso temporal para tomar fourier lo más preciso posible
        for (double j = T_MIN;j<T_MAX; j= j+paso_temporal){
            //printf("Estoy en el segundo FOR");
            //Dudas a como acceder con el paso temporal, preguntar a Pablo
            int indice = (int)((j - T_MIN)/paso_temporal);
            if (indice <= 0 ){
                //Si es menor o igual a 0 suponemos que coge el primer elemento
                Fourier[i] += muestras[0]*cexp(-I*omega*j) * paso_temporal;
            }else if (indice >= n-1){ 
                //Si es mayor o igual al numero de elemnto que hay, se coge el ultimo
                Fourier[i] += muestras[n-1]*cexp(-I*omega*j) * paso_temporal;
            }else{
                //Si no es valido y cogemos el valor calculado
                Fourier[i] += muestras[indice]*cexp(-I*omega*j) * paso_temporal;
            }
            

        }
    }   
    
    MPI_Gather(Fourier + inicio,(fin-inicio),MPI_DOUBLE_COMPLEX,Fourier,(fin-inicio),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD);
}

//Apartador de archivos de entrada y saluda
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const char * MILLON = "txt/MuestraGenerada.txt"; //una muestra de un millo de elementos 🫨🫨
const char * TREINTAMIL = "txt/MuestraGenerada30000.txt"; //una muestra de un 30.000 elementos 🫨🫨
const char * CINCUENTAMIL = "txt/MuestraGenerada50000.txt";
const char * CIENMIL = "txt/MuestraGenerada100000.txt";
const char * CIENTOCINCUENTAMIL = "txt/MuestraGenerada150000.txt";
const char * DOSCIENTOSCINCUENTAMIL = "txt/MuestraGenerada250000.txt";
const char * QUINIENTOSMIL = "txt/MuestraGenerada500000.txt";
const char * DISTINTAS20 = "txt/muestras.txt"; //20.000 muestras de elementos aleatorios 😎😎
const char * FUNCIONA = "txt/funciona.txt"; //muestras de funcionamiento 😎😎
const char * FUNCIONA2 = "txt/funciona2.txt"; //muestras de funcionamiento 😎😎
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const char * SALIDA = "txt/SecuencialDFTMPI.txt"; //salida secuencial (cualquier caso) 😇😇
const char * SALIDA2 = "txt/SecuencialDFT.txt"; //salida secuencial para el caso de 20.000 muestras 😇😇
const char * SALIDA_CONTINUO = "txt/ContinuoDFTMPI.txt"; //salida continuo 🤯🤯
const char * SALIDA_CONTINUO2 = "txt/ContinuoDFT2.txt"; //salida continuo 🤯🤯
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main (){
    
    MPI_Init(NULL,NULL);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    FILE * salida;
    FILE * salida_cont;
    if (rank == 0){
        //printf("Archivo de salida: %s\n",new_salida);
        salida = fopen(SALIDA,"w");
        salida_cont = fopen(SALIDA_CONTINUO,"w");
    
        
        if (!salida){
            printf("Error: No se pudo abrir el archivo de salida\n");
            exit(1);
        }
    
        if (!salida_cont){
            printf("Error: No se pudo abrir el archivo de salida continuo\n");
            exit(1);
        }



    }//Solo quiero que reserve memoria en el 0, este se encargara de escribir, este tendra el fourier

    FILE * entrada = fopen (CIENMIL,"r");
    if (!entrada){
        printf("Error: No se pudo abrir el archivo de entrada\n");
        exit(1);
    } 

    //printf("Soy el procesador %d y voy a leer el archivo de entrada\n",rank);
    
    int Tam_Vector_muestras;
    double *muestras;
    while (fscanf(entrada,"%d",&Tam_Vector_muestras) == 1){
        
        if (rank == 0){muestras = malloc(Tam_Vector_muestras*sizeof(double));}
        double complex *Fourier = malloc(Tam_Vector_muestras*sizeof(double complex));
        double complex *Fourier_cont = malloc(Tam_Vector_muestras*sizeof(double complex));
        double *muestras_local = malloc((Tam_Vector_muestras / size) * sizeof(double));
        
        //printf("Soy el procesador %d y voy a asignar memoria\n",rank);
        if (rank == 0){
            if (!muestras){
                printf("Error: No se pudo asignar memoria muestras\n");
                exit(1);
            }
        }
        
        if (!Fourier){
            printf("Error: No se pudo asignar memoria Fourier\n");
            exit(1);
        }
        if (!Fourier_cont){
            printf("Error: No se pudo asignar memoria Fourier\n");
            exit(1);
        }
        if (rank == 0){
            for (int i=0;i<Tam_Vector_muestras;i++){
                if (fscanf(entrada,"%lf",&muestras[i]) != 1){
                    printf("Error: archivo de entrada de la muestras %d\n",i);
                    exit(1);
                }
            }
        }
        

        MPI_Scatter(muestras, Tam_Vector_muestras / size, MPI_DOUBLE, muestras_local, Tam_Vector_muestras / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //printf("Soy el procesador %d y voy a salir\n",rank);
        //////////////////////////////////////////////////////////////////////////////////////////////////////////
        //Parte de DFT
        //////////////////////////////////////////////////////////////////////////////////////////////////////////
        //clock_t inicio = clock();
        //printf("Soy el procesador %d y voy a iniciar el tiempo\n",rank);
        MPI_Barrier(MPI_COMM_WORLD);
        double primer_tiempo = MPI_Wtime();
        DFTMPI(Fourier,muestras_local,Tam_Vector_muestras/size,rank,size);
        double segundo_tiempo = MPI_Wtime();
        //clock_t fin = clock();
        //double tiempo = (double)(fin-inicio)/CLOCKS_PER_SEC;
        //El tiempo se pasa en ms
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0){fprintf(salida,"%d %lf\n",Tam_Vector_muestras,segundo_tiempo-primer_tiempo);}
        //Para mostrar el vector, esta dentro del archivo 😴😴
        if (rank == 0){
            for (int i=0;i<Tam_Vector_muestras;i++){
                fprintf(salida,"%lf %lf\n",creal(Fourier[i]),cimag(Fourier[i]));
            }
        }
        

        //////////////////////////////////////////////////////////////////////////////////////////////////////////
        //Parte de CFT
        //////////////////////////////////////////////////////////////////////////////////////////////////////////
        //printf("valor tammuestras %d\n",Tam_Vector_muestras);
        //printf("valor tmax - tmin %d\n",T_MAX-T_MIN);
        //printf("valor paso temporal %f\n", (double)(T_MAX-T_MIN) / Tam_Vector_muestras);
        double paso_temporal = (double)(T_MAX-T_MIN) / Tam_Vector_muestras;
        //printf("Valor paso temporal: %f\n",paso_temporal);
        MPI_Barrier(MPI_COMM_WORLD);
        double tercer_tiempo = MPI_Wtime();
        //clock_t inicio_cont = clock();
        CFTMPI(Fourier_cont,muestras_local,Tam_Vector_muestras/size,paso_temporal,rank,size);
        double cuarto_tiempo = MPI_Wtime();
        MPI_Barrier(MPI_COMM_WORLD);
        //clock_t fin_cont = clock();
        //double tiempo_cont = (double)(fin_cont-inicio_cont)/CLOCKS_PER_SEC;
        if (rank == 0){fprintf(salida_cont,"%d %lf\n",Tam_Vector_muestras,cuarto_tiempo-tercer_tiempo);}
        if (rank == 0){
            for (int i=0;i<Tam_Vector_muestras;i++){
                fprintf(salida_cont,"%lf %lf\n",creal(Fourier_cont[i]),cimag(Fourier_cont[i]));
            }
        }
        
        //////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        if (rank == 0){free(muestras);}
        free(Fourier);
        free(Fourier_cont);

    }
    if (rank == 0){
        fclose(salida);
        fclose(salida_cont);
    }

    fclose(entrada);
    
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0){
        printf("Fin del programa, resultados guardados con formato NumeroMuestras Tiempo(Fourier Discreto) en salidaDFT\n");
        printf("Fin del programa, resultados guardados con formato NumeroMuestras Tiempo(Fourier continuo) en salidaCFT\n");
    }
    MPI_Finalize();
    return 0;
    
}