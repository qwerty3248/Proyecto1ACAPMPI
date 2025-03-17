#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <stdlib.h>
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
void DFT(double complex *Fourier,const double *muestras, const int n){
    for (int i=0;i<n;i+=4){
        Fourier[i] = 0.0 + 0.0*I;
        Fourier[i+1] = 0.0 + 0.0*I;
        Fourier[i+2] = 0.0 + 0.0*I;
        Fourier[i+3] = 0.0 + 0.0*I;
        for (int j=0;j<n;j+=4){
            //Primera cuadra
            double angle = -2.0*PI*i*j/n;
            double angle1 = -2.0*PI*(i+1)*j/n;
            double angle2 = -2.0*PI*(i+2)*j/n;
            double angle3 = -2.0*PI*(i+3)*j/n;
            
            //Segunda cuadra
            double angle14 = -2.0*PI*i*(j+1)/n;
            double angle11 = -2.0*PI*(i+1)*(j+1)/n;
            double angle21 = -2.0*PI*(i+2)*(j+1)/n;
            double angle31 = -2.0*PI*(i+3)*(j+1)/n;
            
            //tercerca cuadra
            double angle12 = -2.0*PI*i*(j+2)/n;
            double angle112 = -2.0*PI*(i+1)*(j+2)/n;
            double angle212 = -2.0*PI*(i+2)*(j+2)/n;
            double angle312 = -2.0*PI*(i+3)*(j+2)/n;
            
            
            //cuarta cuadra
            double angle13 = -2.0*PI*i*(j+3)/n;
            double angle113 = -2.0*PI*(i+1)*(j+3)/n;
            double angle213 = -2.0*PI*(i+2)*(j+3)/n;
            double angle313 = -2.0*PI*(i+3)*(j+3)/n;
            
            //Para i
            Fourier[i] += muestras[j]*cexp(I*angle);
            Fourier[i] += muestras[j+1]*cexp(I*angle14);
            Fourier[i] += muestras[j+2]*cexp(I*angle12);
            Fourier[i] += muestras[j+3]*cexp(I*angle13);
            //Para i+1
            Fourier[i+1] += muestras[j]*cexp(I*angle1);
            Fourier[i+1] += muestras[j+1]*cexp(I*angle11);
            Fourier[i+1] += muestras[j+2]*cexp(I*angle112);
            Fourier[i+1] += muestras[j+3]*cexp(I*angle113);
            //Para i+2
            Fourier[i+2] += muestras[j]*cexp(I*angle2);
            Fourier[i+2] += muestras[j+1]*cexp(I*angle21);
            Fourier[i+2] += muestras[j+2]*cexp(I*angle212);
            Fourier[i+2] += muestras[j+3]*cexp(I*angle213);
            //Para i+3
            Fourier[i+3] += muestras[j]*cexp(I*angle3);
            Fourier[i+3] += muestras[j+1]*cexp(I*angle31);
            Fourier[i+3] += muestras[j+2]*cexp(I*angle312);
            Fourier[i+3] += muestras[j+3]*cexp(I*angle313);
            
        }
    }    
}

//Fourier continuo 🤯🤯
void CFT(double complex *Fourier, const double *muestras, const int n, const double paso_temporal){
    for (int i=0;i<n;i+=4){
        Fourier[i] = 0.0 + 0.0*I;
        Fourier[i+1] = 0.0 + 0.0*I;
        Fourier[i+2] = 0.0 + 0.0*I;
        Fourier[i+3] = 0.0 + 0.0*I;
        double omega = 2.0*PI*i/(T_MAX-T_MIN);//La w de la formula que es el omega
        double omega2 = 2.0*PI*(i+1)/(T_MAX-T_MIN);
        double omega3 = 2.0*PI*(i+2)/(T_MAX-T_MIN);
        double omega4 = 2.0*PI*(i+3)/(T_MAX-T_MIN);
        //Aqui ya entra en juego tanto el intervalo del tiempo como los valores que dan la funcion, es decir el tiempo esta entre -1 y 1
        //y los valores de la funcion están en mis muestras
        //Vamos desde el minimo hasta el maximo pero con nuestro paso temporal para tomar fourier lo más preciso posible
        for (double j = T_MIN;j<T_MAX;j+=paso_temporal*4){
            //printf("Estoy en el segundo FOR");
            //Dudas a como acceder con el paso temporal, preguntar a Pablo
            
            int indice = (int)((j - T_MIN)/paso_temporal);
            int indice2 = (int)((j+paso_temporal - T_MIN)/paso_temporal);
            int indice3 = (int)((j+(2*paso_temporal) - T_MIN)/paso_temporal);
            int indice4 = (int)((j+(3*paso_temporal) - T_MIN)/paso_temporal);

            //Para i

            if (indice <= 0 ){
                //Si es menor o igual a 0 suponemos que coge el primer elemento
                Fourier[i] += muestras[0]*cexp(-I*omega*j) * paso_temporal;
            }else if (indice >= n){ 
                //Si es mayor o igual al numero de elemnto que hay, se coge el ultimo
                Fourier[i] += muestras[n-1]*cexp(-I*omega*j) * paso_temporal;
            }else{
                //Si no es valido y cogemos el valor calculado
                Fourier[i] += muestras[indice]*cexp(-I*omega*j) * paso_temporal;
            }

            if (indice2 <= 0 ){
                //Si es menor o igual a 0 suponemos que coge el primer elemento
                Fourier[i] += muestras[0]*cexp(-I*omega*(j+paso_temporal)) * paso_temporal;
            }else if (indice2 >= n){ 
                //Si es mayor o igual al numero de elemnto que hay, se coge el ultimo
                Fourier[i] += muestras[n-1]*cexp(-I*omega*(j+paso_temporal)) * paso_temporal;
            }else{
                //Si no es valido y cogemos el valor calculado
                Fourier[i] += muestras[indice2]*cexp(-I*omega*(j+paso_temporal)) * paso_temporal;
            }

            if (indice3 <= 0 ){
                //Si es menor o igual a 0 suponemos que coge el primer elemento
                Fourier[i] += muestras[0]*cexp(-I*omega*(j+(2*paso_temporal))) * paso_temporal;
            }else if (indice3 >= n){ 
                //Si es mayor o igual al numero de elemnto que hay, se coge el ultimo
                Fourier[i] += muestras[n-1]*cexp(-I*omega*(j+(2*paso_temporal))) * paso_temporal;
            }else{
                //Si no es valido y cogemos el valor calculado
                Fourier[i] += muestras[indice3]*cexp(-I*omega*(j+(2*paso_temporal))) * paso_temporal;
            }

            if (indice4 <= 0 ){
                //Si es menor o igual a 0 suponemos que coge el primer elemento
                Fourier[i] += muestras[0]*cexp(-I*omega*(j+(3*paso_temporal))) * paso_temporal;
            }else if (indice4 >= n){ 
                //Si es mayor o igual al numero de elemnto que hay, se coge el ultimo
                Fourier[i] += muestras[n-1]*cexp(-I*omega*(j+(3*paso_temporal))) * paso_temporal;
            }else{
                //Si no es valido y cogemos el valor calculado
                Fourier[i] += muestras[indice4]*cexp(-I*omega*(j+(3*paso_temporal))) * paso_temporal;
            }


            //Para i+1

            if (indice <= 0 ){
                //Si es menor o igual a 0 suponemos que coge el primer elemento
                Fourier[i+1] += muestras[0]*cexp(-I*omega2*j) * paso_temporal;
            }else if (indice >= n){ 
                //Si es mayor o igual al numero de elemnto que hay, se coge el ultimo
                Fourier[i+1] += muestras[n-1]*cexp(-I*omega2*j) * paso_temporal;
            }else{
                //Si no es valido y cogemos el valor calculado
                Fourier[i+1] += muestras[indice]*cexp(-I*omega2*j) * paso_temporal;
            }

            if (indice2 <= 0 ){
                //Si es menor o igual a 0 suponemos que coge el primer elemento
                Fourier[i+1] += muestras[0]*cexp(-I*omega2*(j+paso_temporal)) * paso_temporal;
            }else if (indice2 >= n){ 
                //Si es mayor o igual al numero de elemnto que hay, se coge el ultimo
                Fourier[i+1] += muestras[n-1]*cexp(-I*omega2*(j+paso_temporal)) * paso_temporal;
            }else{
                //Si no es valido y cogemos el valor calculado
                Fourier[i+1] += muestras[indice2]*cexp(-I*omega2*(j+paso_temporal)) * paso_temporal;
            }

            if (indice3 <= 0 ){
                //Si es menor o igual a 0 suponemos que coge el primer elemento
                Fourier[i+1] += muestras[0]*cexp(-I*omega2*(j+(2*paso_temporal))) * paso_temporal;
            }else if (indice3 >= n){ 
                //Si es mayor o igual al numero de elemnto que hay, se coge el ultimo
                Fourier[i+1] += muestras[n-1]*cexp(-I*omega2*(j+(2*paso_temporal))) * paso_temporal;
            }else{
                //Si no es valido y cogemos el valor calculado
                Fourier[i+1] += muestras[indice3]*cexp(-I*omega2*(j+(2*paso_temporal))) * paso_temporal;
            }

            if (indice4 <= 0 ){
                //Si es menor o igual a 0 suponemos que coge el primer elemento
                Fourier[i+1] += muestras[0]*cexp(-I*omega2*(j+(3*paso_temporal))) * paso_temporal;
            }else if (indice4 >= n){ 
                //Si es mayor o igual al numero de elemnto que hay, se coge el ultimo
                Fourier[i+1] += muestras[n-1]*cexp(-I*omega2*(j+(3*paso_temporal))) * paso_temporal;
            }else{
                //Si no es valido y cogemos el valor calculado
                Fourier[i+1] += muestras[indice4]*cexp(-I*omega2*(j+(3*paso_temporal))) * paso_temporal;
            }

            //Para i+2
            if (indice <= 0 ){
                //Si es menor o igual a 0 suponemos que coge el primer elemento
                Fourier[i+2] += muestras[0]*cexp(-I*omega3*j) * paso_temporal;
            }else if (indice >= n){ 
                //Si es mayor o igual al numero de elemnto que hay, se coge el ultimo
                Fourier[i+2] += muestras[n-1]*cexp(-I*omega3*j) * paso_temporal;
            }else{
                //Si no es valido y cogemos el valor calculado
                Fourier[i+2] += muestras[indice]*cexp(-I*omega3*j) * paso_temporal;
            }

            if (indice2 <= 0 ){
                //Si es menor o igual a 0 suponemos que coge el primer elemento
                Fourier[i+2] += muestras[0]*cexp(-I*omega3*(j+paso_temporal)) * paso_temporal;
            }else if (indice2 >= n){ 
                //Si es mayor o igual al numero de elemnto que hay, se coge el ultimo
                Fourier[i+2] += muestras[n-1]*cexp(-I*omega3*(j+paso_temporal)) * paso_temporal;
            }else{
                //Si no es valido y cogemos el valor calculado
                Fourier[i+2] += muestras[indice2]*cexp(-I*omega3*(j+paso_temporal)) * paso_temporal;
            }

            if (indice3 <= 0 ){
                //Si es menor o igual a 0 suponemos que coge el primer elemento
                Fourier[i+2] += muestras[0]*cexp(-I*omega3*(j+(2*paso_temporal))) * paso_temporal;
            }else if (indice3 >= n){ 
                //Si es mayor o igual al numero de elemnto que hay, se coge el ultimo
                Fourier[i+2] += muestras[n-1]*cexp(-I*omega3*(j+(2*paso_temporal))) * paso_temporal;
            }else{
                //Si no es valido y cogemos el valor calculado
                Fourier[i+2] += muestras[indice3]*cexp(-I*omega3*(j+(2*paso_temporal))) * paso_temporal;
            }

            if (indice4 <= 0 ){
                //Si es menor o igual a 0 suponemos que coge el primer elemento
                Fourier[i+2] += muestras[0]*cexp(-I*omega3*(j+(3*paso_temporal))) * paso_temporal;
            }else if (indice4 >= n){ 
                //Si es mayor o igual al numero de elemnto que hay, se coge el ultimo
                Fourier[i+2] += muestras[n-1]*cexp(-I*omega3*(j+(3*paso_temporal))) * paso_temporal;
            }else{
                //Si no es valido y cogemos el valor calculado
                Fourier[i+2] += muestras[indice4]*cexp(-I*omega3*(j+(3*paso_temporal))) * paso_temporal;
            }


            //Para i+3
            if (indice <= 0 ){
                //Si es menor o igual a 0 suponemos que coge el primer elemento
                Fourier[i+3] += muestras[0]*cexp(-I*omega4*j) * paso_temporal;
            }else if (indice >= n){ 
                //Si es mayor o igual al numero de elemnto que hay, se coge el ultimo
                Fourier[i+3] += muestras[n-1]*cexp(-I*omega4*j) * paso_temporal;
            }else{
                //Si no es valido y cogemos el valor calculado
                Fourier[i+3] += muestras[indice]*cexp(-I*omega4*j) * paso_temporal;
            }

            if (indice2 <= 0 ){
                //Si es menor o igual a 0 suponemos que coge el primer elemento
                Fourier[i+3] += muestras[0]*cexp(-I*omega4*(j+paso_temporal)) * paso_temporal;
            }else if (indice2 >= n){ 
                //Si es mayor o igual al numero de elemnto que hay, se coge el ultimo
                Fourier[i+3] += muestras[n-1]*cexp(-I*omega4*(j+paso_temporal)) * paso_temporal;
            }else{
                //Si no es valido y cogemos el valor calculado
                Fourier[i+3] += muestras[indice2]*cexp(-I*omega4*(j+paso_temporal)) * paso_temporal;
            }

            if (indice3 <= 0 ){
                //Si es menor o igual a 0 suponemos que coge el primer elemento
                Fourier[i+3] += muestras[0]*cexp(-I*omega4*(j+(2*paso_temporal))) * paso_temporal;
            }else if (indice3 >= n){ 
                //Si es mayor o igual al numero de elemnto que hay, se coge el ultimo
                Fourier[i+3] += muestras[n-1]*cexp(-I*omega4*(j+(2*paso_temporal))) * paso_temporal;
            }else{
                //Si no es valido y cogemos el valor calculado
                Fourier[i+3] += muestras[indice3]*cexp(-I*omega4*(j+(2*paso_temporal))) * paso_temporal;
            }

            if (indice4 <= 0 ){
                //Si es menor o igual a 0 suponemos que coge el primer elemento
                Fourier[i+3] += muestras[0]*cexp(-I*omega4*(j+(3*paso_temporal))) * paso_temporal;
            }else if (indice4 >= n){ 
                //Si es mayor o igual al numero de elemnto que hay, se coge el ultimo
                Fourier[i+3] += muestras[n-1]*cexp(-I*omega4*(j+(3*paso_temporal))) * paso_temporal;
            }else{
                //Si no es valido y cogemos el valor calculado
                Fourier[i+3] += muestras[indice4]*cexp(-I*omega4*(j+(3*paso_temporal))) * paso_temporal;
            }

        }
    }    
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
const char * SALIDA = "txt/SecuencialDFTOPT.txt"; //salida secuencial (cualquier caso) 😇😇
const char * SALIDA2 = "txt/SecuencialDFT2OPT.txt"; //salida secuencial para el caso de 20.000 muestras 😇😇
const char * SALIDA_CONTINUO = "txt/ContinuoDFTOPT.txt"; //salida continuo 🤯🤯
const char * SALIDA_CONTINUO2 = "txt/ContinuoDFTOPT2.txt"; //salida continuo 🤯🤯
const char * SALIDA_CONTINUO2 = "txt/ContinuoDFTOPT3.txt"; //salida continuo 🤯🤯
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int main (){
    FILE * entrada = fopen (TREINTAMIL,"r");
    FILE * salida = fopen(SALIDA,"w");
    FILE * salida_cont = fopen(SALIDA_CONTINUO,"w");
    FILE * salida_cont2 = fopen(SALIDA_CONTINUO2,"w");
    FILE * salida_cont3 = fopen(SALIDA_CONTINUO2,"w");

    if (!entrada){
        printf("Error: No se pudo abrir el archivo de entrada\n");
        exit(1);
    }
    if (!salida){
        printf("Error: No se pudo abrir el archivo de salida\n");
        exit(1);
    }

    if (!salida_cont){
        printf("Error: No se pudo abrir el archivo de salida continuo\n");
        exit(1);
    }
    if (!salida_cont2){
        printf("Error: No se pudo abrir el archivo de salida continuo\n");
        exit(1);
    }
    if (!salida_cont3){
        printf("Error: No se pudo abrir el archivo de salida continuo\n");
        exit(1);
    }
    
    int Tam_Vector_muestras;
    while (fscanf(entrada,"%d",&Tam_Vector_muestras) == 1){
        
        double *muestras = malloc(Tam_Vector_muestras*sizeof(double));
        double complex *Fourier = malloc(Tam_Vector_muestras*sizeof(double complex));
        double complex *Fourier_cont = malloc(Tam_Vector_muestras*sizeof(double complex));
        double complex *Fourier_cont2 = malloc(Tam_Vector_muestras*sizeof(double complex));
        double complex *Fourier_cont3 = malloc(Tam_Vector_muestras*sizeof(double complex));
        
        if (!muestras){
            printf("Error: No se pudo asignar memoria muestras\n");
            exit(1);
        }
        if (!Fourier){
            printf("Error: No se pudo asignar memoria Fourier\n");
            exit(1);
        }
        if (!Fourier_cont){
            printf("Error: No se pudo asignar memoria Fourier\n");
            exit(1);
        }
        if (!Fourier_cont2){
            printf("Error: No se pudo asignar memoria Fourier\n");
            exit(1);
        }
        if (!Fourier_cont3){
            printf("Error: No se pudo asignar memoria Fourier\n");
            exit(1);
        }

        for (int i=0;i<Tam_Vector_muestras;i++){
            if (fscanf(entrada,"%lf",&muestras[i]) != 1){
                printf("Error: archivo de entrada de la muestras %d\n",i);
                exit(1);
            }
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////////
        //Parte de DFT
        //////////////////////////////////////////////////////////////////////////////////////////////////////////
        clock_t inicio = clock();
        DFT(Fourier,muestras,Tam_Vector_muestras);
        clock_t fin = clock();
        double tiempo = (double)(fin-inicio)/CLOCKS_PER_SEC;
        //El tiempo se pasa en ms
        fprintf(salida,"%d %lf\n",Tam_Vector_muestras,tiempo);
        //Para mostrar el vector, esta dentro del archivo 😴😴
        /*for (int i=0;i<Tam_Vector_muestras;i++){
            fprintf(salida,"%lf %lf\n",creal(Fourier[i]),cimag(Fourier[i]));
        }*/

        //////////////////////////////////////////////////////////////////////////////////////////////////////////
        //Parte de CFT
        //////////////////////////////////////////////////////////////////////////////////////////////////////////
        //printf("valor tammuestras %d\n",Tam_Vector_muestras);
        //printf("valor tmax - tmin %d\n",T_MAX-T_MIN);
        //printf("valor paso temporal %f\n", (double)(T_MAX-T_MIN) / Tam_Vector_muestras);
        double paso_temporal = (double)(T_MAX-T_MIN) / Tam_Vector_muestras;
        //printf("Valor paso temporal: %f\n",paso_temporal);
        clock_t inicio_cont = clock();
        CFT(Fourier_cont,muestras,Tam_Vector_muestras,paso_temporal);
        clock_t fin_cont = clock();
        double tiempo_cont = (double)(fin_cont-inicio_cont)/CLOCKS_PER_SEC;
        fprintf(salida_cont,"%d %lf\n",Tam_Vector_muestras,tiempo_cont);
        /*for (int i=0;i<Tam_Vector_muestras;i++){
            fprintf(salida_cont,"%lf %lf\n",creal(Fourier_cont[i]),cimag(Fourier_cont[i]));
        }*/
        //////////////////////////////////////////////////////////////////////////////////////////////////////////

        free(muestras);
        free(Fourier);
        free(Fourier_cont);
        free(Fourier_cont2);
        free(Fourier_cont3);

    }

    fclose(entrada);
    fclose(salida);
    fclose(salida_cont);
    fclose(salida_cont2);
    fclose(salida_cont3);

    printf("Fin del programa, resultados guardados con formato NumeroMuestras Tiempo(Fourier Discreto) en salidaDFTOPT\n");
    printf("Fin del programa, resultados guardados con formato NumeroMuestras Tiempo(Fourier continuo) en salidaCFTOPT\n");
    printf("Fin del programa, resultados guardados con formato NumeroMuestras Tiempo(Fourier continuo Simpson) en salidaCFTOPT\n");
    printf("Fin del programa, resultados guardados con formato NumeroMuestras Tiempo(Fourier continuo Trapezio) en salidaCFTOPT\n");

    return 0;
    






}
