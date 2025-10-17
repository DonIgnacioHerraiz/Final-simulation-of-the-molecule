#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <math.h>
#include <errno.h>


#define L_0 1.0

#define FIXED //DEFINIR SI HAY UN EXTREMO FIJO

#ifndef FIXED
void Fuerza_verlet(int N, double x[], double F[], double K);
#else
void Fuerza_verlet(int N, double x[], double F[], double K, double F_cte);
#endif

void Fuerza_euler(int N, double x[], double p[], double F[], double K,double eta, double m);

double Energia_cinetica_instantanea(int N, double v[], double m);

double Energia_potencial_instantanea(int N, double x[], double m, double K);

double Energia_total_instantanea(int N, double x[], double v[], double m, double K);

void escribir_tiempo_en_ultimo_archivo(double tiempo, const char *carpeta, const char *prefijo);

double calcula_radio_giro(int N, double *x);

void procesar_trayectoria(char* archivo_input, int N_start, int N, double K 
    #ifdef FIXED
        , double F_cte
    #endif
                          );

#ifdef FIXED
void procesar_trayectorias_carpeta(double K, int N_start, double F_cte);
#else
void procesar_trayectorias_carpeta(double K, int N_start);
#endif

// Cabecera de la función auxiliar para leer parámetros
int leer_N_desde_parametros(const char *archivo_parametros);

void generar_grafica(double K);