#pragma once

#include "random.h"
#include "funciones_oscilador.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <math.h>


#ifdef FIXED
void un_paso_verlet(double betta[], double b, double a, int N, 
                    double x_antiguo[], double x_nuevo[], 
                    double v_antiguo[], double v_nuevo[], 
                    double F_antiguo[], double F_nuevo[], 
                    double dt, double m, 
                    void (*Fuerza)(int, double[], double[], double, double), 
                    double K, double F_cte);
#else
void un_paso_verlet(double betta[], double b, double a, int N, 
                    double x_antiguo[], double x_nuevo[], 
                    double v_antiguo[], double v_nuevo[], 
                    double F_antiguo[], double F_nuevo[], 
                    double dt, double m, 
                    void (*Fuerza)(int, double[], double[], double), 
                    double K);
#endif

/**
 * Realiza la integración de la trayectoria usando el método de Verlet y guarda los resultados en un archivo.
 */
#ifdef FIXED
void verlet_trayectoria(char* filename_input, double kb, double Temperatura, double alfa, int N,
                        double dt, double m, int pasos,
                        void (*Fuerza)(int, double[], double[], double, double),
                        char* filename_output, double x_0[], double v_0[], double K, double F_cte);
#else
void verlet_trayectoria(char* filename_input, double kb, double Temperatura, double alfa, int N,
                        double dt, double m, int pasos,
                        void (*Fuerza)(int, double[], double[], double),
                        char* filename_output, double x_0[], double v_0[], double K);
#endif

/**
 * Escribe en un fichero los parámetros de la simulación de Verlet.
 */
#ifdef FIXED
void escribe_input_verlet(double kb, double Temperatura, double alfa, int N, double dt, double m, int pasos,
                          double x_0[], double v_0[], char filename[], double K, double F_cte);
#else
void escribe_input_verlet(double kb, double Temperatura, double alfa, int N, double dt, double m, int pasos,
                          double x_0[], double v_0[], char filename[], double K);
#endif

/**
 * Función principal para realizar la simulación completa de Verlet.
 */
#ifdef FIXED
void Verlet(double K, double kb, double Temperatura, double alfa, int N, double dt, double m, int pasos,
            void (*Fuerza)(int, double[], double[], double, double),
            double x_0[], double v_0[], double F_cte);
#else
void Verlet(double K, double kb, double Temperatura, double alfa, int N, double dt, double m, int pasos,
            void (*Fuerza)(int, double[], double[], double),
            double x_0[], double v_0[]);
#endif