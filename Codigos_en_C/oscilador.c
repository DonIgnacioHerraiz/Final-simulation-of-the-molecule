#include <stdio.h>
#include "integracion.h"
#include "funciones_oscilador.h"
#include "random.h"
#include <time.h>



//QUE SOLO ESTE DEFINIDA UNA DE LAS 3
#define SIMULACION
#define ANALISIS
#define GRAFICAS

int main() {
    inicializa_PR(12456); // Inicializa el generador con semilla

    // --- Parámetros físicos y numéricos ---
    double T_fisico = 1500;
    double alfa = 0.5;
    double h = 0.001;
    double kb = 1.0;
    double Temperatura = 1.0;
    double dt = h;
    double m = 1.0;
    double K = 100.0;
    int pasos = (int)(T_fisico / dt);

    // --- Conjuntos de tamaños de cadena ---
    int N_s[5] = {4, 8, 16, 32, 64};

    // --- Conjunto de fuerzas constantes ---
    #ifdef FIXED
    int N_fuerzas = 15;
    double F_cte_vals[15] = {
        0.001, 0.00215443, 0.00464159, 0.01, 0.0215443,
        0.0464159, 0.1, 0.148698, 0.215443, 0.464159,
        1.0, 2.15443, 4.47214, 10.0, 20.0
    };
    #endif

    // --- Ruta base para guardar resultados ---
    char carpeta[256];
    #ifdef FIXED
    sprintf(carpeta, "PARAMETROS/%.1f/FIJOS", K);
    #else
    sprintf(carpeta, "PARAMETROS/%.1f/ESCALA", K);
    #endif

    // --- Bucle principal ---
    #ifdef FIXED
    #ifdef SIMULACION
    int N_actual = 4;  // Fijo si estás en modo FIXED
    printf("Simulando con N = %d\n", N_actual);

    // Inicialización de posiciones y velocidades
    double x_0[3 * N_actual];
    double v_0[3 * N_actual];
    for (int j = 0; j < N_actual; j++) {
        x_0[3*j]   = j;
        x_0[3*j+1] = 0.0;
        x_0[3*j+2] = 0.0;
        v_0[3*j]   = v_0[3*j+1] = v_0[3*j+2] = 0.0;
    }

    // Bucle sobre fuerzas constantes
    for (int f = 0; f < N_fuerzas; f++) {
        double F_cte = F_cte_vals[f];
        printf("  -> F_cte = %.3f\n", F_cte);

        clock_t inicio = clock();
        Verlet(K, kb, Temperatura, alfa, N_actual, dt, m, pasos, Fuerza_verlet, x_0, v_0, F_cte);
        clock_t fin = clock();

        double tiempo_total = (double)(fin - inicio) / CLOCKS_PER_SEC;
        escribir_tiempo_en_ultimo_archivo(tiempo_total, carpeta, "V");
    }
            #endif
    #ifdef ANALISIS
    procesar_trayectorias_carpeta(K,5);
    #endif
    #ifdef GRAFICAS
    generar_grafica(K);
    #endif
    #else
    // --- Caso sin FIXED ---
    
    #ifdef SIMULACION
    for (int i = 0; i < 5; i++) {
        int N_actual = N_s[i];
        printf("Simulando con N = %d\n", N_actual);

        double x_0[3 * N_actual];
        double v_0[3 * N_actual];
        for (int j = 0; j < N_actual; j++) {
            x_0[3*j]   = j;
            x_0[3*j+1] = 0.0;
            x_0[3*j+2] = 0.0;
            v_0[3*j]   = v_0[3*j+1] = v_0[3*j+2] = 0.0;
        }

        clock_t inicio = clock();
        Verlet(K, kb, Temperatura, alfa, N_actual, dt, m, pasos, Fuerza_verlet, x_0, v_0);
        clock_t fin = clock();

        double tiempo_total = (double)(fin - inicio) / CLOCKS_PER_SEC;
        escribir_tiempo_en_ultimo_archivo(tiempo_total, carpeta, "V");
    }
        #endif
        #ifdef ANALISIS
        procesar_trayectorias_carpeta(K,5);
        #endif
        #ifdef GRAFICAS
        generar_grafica(K);
        #endif
    #endif

    return 0;

}
