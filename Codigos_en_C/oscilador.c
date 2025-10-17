#include <stdio.h>
#include "integracion.h"
#include "funciones_oscilador.h"
#include "random.h"
#include <time.h>

int main() {
    inicializa_PR(12456); // Inicializa el generador con semilla

    // --- Parámetros físicos y numéricos ---
    double T_fisico = 150000;
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
    double F_cte_vals[5][6] = {
    {0,0.577,  2.887,  5.773, 11.547, 28.867},  // N = 4
    {0,0.378,  1.890,  3.780,  7.560, 18.900},  // N = 8
    {0,0.258,  1.291,  2.582,  5.164, 12.910},  // N = 16
    {0,0.180,  0.898,  1.796,  3.592,  8.980},  // N = 32
    {0,0.126,  0.630,  1.260,  2.520,  6.300}   // N = 64
};

    #endif

    // --- Bucle principal sobre N ---
    for (int i = 0; i < 5; i++) {
        int N_actual = N_s[i];
        printf("Simulando con N = %d\n", N_actual);

        // --- Inicialización de posiciones y velocidades ---
        double x_0[3 * N_actual];
        double v_0[3 * N_actual];

        for (int j = 0; j < N_actual; j++) {
            x_0[3*j]   = j;    // x
            x_0[3*j+1] = 0.0;  // y
            x_0[3*j+2] = 0.0;  // z

            v_0[3*j]   = 0.0;  // vx
            v_0[3*j+1] = 0.0;  // vy
            v_0[3*j+2] = 0.0;  // vz
        }

        // --- Ruta base para guardar ---
        char carpeta[256];
        #ifdef FIXED
        sprintf(carpeta, "PARAMETROS/%.1f/FIJOS", K);
        #else
        sprintf(carpeta, "PARAMETROS/%.1f/ESCALA", K);
        #endif

        // --- Bucle sobre F_cte si FIXED está definido ---
        #ifdef FIXED
        for (int f = 0; f < 6; f++) {
            double F_cte = F_cte_vals[i][f];

            printf("  -> F_cte = %.2f\n", F_cte);

            clock_t inicio = clock();

            Verlet(K, kb, Temperatura, alfa, N_actual, dt, m, pasos, Fuerza_verlet, x_0, v_0, F_cte);

            clock_t fin = clock();
            double tiempo_total = (double)(fin - inicio) / CLOCKS_PER_SEC;

            escribir_tiempo_en_ultimo_archivo(tiempo_total, carpeta, "V");
        }
        #else
        // --- Caso sin FIXED ---
        clock_t inicio = clock();

        Verlet(K, kb, Temperatura, alfa, N_actual, dt, m, pasos, Fuerza_verlet, x_0, v_0);

        clock_t fin = clock();
        double tiempo_total = (double)(fin - inicio) / CLOCKS_PER_SEC;

        escribir_tiempo_en_ultimo_archivo(tiempo_total, carpeta, "V");
        #endif
    }

    return 0;
}
