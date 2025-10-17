#include "integracion.h"

/**
 * Realiza un paso en la integración del movimiento usando el método de Verlet.
 * @param betta       Array con términos aleatorios para el ruido térmico.
 * @param b           Coeficiente dependiente de alfa y dt.
 * @param a           Coeficiente dependiente de alfa y dt.
 * @param N           Número de partículas o elementos.
 * @param x_antiguo   Array con las posiciones en el paso de tiempo anterior.
 * @param x_nuevo     Array donde se almacenarán las nuevas posiciones.
 * @param v_antiguo   Array con las velocidades en el paso de tiempo anterior.
 * @param v_nuevo     Array donde se almacenarán las nuevas posiciones.
 * @param dt          Paso de tiempo.
 * @param m           Masa de las partículas (asumida igual para todas).
 * @param Fuerza      Puntero a función que calcula las fuerzas; recibe N, las posiciones las velocidades y un real.
 * @param F_antiguo   Array donde se almacenarán las fuerzas en el paso de tiempo anterior.
 * @param F_nuevo     Array donde se almacenarán las nuevas fuerzas.
 * @param K           Constante del oscilador armónico.
 */

 #ifdef FIXED
void un_paso_verlet(double betta[], double b, double a, int N, double x_antiguo[], double x_nuevo[],
                    double v_antiguo[], double v_nuevo[], double F_antiguo[], double F_nuevo[],
                    double dt, double m, void (*Fuerza)(int, double[], double[], double, double),
                    double K, double F_cte)
#else
void un_paso_verlet(double betta[], double b, double a, int N, double x_antiguo[], double x_nuevo[],
                    double v_antiguo[], double v_nuevo[], double F_antiguo[], double F_nuevo[],
                    double dt, double m, void (*Fuerza)(int, double[], double[], double),
                    double K)
#endif
{
    // Actualización de posiciones
    for (int i = 0; i < 3*N; i++) {
        x_nuevo[i] = x_antiguo[i] + v_antiguo[i]*dt*b + F_antiguo[i]*dt*dt*b/(2*m) + b*dt*betta[i];
    }

    // Cálculo de nuevas fuerzas
    #ifdef FIXED
        Fuerza(N, x_nuevo, F_nuevo, K, F_cte);
    #else
        Fuerza(N, x_nuevo, F_nuevo, K);
    #endif

    // Actualización de velocidades
    for (int i = 0; i < 3*N; i++) {
        v_nuevo[i] = a*v_antiguo[i] + (a*F_antiguo[i] + F_nuevo[i])*dt/(2*m) + b*betta[i]/m;
    }
}
/**
    * Realiza la integración del movimiento usando el método de Verlet durante un número dado de pasos.
    * Los resultados se guardan en un archivo.
    * @param kb              Constante de Boltzmann.
    * @param Temperatura     Temperatura del sistema.
    * @param alfa            Coeficiente de fricción.
    * @param N               Número de partículas o elementos.
    * @param dt              Paso de tiempo.  
    * @param m               Masa de las partículas (asumida igual para todas).
    * @param pasos           Número de pasos de tiempo a simular.
    * @param Fuerza          Puntero a función que calcula las fuerzas; recibe N y un array de posiciones.
    * @param filename_output Nombre del archivo donde se guardarán los resultados.
    * @param x_0            Array con las posiciones iniciales.
    * @param v_0            Array con las velocidades iniciales.
 */

#ifdef FIXED
void verlet_trayectoria(char* filename_input, double kb, double Temperatura, double alfa, int N,
                        double dt, double m, int pasos,
                        void (*Fuerza)(int, double[], double[], double, double),
                        char* filename_output, double x_0[], double v_0[], double K, double F_cte)
#else
void verlet_trayectoria(char* filename_input, double kb, double Temperatura, double alfa, int N,
                        double dt, double m, int pasos,
                        void (*Fuerza)(int, double[], double[], double),
                        char* filename_output, double x_0[], double v_0[], double K)
#endif
{
    FILE *archivo = fopen(filename_output, "w");
    if (!archivo) {
        printf("Error al abrir el archivo %s\n", filename_output);
        return;
    }

    double a = (1.0 - alfa * dt / (2.0 * m)) / (1.0 + alfa * dt / (2.0 * m));
    double b = 1.0 / (1.0 + alfa * dt / (2.0 * m));

    fprintf(archivo, "%.6f %d\t%s\n", dt, pasos, filename_input);

    double x_antiguo[3*N];
    double x_nuevo[3*N];
    double v_antiguo[3*N];
    double v_nuevo[3*N];
    double F_antiguo[3*N];
    double F_nuevo[3*N];
    double betta[3*N];
    double Ek, Ep, Et,Rg,Ree;
    double counter = 0;

    for (int i = 0; i < 3*N; i++) {
        x_antiguo[i] = x_0[i];
        v_antiguo[i] = v_0[i];
        F_antiguo[i] = 0.0;
        F_nuevo[i] = 0.0;
    }

    #ifdef FIXED
        Fuerza(N, x_antiguo, F_antiguo, K, F_cte);
    #else
        Fuerza(N, x_antiguo, F_antiguo, K);
    #endif

    for (int paso = 0; paso < pasos; paso++) {
        for (int i = 0; i < 3*N; i++) {
            betta[i] = gaussian() * sqrt(2 * alfa * Temperatura * kb * dt);
        }

        #ifdef FIXED
            un_paso_verlet(betta, b, a, N, x_antiguo, x_nuevo, v_antiguo, v_nuevo,
                           F_antiguo, F_nuevo, dt, m, Fuerza, K, F_cte);
        #else
            un_paso_verlet(betta, b, a, N, x_antiguo, x_nuevo, v_antiguo, v_nuevo,
                           F_antiguo, F_nuevo, dt, m, Fuerza, K);
        #endif

        counter += dt;

        if (counter >= 0.1) {
            fprintf(archivo, "%.6f", paso * dt);
            for (int i = 0; i < 3*N; i++) fprintf(archivo, " %.6f", x_nuevo[i]);
            for (int i = 0; i < 3*N; i++) fprintf(archivo, " %.6f", v_nuevo[i]);
            Ek = Energia_cinetica_instantanea(N, v_nuevo, m);
            Ep = Energia_potencial_instantanea(N, x_nuevo, m, K);
            Et = Energia_total_instantanea(N, x_nuevo, v_nuevo, m, K);
            Rg = calcula_radio_giro(N, x_nuevo);
            Ree=x_nuevo[3*(N-1)+2]-x_nuevo[2];
            fprintf(archivo, " %.6f %.6f %.6f %.6f %.6f\n", Ek, Ep, Et, Rg, Ree);
            counter = 0;
        }

        for (int i = 0; i < 3*N; i++) {
            x_antiguo[i] = x_nuevo[i];
            v_antiguo[i] = v_nuevo[i];
            F_antiguo[i] = F_nuevo[i];
        }
    }
    
    for (int i = 0; i < N; i++) {
            x_0[3*i] = i;
            x_0[3*i+1] = 0;
            x_0[3*i+2] = 0;
            v_0[3*i] = 0;
            v_0[3*i+1] = 0;
            v_0[3*i+2] = 0;
        }


    fclose(archivo);
}


/**
    * Escribe en un fichero los parámetros de la simulación de Verlet. Lo hace en la carpeta PARAMETROS\OSCILADOR con el formato V_i, siendo i el primer número natural tal que no existe un archivo con ese nombre en la carpeta.
    * @param kb              Constante de Boltzmann.
    * @param Temperatura     Temperatura del sistema.
    * @param alfa           Coeficiente de fricción.
    * @param N               Número de partículas o elementos.
    * @param dt              Paso de tiempo.  
    * @param m               Masa de las partículas (asumida igual para todas).
    * @param pasos           Número de pasos de tiempo a simular.
    * @param x_0            Array con las posiciones iniciales.
    * @param v_0            Array con las velocidades iniciales.
    * @param filename       Devuelve el nombre del archivo creado.
 */

#ifdef FIXED
void escribe_input_verlet(double kb, double Temperatura, double alfa, int N, double dt, double m, int pasos,
                          double x_0[], double v_0[], char filename[], double K, double F_cte) {
#else
void escribe_input_verlet(double kb, double Temperatura, double alfa, int N, double dt, double m, int pasos,
                          double x_0[], double v_0[], char filename[], double K) {
#endif
    
    char buffer[256];

    #ifdef FIXED
    #ifdef WLCM
    sprintf(buffer, "PARAMETROS/WLCM/%.1f/FIJOS", K);
    const char* folder = buffer;
    #else
    sprintf(buffer, "PARAMETROS/%.1f/FIJOS", K);
    const char* folder = buffer;
    #endif
    #else 
    #ifdef WLCM
    sprintf(buffer, "PARAMETROS/WLCM/%.1f/ESCALA", K);
    const char* folder = buffer;
    #else
    sprintf(buffer, "PARAMETROS/%.1f/ESCALA", K);
    const char* folder = buffer;
    #endif
    #endif

    FILE* file;
    int k = 0;

    // Buscar el primer archivo V_k.txt que no exista
    while (1) {
        snprintf(filename, 256, "%s/V_%d.txt", folder, k);
        file = fopen(filename, "r");
        if (file) {
            fclose(file);
            k++;
        } else {
            break;
        }
    }

    // Crear archivo nuevo
    file = fopen(filename, "w");
    if (!file) {
        printf("No se pudo crear el archivo %s\n", filename);
        filename[0] = '\0';
        return;
    }

    // --- Cabecera informativa ---
    fprintf(file, "# Archivo de parámetros para simulación de Verlet\n");
    fprintf(file, "# Generado automáticamente por escribe_input_verlet()\n");
    fprintf(file, "# -----------------------------------------------\n\n");

    // --- Parámetros físicos y numéricos ---
    fprintf(file, "K %g\n", K);
    fprintf(file, "kb %g\n", kb);
    fprintf(file, "Temperatura %g\n", Temperatura);
    fprintf(file, "alfa %g\n", alfa);
    fprintf(file, "N %d\n", N);
    fprintf(file, "dt %g\n", dt);
    fprintf(file, "m %g\n", m);
    fprintf(file, "pasos %d\n", pasos);

    // --- Información sobre FIXED ---
    #ifdef FIXED
        fprintf(file, "Modo FIXED: SI\n");
        fprintf(file, "F_cte %g\n", F_cte);
        fprintf(file, "# Nota: La primera partícula está fija.\n");
        fprintf(file, "# Se aplica una fuerza constante F_cte en la dirección z sobre la última partícula.\n");
    #else
        fprintf(file, "Modo FIXED: NO\n");
    #endif

    // --- Condiciones iniciales ---
    fprintf(file, "\n# Posiciones iniciales:\n");
    for (int i = 0; i < 3*N; i++) {
        fprintf(file, "x_0_%d %g\n", i, x_0[i]);
    }

    fprintf(file, "\n# Velocidades iniciales:\n");
    for (int i = 0; i < 3*N; i++) {
        fprintf(file, "v_0_%d %g\n", i, v_0[i]);
    }

    fclose(file);

    printf("Archivo creado: %s\n", filename);
}

/**
    * Lleva a cabo la simulacion completa de Verlet. Los parametros estan en la carpeta PARAMETROS\OSCILADOR, mientras que la trayectoria está en Resultados_simulacion\OSCILADOR\VERLET
    * @param kb              Constante de Boltzmann.
    * @param Temperatura     Temperatura del sistema.
    * @param alfa           Coeficiente de fricción.
    * @param N               Número de partículas o elementos.
    * @param dt              Paso de tiempo.  
    * @param m               Masa de las partículas (asumida igual para todas).
    * @param pasos           Número de pasos de tiempo a simular.
    * @param x_0            Array con las posiciones iniciales.
    * @param v_0            Array con las velocidades iniciales.
    * @param Fuerza       Puntero a función que calcula las fuerzas; recibe N y un array de posiciones.
 */


#ifdef FIXED
void Verlet(double K, double kb, double Temperatura, double alfa, int N, double dt, double m, int pasos,
            void (*Fuerza)(int, double[], double[], double, double),
            double x_0[], double v_0[], double F_cte)
#else
void Verlet(double K, double kb, double Temperatura, double alfa, int N, double dt, double m, int pasos,
            void (*Fuerza)(int, double[], double[], double),
            double x_0[], double v_0[])
#endif
{
    char filename_input[256];

    // --- Crear archivo de parámetros ---
    #ifdef FIXED
        escribe_input_verlet(kb, Temperatura, alfa, N, dt, m, pasos, x_0, v_0, filename_input, K, F_cte);
    #else
        escribe_input_verlet(kb, Temperatura, alfa, N, dt, m, pasos, x_0, v_0, filename_input, K);
    #endif

    // --- Selección de carpeta de salida ---
    char buffer[256];
    #ifdef FIXED
    #ifdef WLCM
    sprintf(buffer, "Resultados_simulacion/WLCM/%.1f/FIJOS", K);
    const char* folder = buffer;
    #else
    sprintf(buffer, "Resultados_simulacion/%.1f/FIJOS", K);
    const char* folder = buffer;
    #endif
    #else 
    #ifdef WLCM
    sprintf(buffer, "Resultados_simulacion/WLCM/%.1f/ESCALA", K);
    const char* folder = buffer;
    #else
    sprintf(buffer, "Resultados_simulacion/%.1f/ESCALA", K);
    const char* folder = buffer;
    #endif
    #endif
    

    char filename_output[256];
    FILE* file;
    int i = 0;

    // Buscar el primer archivo V_i.txt que no exista
    while (1) {
        snprintf(filename_output, sizeof(filename_output), "%s/V_%d.txt", folder, i);
        file = fopen(filename_output, "r");
        if (file) {
            fclose(file);
            i++;
        } else {
            break;
        }
    }

    // Crear el nuevo archivo
    file = fopen(filename_output, "w");
    if (!file) {
        printf("No se pudo crear el archivo %s\n", filename_output);
        return;
    }
    fclose(file);

    // --- Ejecutar simulación Verlet ---
    #ifdef FIXED
        verlet_trayectoria(filename_input, kb, Temperatura, alfa, N, dt, m, pasos,
                           Fuerza, filename_output, x_0, v_0, K, F_cte);
    #else
        verlet_trayectoria(filename_input, kb, Temperatura, alfa, N, dt, m, pasos,
                           Fuerza, filename_output, x_0, v_0, K);
    #endif
}
