#include "funciones_oscilador.h"

#include <sys/stat.h> // mkdir
#include <sys/types.h>

/***********************************************************************
 *  Funcion: void Fuerza_verlet(int N, double x[], double F[], double K)
 *  
 *  Descripción:
 *  Calcula las fuerzas entre partículas consecutivas de un polímero 
 *  unido por resortes (modelo armónico tipo "bond-stretching"), de acuerdo 
 *  con el potencial:
 *      V = (1/2) * K * (|r_{i+1} - r_i| - b)^2
 *  
 *  Si la macro FIXED está definida, se asume que la primera partícula 
 *  está fija (no se calcula su fuerza) y se añade una fuerza constante 
 *  F_cte sobre la última partícula en la dirección z.
 *  
 *  Parámetros:
 *  ------------------------------------------------------------------
 *   N  : número de partículas del polímero.
 *   x[]: vector de posiciones de tamaño 3*N, en el formato:
 *        (x1, y1, z1, x2, y2, z2, ..., xN, yN, zN)
 *   F[]: vector de fuerzas de tamaño 3*N, en el mismo formato que x[].
 *        La función sobrescribe este vector con las fuerzas calculadas.
 *   K  : constante elástica (k_e) del potencial.
 *  
 *  Si FIXED está definida:
 *   F_cte : fuerza constante aplicada sobre la última partícula 
 *            en la componente z (solo sobre z_N).
 *  
 *  Funcionamiento:
 *  ------------------------------------------------------------------
 *  - Se recorren los enlaces entre partículas consecutivas (i, i+1).
 *  - Se calcula la distancia r = |r_{i+1} - r_i|.
 *  - Se aplica la fuerza de tipo Hooke:
 *        F_i   +=  K*(r - b)*(r_{i+1} - r_i)/r
 *        F_{i+1} = -F_i
 *  - Si FIXED está definida, la primera partícula no se modifica 
 *    (su fuerza se mantiene en cero).
 ***********************************************************************/


#ifndef FIXED

void Fuerza_verlet(int N, double x[], double F[], double K)
{
    // Inicializar todas las fuerzas a cero
    for (int i = 0; i < 3*N; i++)
        F[i] = 0.0;

    for (int i = 0; i < N - 1; i++) {
        int i3 = 3*i;
        int j3 = 3*(i+1);

        double dx = x[j3]   - x[i3];
        double dy = x[j3+1] - x[i3+1];
        double dz = x[j3+2] - x[i3+2];

        double r = sqrt(dx*dx + dy*dy + dz*dz);
        double fac = K * (r - L_0) / r;

        double Fx = fac * dx;
        double Fy = fac * dy;
        double Fz = fac * dz;

        F[i3]   +=  Fx;
        F[i3+1] +=  Fy;
        F[i3+2] +=  Fz;

        F[j3]   -=  Fx;
        F[j3+1] -=  Fy;
        F[j3+2] -=  Fz;
    }
}

#else

void Fuerza_verlet(int N, double x[], double F[], double K, double F_cte)
{
    // Dejar F de la primera partícula intacta (está fija)
    for (int i = 3; i < 3*N; i++)
        F[i] = 0.0;

    for (int i = 1; i < N - 1; i++) {
        int i3 = 3*i;
        int j3 = 3*(i+1);

        double dx = x[j3]   - x[i3];
        double dy = x[j3+1] - x[i3+1];
        double dz = x[j3+2] - x[i3+2];

        double r = sqrt(dx*dx + dy*dy + dz*dz);
        double fac = K * (r - L_0) / r;

        double Fx = fac * dx;
        double Fy = fac * dy;
        double Fz = fac * dz;

        F[i3]   +=  Fx;
        F[i3+1] +=  Fy;
        F[i3+2] +=  Fz;

        F[j3]   -=  Fx;
        F[j3+1] -=  Fy;
        F[j3+2] -=  Fz;
    }

    // Añadir fuerza constante sobre la última partícula en dirección z
    F[3*(N-1) + 2] += F_cte;
}

#endif



double Energia_cinetica_instantanea(int N, double v[], double m){
    double K=0;
    for(int i=0;i<3*N;i++){
        K=K+0.5*m*v[i]*v[i];
    }
    return K;
}
double Energia_potencial_instantanea(int N, double x[], double m, double K){
    double V=0;
    for(int i=0;i<N-1;i++){
        double dx=x[3*(i+1)]-x[3*i];
        double dy=x[3*(i+1)+1]-x[3*i+1];
        double dz=x[3*(i+1)+2]-x[3*i+2];
        double r=sqrt(dx*dx+dy*dy+dz*dz);
        V=V+0.5*K*(r-L_0)*(r-L_0);
    }
    return V;
}

double Energia_total_instantanea(int N, double x[], double v[], double m, double K){
    double E=0;
    E=Energia_cinetica_instantanea(N,v,m)+Energia_potencial_instantanea(N,x,m,K);
    return E;
}

void escribir_tiempo_en_ultimo_archivo(double tiempo, const char *carpeta, const char *prefijo) {
    DIR *dir;
    struct dirent *entry;
    int k_max = -1;
    char archivo_objetivo[512];

    dir = opendir(carpeta);
    if (!dir) {
        perror("No se pudo abrir la carpeta");
        return;
    }

    size_t len_prefijo = strlen(prefijo);

    // Recorrer archivos de la carpeta
    while ((entry = readdir(dir)) != NULL) {
        const char *nombre = entry->d_name;

        // Verificar que empieza con prefijo + "_"
        if (strncmp(nombre, prefijo, len_prefijo) == 0 && nombre[len_prefijo] == '_') {
            // Verificar que termina en ".txt"
            const char *ext = strrchr(nombre, '.');
            if (ext && strcmp(ext, ".txt") == 0) {
                // Extraer número entre prefijo_ y .txt
                int k = atoi(nombre + len_prefijo + 1);
                if (k > k_max) {
                    k_max = k;
                    snprintf(archivo_objetivo, sizeof(archivo_objetivo), "%s/%s", carpeta, nombre);
                }
            }
        }
    }
    closedir(dir);

    if (k_max == -1) {
        fprintf(stderr, "No se encontraron archivos con prefijo '%s' en %s\n", prefijo, carpeta);
        return;
    }

    // Abrir archivo en modo append
    FILE *f = fopen(archivo_objetivo, "a");
    if (!f) {
        perror("No se pudo abrir el archivo para escribir");
        return;
    }

    fprintf(f, "\ntiempo de simulacion\t%.6f\n", tiempo);
    fclose(f);
}

double calcula_radio_giro(int N, double *x) {
    double x_cm = 0.0, y_cm = 0.0, z_cm = 0.0;
    double Rg2 = 0.0;
    int i;
    
    // 1. Calcular el centro de masa
    for (i = 0; i < N; i++) {
        x_cm += x[3*i];      // x de partícula i
        y_cm += x[3*i + 1];  // y de partícula i  
        z_cm += x[3*i + 2];  // z de partícula i
    }
    
    x_cm /= N;
    y_cm /= N;
    z_cm /= N;
    
    // 2. Calcular Rg^2
    for (i = 0; i < N; i++) {
        double dx = x[3*i] - x_cm;
        double dy = x[3*i + 1] - y_cm;
        double dz = x[3*i + 2] - z_cm;
        
        Rg2 += dx*dx + dy*dy + dz*dz;
    }
    
    Rg2 /= N;
    
    // 3. Devolver Rg (no Rg^2)
    return sqrt(Rg2);
}
/**
 * Procesa un archivo de trayectoria generado por verlet_trayectoria
 * y escribe promedios y errores de Ek, Ep, Rg, Ree en la carpeta RES_IMPORTANTES
 * Se pasa N para saber cuántas partículas hay y ubicar correctamente las columnas finales.
 */

void procesar_trayectoria(char* archivo_input, int N_start, int N, double K 
    #ifdef FIXED
        ,double F_cte
    #endif
                ) {
    FILE *file = fopen(archivo_input, "r");
    if (!file) {
        printf("No se pudo abrir el archivo %s\n", archivo_input);
        return;
    }

    char linea[32768];
    int linea_actual = 0;
    int n_datos = 0;
    int primera_impresion = 1; // flag para imprimir solo la primera línea

    double sum_Ek=0, sum_Ek2=0;
    double sum_Ep=0, sum_Ep2=0;
    double sum_Rg=0, sum_Rg2=0;
    double sum_Ree=0, sum_Ree2=0;

    while(fgets(linea, sizeof(linea), file)) {
        linea_actual++;
        if(linea_actual <= N_start) continue;

        char *ptr = linea;
        int skip_cols = 1 + 3*N + 3*N; // tiempo + 3*N posiciones + 3*N velocidades
        for(int i=0;i<skip_cols;i++){
            while(*ptr != ' ' && *ptr != '\t' && *ptr != '\0' && *ptr != '\n') ptr++;
            while((*ptr == ' ' || *ptr == '\t') && *ptr != '\0' && *ptr != '\n') ptr++;
        }

        double Ek, Ep, Rg, Ree;
        if(sscanf(ptr, "%lf %lf %*lf %lf %lf", &Ek, &Ep, &Rg, &Ree) != 4) {
            printf("Error leyendo la línea %d\n", linea_actual);
            continue;
        }

        sum_Ek += Ek; sum_Ek2 += Ek*Ek;
        sum_Ep += Ep; sum_Ep2 += Ep*Ep;
        sum_Rg += Rg; sum_Rg2 += Rg*Rg;
        sum_Ree += Ree; sum_Ree2 += Ree*Ree;
        n_datos++;
    }

    fclose(file);

    if(n_datos == 0) {
        printf("No se encontraron datos para procesar.\n");
        return;
    }

    double prom_Ek = sum_Ek / n_datos;
    double prom_Ep = sum_Ep / n_datos;
    double prom_Rg = sum_Rg / n_datos;
    double prom_Ree = sum_Ree / n_datos;

    double err_Ek = sqrt((sum_Ek2 / n_datos - prom_Ek*prom_Ek) / n_datos);
    double err_Ep = sqrt((sum_Ep2 / n_datos - prom_Ep*prom_Ep) / n_datos);
    double err_Rg = sqrt((sum_Rg2 / n_datos - prom_Rg*prom_Rg) / n_datos);
    double err_Ree = sqrt((sum_Ree2 / n_datos - prom_Ree*prom_Ree) / n_datos);

    // Crear carpeta de salida
    char carpeta[512];
    snprintf(carpeta, sizeof(carpeta), "Resultados_simulacion/%.1f/ESCALA/RES_IMPORTANTES", K);
#ifdef _WIN32
    mkdir("Resultados_simulacion");
    mkdir(carpeta);
#else
    mkdir("Resultados_simulacion", 0755);
    mkdir(carpeta, 0755);
#endif

    // Nombre de archivo de salida igual al de entrada
    char *nombre_archivo = strrchr(archivo_input, '/');
    if(!nombre_archivo) nombre_archivo = strrchr(archivo_input, '\\');
    if(nombre_archivo) nombre_archivo++;
    else nombre_archivo = archivo_input;

    char archivo_salida[512];
    snprintf(archivo_salida, sizeof(archivo_salida), "%s/%s", carpeta, nombre_archivo);

    FILE *out = fopen(archivo_salida, "w");
    if(!out) {
        printf("No se pudo crear el archivo %s\n", archivo_salida);
        return;
    }

    fprintf(out, "PROMEDIO_ENERGIA_CINETICA %.6f\n", prom_Ek);
    fprintf(out, "ERROR_ENERGIA_CINETICA %.6f\n", err_Ek);
    fprintf(out, "PROMEDIO_ENERGIA_POTENCIAL %.6f\n", prom_Ep);
    fprintf(out, "ERROR_ENERGIA_POTENCIAL %.6f\n", err_Ep);
    fprintf(out, "PROMEDIO_R_EE %.6f\n", prom_Ree);
    fprintf(out, "ERROR_R_EE %.6f\n", err_Ree);
    fprintf(out, "PROMEDIO_R_G %.6f\n", prom_Rg);
    fprintf(out, "ERROR_R_G %.6f\n", err_Rg);
    fprintf(out, "ERROR_R_G %.6f\n", err_Rg);
    fprintf(out, "N_particulas %d\n", N);
    #ifdef FIXED
    fprintf(out, "F_cte %.6f\n", F_cte);
    #endif

    fclose(out);
    printf("Archivo de resultados creado: %s\n", archivo_salida);
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>

#ifdef FIXED
#define CARPETA_FIJA "FIJOS"
#else
#define CARPETA_FIJA "ESCALA"
#endif


void procesar_trayectorias_carpeta(double K, int N_start
    #ifdef FIXED
        , double F_cte
    #endif
) {
    char carpeta[256];
    snprintf(carpeta, sizeof(carpeta), "Resultados_simulacion/%.1f/%s", K, CARPETA_FIJA);

    DIR *dir = opendir(carpeta);
    if (!dir) {
        printf("No se pudo abrir la carpeta %s\n", carpeta);
        return;
    }

    struct dirent *entry;
    while ((entry = readdir(dir)) != NULL) {
        // Solo procesar archivos que empiecen con 'V_' y terminen con '.txt'
        if (strncmp(entry->d_name, "V_", 2) == 0) {
            char *ext = strrchr(entry->d_name, '.');
            if (ext && strcmp(ext, ".txt") == 0) {
                char nombre_archivo[512];
                snprintf(nombre_archivo, sizeof(nombre_archivo), "%s/%s", carpeta, entry->d_name);
                
                // Construir ruta al archivo de parámetros
                char archivo_parametros[512];
                snprintf(archivo_parametros, sizeof(archivo_parametros), "PARAMETROS/%.1f/%s/%s", K, CARPETA_FIJA, entry->d_name);
                
                // Leer N desde el archivo de parámetros
                int N = leer_N_desde_parametros(archivo_parametros);
                if (N <= 0) {
                    printf("Error: No se pudo leer N del archivo %s\n", archivo_parametros);
                    continue;
                }
                
                #ifdef FIXED
                procesar_trayectoria(nombre_archivo, N_start, N, K, F_cte);
                #else
                procesar_trayectoria(nombre_archivo, N_start, N, K);
                #endif
            }
        }
    }

    closedir(dir);
}

// Función auxiliar para leer N desde el archivo de parámetros
int leer_N_desde_parametros(const char *archivo_parametros) {
    FILE *file = fopen(archivo_parametros, "r");
    if (!file) {
        printf("No se pudo abrir el archivo de parámetros: %s\n", archivo_parametros);
        return -1;
    }
    
    char line[256];
    int N = -1;
    
    while (fgets(line, sizeof(line), file)) {
        // Buscar la línea que contiene "N"
        if (strncmp(line, "N ", 2) == 0) {
            sscanf(line, "N %d", &N);
            break;
        }
    }
    
    fclose(file);
    return N;
}
#ifdef FIXED
#define CARPETA_IMPORTANTE "FIJOS/RES_IMPORTANTES"
#else
#define CARPETA_IMPORTANTE "ESCALA/RES_IMPORTANTES"
#endif

void generar_grafica(double K) {
    char carpeta[256];
    snprintf(carpeta, sizeof(carpeta), "Resultados_simulacion/%.1f/%s", K, CARPETA_IMPORTANTE);

    DIR *dir = opendir(carpeta);
    if (!dir) {
        printf("No se pudo abrir la carpeta %s\n", carpeta);
        return;
    }

    FILE *grafica = NULL;
    char grafica_nombre[512];
    snprintf(grafica_nombre, sizeof(grafica_nombre), "%s/grafica.txt", carpeta);
    grafica = fopen(grafica_nombre, "w");
    if (!grafica) {
        printf("No se pudo crear el archivo %s\n", grafica_nombre);
        closedir(dir);
        return;
    }

    struct dirent *entry;
    while ((entry = readdir(dir)) != NULL) {
        if (strncmp(entry->d_name, "V_", 2) == 0) {
            char *ext = strrchr(entry->d_name, '.');
            if (ext && strcmp(ext, ".txt") == 0) {
                char archivo_nombre[512];
                snprintf(archivo_nombre, sizeof(archivo_nombre), "%s/%s", carpeta, entry->d_name);

                FILE *archivo = fopen(archivo_nombre, "r");
                if (!archivo) {
                    printf("No se pudo abrir el archivo %s\n", archivo_nombre);
                    continue;
                }

                char linea[256];
                double Prom_Rg, Error_Rg;
                double Prom_Ree, Error_Ree;
                double F_cte;
                int N_particulas;

                // Leer línea por línea y extraer datos
                while (fgets(linea, sizeof(linea), archivo)) {
                    if (strstr(linea, "PROMEDIO_R_G")) sscanf(linea, "PROMEDIO_R_G %lf", &Prom_Rg);
                    else if (strstr(linea, "ERROR_R_G")) sscanf(linea, "ERROR_R_G %lf", &Error_Rg);
                    else if (strstr(linea, "PROMEDIO_R_EE")) sscanf(linea, "PROMEDIO_R_EE %lf", &Prom_Ree);
                    else if (strstr(linea, "ERROR_R_EE")) sscanf(linea, "ERROR_R_EE %lf", &Error_Ree);
                    else if (strstr(linea, "F_cte")) sscanf(linea, "F_cte %lf", &F_cte);
                    else if (strstr(linea, "N_particulas")) sscanf(linea, "N_particulas %d", &N_particulas);
                }
                fclose(archivo);

                // Escribir en grafica.txt según FIXED
#ifdef FIXED
                fprintf(grafica, "%.6f %.6f %.6f\n", F_cte, Prom_Ree, Error_Ree);
#else
                fprintf(grafica, "%d %.6f %.6f\n", N_particulas, Prom_Rg, Error_Rg);
#endif
            }
        }
    }

    fclose(grafica);
    closedir(dir);

    printf("Archivo grafica.txt creado en %s\n", carpeta);
}