import os
import matplotlib.pyplot as plt
import numpy as np

# --- Uso ---
prefijo = "V"  # "V", "E-M", "R-K"
A = "0.5"      # Parámetro de carpeta
mode = "E"     # Solo estancias

def obtener_archivos(carpeta_in, prefijo):
    """Obtiene todos los archivos que contienen el prefijo en su nombre y terminan en .txt"""
    archivos = sorted([f for f in os.listdir(carpeta_in) if prefijo in f and f.endswith('.txt')])
    if not archivos:
        print(f"No se encontraron archivos con prefijo '{prefijo}' en {carpeta_in}")
    return archivos

def extraer_numero_archivo(nombre_archivo):
    """Extrae el primer número que aparece en el nombre del archivo"""
    return ''.join(filter(str.isdigit, nombre_archivo))

def plot_histograms_with_theory_estancias(prefijo, carpeta_entrada, carpeta_salida, carpeta_parametros):
    os.makedirs(carpeta_salida, exist_ok=True)
    archivos = obtener_archivos(carpeta_entrada, prefijo)

    for archivo in archivos:
        ruta_in = os.path.join(carpeta_entrada, archivo)

        x_vals, y_vals = [], []
        t_med = None

        # Leer datos del archivo
        with open(ruta_in, 'r') as f:
            for linea in f:
                linea = linea.strip()
                if linea.startswith("# Tiempo medio:"):
                    t_med = float(linea.split(":")[1].strip())
                    continue
                if linea.startswith("#") or not linea:
                    continue
                parts = linea.split()
                if len(parts) == 2:
                    try:
                        x, y = map(float, parts)
                        x_vals.append(x)
                        y_vals.append(y)
                    except ValueError:
                        continue

        if t_med is None:
            print(f"No se encontró tiempo medio en {archivo}, se omite.")
            continue

        x_vals = np.array(x_vals)
        y_vals = np.array(y_vals)
        if len(x_vals) == 0:
            print(f"No se encontraron datos en {archivo}, se omite.")
            continue

        widths = np.diff(x_vals)
        widths = np.append(widths, widths[-1]) if len(widths) > 0 else np.array([1.0])

        # Normalizar histograma
        area_histograma = np.sum(y_vals * widths)
        y_vals_normalized = y_vals / area_histograma

        # Leer parámetros (opcional, solo para etiqueta teórica)
        num_archivo = extraer_numero_archivo(archivo)
        archivo_param = f"{prefijo}_{num_archivo}.txt"
        ruta_param = os.path.join(carpeta_parametros, archivo_param)
        eta = None
        try:
            with open(ruta_param, 'r') as f:
                lineas = f.readlines()
                eta = float(lineas[3].strip().split()[1])
        except FileNotFoundError:
            print(f"Archivo de parámetros no encontrado: {archivo_param}, se omite función teórica")
        except Exception as e:
            print(f"Error leyendo parámetros: {e}")

        # Crear figura
        plt.figure(figsize=(8,5))
        plt.bar(x_vals,np.log( y_vals), width=widths, align='center', 
                color='skyblue', edgecolor='black', label='Histograma')

        # Función teórica
        xx = np.linspace(min(x_vals), max(x_vals), 500)
        if eta is not None and eta != 0:
            yy_unnormalized = np.exp(-xx / t_med) / t_med
            Z = np.trapz(yy_unnormalized, xx)
            yy = yy_unnormalized / Z
            ## plt.plot(xx, np.log(yy), 'r-', label=f'Función teórica (eta={eta})')

        plt.xlabel('Tiempo de estancia')
        plt.ylabel('Densidad de probabilidad')
        plt.title(f'Histograma tiempos de estancia - {archivo}')
        plt.legend()
        plt.tight_layout()

        ruta_out = os.path.join(carpeta_salida, archivo.replace('.txt', '.png'))
        plt.savefig(ruta_out, dpi=300)
        plt.close()
        print(f"Guardado: {ruta_out}")

# --- Configuración de rutas ---
if prefijo == "V":
    carpeta_entrada    = fr"Resultados_simulacion\DOBLE_POZO\VERLET\{A}\ESTANCIAS\HISTOGRAMA"
    carpeta_salida     = fr"Graficas\DOBLE_POZO\VERLET\{A}\ESTANCIAS"
    carpeta_parametros = fr"PARAMETROS\DOBLE_POZO\VERLET\{A}"
elif prefijo == "E-M":
    carpeta_entrada    = fr"Resultados_simulacion\DOBLE_POZO\EULER-MARUYAMA\{A}\ESTANCIAS\HISTOGRAMA"
    carpeta_salida     = fr"Graficas\DOBLE_POZO\EULER-MARUYAMA\{A}\ESTANCIAS"
    carpeta_parametros = fr"PARAMETROS\DOBLE_POZO\EULER-MARUYAMA\{A}"
elif prefijo == "R-K":
    carpeta_entrada    = fr"Resultados_simulacion\DOBLE_POZO\RUNGE-KUTTA\{A}\ESTANCIAS\HISTOGRAMA"
    carpeta_salida     = fr"Graficas\DOBLE_POZO\RUNGE-KUTTA\{A}\ESTANCIAS"
    carpeta_parametros = fr"PARAMETROS\DOBLE_POZO\RUNGE-KUTTA\{A}"
else:
    raise ValueError("FLAG no válido")

plot_histograms_with_theory_estancias(prefijo, carpeta_entrada, carpeta_salida, carpeta_parametros)

