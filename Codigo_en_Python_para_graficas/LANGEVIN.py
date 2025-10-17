import numpy as np
import matplotlib.pyplot as plt
import os

# --- Configuración ---
# Ruta al archivo de datos.
# Asegúrate de que esta ruta sea correcta desde donde ejecutes el script.
file_path = r'Resultados_simulacion/100.0/FIJOS/RES_IMPORTANTES/grafica.txt'
N = 4  # Número de partículas (monómeros)

# Ruta personalizada para guardar la gráfica
ruta_guardado = r'Graficas/100.0/ree_vs_fuerza.png'


# --- Función para el modelo teórico ---
def theoretical_ree(f):
    """Calcula el valor teórico de R_ee para una fuerza dada."""
    # Se añade un valor muy pequeño a f para evitar la división por cero si f=0.
    return (N-1) * (np.cosh(f) / np.sinh(f) - 1 / (f + 1e-9))


# --- Lógica Principal ---
try:
    # 1. Cargar los datos del archivo de texto.
    #    La opción unpack=True asigna cada columna a una variable diferente.
    f_cte, r_ee_exp, error_r_ee = np.loadtxt(file_path, unpack=True)

    # 2. Preparar los datos para la curva teórica.
    #    Creamos un array de fuerzas mucho más denso (500 puntos) para que la
    #    curva se dibuje de forma suave y continua.
    #    El rango va desde el mínimo al máximo de las fuerzas experimentales.
    f_teorico = np.linspace(f_cte.min(), f_cte.max(), 500)
    
    #    Calculamos los valores de R_ee para la curva teórica.
    r_ee_teorico = theoretical_ree(f_teorico)

    # 3. Crear la gráfica.
    plt.style.use('seaborn-v0_8-whitegrid') # Estilo visual agradable
    fig, ax = plt.subplots(figsize=(10, 6))

    #    Dibujar los datos experimentales con sus barras de error.
    ax.errorbar(f_cte, r_ee_exp, yerr=error_r_ee,
                fmt='o',             # 'o' para marcadores circulares
                capsize=5,           # Añade "tapones" a las barras de error
                color='royalblue',
                label='Datos de Simulación')

    #    Dibujar la curva teórica.
    ax.plot(f_teorico, r_ee_teorico,
            color='red',
            linestyle='-',       # Línea sólida
            linewidth=2,
            label='Curva Teórica: $3(\cosh(F) - 1/F)$')

    # 4. Añadir títulos y etiquetas para mayor claridad.
    ax.set_title('R_ee frente a fuerza constante aplicada', fontsize=16)
    ax.set_xlabel('Fuerza Constante Aplicada ($F_{cte}$)', fontsize=12)
    ax.set_ylabel('Distancia Extremo a Extremo ($R_{ee}$)', fontsize=12)
    
    #    Añadir la leyenda para identificar cada serie.
    ax.legend()

    # 5. Crear el directorio si no existe
    os.makedirs(os.path.dirname(ruta_guardado), exist_ok=True)
    
    # 6. Guardar la gráfica antes de mostrarla
    plt.savefig(ruta_guardado, dpi=300, bbox_inches='tight')
    print(f"✅ Gráfica guardada en: {ruta_guardado}")
    
    # 7. Mostrar la gráfica.
    plt.show()

except FileNotFoundError:
    print(f"❌ Error: No se pudo encontrar el archivo en la ruta especificada.")
    print(f"   Ruta buscada: '{file_path}'")
    print("   Por favor, verifica que la ruta es correcta y que el archivo existe.")
except Exception as e:
    print(f"Ha ocurrido un error inesperado: {e}")