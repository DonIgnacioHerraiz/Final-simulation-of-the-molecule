import matplotlib.pyplot as plt
import numpy as np
import os

def plot_grafica_txt(ruta_archivo, guardar_imagen=False, titulo=None, b=1.0, usar_log=False, ruta_guardado_personalizada=None):
    """
    Lee un archivo .txt con formato (N, Rg, error) y crea una gráfica comparando
    con la ley teórica Rg = b * sqrt((N^2 - 1) / (6N)).

    Parámetros:
    - ruta_archivo: ruta al archivo .txt
    - guardar_imagen: si True, guarda la imagen en el mismo directorio
    - titulo: título personalizado para la gráfica
    - b: longitud media de enlace (por defecto 1.0)
    - usar_log: si True, usa escala log-log
    - ruta_guardado_personalizada: ruta personalizada para guardar la imagen
    """
    try:
        # Leer datos
        datos = np.loadtxt(ruta_archivo)
        if datos.ndim == 1:
            N = np.array([datos[0]])
            Rg = np.array([datos[1]])
            error = np.array([datos[2]]) if len(datos) >= 3 else np.zeros(1)
        else:
            N = datos[:, 0]
            Rg = datos[:, 1]
            error = datos[:, 2] if datos.shape[1] >= 3 else np.zeros_like(N)

        # Curva teórica
        Rg_teo = b * np.sqrt((N**2 - 1) / (6 * N))

        # Crear figura
        plt.figure(figsize=(9, 6))

        # Datos simulados con barras de error
        plt.errorbar(
            N, Rg, yerr=error, fmt='o', capsize=5, elinewidth=1.5,
            markersize=6, label='Simulación', color='royalblue'
        )

        # Curva teórica
        plt.plot(
            N, Rg_teo, '--', color='darkorange', linewidth=2,
            label=r'Teoría $R_g = b \sqrt{\frac{N^2-1}{6N}}$'
        )

        # Personalización de la gráfica
        if titulo is None:
            titulo = "Radio de giro vs número de monómeros"

        plt.title(titulo, fontsize=14)
        plt.xlabel("Número de monómeros (N)", fontsize=12)
        plt.ylabel(r"$R_g$ (unidades reducidas)", fontsize=12)

        if usar_log:
            plt.xscale("log")
            plt.yscale("log")
            plt.grid(True, which="both", ls="--", alpha=0.3)
        else:
            plt.grid(True, ls="--", alpha=0.4)

        plt.legend(fontsize=11)
        plt.tight_layout()

        # Mostrar
        plt.show()

        # Guardar imagen
        if guardar_imagen:
            if ruta_guardado_personalizada:
                # Usar ruta personalizada si se proporciona
                ruta_guardado = ruta_guardado_personalizada
            else:
                # Guardar en el mismo directorio del archivo original
                directorio = os.path.dirname(ruta_archivo)
                nombre_base = os.path.splitext(os.path.basename(ruta_archivo))[0]
                ruta_guardado = os.path.join(directorio, f"{nombre_base}_plot.png")
            
            # Crear directorio si no existe
            os.makedirs(os.path.dirname(ruta_guardado), exist_ok=True)
            
            plt.savefig(ruta_guardado, dpi=300, bbox_inches="tight")
            print(f"✅ Gráfica guardada en: {ruta_guardado}")

    except Exception as e:
        print(f"⚠️ Error al procesar el archivo: {e}")


# Ejemplo de uso - GUARDANDO EN Graficas\100.0
ruta = r"Resultados_simulacion\100.0\ESCALA\RES_IMPORTANTES\grafica.txt"

# Crear la ruta personalizada para guardar en Graficas\100.0
ruta_guardado = r"Graficas\100.0\grafica_Rg_vs_N.png"

plot_grafica_txt(
    ruta, 
    guardar_imagen=True, 
    titulo="Dependencia de Rg con N (K=100, L₀=1)", 
    b=1.0,
    ruta_guardado_personalizada=ruta_guardado  # <-- NUEVO PARÁMETRO
)
