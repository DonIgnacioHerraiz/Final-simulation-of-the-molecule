import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def animar_molecula_con_lineas(archivo, N, dt):
    """
    Animación de la molécula con líneas entre partículas y tiempo.

    :param archivo: ruta del archivo de simulación
    :param N: número de partículas
    :param dt: paso temporal entre frames (segundos)
    """

    # --- Leer archivo ---
    with open(archivo, 'r') as f:
        cabecera = f.readline()  # saltar cabecera
        data = f.readlines()

    # --- Extraer posiciones ---
    pos = []
    for linea in data:
        tokens = linea.strip().split()
        if len(tokens) < 1 + 3*N:
            continue
        x = [float(tokens[i]) for i in range(1, 3*N+1)]
        pos.append(x)
    pos = np.array(pos)  # shape: (frames, 3*N)

    # --- Configurar figura ---
    fig = plt.figure(figsize=(7,7))
    ax = plt.axes(projection='3d')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # Limites dinámicos según posiciones
    all_x = pos[:, 0::3].flatten()
    all_y = pos[:, 1::3].flatten()
    all_z = pos[:, 2::3].flatten()
    margin = 0.5
    ax.set_xlim(all_x.min()-margin, all_x.max()+margin)
    ax.set_ylim(all_y.min()-margin, all_y.max()+margin)
    ax.set_zlim(all_z.min()-margin, all_z.max()+margin)

    # --- Inicializar puntos y líneas ---
    scatter = ax.scatter(pos[0,0::3], pos[0,1::3], pos[0,2::3], c='blue', s=50)
    line, = ax.plot(pos[0,0::3], pos[0,1::3], pos[0,2::3], c='red')
    tiempo_text = ax.text2D(0.05, 0.95, "", transform=ax.transAxes)

    # --- Función de actualización ---
    def update(frame):
        x = pos[frame, 0::3]
        y = pos[frame, 1::3]
        z = pos[frame, 2::3]

        scatter._offsets3d = (x, y, z)
        line.set_data(x, y)
        line.set_3d_properties(z)

        tiempo_text.set_text(f"Tiempo = {frame*dt:.2f}")
        return scatter, line, tiempo_text

    # --- Crear animación ---
    anim = FuncAnimation(fig, update, frames=len(pos), interval=50, blit=False)

    plt.show()


# --- Uso ---
archivo = r"Resultados_simulacion\100.0\ESCALA\V_0.txt"
N = 16           # número de partículas
dt = 10        # tiempo entre frames (coherente con counter >= 0.1 en tu simulación)
animar_molecula_con_lineas(archivo, N, dt)
