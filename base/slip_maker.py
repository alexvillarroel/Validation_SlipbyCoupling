import geostochpy
import random
from scipy.interpolate import griddata
import numpy as np
def generate_scenarios(M, n_slips, sigma_L=0.18**2, alpha_W=0.17**2):
    """
    Genera escenarios de longitud (L) y ancho (W) para un evento sísmico dado.

    Parámetros:
    M (float): Magnitud del evento sísmico.
    num_scenarios (int): Número de escenarios a generar.
    sigma_L (float): Desviación estándar para log10(L).
    alpha_W (float): Desviación estándar para log10(W).
    Paper:"Scaling Relations of Earthquake Source Parameter Estimates
        with Special Focus on Subduction Environment
        by Lilian Blaser, Frank Krüger, Matthias Ohrnberger, and Frank Scherbaum"
    Retorna:
    list: Lista de tuplas (L, W) con los escenarios generados.
    """
    # Ecuaciones para log10(L) y log10(W)
    mu_L = -2.37 + 0.57 * M
    mu_W = -1.86 + 0.46 * M

    # Generar valores de log10(L) y log10(W)
    log10_L = np.random.normal(mu_L, sigma_L, n_slips)
    log10_W = np.random.normal(mu_W, alpha_W, n_slips)

    # Convertir log10(L) y log10(W) a L y W
    L = 10 ** log10_L
    W = 10 ** log10_W
    W[W > 180] = 180
    # Combinar L y W en una lista de escenarios
    scenarios = list(zip(L, W))
    return scenarios
# Adjust the scenarios inside of the mesh
#ex with meshgrid_Iquique.npz
data='meshgrid_Iquique.npz'
mesh=np.load(data)
median_lock=np.loadtxt('Chile_locking.txt')
median_lock = median_lock[~np.isnan(median_lock).any(axis=1)]
x=median_lock[:,0]
y=median_lock[:,1]
z=median_lock[:,2]
slabcoupling=np.column_stack((x, y, z))

#save the data in variables as the same name of mesh.files
X_grid=mesh['X_grid']
Y_grid=mesh['Y_grid']
dep=mesh['dep']
dip=mesh['dip']
strike=mesh['strike']
rake=mesh['rake']
length=mesh['length']
width=mesh['width']
dx=mesh['dx']
dy=mesh['dy']
nx=mesh['nx']
ny=mesh['ny']
print("length",length)
taper_coupling= griddata((slabcoupling[:,0], slabcoupling[:,1]), slabcoupling[:,2], (X_grid, Y_grid), method='linear', fill_value=0, rescale=True)

#generate the scenarios
import time
tic = time.time()
Mw=8.8
scenarios = generate_scenarios(Mw, 1000)
for i in range(len(scenarios)):
    print("Scenario",i+1,":",scenarios[i])
    length2 = scenarios[i][0]
    if length2>length:
        length2=length
    width2 = scenarios[i][1]
    northlat = np.max(Y_grid)
    southlat = np.min(Y_grid)
    new_ny=int(length2/dy)
    new_nx=int(width2/dx)
    # Selecciona aleatoriamente el índice inicial del rango
    start_index_nx = random.randint(0, nx-new_nx)
    start_index_ny = random.randint(0, ny-new_ny)
    slice_x = slice(start_index_nx, start_index_nx + new_nx)
    slice_y = slice(start_index_ny, start_index_ny + new_ny)

    # Aplicar los slices a todas las matrices de manera consistente
    X_grid_new = X_grid[slice_y, slice_x]
    Y_grid_new = Y_grid[slice_y, slice_x]
    dep_new = dep[slice_y, slice_x]
    dip_new = dip[slice_y, slice_x]
    strike_new = strike[slice_y, slice_x]
    rake_new = rake[slice_y, slice_x]
    taper_coupling_new = taper_coupling[slice_y, slice_x]
    # tuple of files of the folder base
    villarroel_taper = geostochpy.taper_except_trench_tukey(dep_new, alpha_dip=0.3, alpha_strike=0.3)
    taper_coupling_new = taper_coupling_new * villarroel_taper
    media, rigidez = geostochpy.media_slip(Mw, dx * 1000, dy * 1000, dep)
    mu = geostochpy.matriz_medias_villarroel(media, taper_coupling_new)
    C = geostochpy.matriz_covarianza_von_karman(dip_new, dep_new, X_grid_new, Y_grid_new, length2, width2)
    Slip = geostochpy.distribucion_slip_optimizada(C, mu, new_nx*new_ny-1)
    Slip, rigidez, Mo_original, Mo_deseado = geostochpy.escalar_magnitud_momento(Mw, Slip, dep_new, dy * 1000, dx * 1000, prem=True)
toc=time.time()
print("Time:",toc-tic)