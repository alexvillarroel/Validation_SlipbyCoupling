import geostochpy
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import pandas as pd
# read csv file
data = pd.read_csv('trench_sam.csv')
# get the data of the trench
# Extraer la columna 1 como longitud
longitud = data.iloc[:, 0]
# Extraer la columna 2 como latitud
latitud = data.iloc[:, 1]
slabdep,slabdip,slabstrike,slabrake=geostochpy.load_files_slab2(zone='south_america',rake=True)
# interp slabstrike
slabstrike[:,0] += 360
# interp slabstrike
strike=slabstrike[:,2]
strike[strike>180]-=360
slab_fosa=griddata((slabstrike[:,0],slabstrike[:,1]),strike,(longitud+1,latitud),method='linear')
data=np.column_stack((longitud-360,latitud,slab_fosa))
# sort data by latitude, from minor to mayor
data=data[data[:,1].argsort()]
# save as txt longitud,latitud,slab_fosa
np.savetxt('trench_fosa.txt',data,fmt='%f %f %f')