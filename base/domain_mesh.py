import geostochpy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
#
# Let's start defining latitudes and longitudes of the domain of each segment
# for this, we will use the coordinates of data_complilation.xlsx
#
data_segments = pd.read_excel('data_compilation.xlsx')
data_iquique = [data_segments.query('Segment ==  "Iquique Segment" ')['Start'].max(), data_segments.query('Segment ==  "Iquique Segment" ')['End'].min()]
data_vallenar = [data_segments.query('Segment ==  "Vallenar Segment" ')['Start'].max(), data_segments.query('Segment ==  "Vallenar Segment" ')['End'].min()]
data_valparaiso = [data_segments.query('Segment ==  "Valparaiso Segment" ')['Start'].max(), data_segments.query('Segment ==  "Valparaiso Segment" ')['End'].min()]
data_concepcion = [data_segments.query('Segment ==  "Concepcion Segment" ')['Start'].max(), data_segments.query('Segment ==  "Concepcion Segment" ')['End'].min()]
dict={'Iquique':data_iquique,'Vallenar':data_vallenar,'Valparaiso':data_valparaiso,'Concepcion':data_concepcion}
# make a meshgrid with superior and inferior limits of each segment
#
# first, we need load the trench coordinates and slab geometry of the segment
# for this, we can use the geostochpy library
lonsfosa, latsfosa,strikefosa  = np.loadtxt('trench_fosa.txt',unpack=True)
slabdep,slabdip,slabstrike,slabrake=geostochpy.load_files_slab2(zone='south_america',rake=True)
#Next, we need to obtain the depths at each subfault. We can achieve this by
#  interpolating the Slab2 data using geostochpy.interp_slabtofault.
# lon,lat,lon_flat,lat_flat=geostochpy.make_fault_alongstriketrench(lonsfosa, latsfosa,strikefosa,north, nx, ny, width, length)

for key,value in dict.items():
    print(key,value)
    #we need compute the length of the segment, that its the absolute diference between value[1] and value[0] and translate to km
    length = float(int(abs(value[0]-value[1])*111.111))
    # width, will be of 180 km
    width = 180
    # nx, ny, will be the number of subfaults in the x and y direction
    # in all the segments we will use a number of subfault of form that each subfault will have a dx=dy=10 km
    nx = int(width/10)
    ny = int(length/10)

    lon1,lat1,lon_flat,lat_flat=geostochpy.make_fault_alongstriketrench(lonsfosa, latsfosa,strikefosa,value[0], nx, ny, width, length)
    X_grid,Y_grid,dep,dip,strike,rake=geostochpy.interp_slabtofault(lon_flat,lat_flat,nx,ny,slabdep,slabdip,slabstrike,slabrake)
    #let's save the data in a npz file
    np.savez('meshgrid_'+key+'.npz',X_grid=lon1,Y_grid=lat1,dep=dep,dip=dip,strike=strike,rake=rake,length=length,width=width,nx=nx,ny=ny,dx=width/nx,dy=length/ny)
 