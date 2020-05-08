# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 20:57:07 2020

@author: caitl
"""
from obspy.geodetics import locations2degrees
from obspy.geodetics import degrees2kilometers
from obspy.taup import plot_ray_paths
import matplotlib.pyplot as plt
from obspy.taup.tau import plot_ray_paths
from obspy.taup import plot_travel_times
import matplotlib.pyplot as plt
from obspy import read
from obspy.clients.arclink import Client
from obspy import UTCDateTime
import numpy as np
import matplotlib.pyplot as plt
from obspy import read
from obspy.taup import TauPyModel
import os
import cartopy.crs as ccrs
from cartopy.geodesic import Geodesic
import shapely
from obspy.core import read
from obspy.signal.trigger import ar_pick

directory =r'C:\Users\caitl\Documents\390\seismic'
os.chdir(directory)
#%%

plot_ray_paths(source_depth=100,phase_list=['S'],min_degrees=0, max_degrees=360,npoints=25)

fig, ax = plt.subplots(subplot_kw=dict(polar=True))
ax = plot_ray_paths(source_depth=100, ax=ax, fig=fig, phase_list=['P', 'PKP'],
                    npoints=1000, legend=True, verbose=True)

#%%

plot_ray_paths(10, min_degrees=0, max_degrees=35, npoints=10, plot_type='spherical', 
               phase_list=['S', "SKS"], model='iasp91', plot_all=True, legend=True, 
               label_arrivals=False, verbose=False, fig=None, show=True, ax=None)


#toggle between cartsian and spherical
#%%
fig,ax=plt.subplots()
ax=plot_travel_times(source_depth=33, min_degrees=15, max_degrees=40, phase_list=["S","P"], plot_all=True,npoints=20, ax=ax, fig=fig, verbose=False)
ax.grid()
#%%

st = read("quake.sac_.gz")
starttime = UTCDateTime("2014-04-13T12:40:00")
st.trim(starttime,starttime+1500)
st.plot()

#%%

def main():
    # estimate the epicentral distance. This is one of your free parameters:
    Delta = 32 # in degrees
    # estimate the origin time of the earthquake; your other free parameter:
    t0=UTCDateTime("2014-04-13T12:36:20")
    maxamp = readandplotseismogram("quake.sac_.gz")
    computeandplottts(Delta,t0,maxamp)
    # tighten up the axes, and sh1w:
    plt.axis('tight')
    plt.show()

def readandplotseismogram(seismogram):
    '''read and plot the seismogram'''
    # read the stream
    st = read(seismogram)
    starttime = UTCDateTime("2014-04-13T12:35:00")
    st.trim(starttime,starttime+1500)
    # the trace in the stream is
    tr = st[0]
    # Make a time array and an amps array:
    t1= tr.stats.starttime.timestamp
    t2= tr.stats.endtime.timestamp
    npts = tr.stats.npts
    times = np.linspace(t1,t2,npts)
    amps = tr.data
    maxamp = max(amps)
    #Plot the seismogram against a grid
    plt.figure(figsize=(12,3))
    plt.plot(times,amps,color='blue')
    #Converting tick marks to actual times:
    locs, labels = plt.xticks()
    new_label=[]

    for loc in locs:
        new_label.append(UTCDateTime(loc).strftime("%H-%M-%S"))
        plt.xticks(locs,new_label,rotation=10)
        plt.xlabel('Time (HH-MM-SS UTC)')
        plt.ylabel('Displacement')
        plt.grid()
    return maxamp


def computeandplottts(Delta,t0,maxamp):
    # compute arrival times based on a epicentral distance (and depth, but
    # I fixed depth to a constant value)
    model = TauPyModel(model="iasp91")     
    arrivals = model.get_travel_times(distance_in_degree=Delta,source_depth_in_km=32)
    #Construct vertical lines to show arrival times of predicted seismic waves
    for arrival in arrivals:
        dummy = t0+ arrival.time
        if arrival.name == "P":
            plt.vlines(dummy,-maxamp/2,maxamp/2)
            plt.text(dummy,maxamp/2+0.05*maxamp,arrival.name)
        if arrival.name == "PP":
            plt.vlines(dummy,-maxamp/2,maxamp/2)
            plt.text(dummy,maxamp/2+0.05*maxamp,arrival.name)
        if arrival.name == "S":
            plt.vlines(dummy,-maxamp/2,maxamp/2)
            plt.text(dummy,maxamp/2+0.05*maxamp,arrival.name)
        if arrival.name == "SS":
            plt.vlines(dummy,-maxamp/2,maxamp/2)
            plt.text(dummy,maxamp/2+0.05*maxamp,arrival.name)
    return()

# this will actually run the code if called stand-alone:
if __name__ == '__main__':
    main()


#%%


# compute points on a circle on the Earth, centered on (lon,lat): 
lon=174.704
eqlon=159.61 
lat= -41.309
eqlat= -4.23
# a map:
plt.figure()
ax = plt.axes(projection=ccrs.Robinson(central_longitude=180))
ax.set_extent([120,230,5,-80])
ax.stock_img()

# create a polygon of points on a circle:
circle_points = Geodesic().circle(lon=lon, lat=lat, radius= 1000*rad, n_samples=100, endpoint=False) 
geom = shapely.geometry.Polygon(circle_points) 

# add the polygon to the map with the proper projection: 
ax.add_geometries((geom,), crs=ccrs.Geodetic(), facecolor='None' , edgecolor="navy", linewidth=2)
plt.scatter(lon,lat,marker='o',color='darkorange',transform=ccrs.Geodetic(),label="SNZO")
#plt.plot(eqlon,eqlat,marker='o',color='crimson',transform=ccrs.Geodetic(), label="Earthquake", ms=10)
plt.grid()
plt.legend()
plt.show()



