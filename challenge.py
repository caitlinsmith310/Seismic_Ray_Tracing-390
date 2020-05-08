# -*- coding: utf-8 -*-
"""
Created on Fri May  8 14:46:19 2020

@author: caitl
"""


from obspy.geodetics import locations2degrees
from obspy.geodetics import degrees2kilometers
from obspy import UTCDateTime
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
from obspy.clients.fdsn import Client
client = Client("IRIS")

lon=106.65
eqlon=159.61 
lat= 35.08
eqlat= -4.23


rad=locations2degrees(lat, lon, eqlat, eqlon)


t = UTCDateTime("2014-04-13T12:36:20")
st = client.get_waveforms("HV", "UWB", "*", "*Z", t, t + 60 * 45)
st.plot()
#st.write("japan.mseed", format="MSEED")  

#%%



def main():
    # estimate the epicentral distance. This is one of your free parameters:
    Delta = 42 # in degrees
    # estimate the origin time of the earthquake; your other free parameter:
    t0=UTCDateTime("2014-04-13T12:36:20")
    maxamp = readandplotseismogram("japan.mseed")
    computeandplottts(Delta,t0,maxamp)
    # tighten up the axes, and sh1w:
    plt.axis('tight')
    plt.show()

def readandplotseismogram(seismogram):
    '''read and plot the seismogram'''
    # read the stream
    st = read(seismogram)
    starttime = UTCDateTime("2014-04-13T12:35:00")
    st.trim(starttime,starttime+2000)
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

eqlon=162.05 
eqlat= -11.463

#SNZO
lat1= -41.309
lon1=174.704
rad1=degrees2kilometers(32)


#JCP
lat2=27.0957
lon2=142.1849
rad2=degrees2kilometers(42)

#UWB
lat3=19.425
lon3=-155.226
rad3=degrees2kilometers(50)



# a map:
plt.figure()
ax = plt.axes(projection=ccrs.Robinson(central_longitude=180))
ax.set_extent([80,270,35,-80])
ax.stock_img()


# create a polygon of points on a circle:
circle_points1 = Geodesic().circle(lon=lon1, lat=lat1, radius= 1000*rad1, n_samples=100, endpoint=False) 
geom1 = shapely.geometry.Polygon(circle_points1) 

circle_points2 = Geodesic().circle(lon=lon2, lat=lat2, radius= 1000*rad2, n_samples=100, endpoint=False) 
geom2 = shapely.geometry.Polygon(circle_points2) 

circle_points3 = Geodesic().circle(lon=lon3, lat=lat3, radius= 1000*rad3, n_samples=100, endpoint=False) 
geom3 = shapely.geometry.Polygon(circle_points3) 

# add the polygon to the map with the proper projection: 
ax.add_geometries((geom1,), crs=ccrs.Geodetic(), facecolor='None' , edgecolor="navy", linewidth=2)
plt.scatter(lon1,lat1,marker='o',color='navy',transform=ccrs.Geodetic(),label="SNZO")

ax.add_geometries((geom2,), crs=ccrs.Geodetic(), facecolor='None' , edgecolor="darkorange", linewidth=2)
plt.scatter(lon2,lat2,marker='o',color='darkorange',transform=ccrs.Geodetic(),label="JCP")

ax.add_geometries((geom3,), crs=ccrs.Geodetic(), facecolor='None' , edgecolor="forestgreen", linewidth=2)
plt.scatter(lon3,lat3,marker='o',color='forestgreen',transform=ccrs.Geodetic(),label="UWB")

#plt.plot(eqlon,eqlat,marker='o',color='crimson',transform=ccrs.Geodetic(), label="Earthquake", ms=10)
plt.grid()
plt.legend()
plt.show()





