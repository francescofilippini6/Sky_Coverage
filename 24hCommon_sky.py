import math as m
import numpy as np
import csv
import matplotlib.pyplot as plt
import pandas as pd
import random as rnn
import astropy.coordinates as coord
import astropy.units as u
from astropy.io import ascii
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, concatenate
from astropy import coordinates as coord
from astropy.coordinates.tests.utils import randomly_sample_sphere
from astropy.time import Time
from astropy import units as u


#...........LHASOO & ARCA...............
locationlh = locationlh = EarthLocation(lon=-100.13 * u.deg, lat=29.35 * u.deg, height=4410 * u.m)
locationAr = locationAr = EarthLocation(lon=(-15.89+180) * u.deg, lat=-36.41 * u.deg, height=-12000 * u.m)
totallh=[]
totalArca=[]
sim_points=10000
x_bins = np.linspace(0,360,181)
y_bins = np.linspace(-90,90,91)

for t in range(24):
    obstime = Time('2020-12-27T20:00') + (1*t) *u.hour
    frame_Lhasoo = AltAz(obstime=obstime, location=locationlh)
    frame_Arca = AltAz(obstime=obstime, location=locationAr)
    #..............Alt-Az coordinate generation for Lhasoo sampling................
    Lhazz=[]
    Lhalt=[]
    for a in range(0,sim_points):
        Lhazz.append(rnn.uniform(0,360))
        Lhalt.append(rnn.uniform(0,90))
    #..............Alt-Az coordinate generation for ARCA sampling................
    Arazz=[]
    Aralt=[]
    for a in range(0,sim_points):
        Arazz.append(rnn.uniform(0,360))
        Aralt.append(rnn.uniform(0,90))
    #..............Alt-Az SkyCoord object for Lhasoo and conversion to galactic................
    Lhasoo_local_sky=SkyCoord(Lhazz,Lhalt,frame=frame_Lhasoo, unit=u.deg)
    lhasoo_galactic=Lhasoo_local_sky.transform_to('galactic')
    totallh.append(lhasoo_galactic)
    #    print(lhasoo_galactic.l.radian)
    #..............Alt-Az SkyCoord object for Lhasoo and conversion to galactic................
    ARCA_local_sky=SkyCoord(Arazz,Aralt,frame=frame_Arca, unit=u.deg)
    ARCA_galactic=ARCA_local_sky.transform_to('galactic')
    totalArca.append(ARCA_galactic)

#making a second cycle , dividing the generation and the plotting phase
f2=plt.figure()
f2.suptitle("Cumulative simulation with 24 x "+str(sim_points)+"points")
cumulative=f2.add_subplot(111,projection='aitoff')
cumulative.grid(True)
pixelcounter=0

for hou in range(24):
    arca2D=np.histogram2d(totalArca[hou].l.deg,totalArca[hou].b.deg,bins=(x_bins,y_bins))
    Lhasoo2D=np.histogram2d(totallh[hou].l.deg,totallh[hou].b.deg,bins=(x_bins,y_bins))
    #print(arca2D)
    intersect_col=[]
    intersect_row=[]
    for col in range(180):
        for row in range(90):
            if arca2D[0][col][row] > 0 and Lhasoo2D[0][col][row] > 0:
                #print("ciccio")
                #uso Arca per estrapolare il binning, tanto Lhasso e arca hanno lo stesso binning
                pixelcounter+=1
                intersect_col.append(Lhasoo2D[1][col])
                intersect_row.append(Lhasoo2D[2][row])
            else:
                continue
    
    #------create the object SkyCoord with the intersection
    Common_sky=SkyCoord(intersect_col,intersect_row,frame='galactic', unit=u.deg)
    print(len(Common_sky.l.deg),len(Common_sky.b.deg))
    f=plt.figure()
    f.suptitle("simulation with"+str(sim_points)+"points")
    sca=f.add_subplot(311, projection='aitoff')
    sca.grid(True)
    sca.plot(totalArca[hou].l.wrap_at('180d').radian,totalArca[hou].b.radian,'+',markersize=2, color='b')
    sca.set_title(str(hou)+'_hour_Arca')

    

    sca=f.add_subplot(312, projection='aitoff')
    sca.grid(True)
    sca.plot(totallh[hou].l.wrap_at('180d').radian,totallh[hou].b.radian,'+',markersize=2, color='r')
    sca.set_title(str(hou)+'_hour_Lhasoo')
    
    hh=f.add_subplot(313, projection='aitoff')
    hh.grid(True)
    hh.plot(Common_sky.l.wrap_at('180d').radian,Common_sky.b.radian,'+',markersize=2, color='g')
    hh.set_title(str(hou)+'_hour_intersection')
    #plt.show()
    plt.tight_layout()
    f.savefig('/Users/francescofilippini/Desktop/Sky_coverage/Common_sky_plot/'+str(hou)+'_hour_common_sky.png')
    print("saved figure n."+str(hou))
    #plt.show()

    cumulative.plot(Common_sky.l.wrap_at('180d').radian,Common_sky.b.radian,'+',markersize=2, color='g')
f2.savefig('/Users/francescofilippini/Desktop/Sky_coverage/Common_sky_plot/Cumulative_intersection.png')

    
