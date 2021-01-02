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
sim_points=1000
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
    #..............Alt-Az SkyCoord object for Lhasoo and conversion to galactic................
    ARCA_local_sky=SkyCoord(Arazz,Aralt,frame=frame_Arca, unit=u.deg)
    ARCA_galactic=ARCA_local_sky.transform_to('galactic')
    totalArca.append(ARCA_galactic)
    

x_bins = np.linspace(-180,180,181)
y_bins = np.linspace(-90,90,91)
arca2D=np.histogram2d(totalArca[1].l.wrap_at('180d').deg,totalArca[1].b.deg,bins=(x_bins,y_bins))
Lhasoo2D=np.histogram2d(totallh[1].l.wrap_at('180d').deg,totallh[1].b.deg,bins=(x_bins,y_bins))
print(arca2D)
intersect_col=[]
intersect_row=[]
for col in range(180):
    for row in range(90):
        if arca2D[0][col][row] > 0 and Lhasoo2D[0][col][row] > 0:
            print("ciccio")
            intersect_col.append(col)
            intersect_row.append(row)
        else:
            continue


#------create the object SkyCoord with the intersection
Common_sky=SkyCoord(intersect_col,intersect_row,frame='galactic', unit=u.deg)
f=plt.figure()

sca=f.add_subplot(131, projection='aitoff')
sca.grid(True)
sca.plot(totalArca[1].l.wrap_at('180d').radian,totalArca[1].b.radian,'o',markersize=2, color='b')
sca.set_title('first_hour_Arca')

sca=f.add_subplot(132, projection='aitoff')
sca.grid(True)
sca.plot(totallh[1].l.wrap_at('180d').radian,totallh[1].b.radian,'o',markersize=2, color='r')
sca.set_title('first_hour_Lhasoo')

hh=f.add_subplot(133, projection='aitoff')
#hh.plot(intersect_col.radian,intersect_row.radian,'o',markersize=2, color='g')
hh.plot(Common_sky.l.wrap_at('180d').radian,Common_sky.b.radian,'o',markersize=2, color='g')
#hh.set_xlim(-180, 180)
#hh.set_ylim(-90, 90)
hh.set_title('first_hour_intersection')
plt.show()
f.savefig('intersection.png')
#intersection=np.histogram2d(intersect_col, intersect_row,bins=(x_bins,y_bins))
#print(intersection[0])

#plt.figure()
#myextent  =[x_bins[0],x_bins[-1],y_bins[0],y_bins[-1]]
#plt.imshow(intersection,origin='low',extent=myextent,interpolation='nearest',aspect='auto')
#plt.colorbar()

#plt.imshow(intersection,cmap="viridis",aspect='auto',origin='lower')
#plt.imshow(intersection, interpolation='none', origin='low',extent=[x_bins[0], x_bins[-1], y_bins[0], y_bins[-1]])

