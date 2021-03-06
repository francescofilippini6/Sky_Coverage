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


#...........LHASOO...............
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
    

    
#..............Hour by Hour merged sky coverage (24 subplots)................
figure, axes = plt.subplots(nrows=8, ncols=3, subplot_kw={'projection':  'aitoff'})
counter=0
for r in range(8):
    for c in range(3):
        axes[r,c].plot(totallh[counter].l.wrap_at('180d').radian,totallh[counter].b.radian,'o',markersize=0.5, alpha=0.1, color='r', label='LHASOO')
        axes[r,c].plot(totalArca[counter].l.wrap_at('180d').radian,totalArca[counter].b.radian,'o',markersize=0.5, alpha=0.1, color='b',  label='ARCA')
        axes[r,c].grid(True)
        counter+=1


#..............24 Hour cumulative (1 subplots)................
f2=plt.figure()
cumulative1=f2.add_subplot(111, projection='aitoff')
#cumulative2=f2.add_subplot(122, projection='aitoff')
cumulative1.grid(True)
#cumulative1.grid(True)
for a in range(24):
    cumulative1.plot(totallh[a].l.wrap_at('180d').radian,totallh[a].b.radian,'+',markersize=0.5, alpha=0.1, color='r',label='LHASOO')
    cumulative1.plot(totalArca[a].l.wrap_at('180d').radian,totalArca[a].b.radian,'+',markersize=0.5, alpha=0.1, color='b',label='ARCA')
#f2.legend(axes[1,1],label)
plt.show()
f2.savefig('Merged_Coverage.png')
figure.savefig('Merged_Coverage_HourbyHour.png')
