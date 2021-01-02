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



#times = Time.now() + np.linspace(-5, 5, 300)*u.hour
#lapalma = EarthLocation.from_geocentric(5327448.9957829, -1718665.73869569, 3051566.90295403, unit='m')
#aa_frame = coord.AltAz(obstime=times[:, np.newaxis], location=lapalma)
#aa_coos = coos.transform_to(aa_frame)
#crab = SkyCoord(ra='05h34m31.94s', dec='22d00m52.2s')

location = location = EarthLocation(lon=(-15.89+180) * u.deg, lat=-36.41 * u.deg, height=-12000 * u.m)
total=[]
for t in range(24):
    obstime = Time('2020-12-27T20:00') + (1*t) *u.hour
    frame_hor = AltAz(obstime=obstime, location=location)
    azz=[]
    alt=[]
    for a in range(0,10000):
        azz.append(rnn.uniform(0,360))
        alt.append(rnn.uniform(0,90))

    local_sky=SkyCoord(azz,alt,frame=frame_hor, unit=u.deg)
    gal=local_sky.transform_to('galactic')
    total.append(gal)
#print(len(gal))
f=plt.figure()
#sc=f.add_subplot(131)#, projection='aitoff')
#sc.grid(True)
#for a in range(24):
#    sc.plot(total[a].l.wrap_at('180d').deg,total[a].b.deg,'o',markersize=1, alpha=0.1, color='b')
#sc.set_title('cumulative on 24h')

#sca=f.add_subplot(132)
#sca.grid(True)
#sca.plot(total[1].l.wrap_at('180d').deg,total[1].b.deg,'o',markersize=2, color='g')
#sca.set_title('first_hour')
cc=[]
dd=[]
for a in range(24):
    cc.append(total[a].l)
    dd.append(total[a].b)
print("Len di cc",len(cc))
print("Len di dd",len(dd))

coords = SkyCoord(cc, dd, frame='galactic',unit='degree')
color_map = plt.cm.Spectral_r
f.add_subplot(111, projection='aitoff')
#for a in range(24):
image=plt.hexbin(coords.l.wrap_at('180d'),coords.b, cmap=color_map,gridsize=1000, mincnt=1, bins='log')

plt.xlabel('Galactic latitude (l)')
plt.ylabel('Galactic longitude (b)')
plt.grid(True)
plt.colorbar(image, spacing='uniform', extend='max')
plt.show()    
#x_bins = np.linspace(-180,180,180)
#x_bins = np.linspace(-np.pi,np.pi,180)
#y_bins = np.linspace(-np.pi/2,np.pi/2,180) 
#y_bins = np.linspace(-90,90,90)

#hh=f.add_subplot(133, projection='aitoff')
#hh.hist2d(total[1].l.deg,total[1].b.deg, bins =[x_bins,y_bins])
#hh.set_title('first_hour_2D')
#f.savefig('ARCA_Coverage_histo.png')
plt.show()
