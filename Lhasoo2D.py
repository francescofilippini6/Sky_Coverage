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

location = location = EarthLocation(lon=-100.13 * u.deg, lat=29.35 * u.deg, height=4410 * u.m)
total=[]
for t in range(24):
    obstime = Time('2020-12-27T20:00') + (1*t) *u.hour
    frame_hor = AltAz(obstime=obstime, location=location)
    azz=[]
    alt=[]
    for a in range(0,1000):
        azz.append(rnn.uniform(0,360))
        alt.append(rnn.uniform(0,90))

    local_sky=SkyCoord(azz,alt,frame=frame_hor, unit=u.deg)
    gal=local_sky.transform_to('galactic')
    total.append(gal)
#print(len(gal))
f=plt.figure()
sc=f.add_subplot(131)#, projection='aitoff')
sc.grid(True)
for a in range(24):
    sc.plot(total[a].l.wrap_at('180d').deg,total[a].b.deg,'o',markersize=1, alpha=0.1, color='b')
sc.set_title('cumulative on 24h')

sca=f.add_subplot(132)#, projection='aitoff')
sca.grid(True)
sca.plot(total[1].l.wrap_at('180d').deg,total[1].b.deg,'o',markersize=2, color='g')
sca.set_title('first_hour')



x_bins = np.linspace(-180,180,180)
y_bins = np.linspace(-90,90,90)

hh=f.add_subplot(133)
hh.hist2d(total[1].l.wrap_at('180d').deg,total[1].b.deg, bins =[x_bins,y_bins])
hh.set_title('first_hour_2D')
f.savefig('Lhasoo_2D.png')
plt.show()
