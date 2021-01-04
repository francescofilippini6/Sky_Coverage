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
    for a in range(0,10000):
        azz.append(rnn.uniform(0,360))
        alt.append(rnn.uniform(0,90))

    local_sky=SkyCoord(azz,alt,frame=frame_hor, unit=u.deg)
    gal=local_sky.transform_to('galactic')
    total.append(gal)
#print(len(gal))

plt.subplot(111, projection='aitoff')
plt.grid(True)
for a in range(24):
    plt.plot(total[a].l.wrap_at('180d').radian,total[a].b.radian,'+',markersize=0.5, alpha=0.1, color='r')
#plt.plot(ra_rad, dec_rad, 'o', markersize=2, alpha=0.3)
#plt.scatter(eq.ra.wrap_at('180d').radian, eq.dec.radian)
plt.savefig('Lhasoo_Coverage.png')
plt.show()