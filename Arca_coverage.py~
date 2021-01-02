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

plt.subplot(111, projection='aitoff')
plt.grid(True)
for a in range(24):
    plt.plot(total[a].l.wrap_at('180d').radian,total[a].b.radian,'+',markersize=0.5, alpha=0.1, color='b')
#plt.plot(ra_rad, dec_rad, 'o', markersize=2, alpha=0.3)
#plt.scatter(eq.ra.wrap_at('180d').radian, eq.dec.radian)
plt.savefig('ARCA_Coverage_several.png')
plt.show()
