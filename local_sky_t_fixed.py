import math as m
import numpy as np
import csv
import matplotlib.pyplot as plt
import pandas as pd
import random as rnn
import astropy.coordinates as coord
import astropy.units as u
from astropy.io import ascii
from astropy.coordinates import SkyCoord

from astropy.coordinates import SkyCoord, EarthLocation
from astropy import coordinates as coord
from astropy.coordinates.tests.utils import randomly_sample_sphere
from astropy.time import Time
from astropy import units as u
import numpy as np
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

times = Time.now() + np.linspace(-5, 5, 300)*u.hour
lapalma = EarthLocation.from_geocentric(5327448.9957829, -1718665.73869569, 3051566.90295403, unit='m')
aa_frame = coord.AltAz(obstime=times[:, np.newaxis], location=lapalma)
#aa_coos = coos.transform_to(aa_frame)
obstime = Time('2020-12-27T20:00') +24*u.hour# np.linspace(0, 6, 10) * u.hour
location = location = EarthLocation(lon=(-15.89+180) * u.deg, lat=-36.41 * u.deg, height=20 * u.m)
frame_hor = AltAz(obstime=obstime, location=location)
#crab = SkyCoord(ra='05h34m31.94s', dec='22d00m52.2s')
azz=[]
alt=[]
for a in range(0,300):
    azz.append(rnn.uniform(0,360))
    alt.append(rnn.uniform(0,90))
#local_sky=SkyCoord(np.linspace(0, 360, 128),np.zeros(128),frame=frame_hor, unit=u.deg)
#local_sky=SkyCoord(np.linspace(0, 360, 128),np.linspace(-90,90,128),frame=frame_hor, unit=u.deg)
local_sky=SkyCoord(azz,alt,frame=frame_hor, unit=u.deg)

gal = local_sky.transform_to('galactic')
#print(len(crab_altaz))

plt.subplot(111, projection='aitoff')
plt.grid(True)
plt.scatter(gal.l.wrap_at('180d').radian, gal.b.radian)
#plt.scatter(eq.ra.wrap_at('180d').radian, eq.dec.radian)
plt.savefig('galactic_&_equatorial.png')
plt.show()
