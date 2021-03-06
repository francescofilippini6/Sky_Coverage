import math as m
import numpy as np
import csv
import matplotlib.pyplot as plt
import pandas as pd
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
obstime = Time('2010-01-01T20:00') + np.linspace(0, 6, 10000) * u.hour
location = EarthLocation(lon=-17.89 * u.deg, lat=28.76 * u.deg, height=2200 * u.m)
frame = AltAz(obstime=obstime, location=location)
crab = SkyCoord(ra='05h34m31.94s', dec='22d00m52.2s')
crab_altaz = crab.transform_to(frame)
print(len(crab_altaz))

plt.subplot(111, projection='aitoff')
plt.grid(True)
plt.scatter(crab_altaz.az.wrap_at('180d').radian, crab_altaz.alt.radian)
plt.savefig('crab_movement.png')
plt.show()
