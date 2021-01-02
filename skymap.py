import math as m
import numpy as np
import csv
import matplotlib.pyplot as plt
import pandas as pd
import astropy.coordinates as coord
import astropy.units as u
from astropy.io import ascii
from astropy.coordinates import SkyCoord

l_ncp=122.93192
alpha_g=192.85948
delta_g=27.12825

arca_lat=36,41
arca_lon=15,88
#def equatorial_to_galactic(alpha,delta):
#    np.cos(l_ncp-l)*np.cos(b)=np.sin(delta)*np.cos(m.radians(delta_g))-np.cos(delta)*np.sin(m.radians(delta_g))*np.cos(alpha-m.radians(alpha_g))
#    np.sin(l_ncp-l)*np.cos(b)=np.cos(delta)*np.sin(alpha-m.radians(alpha_g))
#    np.sin(b)=np.sin(delta)*np.sin(delta_g)+np.cos(delta)*np.cos(m.radians(delta_g))*np.cos(alpha-m.radians(alpha_g))
    
eq = SkyCoord(np.linspace(0, 360, 128),np.zeros(128),frame='icrs', unit=u.deg)
gal = eq.galactic
#gal = SkyCoord(xarr[:], yarr[:], frame='galactic', unit=u.deg)


plt.subplot(111, projection='aitoff')
plt.grid(True)
plt.scatter(gal.l.wrap_at('180d').radian, gal.b.radian)
plt.scatter(eq.ra.wrap_at('180d').radian, eq.dec.radian)
plt.savefig('galactic_&_equatorial.png')
plt.show()
