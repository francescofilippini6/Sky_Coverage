from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.coordinates.erfa_astrom import erfa_astrom, ErfaAstromInterpolator
from astropy.time import Time
from time import perf_counter
import numpy as np
import astropy.units as u


# array with 10000 obstimes
obstime = Time('2010-01-01T20:00') + np.linspace(0, 6, 10000) * u.hour
location = location = EarthLocation(lon=-17.89 * u.deg, lat=28.76 * u.deg, height=2200 * u.m)
frame = AltAz(obstime=obstime, location=location)
crab = SkyCoord(ra='05h34m31.94s', dec='22d00m52.2s')

# transform with default transformation and print duration
t0 = perf_counter()
crab_altaz = crab.transform_to(frame)  
print(f'Transformation took {perf_counter() - t0:.2f} s')  


# transform with interpolating astrometric values
t0 = perf_counter()
with erfa_astrom.set(ErfaAstromInterpolator(300 * u.s)):
    crab_altaz_interpolated = crab.transform_to(frame)
print(f'Transformation took {perf_counter() - t0:.2f} s')  


err = crab_altaz.separation(crab_altaz_interpolated)
print(f'Mean error of interpolation: {err.to(u.microarcsecond).mean():.4f}')  


# To set erfa_astrom for a whole session, use it without context manager:
erfa_astrom.set(ErfaAstromInterpolator(300 * u.s)) 
