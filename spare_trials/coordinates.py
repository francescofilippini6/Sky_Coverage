from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation

c = SkyCoord(10, 20, unit="deg")
c1 = SkyCoord([1, 2, 3], [-30, 45, 8], frame="icrs", unit="deg")
c2 = SkyCoord(frame="galactic", l="1h12m43.2s", b="+1d12m43s")
#t = Time(56370.5, format='mjd', scale='utc')
loc = EarthLocation('149d33m00.5s','-30d18m46.385s',236.87*u.m)
sc = SkyCoord(1*u.deg, 2*u.deg)

axa = (projection="mollweide")
axa.scatter(np.array(np.pi/2,np.pi/2))
plt.show()
