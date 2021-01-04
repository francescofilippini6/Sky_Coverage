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
import healpy as hp

#...........LHASOO & ARCA...............
locationlh = locationlh = EarthLocation(lon=-100.13 * u.deg, lat=29.35 * u.deg, height=4410 * u.m)
locationAr = locationAr = EarthLocation(lon=(-15.89+180) * u.deg, lat=-36.41 * u.deg, height=-12000 * u.m)
totallh=[]
totalArca=[]
sim_points=10000
x_bins = np.linspace(0,360,181)
y_bins = np.linspace(-90,90,91)

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
    #    print(lhasoo_galactic.l.radian)
    #..............Alt-Az SkyCoord object for Lhasoo and conversion to galactic................
    ARCA_local_sky=SkyCoord(Arazz,Aralt,frame=frame_Arca, unit=u.deg)
    ARCA_galactic=ARCA_local_sky.transform_to('galactic')
    totalArca.append(ARCA_galactic)


NSIDE = 32 # this sets the side resolution of each one of the 12 basic-macro pixels
print(
    "Approximate resolution at NSIDE {} is {:.2} deg".format(
        NSIDE, hp.nside2resol(NSIDE, arcmin=True) / 60
    )
)
NPIX = hp.nside2npix(NSIDE)
print("number of pixels in the map",NPIX)



#--------------------------------------------------------------
# Setting of healpy resolution map and n. of pixels
#--------------------------------------------------------------
def cat2hpx(lon, lat, nside, radec=True):
    """
    Convert a catalogue to a HEALPix map of number counts per resolution
    element.

    Parameters
    ----------
    lon, lat : (ndarray, ndarray)
        Coordinates of the sources in degree. If radec=True, assume input is in the icrs
        coordinate system. Otherwise assume input is glon, glat

    nside : int
        HEALPix nside of the target map

    radec : bool
        Switch between R.A./Dec and glon/glat as input coordinate system.

    Return
    ------
    hpx_map : ndarray
        HEALPix map of the catalogue number counts in Galactic coordinates

    """

    npix = hp.nside2npix(nside)

    if radec:
        eq = SkyCoord(lon, lat, 'icrs', unit='deg')
        l, b = eq.galactic.l.value, eq.galactic.b.value
    else:
        l, b = lon, lat

    # conver to theta, phi
    theta = np.radians(90. - b)
    phi = np.radians(l)

    # convert to HEALPix indices
    indices = hp.ang2pix(nside, theta, phi)

    idx, counts = np.unique(indices, return_counts=True)

    # fill the fullsky map
    hpx_map = np.zeros(npix, dtype=int)
    hpx_map[idx] = counts

    return hpx_map

#--------------------------------------------------------------
#setting to zero all the map m
#--------------------------------------------------------------
#m = np.zeros(hp.nside2npix(NSIDE))

#hpx_map = cat2hpx(total[1].l.wrap_at('180d').deg, total[1].b.deg, nside=32, radec=False)
#pixel_indices = hp.ang2pix(NSIDE, theta, phi)
#m[pixel_indices] = data

#hp.mollview(hpx_map, title="Mollview image RING")
#hp.mollview(np.log10(cat2hpx(total[1].l.wrap_at('180d').deg, total[1].b.deg, nside=32, radec=False)+1))
#hp.graticule()
#plt.show()

#---------------------
#producing the galactic equator
#---------------------
equator = SkyCoord(np.linspace(0,360,256), np.zeros(256), frame='icrs', unit=u.deg)
equator_in_gal = equator.galactic


#------------------------------------------------
#making a second cycle , dividing the generation and the plotting phase
#------------------------------------------------

f2=plt.figure()
f2.suptitle("Cumulative simulation with 24 x "+str(sim_points)+"points")
cumulative=f2.add_subplot(111)
cumulative.grid(True)

sky_coverage=[]
for hou in range(24):
    arca2D=np.histogram2d(totalArca[hou].l.deg,totalArca[hou].b.deg,bins=(x_bins,y_bins))
    Lhasoo2D=np.histogram2d(totallh[hou].l.deg,totallh[hou].b.deg,bins=(x_bins,y_bins))
    #print(arca2D)
    intersect_col=[]
    intersect_row=[]
    pixelcounter=0
    for col in range(180):
        for row in range(90):
            if arca2D[0][col][row] > 0 and Lhasoo2D[0][col][row] > 0:
                #print("ciccio")
                #uso Arca per estrapolare il binning, tanto Lhasso e arca hanno lo stesso binning
                pixelcounter+=1
                intersect_col.append(Lhasoo2D[1][col])
                intersect_row.append(Lhasoo2D[2][row])
            else:
                continue
    
    #------create the object SkyCoord with the intersection
    Common_sky=SkyCoord(intersect_col,intersect_row,frame='galactic', unit=u.deg)
    #print(len(Common_sky.l.deg),len(Common_sky.b.deg))
    #-------------------------------------------------------
    #          starting the plotting functions
    #-------------------------------------------------------
    f=plt.figure(figsize=(20, 10))
    f.suptitle("simulation with"+str(sim_points)+"points")
    #-------------------------------------------------------
    sca=f.add_subplot(231)
    sca.grid(True)
    sca.plot(totalArca[hou].l.wrap_at('180d').deg,totalArca[hou].b.deg,'+',markersize=2, color='b')
    sca.plot(equator_in_gal.l.wrap_at('180d').deg, equator_in_gal.b.deg,'*',markersize=2, color='k')
    sca.set_title(str(hou)+'_hour_Arca_Plane')
    sca.set_xlim(180,-180)
    
    #-------------------------------------------------------
    sc=f.add_subplot(232)
    sc.grid(True)
    sc.plot(totallh[hou].l.wrap_at('180d').deg,totallh[hou].b.deg,'+',markersize=2, color='r')
    sc.plot(equator_in_gal.l.wrap_at('180d').deg, equator_in_gal.b.deg,'*',markersize=2, color='k')
    sc.set_title(str(hou)+'_hour_Lhasoo_Plane')
    sc.set_xlim(180,-180)
    #-------------------------------------------------------
    hh=f.add_subplot(233)
    hh.grid(True)
    hh.plot(Common_sky.l.wrap_at('180d').deg,Common_sky.b.deg,'+',markersize=2, color='g')
    hh.plot(equator_in_gal.l.wrap_at('180d').deg, equator_in_gal.b.deg,'*',markersize=2, color='k')
    hh.set_title(str(hou)+'_hour_intersection_Plane')
    hh.set_xlim(180,-180)
    #-------------------------------------------------------

    m1=f.add_subplot(234)
    hp.visufunc.mollview(np.log10(cat2hpx(totalArca[hou].l.wrap_at('180d').deg, totalArca[hou].b.deg, nside=32, radec=False)+1),hold=True,title=str(hou)+'_hour_ARCA_Moll')
    hp.graticule()
    #-------------------------------------------------------

    m2=f.add_subplot(235)
    hp.visufunc.mollview(np.log10(cat2hpx(totallh[hou].l.wrap_at('180d').deg, totallh[hou].b.deg, nside=32, radec=False)+1),hold=True,title=str(hou)+'_hour_Lhasoo_Moll')
    hp.graticule()
    #-------------------------------------------------------

    m3=f.add_subplot(236)
    hp.mollview(np.log10(cat2hpx(Common_sky.l.wrap_at('180d').deg, Common_sky.b.deg, nside=32, radec=False)+1),hold=True,title=str(hou)+'_hour_intersection_Moll')
    hp.graticule()
    #-------------------------------------------------------
    #plt.show()

    #f.canvas.manager.full_screen_toggle() # toggle fullscreen mode
    f.savefig('/Users/francescofilippini/Desktop/Sky_coverage/Common_sky_plot/'+str(hou)+'_hour_common_sky_A.png')
    print("saved figure n."+str(hou))
    #plt.show()
    #-------------------------------------------------------
    #building the sky coverage considering the pixels in commmon,
    #each is a square of 2 deg x 2 deg = 4 deg square / the solid angle in degree (41253) x 100 to have the percentage
    #-------------------------------------------------------
    intersection_solid_angle = (pixelcounter*4*100)/41253  
    sky_coverage.append(intersection_solid_angle)
    #-----------------------------------------------------------------------
    cumulative.plot(Common_sky.l.wrap_at('180d').deg,Common_sky.b.deg,'+',markersize=2, color='g')
    #-------------------------------------------------------
    
f2.savefig('/Users/francescofilippini/Desktop/Sky_coverage/Common_sky_plot/Cumulative_intersection_A.png')

#-------------------------------------------------------
#plotting the skycoverage function over time
#-------------------------------------------------------
t = np.arange(0, 24, 1)
f3=plt.figure()
ax=f3.add_subplot(111)
ax.plot(t,sky_coverage,'bo',t,sky_coverage,'r--')
ax.set_title("Sky coverage vs time")
ax.set_xlabel("Daily time (h)")
ax.set_ylabel("Sky coverage (%)")
f3.savefig('/Users/francescofilippini/Desktop/Sky_coverage/Common_sky_plot/Sky_coverage_vs_time_A.png')
    
