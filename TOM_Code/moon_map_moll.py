from skyfield.api import load
from skyfield.api import Topos
import math

from matplotlib import pyplot as plt
import numpy as np
import healpy as hp

NSIDE = 64 # this sets the side resolution of each one of the 12 basic-macro pixels
m = np.arange(hp.nside2npix(NSIDE))
print ("number of pixels:", len(m))


def make_map(): # function used to generate the pixel maps

#load the astronomical objects
  planets = load('de421.bsp')
  earth, moon = planets['earth'], planets['moon']

#set the ANTARES coordinates (Bologna is for fun)
  antares = earth + Topos('42.08 N', '6.1666666667 E')
  bologna = earth + Topos('43.1479 N', '12.1097 E')


#define the time-structure
  ts = load.timescale()

#set the time range
  t_start = ts.utc(2007,1,1,0,0,0)
  t_stop  = ts.utc(2017,1,1,0,0,0)

#getting total time to process (in adu - to be investigated)
  total_time_frac, total_time =math.modf(t_stop.tt-t_start.tt)

#defining the time steps
  step =3600/86400 #note: this is the moon-paper and analysis setting
 #step = 10 # this is for quick checks

  t_cursor = t_start # at the beginning the cursor equals the start


  


  import progressbar # this 
  from time import sleep
  bar = progressbar.ProgressBar(maxval=total_time+10, \
    widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])


  bar_index=0;
  bar.start()

# note on the angular position of objects: 
# with Skyfield, altaz() is only available  with apparent coordinates.
# Beware:  apparent() put into game optical aberration and light deflection in atmosphere. 
# In principle this is not needed.


  myPiOver2=math.pi/2.


# =========== HEALPY ===============
# HEALPY is the Python module derived from the HEALPIX library 
# for plotting the MOON  coordinates 
# in a Mollweide projection  (equal area) map.
#
# Useful reference can be found here:
#
# http://healpy.readthedocs.io/en/latest/index.html
# https://workshops.ift.uam-csic.es/uploads/charla/150/Healpix.pdf
# http://healpix.sourceforge.net/
# http://healpy.readthedocs.io/en/latest/generated/healpy.visufunc.mollview.html?highlight=mollview
#


# definition of the pixel arrays for the Mollweide projection map
  pixel_array_ring =np.ndarray(len(m))  # RING MODE
  pixel_array_nest =np.ndarray(len(m))	# NEST MODE


  time_lap = 0.0
  
  array_index=0; # an utility array

 # outfile=open("moon_positions.txt","w")

  while t_cursor.tt < t_stop.tt:
      astrometric = antares.at(t_cursor).observe(moon)
      alt, az, d = astrometric.apparent().altaz()

#      outfile.write("%d\n" %az.degrees)

      if (alt.degrees>=0) :
#      if (alt.degrees<0) :
#      if (1>0) :
          zen_rad=(myPiOver2-alt.radians)
#          pixind= hp.ang2pix(NSIDE,alt.radians,az.radians,nest=False) #inverted
          pixind= hp.ang2pix(NSIDE,zen_rad,az.radians,nest=False)
          pixel_array_ring[pixind]=pixel_array_ring[pixind]+1
          
#          pixind= hp.ang2pix(NSIDE,alt.radians,az.radians,nest=True) #inverted
          pixind= hp.ang2pix(NSIDE,zen_rad,az.radians,nest=True)
          pixel_array_nest[pixind]=pixel_array_nest[pixind]+1
        #print("pixel:", pixind)
		


      t_cursor = ts.tt(jd=(t_cursor.tt + step))
      array_index+=1
      lap_frac,lap=math.modf(t_cursor.tt-t_start.tt)
      if(time_lap<lap):
          time_lap=lap
          bar.update(time_lap)
 # 
  bar.finish()
  
# saving the pixel arrays into files
  np.savetxt("pixel_map_ring",pixel_array_ring)
  np.savetxt("pixel_map_nest",pixel_array_nest)

def read_map(): # function to read and draw the pixel maps
	
# RING	



  pixel_array_ring= np.fromfile("pixel_map_ring",sep=" ")	
  print(len(pixel_array_ring))
  hp.mollview(pixel_array_ring, nest=False, title="Mollview image RING",unit=r'Number of occurrences per pixel', min=0, max=16, norm='none',xsize=2000,fig=1)
  #hp.write_map("my_map_ring.fits",pixel_array_ring)
  hp.visufunc.graticule(dpar=10,dmer=10,local=True)
## NEST  
#  pixel_array_nest= np.fromfile("pixel_map_nest",sep=" ")	
#  print(len(pixel_array_nest))
#  hp.mollview(pixel_array_nest, nest=True, title="Mollview image NEST", fig=2)
#  hp.visufunc.graticule(dpar=10,dmer=10,local=True)

  fig1 = plt.figure(1, figsize=(5, 3.75))
#  fig1=plt.figure(1)
#  ax = fig1.add_subplot(2,1,1,projection='mollweide')

#  fig2 = plt.figure(2, figsize=(5, 3.75))

  
#fig = plt.figure(2, figsize=(5, 3.75))
#hp.mollview(mask, min=-1, max=1, title='Raw WMAP data',
#            unit=r'$\Delta$T (mK)', fig=2)
            
#  fig1.axes[0].set_title(r'Azimuth') 
#  fig2.axes[1].set_title(r'Zenith')          
#  fig1.axes[1].texts[0].set_fontsize(8)
  plt.show()


options = {0:make_map,
           1:read_map,
}
# you need to change the [<option>] below:
options[1]()



