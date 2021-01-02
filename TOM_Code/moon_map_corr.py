from skyfield.api import load
from skyfield.api import Topos
import math

#load the astronomical objects
planets = load('de421.bsp')
earth, moon = planets['earth'], planets['moon']

#set the ANTARES coordinates (Bologna is for fun)
antares = earth + Topos('42.08 N', '6.1666666667 E')
bologna = earth + Topos('43.1479 N', '12.1097 E')

##set the time interval (not used) ---> to be deleted
#start_year=2007;
#stop_year=2015;
#months = 12*(stop_year-start_year)



#define the time-structure
ts = load.timescale()

#set the time range
t_start = ts.utc(2009,1,1,0,0,0)
t_stop  = ts.utc(2015,1,1,0,0,0)

#getting total time to process (in adu - to be investigated)
total_time_frac, total_time =math.modf(t_stop.tt-t_start.tt)




array_index=0;
step =3600/86400 #note: this is the moon-paper and analysis setting
#step = 10
t_cursor = t_start



from array import *
zen_array=array('d')
azi_array=array('d')
azicorr_array=array('d')

xstereo_array=array('d')
ystereo_array=array('d')

xlambert_array=array('d')
ylambert_array=array('d')



import progressbar # this 
from time import sleep
bar = progressbar.ProgressBar(maxval=total_time+10, \
    widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])


bar_index=0;
bar.start()

#note on the angular position of objects: 
# with Skyfield, altaz() is only available  with apparent coordinates.
# Beware:  apparent() put into game optical aberration and light deflection in atmosphere. 
# In principle this is not needed.


myPiOver2=math.acos(-1.)/2

xtmp_r=[1000000,0]
ytmp_r=[1000000,0]

time_lap = 0.0
while t_cursor.tt < t_stop.tt:
    astrometric = antares.at(t_cursor).observe(moon)
    alt, az, d = astrometric.apparent().altaz()


    if (alt.degrees>=0) :
        zen = 90-alt.degrees
        zen_rad_over_2=(myPiOver2-alt.radians)/2
        zen_array.insert(array_index,zen)
        azi_array.insert(array_index,az.degrees)
        azicorr_array.insert(array_index,az.degrees*math.cos(alt.radians))
# Steroscopic
#        xtmp=math.cos(az.radians)/math.tan(zen_rad_over_2)
#        ytmp=math.sin(az.radians)/math.tan(zen_rad_over_2)
# Lambert
        xtmp=math.cos(az.radians)*2*math.cos(zen_rad_over_2)
        ytmp=math.sin(az.radians)*2*math.cos(zen_rad_over_2)


        
        if (xtmp<xtmp_r[0]):
            xtmp_r[0]=xtmp
			
        if (xtmp>xtmp_r[1]):
           xtmp_r[1]=xtmp

        if(ytmp<ytmp_r[0]):
           ytmp_r[0]=ytmp
			
        if(ytmp>ytmp_r[1]):
           ytmp_r[1]=ytmp
			
			
			
        xstereo_array.insert(array_index, xtmp)
        ystereo_array.insert(array_index, ytmp)
#        print (array_index,xtmp,ytmp)
#        print(array_index)
#    lzen = zen
#    laz = az.degrees
#    print (t_cursor,lzen,laz)



    t_cursor = ts.tt(jd=(t_cursor.tt + step))
    array_index+=1
    lap_frac,lap=math.modf(t_cursor.tt-t_start.tt)
    if(time_lap<lap):
        time_lap=lap
        bar.update(time_lap)
 # 
bar.finish()

print('x in: [',xtmp_r[0],',',xtmp_r[1],']')
print('y in: [',ytmp_r[0],',',ytmp_r[1],']')



import matplotlib.pyplot as plt
import matplotlib.dates as md
from matplotlib.colors import LogNorm
from matplotlib.ticker import NullFormatter

nullfmt = NullFormatter()

#plt.plot(azi_array,zen_array,'ro')
#plt.ylabel('zenith (deg)')
#plt.xlabel('azimut (deg)')
#plt.show()

left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
bottom_h = left_h = left + width + 0.02


rect_H2D = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.2]
rect_histy = [left_h, bottom, 0.2, height]



# start with a rectangular Figure
plt.figure(1, figsize=(8, 8))

# define the axes for the three histos: one 2D histo + two 1D projection histos
axH2D = plt.axes(rect_H2D)
axHistx = plt.axes(rect_histx)
axHisty = plt.axes(rect_histy)

# no labels on projection histos
axHistx.xaxis.set_major_formatter(nullfmt)
axHisty.yaxis.set_major_formatter(nullfmt)


import numpy as np
nbins=80 # this is the moon-paper settings for both nbinx and nbiny
nbinx=nbins #13 # setting useful for MC studies
nbiny=nbins #4 # setting useful for MC studies



# Histo ranges

#rangex=[0,360]
rangex=[40,360]
#rangey=[0,90]
rangey=[10,90]


#creation of the 2D histo with adding the entries given the arrays for the angles

myX_array=azicorr_array
myY_array=zen_array


#myX_array=xstereo_array
#myY_array=ystereo_array


(n, mybinx,mybiny, patches) = axH2D.hist2d(myX_array, myY_array, bins=(nbinx,nbiny), range=[rangex,rangey], norm=LogNorm())

# dump a file with the bin-by-bin histo entries 
file_grid=open("moon_grid","w")
for i in reversed(range(nbiny)):
 for j in reversed(range(nbinx)):
     file_grid.write("%d\n" %n[j,i])
#     if n[j,i]>0 :
#	     print("(",mybinx[j],",",mybiny[i],")\t->\t",n[j,i],"\n", end='', flush=True)
#print("")
# file_grid.write("\n")
#
#file_grid.close()


axH2D.set_ylabel('zenith (deg)')
axH2D.set_xlabel('azimuth * cos(altitude) (deg)')

#axH2D.set_ylabel('Y stereo-zenith (deg)')
#axH2D.set_xlabel('X stereo-azimuth (deg)')

#creation of the two 1D histos 
axHistx.hist(myX_array, bins=nbinx, range=rangex)
axHisty.hist(myY_array, bins=nbiny, range=rangey, orientation='horizontal')

#axHistx.set_xlim(axH2D.get_xlim())
axHistx.set_xlim(rangex)
#axHisty.set_ylim(axH2D.get_ylim())
axHisty.set_ylim(rangey)


print("Range of H2D: ",axH2D.get_xlim(),axH2D.get_ylim())
print("Range of H1D-x: ",axHistx.get_xlim())
print("Range of H1D-y: ",axHisty.get_ylim())
 
 
 
plt.show()

# ----- attempt of new figure... to be better finalised

plt.figure(1, figsize=(6, 6))
axH2D_2 = plt.axes(rect_H2D)
nbins2=100
nbinx2=nbins2
nbiny2=nbins2
rangex2=[-2,2]
rangey2=[-2,2]
(n2, mybinx2,mybiny2, patches2) = axH2D_2.hist2d(xstereo_array, ystereo_array, bins=(nbinx2,nbiny2), range=[rangex2,rangey2], norm=LogNorm())
#(n2, mybinx2,mybiny2, patches2) = axH2D_2.hist2d(xstereo_array, ystereo_array, bins=(nbinx2,nbiny2), range=[rangex2,rangey2])
axH2D_2.set_ylabel('y-lambert')
axH2D_2.set_xlabel('x-lambert')

plt.show()


print ("done!")
    
