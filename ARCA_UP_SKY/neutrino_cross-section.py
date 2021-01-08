import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math as m
import random as ran
from scipy.optimize import curve_fit
from scipy.integrate import quad
#-----------------------------------------------
# dataframe with the cross section points as in
# "Ultrahigh-Energy Neutrino Interactions" Raj Gandhi
#-----------------------------------------------
ncs={
    'Energy (GeV)':[10**1,10**2,10**3,10**4,10**5,10**6,10**7,10**8,10**9,10**10,10**11,10**12],
    'CrossSectionCC':[0.777*10**(-37),0.697*10**(-36),0.625*10**(-35),0.454*10**(-34),0.196*10**(-33),0.611*10**(-33),0.176*10**(-32),0.478*10**(-32),0.123*10**(-31),0.301*10**(-31),0.706*10**(-31),0.159*10**(-30)],
    'minelasticityCC':[0.483,0.477,0.472, 0.426,0.332,0.273,0.250,0.237,0.225,0.216,0.208,0.205]}

df = pd.DataFrame(data=ncs)
print(df)

#-----------------------------------------------
# Fit of the linear part till 10^12 eV
#-----------------------------------------------
def fitfuncUND(x, a, b):
	return a * x + b

popt, _ = curve_fit(fitfuncUND,df['Energy (GeV)'][:3],df['CrossSectionCC'][:3])
a,b=popt
print("LINEAR PARAM. FIT",a,b)

#-----------------------------------------------
#extrapolate the function for 10^16<E<10^21 eV as in article
#-----------------------------------------------

const=2.69*10**-36
index=0.402

def fitfuncOVE(x):
        ss=[]
        if hasattr(x, "__len__"):
                for a in x:
                        ss.append(const*np.power(a,index))
                return ss         
        else:
                return const*np.power(x,index)
       
#----------------------------------------------
# defining the interaction length in mwe
#-----------------------------------------------
Na=6.022*10**23

def buildingLint(x):
        if x <= 10**3:
                return 1/(Na*fitfuncUND(x, a, b))
        elif x>= 10**6:
                return 1/(Na*fitfuncOVE(x))
        else:
                print("region between 10**12 and 10**16, still problem")
                return 0
def Lint(x):
        if hasattr(x, "__len__"):
                ss=[]
                for a in x:
                        ss.append(buildingLint(a))
                return ss
        else:
                return buildingLint(x)        



#-------------------------------------------------------------
# Earth density profile wrt x=r/EarthRadius
#-------------------------------------------------------------
EarthRadius=6371   #km
radia=[1221.5,3480,5701,5771,5971,6151,6346.6,6356,6368,6371]
def buildingDensity(r):
        ratio=r/EarthRadius
        status=0
        if r< radia[0]:
                status=10
                return (13.0885-8.8381*ratio**2,status)
        elif r>radia[0] and r< radia[1]:
                status=9
                return (12.5815-1.2638*ratio-3.6426*ratio**2-5.5281*ratio**3,status)
        elif r> radia[1] and r<radia[2]:
                status=8
                return (7.9565 - 6.4761*ratio+5.5283*ratio**2-3.0807*ratio**3,status)
        elif r>radia[2] and r<radia[3]:
                status=7
                return (5.3197-1.4836*ratio,status)
        elif r>radia[3] and r<radia[4]:
                status=6
                return (11.2494-8.0298*ratio,status)
        elif r> radia[4] and r<radia[5]:
                status=5
                return (7.1089 - 3.8045*ratio,status)
        elif r> radia[5] and r< radia[6]:
                status=4
                return (2.691 + 0.6924 *ratio,status)
        elif r> radia[6] and r<radia[7]:
                status=3
                return (2.9,status)
        elif r> radia[7] and r< radia[8]:
                status=2
                return (2.6,status)
        elif r> radia[8] and r<=EarthRadius:
                status=1
                return (1.02,status)
        else:
                print ("Radius out of limits")
                return 0

def density_profile(r):
        """
        call the buildingDensity function for r = scalar or list
        return the density value for the given r
        """

        ss=[]
        if hasattr(r, "__len__"):
                for a in r:
                        ss.append(buildingDensity(a)[0])
                return ss
        else:
                return buildingDensity(r)[0]



def int_lentgth_earth(theta):
        """
        integral of rho x dl (g/cm^2)
        each integral calculated in angle x from cos-1(D/Ri) to cos-1(D/Ri+1)
        where Ri is the radius of the i-th concentric Earth shell
        for i = 0 Ri=0
        """
        distance_from_center = EarthRadius*np.sin(np.radians(180-theta))
        print("D",distance_from_center)
        status=buildingDensity(distance_from_center)[1]
        print("Numbers of traversed layers",status)
        print("radia",radia[-status:])
        
        integral_extrema=[0]
        for rr in radia[-status:]:
                integral_extrema.append(np.arccos(distance_from_center/rr))
        print("Extremis",integral_extrema)
        final=[]
        for a in range(status):
                integrand = lambda x: 2*buildingDensity(distance_from_center/np.cos(x))[0] * distance_from_center*10**5/(np.cos(x)*np.cos(x))
                Integral=quad(integrand,integral_extrema[a],integral_extrema[a+1])
                final.append(Integral[0])
        print("FINAL value of the integral contributions",final)
        int_length=sum(final)
        print("sum of contributions:",int_length)
        return(int_length)

atm_sup_limit=455
def atmosphere(h):
        """
        model of the density profile of the atmosphere
        h in km
        """
        ss=[]
        if hasattr(h, "__len__"):
                for a in h:
                        if a > EarthRadius+5 and  a < EarthRadius+15:  # h in km
                                a-=(EarthRadius+5)
                                ss.append(1.225 *10**-3 * m.exp(-a/9.192))
                        elif a > EarthRadius+15 and a< EarthRadius+atm_sup_limit:
                                a-=(EarthRadius+5)
                                ss.append(1.944 *10**-3 * m.exp(-a/6.452))
                        else:
                                print("Too high!!")
                                ss.append(0)
                return(ss)
        else:
                if h > EarthRadius+5 and  h < EarthRadius+15:  # h in km
                        h-=(EarthRadius+5)
                        return 1.225 *10**-3 * m.exp(-h/9.192)
                elif h > EarthRadius+15 and h< EarthRadius+atm_sup_limit:
                        h-=(EarthRadius+5)
                        return 1.944 *10**-3 * m.exp(-h/6.452)
                else:
                        print("Too high!!")
                        return 0

def int_length_atm(theta):
        """
        integral of rho x dl (g/cm^2)
        each integral calculated in angle x from cos-1(D/Ri) to cos-1(D/Ri+1)
        where Ri is: 
        EarthRadius+5 km
        EarthRadius+15 km
        EarthRadius+atm_sup_limit km
        """
        distance_from_center = EarthRadius*np.sin(np.radians(180-theta))
        print("D",distance_from_center)
        integral_extrema=[]
        integral_extrema.append(np.arccos(distance_from_center/EarthRadius))
        integral_extrema.append(np.arccos(distance_from_center/(EarthRadius+15)))
        integral_extrema.append(np.arccos(distance_from_center/(EarthRadius+atm_sup_limit)))
        print("Extremis",integral_extrema)
        final=[]
        for a in range(2):
                integrand = lambda x: atmosphere(distance_from_center/np.cos(x)) * distance_from_center*10**5/(np.cos(x)*np.cos(x))
                Integral=quad(integrand,integral_extrema[a],integral_extrema[a+1])
                final.append(Integral[0])
        print("FINAL values for ATMOSPHERE",final)
        int_length=sum(final)
        print("ATMOSPHERE sum of contributions:",int_length)
        return(int_length)


def sea(h):
        """
        to be deleted or implemented in case of a more sofisticated 
        sea-water density profile
        By now assumed const = 1.040 g/cm^2
        """
        ss=[]
        if hasattr(h, "__len__"):
                for a in h:
                        if a > EarthRadius*10**5 and  a < EarthRadius*10**5+2000:  # h in cm
                                a-=EarthRadius*10**5
                                ss.append(1.040)
                        elif a > EarthRadius*10**5+2000 and a< EarthRadius*10**5+5000:
                                a-=EarthRadius*10**5
                                ss.append(1.040)
                        else:
                                print("Too high!!")
                                ss.append(0)
                return(ss)
        else:
                if h > EarthRadius*10**5 and  h < EarthRadius*10**5+2000:  # h in cm
                        h-=EarthRadius*10**5
                        return 1.040
                elif h > EarthRadius*10**5+2000 and h< EarthRadius*10**5+5000:
                        h-=EarthRadius*10**5
                        return 1.040
                else:
                        print("Too high!!")
                        return 0
        

def sea_profile(theta):
        """
        product of rho x l
        with l = h/cos(theta)
        Should be an integral if sea(h) implemented in a different way
        """
        h=3.5*10**5
        return 1.040*h/abs(np.cos(np.radians(theta)))


def total_slant(theta):
        ss=[]
        if hasattr(theta, "__len__"):
                for a in theta:
                        if a < 90:
                                ss.append(int_length_atm(a)+sea_profile(a))
                        else:
                                ss.append(int_length_atm(a)+sea_profile(a)+int_lentgth_earth(a))
                return(ss)
        else:
                if theta < 90:
                        return int_length_atm(theta)+sea_profile(theta)
                else:
                        return int_length_atm(theta)+sea_profile(theta)+int_lentgth_earth(theta)
                
#-------------------------------------------------------------
# Plotting functions
#-------------------------------------------------------------

x_line = np.linspace(10**1,10**3, 100)
y_line = fitfuncUND(x_line, a, b)

x_lineO = np.linspace(10**6,10**12, 100)
y_lineO=fitfuncOVE(x_lineO)
radius=np.linspace(0,EarthRadius,100)
hatm=np.linspace(EarthRadius+5,EarthRadius+atm_sup_limit,100)
hsea=np.linspace(EarthRadius*10**5+1,EarthRadius*10**5+4999,100)
ttheta=np.linspace(0.1,89.9,100)
stheta=np.linspace(90.01,179.9,100)
cc=[]
for aa in stheta:
     cc.append(int_lentgth_earth(aa))   
#print(stheta)
dd=[]
for bb in stheta:
        dd.append(int_length_atm(bb))



zenithy=np.linspace(0,179.9,1000)
#-------------------------------------------------------------
# Neutrino functions: Cross Section CC, Bjorken y, Interaction Length in (g/cm^2)
#-------------------------------------------------------------
f=plt.figure()
ax1=f.add_subplot(311)
df.plot(kind='scatter',x='Energy (GeV)',y='CrossSectionCC',color='red',ax=ax1, logx=True, logy=True)
ax1.plot(x_line,y_line,'b--')
ax1.plot(x_lineO,y_lineO,'g--')
ax1.set_title('nu N CC Cross Section')
ax2=f.add_subplot(312)
df.plot(kind='scatter',x='Energy (GeV)',y='minelasticityCC',color='g',ax=ax2,logx=True,logy=True)
ax2.set_title('Mean Bjorken y nu N CC')
ax3=f.add_subplot(313)
ax3.plot(x_line,Lint(x_line),'r--')
ax3.plot(x_lineO,Lint(x_lineO),'r--')
ax3.set_yscale('log')
ax3.set_xscale('log')
ax3.set_ylabel('Lint (g/cm^2)')
ax3.set_xlabel('E_nu (GeV)')
ax3.set_title('Interaction Length')
f.suptitle('Neutrino property functions')
#-------------------------------------------------------------
# Earth: 0) density profile vs radius; 4) slant depth (g/cm^2) vs angle
#-------------------------------------------------------------

f2=plt.figure()
ax=f2.add_subplot(121)
ax.plot(radius, density_profile(radius), 'b-')
ax.set_ylabel('Density')
ax.set_xlabel('Earth radius (km)')
ax.set_title('Earth density profile')
ax4=f2.add_subplot(122)
#ax4.plot(stheta,int_lentgth_earth(stheta),'b--')
ax4.plot(stheta,cc,'b--')
ax4.set_title('slant depth')
ax4.set_ylabel('x rho (g/cm^2)')
ax4.set_xlabel('zenith angle (deg)')
f2.suptitle('Earth')
#-------------------------------------------------------------
# Atmosphere: 5) density profile vs height ; 6) slanth depth (g/cm^2) vs angle
#-------------------------------------------------------------

f3=plt.figure()
ax5=f3.add_subplot(121)
ax5.plot(hatm, atmosphere(hatm),'+',markersize=2, color='g')
ax5.set_ylabel('Density')
ax5.set_xlabel('Atmoosphere height (km)')
ax5.set_yscale('log')
ax6=f3.add_subplot(122)
#ax6.plot(stheta,int_length_atm(stheta),'b--')
ax6.plot(stheta,dd,'g--')
ax6.set_ylabel('x rho (g/cm^2)')
ax6.set_xlabel('zenith angle (deg)')
ax6.set_yscale('log')
f3.suptitle('Atmosphere')
#-------------------------------------------------------------
# Sea: 7) density profile vs height (0 = sea bed, 3.5*10**5 sea level for ARCA;
#      8) slanth depth (g/cm^2) vs angle
#-------------------------------------------------------------

f4=plt.figure()
ax7=f4.add_subplot(121)
ax7.plot(hsea, sea(hsea),'+',markersize=2, color='b')
ax7.set_ylabel('Sea Density')
ax7.set_xlabel('Sea height (km)')
ax8=f4.add_subplot(122)
#ax6.plot(stheta,int_length_atm(stheta),'b--')
#ax8.plot(ttheta,sea_profile((3.5+EarthRadius)*10**5,ttheta),'g--')
ax8.plot(stheta,sea_profile(stheta),'g--')
#ax8.plot(stheta,sea_profile(3.5+EarthRadius),'g--')
ax8.set_ylabel('x rho (g/cm^2)')
ax8.set_xlabel('zenith angle (deg)')
ax8.set_yscale('log')
f4.suptitle('Sea')

#-------------------------------------------------------------
# TOTAL SLANT DEPTH vs ZENITH ANGLE
#-------------------------------------------------------------
f5=plt.figure()
ax9=f5.add_subplot(111)
ax9.plot(zenithy,total_slant(zenithy),'+',markersize=2, color='r')
ax9.set_ylabel('x rho (g/cm^2)')
ax9.set_xlabel('zenith angle (deg)')
ax9.set_yscale('log')
f5.suptitle('TOTAL SLANT DEPTH vs ZENITH')
plt.show()
