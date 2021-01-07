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
        ss=[]
        if hasattr(r, "__len__"):
                for a in r:
                        ss.append(buildingDensity(a)[0])
                return ss
        else:
                return buildingDensity(r)[0]



def int_lentgth_earth(theta):
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


def atmosphere(h):
        if h < 10:  # h in km
                return 1.225 *10**-3 * m.exp(h/9.192)
        elif h > 10 and h< 180:
                return 1.944 *10**-3 * m.exp(h/6.452)
        else:
                print("Too high!!")
                return 0

def int_length_atm(theta):
        distance_from_center = EarthRadius*np.sin(np.radians(180-theta))
        print("D",distance_from_center)
        integral_extrema=[]
        integral_extrema.append(np.arccos(distance_from_center/EarthRadius))
        integral_extrema.append(np.arccos(distance_from_center/EarthRadius+10))
        integral_extrema.append(np.arccos(distance_from_center/EarthRadius+10))
        print("Extremis",integral_extrema)
        final=[]
        for a in range(2):
                integrand = lambda x: atmosphere(distance_from_center/np.cos(x)) * distance_from_center*10**5/(np.cos(x)*np.cos(x))
                Integral=quad(integrand,integral_extrema[a],integral_extrema[a+1])
                final.append(Integral[0])
        print("FINAL value of the integral contributions",final)
        int_length=sum(final)
        print("sum of contributions:",int_length)
        return(int_length)

#print(91)
#int_lentgth_earth(91)
int_length_atm(179.9)
        
#-------------------------------------------------------------
# Plotting functions
#-------------------------------------------------------------

x_line = np.linspace(10**1,10**3, 100)
y_line = fitfuncUND(x_line, a, b)

x_lineO = np.linspace(10**6,10**12, 100)
y_lineO=fitfuncOVE(x_lineO)
radius=np.linspace(0,EarthRadius,100)
tttheta=np.linspace(90.01,179.9,100)
cc=[]
for aa in tttheta:
     cc.append(int_lentgth_earth(aa))   
print(tttheta)
f=plt.figure()
ax1=f.add_subplot(311)
df.plot(kind='scatter',x='Energy (GeV)',y='CrossSectionCC',color='red',ax=ax1, logx=True, logy=True)
ax1.plot(x_line,y_line,'b--')
ax1.plot(x_lineO,y_lineO,'g--')
ax2=f.add_subplot(312)
df.plot(kind='scatter',x='Energy (GeV)',y='minelasticityCC',color='g',ax=ax2,logx=True,logy=True)
ax3=f.add_subplot(313)
ax3.plot(x_line,Lint(x_line),'r--')
ax3.plot(x_lineO,Lint(x_lineO),'r--')
ax3.set_yscale('log')
ax3.set_xscale('log')
ax3.set_ylabel('Lint (mwe)')
ax3.set_xlabel('E_nu (GeV)')

f2=plt.figure()
ax=f2.add_subplot(121)
ax.plot(radius, density_profile(radius), 'b-')
ax.set_ylabel('Density')
ax.set_xlabel('Earth radius (km)')
ax4=f2.add_subplot(122)
#ax4.plot(tttheta,int_lentgth_earth(tttheta),'b--')
ax4.plot(tttheta,cc,'b--')
ax4.set_ylabel('xrho (g/cm^2')
ax4.set_xlabel('zenith angle (deg)')
plt.show()
