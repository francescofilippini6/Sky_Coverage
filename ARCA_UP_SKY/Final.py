import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
import math as m
import random as ran
from scipy.optimize import curve_fit
from scipy.integrate import quad
from mpl_toolkits.mplot3d import Axes3D
from decimal import Decimal
#-----------------------------------------------
# dataframe with the cross section points as in
# "Ultrahigh-Energy Neutrino Interactions" Raj Gandhi
#-----------------------------------------------
ncs={
        'Energy (GeV)':[10**1,10**2,10**3,10**4,10**5,10**6,10**7,10**8,10**9,10**10,10**11,10**12],
        'CrossSectionCC':[0.777*10**(-37),0.697*10**(-36),0.625*10**(-35),0.454*10**(-34),0.196*10**(-33),0.611*10**(-33),0.176*10**(-32),0.478*10**(-32),0.123*10**(-31),0.301*10**(-31),0.706*10**(-31),0.159*10**(-30)],
        'minelasticityCC':[0.483,0.477,0.472, 0.426,0.332,0.273,0.250,0.237,0.225,0.216,0.208,0.205]}

simulatedncs={
        'Energy (GeV)':[1.00*10,1.26*10,1.58*10,2.00*10,2.51*10,3.16*10,3.98*10,5.01*10,6.31*10,7.94*10,1*10**2,1.26*10**2,1.58*10**2,2.00*10**2,2.51*10**2,3.16*10**2,3.98*10**2,5.01*10**2,6.31*10**2,7.94*10**2,1*10**3,1.26*10**3,1.58*10**3,2.00*10**3,2.51*10**3,3.16*10**3,3.98*10**3,5.01*10**3,6.31*10**3,7.94*10**3,1*10**4,1.26*10**4,1.58*10**4,2.00*10**4,2.51*10**4,3.16*10**4,3.98*10**4,5.01*10**4,6.31*10**4,7.94*10**4,1*10**5,1.26*10**5,1.58*10**5,2.00*10**5,2.51*10**5,3.16*10**5,3.98*10**5,5.01*10**5,6.31*10**5,7.94*10**5,1*10**6,1.26*10**6,1.58*10**6,2.00*10**6,2.51*10**6,3.16*10**6,3.98*10**6,5.01*10**6,6.31*10**6,7.94*10**6,1*10**7],
        'CrossSCC': [ 8.08*10**-38, 1*10**-37, 1.24*10**-37, 1.55*10**-37, 1.92*10**-37, 2.39*10**-37, 2.96*10**-37,3.69*10**-37,4.58*10**-37,5.71*10**-37,7.1*10**-37,8.85*10**-37,1.1*10**-36,1.38*10**-36,1.71*10**-36,2.14*10**-36,2.67*10**-36,3.32*10**-36,4.14*10**-36,5.15*10**-36,6.39*10**-36,7.93*10**-36,9.81*10**-36,1.21*10**-35,1.49*10**-35,1.83*10**-35,2.23*10**-35,2.71*10**-35,3.28*10**-35,3.95*10**-35,4.72*10**-35,5.62*10**-35,6.65*10**-35,7.83*10**-35,9.16*10**-35,1.07*10**-34,1.24*10**-34,1.43*10**-34,1.64*10**-34,1.88*10**-34,2.14*10**-34,2.44*10**-34,2.76*10**-34,3.13*10**-34,3.53*10**-34,3.97*10**-34,4.47*10**-34,5.01*10**-34,5.62*10**-34,6.29*10**-34,7.02*10**-34,7.84*10**-34,8.73*10**-34,9.72*10**-34,1.08*10**-33,1.2*10**-33,1.33*10**-33,1.47*10**-33,1.63*10**-33,1.8*10**-33,1.99*10**-33],
        'CrossNC':[2.49*10**-38,3.09*10**-38, 3.84*10**-38, 4.78*10**-38, 5.93*10**-38, 7.37*10**-38, 9.17*10**-38, 1.14*10**-37, 1.42*10**-37, 1.77*10**-37, 2.2*10**-37, 2.74*10**-37, 3.42*10**-37, 4.27*10**-37, 5.33*10**-37, 6.65*10**-37, 8.3*10**-37, 1.04*10**-36, 1.29*10**-36, 1.61*10**-36, 2*10**-36, 2.49*10**-36, 3.09*10**-36, 3.83*10**-36, 4.73*10**-36, 5.83*10**-36, 7.16*10**-36, 8.75*10**-36, 1.07*10**-35, 1.29*10**-35, 1.56*10**-35, 1.87*10**-35, 2.23*10**-35, 2.64*10**-35, 3.11*10**-35, 3.65*10**-35, 4.27*10**-35, 4.96*10**-35, 5.74*10**-35, 6.61*10**-35, 7.59*10**-35, 8.68*10**-35, 9.89*10**-35, 1.12*10**-34, 1.28*10**-34, 1.44*10**-34, 1.63*10**-34, 1.83*10**-34, 2.06*10**-34, 2.32*10**-34, 2.6*10**-34, 2.91*10**-34, 3.25*10**-34, 3.62*10**-34, 4.04*10**-34, 4.5*10**-34, 5*10**-34, 5.56*10**-34, 6.16*10**-34, 6.83*10**-34, 7.56*10**-34],

        'Total_CS':[]

}

for a in range(len(simulatedncs['CrossSCC'])):
        simulatedncs['Total_CS'].append(simulatedncs['CrossSCC'][a]+simulatedncs['CrossNC'][a])

df = pd.DataFrame(data=ncs)
simulated = pd.DataFrame(data=simulatedncs)
#print(simulated)
#print(ncs['Energy (GeV)'][1])
#-----------------------------------------------
# Fit of the linear part till 10^12 eV
#-----------------------------------------------
def fitfuncUND(x, a, b):
	return a * x + b

#popt, _ = curve_fit(fitfuncUND,df['Energy (GeV)'][:3],df['CrossSectionCC'][:3])
#a,b=popt

allm=[]
allq=[]
for ele in range(len(ncs['minelasticityCC'])-1):
        params, _ =curve_fit(fitfuncUND,ncs['Energy (GeV)'][ele:ele+2],ncs['minelasticityCC'][ele:ele+2])
        allmb,allqb=params
        allm.append(allmb)
        allq.append(allqb)
#print("allm",allm)
#print("allq",allq)

def Bjorken(x):
        ss=[]
        if hasattr(x, "__len__"):
                for b in x:
                        for a in range(len(ncs['minelasticityCC'])-1):
                                if b>=ncs['Energy (GeV)'][a] and b<ncs['Energy (GeV)'][a+1]:
                                        ss.append(fitfuncUND(b,allm[a], allq[a]))
                return ss
        else:
                for a in range(len(ncs['minelasticityCC'])-1):
                        if x>=ncs['Energy (GeV)'][a] and x<ncs['Energy (GeV)'][a+1]:
                                #print("Single Energy Value Bjorken")
                                return fitfuncUND(x,allm[a], allq[a])
#def fitBjorken(x, a, b, c, d, e, f, g):
#	return a * x + b

#fittedBjorken=np.polyfit(df['Energy (GeV)'],df['minelasticityCC'],10)
#funcY=np.poly1d(fittedBjorken)
#-----------------------------------------------
#extrapolate the function for 10^16<E<10^21 eV as in article
#-----------------------------------------------

constCC=2.69*10**-36
indexCC=0.402

constNC=1.06*10**-36
indexNC=0.408

def fitfuncOVE(const,index,x):
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


def closest_energy_value(x):
        ss=[]
        for energy in simulated['Energy (GeV)']:
                ss.append(abs(energy-x))
        #minimum=min(ss)
        min_index = ss.index(min(ss))
        #print(simulated['Energy (GeV)'][min_index])
        return (simulated['Total_CS'][min_index],simulated['CrossSCC'][min_index],simulated['CrossNC'][min_index])

Na=6.022*10**23
def buildingLint(x,interaction):
        cross_interaction=0
        if x <= 1*10**7:
                if interaction=='NC':
                        cross_interaction=closest_energy_value(x)[2]
                elif interaction=='CC':
                        cross_interaction=closest_energy_value(x)[1]
                elif interaction=='Total':
                        cross_interaction=closest_energy_value(x)[0]
                else:
                        print("No other known type of interactions")
                return 1/(Na*cross_interaction)

        elif x > 1* 10**7:
                if interaction=='NC':
                        cross_interaction=fitfuncOVE(constNC,indexNC,x)
                elif interaction=='CC':
                        cross_interaction=fitfuncOVE(constCC,indexCC,x)
                elif interaction=='Total':
                        cross_interaction=fitfuncOVE(constNC,indexNC,x)+fitfuncOVE(constCC,indexCC,x)
                else:
                        print("No other known type of interactions")
                return 1/(Na*cross_interaction)
        else:
                print("Energy too high > 10^21 eV or too low < 10 GeV")
                return 0

def Lint(x,interaction):
        if hasattr(x, "__len__"):
                ss=[]
                for a in x:
                        ss.append(buildingLint(a,interaction))
                return ss
        else:
                return buildingLint(x,interaction)        



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
        #print("D",distance_from_center)
        status=buildingDensity(distance_from_center)[1]
        #print("Numbers of traversed layers",status)
        #print("radia",radia[-status:])
        
        integral_extrema=[0]
        for rr in radia[-status:]:
                integral_extrema.append(np.arccos(distance_from_center/rr))
        #print("Extremis",integral_extrema)
        final=[]
        for a in range(status):
                integrand = lambda x: 2*buildingDensity(distance_from_center/np.cos(x))[0] * distance_from_center*10**5/(np.cos(x)*np.cos(x))
                Integral=quad(integrand,integral_extrema[a],integral_extrema[a+1])
                final.append(Integral[0])
        #print("FINAL value of the integral contributions",final)
        int_length=sum(final)
        #print("sum of contributions:",int_length)
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
                                #print("Too high!!")
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
                        #print("Too high!!")
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
        #print("D",distance_from_center)
        integral_extrema=[]
        #added a +3.5 on first integral
        integral_extrema.append(np.arccos(distance_from_center/(EarthRadius+3.5)))
        integral_extrema.append(np.arccos(distance_from_center/(EarthRadius+15)))
        integral_extrema.append(np.arccos(distance_from_center/(EarthRadius+atm_sup_limit)))
        #print("Extremis",integral_extrema)
        final=[]
        for a in range(2):
                integrand = lambda x: atmosphere(distance_from_center/np.cos(x)) * distance_from_center*10**5/(np.cos(x)*np.cos(x))
                Integral=quad(integrand,integral_extrema[a],integral_extrema[a+1])
                final.append(Integral[0])
        #print("FINAL values for ATMOSPHERE",final)
        int_length=sum(final)
        #print("ATMOSPHERE sum of contributions:",int_length)
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
                        if a > EarthRadius and  a < EarthRadius+2:  # h in km
                                #a-=EarthRadius
                                ss.append(1.040)
                        elif a > EarthRadius+2 and a< EarthRadius+5:
                                #a-=EarthRadius
                                ss.append(1.040)
                        else:
                                ss.append(0)
                return(ss)
        else:
                if h > EarthRadius and  h < EarthRadius+2:  # h in cm
                        # h-=EarthRadius
                        return 1.040
                elif h > EarthRadius+2 and h< EarthRadius+5:
                        #h-=EarthRadius*10**5
                        return 1.040
                else:
                        #print("Too high!!")
                        return 0
        

#def sea_profile(theta):
#        """
#        product of rho x l
#        with l = h/cos(theta)
#        Should be an integral if sea(h) implemented in a different way
#        """
#        h=3.5*10**5
#       return 1.040*h/abs(np.cos(np.radians(theta)))

#--------------------------------------------------
# sea profile integrate to account for the Earth curvature
#--------------------------------------------------
def sea_profile(theta):
        """
        """
        distance_from_center = EarthRadius*np.sin(np.radians(180-theta))
        integrand = lambda x: sea(distance_from_center/np.cos(x)) * distance_from_center*10**5/(np.cos(x)*np.cos(x))
        Integral=quad(integrand,np.arccos(distance_from_center/EarthRadius),np.arccos(distance_from_center/(EarthRadius+3.5)))
        return(Integral[0])
        
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
                


#def number_of_Lint(energy,theta):
#        final=np.zeros((len(energy),len(theta)))
#        encounter=0
#        for x in energy:
#                lint=Lint(x,'Total')
#                angcounter=0
#                for ang in theta:
#                        final[encounter][angcounter]=total_slant(ang)/lint
#                        angcounter+=1
#                encounter+=1
#        return final
#------------------------------------------------------
#improved performances
#------------------------------------------------------
def number_of_Lint(energy,theta):
        status=0
        if hasattr(energy, "__len__") and  hasattr(theta, "__len__"):
                final=np.zeros((len(energy),len(theta)))
                encounter=0
                lint=Lint(energy,'Total')
                total=total_slant(theta)
                for x in lint:
                        angcounter=0
                        for y in total:
                                final[encounter][angcounter]=y/x
                                angcounter+=1
                        encounter+=1
                status=1
                return final,status
        elif hasattr(energy, "__len__") and np.isscalar(theta):
                final=np.zeros(len(energy))
                status=2
                print("STATO",status)
                encounter=0
                lint=Lint(energy,'Total')
                total=total_slant(theta)
                for x in lint:
                        final[encounter]=total/x
                        encounter+=1
                return final,status
        elif hasattr(theta, "__len__") and np.isscalar(energy):
                final=np.zeros(len(theta))
                status=3
                print("STATO",status)
                angcounter=0
                lint=Lint(energy,'Total')
                total=total_slant(theta)
                for y in total:
                        final[angcounter]=y/lint
                        angcounter+=1
                return final,status
        
        else:
                status=4
                print("STATO",status)
                return total_slant(theta)/Lint(energy,'Total'),status
#------------------------------------------------------
# neutrino survival probability
#------------------------------------------------------
def Survival_probability_mean(energy,theta):
        Lint,status=number_of_Lint(energy,theta)
        return np.exp(-Lint)

def muon_range(energy):
        alpha=0.274   #GeV/m
        beta=3.49*10**-4  #1/m
        Ec=alpha/beta
        ss=[]
        rho=1 #g/cm^3 
        if hasattr(energy, "__len__"):
                for x in energy:
                        E_muon=x*(1-Bjorken(x))
                        if E_muon<Ec:
                                ss.append(E_muon*rho/alpha*10**2)   #g/cm^2
                        else:
                                ss.append(rho/beta*m.log(1+beta*E_muon/alpha)*10**2)  #g/cm^2
                return ss
        else:
                E_muon=energy*(1-Bjorken(energy))
                if E_muon<Ec:
                        return E_muon*rho/alpha*10**2
                else:
                        return rho/beta* m.log(1+beta*E_muon/alpha) *10**2

lambda_ass=60*10**2   #cm  
def radius_can(energy):
        radius=[]
        if hasattr(energy, "__len__"):
                for a in energy:
                        radius.append(muon_range(a)+np.cos(np.radians(42))*lambda_ass)
                return radius
        else:
                return muon_range(energy)+np.cos(np.radians(42))*lambda_ass
 



def nu_survival_can_level(energy,theta):
        numb_Lint,status=number_of_Lint(energy,theta)
        print("$P_S$")
        lint=Lint(energy,'Total')
        #print("LINT",Lint)
        print("status",status)
        if status==1:
                probability=np.zeros((len(energy),len(theta)))
                Rmu=radius_can(energy)
                #print("Radius can",Rmu)
                slant=total_slant(theta)
                print("slant",len(slant),len(theta))
                for ene in range(len(Rmu)):
                        for ang in range(len(theta)):
                                if slant[ang]<=Rmu[ene]:
                                        probability[ene][ang]=1
                                else:
                                        index=numb_Lint[ene][ang]-Rmu[ene]*np.power(lint[ene],-1)
                                        probability[ene][ang]=np.exp(-index)
                return probability
        elif status==2:
                probability=[]
                Rmu=radius_can(energy)
                slant=total_slant(theta)
                for ene in range(len(energy)):
                        if slant<=Rmu[ene]:
                                probability.append(1)
                        else:
                                index=numb_Lint[ene]-Rmu[ene]*np.power(lint[ene],-1)
                                probability.append(np.exp(-index))
                return probability
        elif status==3:
                probability=[]
                Rmu=radius_can(energy)
                slant=total_slant(theta)
                for ang in range(len(theta)):
                        if slant[ang]<=Rmu:
                                #print(theta[ang])
                                probability.append(1)
                        else:
                                index=numb_Lint[ang]-Rmu*np.power(lint,-1)
                                probability.append(np.exp(-index))
                return probability
        
        elif status==4:
                Rmu=radius_can(energy)
                slant=total_slant(theta)
                if slant<=Rmu:
                        return(1)
                else:
                        index=numb_Lint-Rmu*np.power(lint,-1)
                        return(np.exp(-index))


def nu_interaction_inside_can(energy,theta):
        numb_Lint,status=number_of_Lint(energy,theta)
        lint=Lint(energy,'Total')  #possibility to set the CC interaction..
        surv=Survival_probability_mean(energy,theta)
        print("$P_i$")
        print("status",status)
        if status==1:
                probability=np.zeros((len(energy),len(theta)))
                Rmu=radius_can(energy)
                #print("Radius can",Rmu)
                slant=total_slant(theta)
                #print("slant",slant,theta)
                for ene in range(len(Rmu)):
                        for ang in range(len(theta)):
                                if slant[ang]<=Rmu[ene]:
                                        probability[ene][ang]=1-surv[ene][ang]
                                else:
                                        index=Rmu[ene]*np.power(lint[ene],-1)
                                        probability[ene][ang]=1-np.exp(-index)
                return probability
       
        elif status==2:
                probability=[]
                Rmu=radius_can(energy)
                slant=total_slant(theta)
                for ene in range(len(energy)):
                        if slant<=Rmu[ene]:
                                probability.append(1-surv[ene])
                        else:
                                index=Rmu[ene]*np.power(lint[ene],-1)
                                probability.append(1-np.exp(-index))
                return probability

        elif status==3:
                probability=[]
                Rmu=radius_can(energy)
                slant=total_slant(theta)
                for ang in range(len(theta)):
                        if slant[ang]<=Rmu:
                                probability.append(1-surv[ang])
                                #print(theta[ang])
                        else:
                                index=Rmu*np.power(lint,-1)
                                probability.append(1-np.exp(-index))
                return probability
        
        elif status==4:
                Rmu=radius_can(energy)
                slant=total_slant(theta)
                if slant<=Rmu:
                        return(1-surv)
                else:
                        index=Rmu*np.power(lint,-1)
                        return(1-np.exp(-index))


#-------------------------------------------------------------

        
#if total_slant(theta)<=muon_range(energy):
        #          return 1-Survival_probability_mean(energy,theta)
        # else:
        #         fraction=radius_can(energy)*np.power(Lint(energy,'Total'),-1)
        #         return np.ones(len(energy))-np.exp(-fraction)
       
def final_prbability(energy,theta):
        print("final")
        numb_Lint,status=number_of_Lint(energy,theta)
        interaction=nu_interaction_inside_can(energy,theta)
        survival=nu_survival_can_level(energy,theta)
        if status==1:
                probability=np.zeros((len(energy),len(theta)))
                for ene in range(len(energy)):
                        for ang in range(len(theta)):
                                probability[ene][ang]=survival[ene][ang]*interaction[ene][ang]
                return probability

        elif status==2:
                probability=[]
                for ene in range(len(energy)):
                        probability.append(survival[ene]*interaction[ene])
                return probability

        elif status==3:
                #print("Final",status)
                probability=[]
                for ang in range(len(theta)):
                        probability.append(survival[ang]*interaction[ang])
                return probability
        elif status==4:
                return survival*interaction
        
year=3.15*10**7 #seconds
A_geom=1*10**10 #cm^2

def A_eff(energy,theta):
        over=[]
        under=[]
        for ang in theta:
                if ang >90:
                        under.append(ang)
                else:
                        over.append(ang)
        AeffOver=[]
        AeffUnder=[]
        reco_eff=0.1
        A_times_recoeff=A_geom*reco_eff*10**-4 #area in m^2 
        
        #working configuration
        #for ene in energy:
        #AeffOver.append(sum(final_prbability(ene,over)))
                #        AeffOver.append(final_prbability(ene,over[1]))
                #AeffUnder.append(sum(final_prbability(ene,under)))
                #       AeffUnder.append(final_prbability(ene,under[-1]))
        for ang in theta:
                ciccio1=[]
                ciccio2=[]
                if ang < 90:
                        for ene in energy:
                                ciccio1.append(A_times_recoeff*final_prbability(ene,ang))
                        AeffOver.append(ciccio1)
                else:
                        for ene in energy:
                                ciccio2.append(A_times_recoeff*final_prbability(ene,ang))
                        AeffUnder.append(ciccio2)
        print(len(AeffOver))
        #AeffOver=np.array(AeffOver)#*A_times_recoeff
        #AeffUnder=np.array(AeffUnder)#*A_times_recoeff
        return (AeffOver,AeffUnder)
                
def nu_flux(energy):
        ss=[]
        if hasattr(energy, "__len__"):
                for ene  in energy:
                        ss.append(10**-8*ene**-2)
                        
                return ss
        else:
                return 10**-8*energy**-2
 
        
def Numb_events(energy,theta):
        numb_Lint,status=number_of_Lint(energy,theta)
        factor=year*A_geom*2*np.pi
        proba=final_prbability(energy,theta)
        flux=nu_flux(energy)
        if status==1:
                events=np.zeros((len(energy),len(theta)))
                for ene in range(len(energy)):
                        for ang in range(len(theta)):
                                events[ene][ang]=proba[ene][ang]*flux[ene]*factor#*np.sin(theta[ang])
                return events
        elif status==2:
                events=[]
                for ene in range(len(energy)):
                        events.append(proba[ene]*flux[ene]*factor)#*np.sin(theta))
                return events
        
        elif status==3:
                events=[]
                for ang in range(len(theta)):
                        events.append(proba[ang]*flux*factor)#*np.sin(theta[ang]))
                return events
        elif status==4:
                return proba*flux*factor*np.sin(theta)


def sum_of_events(energy,theta):
        #limit1=a
        #limit2=b
        events=Numb_events(energy,theta)
        sums=[]
        for ang in theta:
                sums.append(sum(Numb_events(energy,ang)))
        return sums

#-------------------------------------------------------------
# Plotting functions
#-------------------------------------------------------------

x_lint = np.logspace(1,11, 100)#np.linspace(10,0.99*10**12, 100000)
x_lineO = np.logspace(6,12, 100)#x_lineO = np.linspace(10**6,10**12, 100)
y_lineO=fitfuncOVE(constCC,indexCC,x_lineO)
y_lineNC=fitfuncOVE(constNC,indexNC,x_lineO)
y_line_total=[]
for a in range(len(y_lineO)):
        y_line_total.append(y_lineO[a]+y_lineNC[a])
        
radius=np.linspace(0,EarthRadius,100)
hatm=np.linspace(EarthRadius+5,EarthRadius+atm_sup_limit,100)
hsea=np.linspace(EarthRadius+0.1,EarthRadius+4.9999,100)
ttheta=np.linspace(0.1,89.9,100)
stheta=np.linspace(90.01,179.9,100)
cc=[]
ff=[]
for aa in stheta:
     cc.append(int_lentgth_earth(aa))   
     ff.append(sea_profile(aa))
     #print(stheta)
dd=[]
for bb in stheta:
        dd.append(int_length_atm(bb))

#-------------------------------------------------------------
#surface plot number generators
#-------------------------------------------------------------
Number_of_points=100
#x_energy = np.linspace(10,10**11, Number_of_points)
x_energy = np.logspace(1,7, Number_of_points)
zenithy=np.linspace(0,179.9,Number_of_points)
#zenithy=np.logspace(0,2,Number_of_points)

#-------------------------------------------------------------
# Neutrino functions: Cross Section CC, Bjorken y, Interaction Length in (g/cm^2)
#-------------------------------------------------------------
f=plt.figure()
ax1=f.add_subplot(311)
df.plot(kind='scatter',x='Energy (GeV)',y='CrossSectionCC',color='green',ax=ax1, logx=True, logy=True)
simulated.plot(kind='scatter',x='Energy (GeV)',y='Total_CS',color='red',ax=ax1, logx=True, logy=True)
simulated.plot(kind='scatter',x='Energy (GeV)',y='CrossSCC',color='green',ax=ax1, logx=True, logy=True)
simulated.plot(kind='scatter',x='Energy (GeV)',y='CrossNC',color='blue',ax=ax1, logx=True, logy=True)
#ax1.plot(x_line,y_line,'b--')
ax1.plot(x_lineO,y_lineO,'+',markersize=2, color='g')#'g--')
ax1.plot(x_lineO,y_lineNC,'+',markersize=2, color='b')#'b--')
ax1.plot(x_lineO,y_line_total,'+',markersize=2, color='r')
ax1.set_title('nu N CC Cross Section')
ax2=f.add_subplot(312)
df.plot(kind='scatter',x='Energy (GeV)',y='minelasticityCC',color='g',ax=ax2,logx=True,logy=True)
ax2.plot(x_lint,Bjorken(x_lint),'+',markersize=2, color='r')
ax2.set_title('Mean Bjorken y nu N CC')
ax3=f.add_subplot(313)
#ax3.plot(x_line,Lint(x_line),'r--')
ax3.plot(x_lint,Lint(x_lint,'CC'),'+',markersize=2, color='g')#'g--')
ax3.plot(x_lint,Lint(x_lint,'Total'),'+',markersize=2, color='r')#'r--')
ax3.plot(x_lint,Lint(x_lint,'NC'),'+',markersize=2, color='b')#'b--')
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
ax.plot(radius, density_profile(radius),'b-')
ax.set_ylabel('Density')
ax.set_xlabel('Earth radius (km)')
ax.set_title('Earth density profile')
ax4=f2.add_subplot(122)
#ax4.plot(stheta,int_lentgth_earth(stheta),'b--')
ax4.plot(stheta,cc,'+',markersize=2, color='b')
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
ax6.plot(stheta,dd,'+',markersize=2, color='b')#'g--')
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
ax8.plot(stheta,ff,'g--')
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
#ax9.plot(stheta,cc,'b--')
#ax9.plot(stheta,dd,'g--')
#ax9.plot(stheta,sea_profile(stheta),'g--')
ax9.set_ylabel('x rho (g/cm^2)')
ax9.set_xlabel('zenith angle (deg)')
ax9.set_yscale('log')
f5.suptitle('TOTAL SLANT DEPTH vs ZENITH')



#-------------------------------------------------------------
# SURFACE PLOTS:
#   10) Number on Lint; 11) survival probability; 12) survival at can level
#   13) Interaction inside can;  14) 12*13= final probability
#-------------------------------------------------------------
fig = plt.figure()
ax10 = fig.add_subplot(121, projection='3d')
ax10.set_title('zenith energy # of Lint')
X1=x_energy
Y1=zenithy
Z1=number_of_Lint(X1,Y1)[0]
X1, Y1 = np.meshgrid(X1, Y1)
surf = ax10.plot_surface(X1, Y1, Z1, cmap=cm.coolwarm,linewidth=0, antialiased=False)
ax10.set_xlabel('energy (GeV)')
ax10.set_ylabel('zenith angle (deg)')
#ax10.colorbar(surf, shrink=0.5, aspect=5)
ax11 = fig.add_subplot(122, projection='3d')
Z2=Survival_probability_mean(x_energy,zenithy)
ax11.set_xlabel('energy (GeV)')
ax11.set_ylabel('zenith angle (deg)')
surf1=ax11.plot_surface(X1, Y1, Z2, cmap=cm.coolwarm,linewidth=0, antialiased=False)
#plt.imshow(number_of_Lint(x_lint,zenithy),origin='lower',interpolation='none')
#fig.suptitle("distribution")
fig.colorbar(surf1, shrink=0.5, aspect=5)

figg= plt.figure()
axmuon = figg.add_subplot(111)
axmuon.set_title('muon range')
axmuon.plot(x_energy,muon_range(x_energy),'+',markersize=2, color='r')
axmuon.plot(x_energy,radius_can(x_energy),'+',markersize=2, color='g')
axmuon.set_xlabel('energy (GeV)')
axmuon.set_ylabel('Range (g/cm^2)')
axmuon.set_xscale('log')


figmuon = plt.figure()
ax12 = figmuon.add_subplot(231, projection='3d')
X2=x_energy
Y2=zenithy
Z3=nu_survival_can_level(X2,Y2)
Z4=nu_interaction_inside_can(X2,Y2)
Z5=final_prbability(X2,Y2)
X2, Y2 = np.meshgrid(X2, Y2)
ax12.set_xlabel('energy (GeV)')
ax12.set_ylabel('zenith angle (deg)')
ax12.set_title("survival at can level")
surf2=ax12.plot_surface(X2, Y2, Z3, cmap=cm.coolwarm,linewidth=0, antialiased=False)


ax13 = figmuon.add_subplot(232, projection='3d')
ax13.set_xlabel('energy (GeV)')
ax13.set_ylabel('zenith angle (deg)')
surf3=ax13.plot_surface(X2, Y2, Z4, cmap=cm.coolwarm,linewidth=0, antialiased=False)
ax13.set_title("Interaction inside can")
#plt.imshow(number_of_Lint(x_lint,zenithy),origin='lower',interpolation='none')
ax14 = figmuon.add_subplot(233, projection='3d')
ax14.set_xlabel('energy (GeV)')
ax14.set_ylabel('zenith angle (deg)')
ax14.set_title("Final Probability")
surf2=ax14.plot_surface(X2, Y2, Z5, cmap=cm.coolwarm,linewidth=0, antialiased=False)


#print("_---------------______________------___DEBUG_______------_______---_--_-__-_-_-_--_-__")
E_fixed=10**7
ax15 = figmuon.add_subplot(234)
ax15.plot(zenithy,nu_survival_can_level(E_fixed,zenithy),'+',markersize=2)

ax16 = figmuon.add_subplot(235)
ax16.plot(zenithy,nu_interaction_inside_can(E_fixed,zenithy),'+',markersize=2)

ax17 = figmuon.add_subplot(236)
ax17.plot(zenithy,final_prbability(E_fixed,zenithy),'+',markersize=2)

TeVSky = np.logspace(3,6, Number_of_points)
PeVSky = np.logspace(6,9, Number_of_points)
EeVSky = np.logspace(9,11, Number_of_points)

Nevents=plt.figure()
ax18=Nevents.add_subplot(231,projection='3d')
ax18.set_title('Number of events')
X3=x_energy
Y3=zenithy
Z6=Numb_events(X3,Y3)
X3, Y3 = np.meshgrid(X3, Y3)
surf3=ax18.plot_surface(X3, Y3, Z6, cmap=cm.coolwarm,linewidth=0, antialiased=False)
ax19=Nevents.add_subplot(232)
ax19.set_title(Decimal(str(E_fixed)))
ax19.plot(zenithy,Numb_events(E_fixed,zenithy),'+',markersize=2)
ax20=Nevents.add_subplot(233)
ax20.set_title('TeV Sky (1 TeV-1PeV)')
ax20.plot(zenithy,sum_of_events(TeVSky,zenithy),'+',markersize=2,color='r')
ax21=Nevents.add_subplot(234)
ax21.set_title('PeV Sky (1PeV-1EeV)')
ax21.plot(zenithy,sum_of_events(PeVSky,zenithy),'+',markersize=2,color='b')
ax22=Nevents.add_subplot(235)
ax22.set_title('EeV Sky (1EeV-$10^21$ eV)')
ax22.plot(zenithy,sum_of_events(EeVSky,zenithy),'+',markersize=2,color='g')


aeff=plt.figure()
ax23=aeff.add_subplot(121)
Over,Under=A_eff(TeVSky,zenithy)
#print(Over)
ax23.set_title('$A_{eff}$ Over')
sumsOver=np.zeros(len(Over[0]))
for aefff in range(len(Over)):
        ax23.plot(TeVSky,Over[aefff],'+',markersize=2)#,color='g')
        sumsOver=sumsOver+np.array(Over[aefff])
ax23.plot(TeVSky,sumsOver,'g--')
ax23.set_xscale('log')
ax23.set_yscale('log')

ax24=aeff.add_subplot(122)
ax24.set_title('$A_{eff}$ Under')
sumsUnder=np.zeros(len(Under[0]))
for ae in range(len(Under)):
        ax24.plot(TeVSky,Under[ae],'+',markersize=2)#,color='g')
        sumsUnder=sumsUnder+np.array(Under[ae])
ax24.plot(TeVSky,sumsUnder,'g--')
ax24.set_xscale('log')
ax24.set_yscale('log')

plt.show()
 