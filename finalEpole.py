import numpy as np
import math
import csv
import matplotlib.pyplot as plt
import scipy
from scipy.fftpack import fftshift
from scipy.fftpack import diff
from scipy import signal
from datetime import datetime
import math
from math import sqrt
from statistics import mode
from datetime import date
import datetime
import jdcal
from jdcal import jd2gcal
from PyAstronomy import pyasl
import julian


#f=open ('magnetogram.txt','r')
#f=open ('2005-01-21/2005-01-21.txt','r')
#f=open ('haloween_CLF.txt','r')
f=open ('/Users/Taninka/Documents/skola,magister/diplomka/python/modelovanie GIC/haloween.txt','r')
#f=open ('magnetogram.txt','r')
lines=f.readlines()[1:]
f.close()

x=[]
y=[]
t=[]
d=[]

for line in lines:
    p=line.split()
    d.append((p[0]))
    x.append(float(p[3]))
    y.append(float(p[4]))
    t.append((p[1]))
print(d[0],d[-1])
    
t=np.array(t)
d=np.array(d)
jds=[]
for i in range(len(t)):
#    print(i)
    (h,m,s)=t[i].split(':')
    (year,month,day)=d[i].split('-')
    #(day,month,year)=d[i].split('.')
    jds.append(sum(jdcal.gcal2jd(year, month, day))+(float(h)/24+float(m)/60/24+float(s)/3600/24))
    #result.append(int(h)*60*60+int(m)*60+float(s))

#pocitej diference mezi jednotlivymi casovymi index
dts=np.roll(jds,-1)-jds
# a najdi modus (nejcasteji hodnota). To bere v uvahu moznost chybejicich dat
dt=mode(dts)  # in days
# generuj novou casovou osu, ktera bude spojita bez der, s krokem z prechoziho a rozsahem od minima do maxima puvodni osy
time_axis=(range(0,int((max(jds)-min(jds))/dt)+1))*dt+min(jds)
# interpoluj data na novou casovou osu a konvertuj do T
Bx=np.interp(time_axis,jds,np.array(x),right=0.0, left=0.0)*1E-9
By=np.interp(time_axis,jds,np.array(y),right=0.0, left=0.0)*1E-9

#kontrola, ze se data nezmenila
#plt.plot(jds, np.array(y)*1E-9, 'r', time_axis, By, 'b')
#plt.show()

#plt.plot(time_axis, (Bx-np.mean(Bx))*1e9, 'b', time_axis, (By-np.mean(By))*1e9, 'r')
#plt.xlabel('Time [seconds from beginning]')
#plt.ylabel('B-mean(B) [nT]')
#plt.legend(['$B_x$', '$B_y$'])
#plt.show()

#stop
#vypocitej derivace. Diky ekvidistantni ose je mozne pouzit numpy.diff a proste delit casovym krokem
# melo by tedy byl v T/s
#dBx=np.diff(Bx,1)/(dt*3600*24)
#dBy=np.diff(By,1)/(dt*3600*24)
#dBx=np.append(dBx,[0])
#dBy=np.append(dBy,[0])

# nebo pouzit diff ze scipy, je treba nastavit periodu, cili delku intervalu ve fyzikalnich jednotkach
dBx=diff(Bx,order=1, period=len(Bx)*(dt*3600*24))
dBy=diff(By,order=1, period=len(Bx)*(dt*3600*24))

#plt.plot(time_axis, dBx*1e9, 'b', time_axis, dBy*1e9, 'r')
#plt.xlabel('Time [seconds from beginning]')
#plt.ylabel('dB/dt [nT/s]')
#plt.legend(['$dB_x/dt$', '$dB_y/dt$'])
#plt.show()


sigma=0.001 # Ohm/m
mu=4*math.pi*1E-7 # SI units

# VERSION 1
def g(t):
    if (t)<0:
        return 0
    else:
        return 1.0/pow(t,0.5)
    
G=np.vectorize(g)

#tshift=(max(time_axis)+min(time_axis))/2.0 # konvolucni jadro musi mit nulovy bod  na stredu osy
#Ex1=1/(sqrt(math.pi*sigma*mu))*signal.fftconvolve(dBy, G(time_axis-tshift),mode='same')*(dt)
#Ey1=-1/(sqrt(math.pi*sigma*mu))*signal.fftconvolve(dBx, G(time_axis-tshift),mode='same')*(dt)

#plt.plot(time_axis, Ex1*1000, 'b', time_axis, Ey1*1000, 'r')
#plt.xlabel('Time [Julian date]')
#plt.ylabel('E [V/km]')
#plt.legend(['$E_x$', '$E_y$'])
#plt.savefig('epole.pdf', format='PDF')
#plt.show()

# VERSION 2
def H(t):
    if (t)==0.0:
        return 0.5
    if (t)<0:
        return 0
    if (t)>0:
        return 1

def chi(t,tau):
    temp=2.0/(pow(math.pi,0.5))*(pow(t,0.5)*H(t)-pow(t-tau,0.5)*H(t-tau))
    return(temp.real)

Chi=np.vectorize(chi,excluded=['tau'])

tshift=(max(time_axis)+min(time_axis))/2.0
#plt.plot(time_axis-tshift,Chi(time_axis-tshift,dt))
#plt.show()
#print(dt)
#stop

print(len(time_axis))
print(len(dBy))
sigma=0.001 # Ohm/m
mu=4*math.pi*1E-7 # SI units
#sh=(max(time_axis)+min(time_axis))/2
#sh=max(time_axis)
day2s=3600*24
Ex2=1/(sqrt(sigma*mu))*signal.fftconvolve(dBy,Chi((time_axis-tshift)*day2s,(dt*day2s)), mode='same')
Ey2=-1/(sqrt(sigma*mu))*signal.fftconvolve(dBx, Chi((time_axis-tshift)*day2s,(dt*day2s)), mode='same')
#Ex=1/(sqrt(sigma*mu))*signal.fftconvolve(Chi(time_axis-tshift,dt),dBy, mode='same')
#Ey=-1/(sqrt(sigma*mu))*signal.fftconvolve(Chi(time_axis-tshift,dt),dBx, mode='same')
print(time_axis[0])
print(min(time_axis))
print(len(time_axis))

date_axis=[]
d=[]
for i in range(len(time_axis)):
    datee=pyasl.daycnv(time_axis[i],mode='dt')
    date_axis.append(datee)
print(date_axis[0])   

plt.figure(figsize=(13,6))
plt.plot(date_axis, Ex2*1e3, 'g',label= 'X component')


plt.xlabel('date',size=17)
plt.ylabel('geoelectric field [V/km]',size=17)
plt.tight_layout()
#plt.legend(['$E_x$', '$E_y$'],fontsize = 'xx-large')
plt.xticks(size = 12)
plt.yticks(size = 12)
plt.legend(loc='upper left',fontsize = 'xx-large')
#plt.savefig("/Users/Taninka/Documents/skola,magister/diplomka/python/modelovanie GIC/EPxfinal.pdf", format="PDF")
plt.show()
#print(np.min(Ex2*1e3), np.max(Ex2*1e3), np.min(Ey2*1e3), np.max(Ey2*1e3))

plt.figure(figsize=(13,6))
plt.plot(date_axis, Ey2*1e3, 'b',label= 'Y component')
plt.xlabel('date',size=17)
plt.ylabel('geoelectric field [V/km]',size=17)
plt.tight_layout()
plt.xticks(size = 12)
plt.yticks(size = 12)
plt.legend(loc='upper left',fontsize = 'xx-large')
#plt.savefig("/Users/Taninka/Documents/skola,magister/diplomka/python/modelovanie GIC/EPyfinal.pdf", format="PDF")
plt.show()


#==========================================================================================
elekpole=open('/Users/Taninka/Documents/skola,magister/diplomka/python/modelovanie GIC/elpole_new.txt', 'w')
for i in range(len(time_axis)):
    (y,m,d,dpart)=jdcal.jd2gcal(2400000.5,time_axis[i]-2400000.5)
    # rok.mesic.den cast_dne E_(W->E) E(S->N)
    elekpole.write(str(y)+'.'+str(m)+'.'+str(d)+' '+str(dpart)+' '+str(Ex2[i]*1e3)+' '+str(Ey2[i]*1e3)+'\n')
elekpole.close()
#==========================================================================================


dlzka_vedeniY=284 #km
Ry=0.0289 #ohm/km
Rx=0.0211 #ohm/km
dlzka_vedeniX=138 #km

Vx=max(Ex2*1e3)*dlzka_vedeniX
Vy=max(Ey2*1e3)*dlzka_vedeniY

Jx=Vx/(Rx*dlzka_vedeniX)
Jy=Vy/(Ry*dlzka_vedeniY)
##
##plt.figure(figsize=(13,6))
##plt.plot(date_axis,El_pole_x*1e3,'r',date_axis,Jx,'b')
##plt.show()

Exmax=max(Ex2*1e3)
Eymax=max(Ey2*1e3)

#N=['2005-1-21','2006-12-14','2011-8-5','2011-9-26','2014-2-27','2015-9-5','2015-12-20','2016-3-6','2016-5-8','2017-3-25']

print('Vysledne hodnoty GEP a GIC pre')
print('Maximum Ex: '+ str(Exmax))
print('x-ova komponenta GIC: '+str(Jx))
print('Maximum Ey: '+ str(Eymax))
print('y-ova komponenta GIC: '+str(Jy))




