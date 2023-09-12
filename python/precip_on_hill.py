import numpy as np
import scipy.interpolate as sci
from netCDF4 import Dataset
import os
import getpass

import matplotlib.pyplot as plt

username=getpass.getuser()


# assume 2.5 m/s e.g. 
u=2.5 # m/s


# the hill
tau=6000.
rt=tau
t=np.linspace(0,rt,100)
hill=tau*0.6/(2.*np.pi)*(1.-np.cos(2.*np.pi/tau*t))/1000.
t=(t+10)*u/1000.
t=np.append(t,t[::-1])
hill=np.append(hill,np.zeros((len(hill),1)))
l1=len(t)

# note the top of the hill is t[0:l1/2] and hill[0:l1/2]




# read data needed
outputDir='/tmp/' + username + '/'
fileName=outputDir + 'output.nc'

nc=Dataset(fileName)

time=nc['time'][:]*u/1000. # although called "time" this is actually distance in km
dt = nc['time'][1]-nc['time'][0]
z=nc['z'][:]/1000. # height in km 
q=nc['precip'][:,:,0] # precipitation rate defined on time x height grid 

nc.close()




# we need to know the height of the hill on the time grid
# create the interpolant
f_hill_interp = sci.interp1d(t[0:int(l1/2)],hill[0:int(l1/2)],fill_value='extrapolate')

# use interpolant
hill1 = np.maximum(f_hill_interp(time), 0.) # limit so can't go below ground as extrapolation sometimes does


# create a 1-d interpolant for precipitation field for every time level
f_interp_precip = []
for i in range(len(time)):
    # this is z vs precip for every time level
    f_interp_precip.append(sci.interp1d(z,q[i,:],fill_value='extrapolate')) 
    
# now we have the interpolant, use it to get the precipitation on every time level on 
# the contour of the hill
rain1=np.zeros(len(time))
for i in range(len(time)):
    rain1[i]=f_interp_precip[i](hill1[i])
    
plt.ion()
plt.subplot(211)
plt.plot(time,rain1)
plt.xlabel('distance (km)')
plt.ylabel('Precipitation rate (mm/hr)')

plt.subplot(212)
plt.plot(time,np.cumsum(rain1*dt/3600))
plt.xlabel('distance (km)')
plt.ylabel('Precipitation total (mm)')

plt.show()

