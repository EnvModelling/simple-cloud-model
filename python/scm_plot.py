# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 22:35:33 2020

@author: mccikpc2
"""

from netCDF4 import Dataset
import numpy as np
from matplotlib import rc
#rc('font',family='serif')
#rc('text',usetex = True)
from svp import svp as svp
import matplotlib as mpl
mpl.use('agg')

import os
import getpass

import matplotlib.pyplot as plt
#import pylab as plt

username=getpass.getuser()

def plot_model_run():
    u=2.5 # m/s
    # create a patch
    tau=6000.
    rt=tau
    t=np.linspace(0,rt,100)
    hill=tau*0.6/(2.*np.pi)*(1.-np.cos(2.*np.pi/tau*t))/1000.
    t=(t+10)*u/1000.
    t=np.append(t,t[::-1])
    hill=np.append(hill,np.zeros((len(hill),1)))
    l1=len(t)
    pgon=np.zeros((l1,2))
    pgon[:,0]=t
    pgon[:,1]=hill    
    pgon1=plt.Polygon(pgon,color='r',alpha=1)
    pgon2=plt.Polygon(pgon,color='r',alpha=1)
    pgon3=plt.Polygon(pgon,color='r',alpha=1)
    pgon4=plt.Polygon(pgon,color='r',alpha=1)
    
    outputDir='/tmp/' + username + '/'
    fileName=outputDir + 'output.nc'
    
    nc=Dataset(fileName)
    
    hillgrid=tau*0.6/(2.*np.pi)*(1.-np.cos(2.*np.pi/tau*nc['time'][:]))/1000.
    time=nc['time'][:]*u/1000.*60.
    
    z=nc['z'][:]
    q=nc['q'][:,:,:]

    max1=0.
    for i in range(len(time)):
        ind,=np.where(z/1000>=hillgrid[i])
        max1=max(max1,max(q[i,ind,23]*1000))


    
    m1=np.max(q[0,:,14]/1.e6)
    #plt.ion()
    fig=plt.figure(figsize=(12,6))     
    ax=plt.subplot(221)
    plt.pcolor(time/60,z/1000.,q[:,:,14].T/1.e6)
    plt.xlabel('distance (km)')
    plt.ylabel('z (km)')
    plt.clim((0,m1))
    plt.text(0.1,0.9,'(a) CDNC',color='white',transform=ax.transAxes)
    cbar=plt.colorbar()
    cbar.set_label('number of cloud drops (cm$^{-3}$)')
    ax.add_patch(pgon1)
    
    ax=plt.subplot(222)
    plt.pcolor(time/60,z/1000.,q[:,:,15].T*1000.)
    plt.xlabel('distance (km)')
    plt.ylabel('z (km)')
    plt.text(0.1,0.9,'(b) $q_c$',color='white',transform=ax.transAxes)
    cbar=plt.colorbar()
    cbar.set_label('mass of cloud drops (g kg$^{-1}$)')
    ax.add_patch(pgon2)
    
    ax=plt.subplot(223)
    plt.pcolor(time/60,z/1000.,q[:,:,23].T*1000.)
    plt.clim((0,max1))
    plt.xlabel('distance (km)')
    plt.ylabel('z (km)')
    plt.text(0.1,0.9,'(c) $q_r$',color='white',transform=ax.transAxes)
    cbar=plt.colorbar()
    cbar.set_label('mass of rain drops (g kg$^{-1}$)')
    ax.add_patch(pgon3)
    
    ax=plt.subplot(224)
    plt.pcolor(time/60,z/1000.,q[:,:,30].T/1000.)
#     plt.plot(time/60.,hillgrid)
    plt.xlabel('distance (km)')
    plt.ylabel('z k(m)')
    plt.text(0.1,0.9,'(d) $N_{ice}$',color='white',transform=ax.transAxes)
    cbar=plt.colorbar()
    cbar.set_label('number of ice crystals (L$^{-1}$)')
    ax.add_patch(pgon4)
    
    
#     fig2=plt.figure()
#     plt.plot(time/60,np.sum(nc['precip'][:,0,0:2],axis=1))
#     plt.xlabel('time (mins)')
#     plt.ylabel('Precipitation (mm hr$^{-1}$)')
    
    
    nc.close()
    #plt.show()
    
    fig.savefig('/tmp/' + username + '/scm_plot.png',format='png')
    
if __name__=="__main__":
    plot_model_run()
