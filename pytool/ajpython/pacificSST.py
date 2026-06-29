#!/usr/bin/env python
# coding: utf-8

# In[1]:


import glob
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import cartopy as cart
import cartopy.crs as ccrs
from netCDF4 import Dataset
# from mpl_toolkits.basemap import Basemap
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter

def loadgrid(gridname, varname=None, pacific=False):
    """
    Load grid from each model/observation data for the Pacific ocean.(LON[99.5~290.5],LAT[-15.5~66.5])
    gridname = ['ARGO_IPRC', 'ARGO_MOAA', 'HadelySST', 'NEMO_ORCA2_1m', 'glosea5', 'Levitus', 'SIO']
    """
    if gridname=='ARGO_IPRC':
        dirGrid='/home/share/ARGO/IPRC/'
    elif gridname == 'ARGO_MOAA':
        dirGrid='/home/share/ARGO/MOAA/TS_'
    elif gridname == 'HadleySST':
        dirGrid='/home/share/HadleySST/'
    elif gridname == 'NEMO_ORCA2_1m':
        dirGrid='/home/ajin05/NEMO/release-4.0/cfgs/global2deg/EXP00/ORCA2_1m/100yrs/'
    elif gridname == 'Levitus':
        dirGrid = '/home/share/EN4.2.1/Levitus/'
    elif gridname == 'SIO':
        dirGrid = '/home/share/ARGO/SIO/'
    elif gridname=='glosea5':
        dirGrid='/home/ajin05/NEMO/gs5_hcst_gc20/data/19910101/'
        
    if varname is None:
        varname=['XC','YC'];
    
    if  gridname == 'glosea5':
        files = sorted(glob.glob(dirGrid+'*seatmp.nc'))
    else:
        files = sorted(glob.glob(dirGrid+'*.nc'))
    
    f = Dataset(files[0], 'r')
        
    class grd(object):
        if pacific==False:
            if gridname=='ARGO_IPRC':
                XC = f.variables['LONGITUDE'][:].data 
                YC = f.variables['LATITUDE'][:].data
                ZC = f.variables['LEVEL'][:].data
            elif gridname=='ARGO_MOAA':
                XC = f.variables['LONGITUDE'][:]
                YC = f.variables['LATITUDE'][:]
            elif gridname=='HadleySST':
                XC = f.variables['longitude'][:].data
                YC = f.variables['latitude'][:].data 
            elif gridname=='NEMO_ORCA2_1m':
                XC=f.variables['nav_lon'][:]
                YC=f.variables['nav_lat'][:]
                ZC=f.variables['deptht'][:].data
            elif gridname == 'Levitus':
                XC = f.variables['lon'][:].data
                YC = f.variables['lat'][:].data
                ZC = f.variables['depth_bnds'][:][:,0].data
            elif gridname == 'SIO':
                XC = f.variables['LONGITUDE'][:].data
                YC = f.variables['LATITUDE'][:].data
                ZC = f.variables['PRESSURE'][:].data
            elif gridname=='glosea5':
                XC = f.variables['x'][:].data
                YC = f.variables['y'][:].data
                ZC = f.variables['deptht'][:].data
        elif pacific==True:
            if gridname=='ARGO_IPRC':
                XC = f.variables['LONGITUDE'][99:291].data # lon = 99.5~290.5 deg
                YC = f.variables['LATITUDE'][74:157].data # len = -15.5~66.5 deg
                ZC = f.variables['LEVEL'][:].data
            elif gridname=='ARGO_MOAA':
                XC = f.variables['LONGITUDE'][:]
                YC = f.variables['LATITUDE'][:]
            elif gridname=='HadleySST':
                XC = np.concatenate([f.variables['longitude'][279:].data,360+f.variables['longitude'][:111].data]) # lon = 99.5~290.5 deg
                YC = f.variables['latitude'][23:-74].data # len = -15.5~66.5 deg
            elif gridname=='NEMO_ORCA2_1m':
                XC=f.variables['nav_lon'][:]
                YC=f.variables['nav_lat'][:]
                ZC=f.variables['deptht'][:].data
            elif gridname == 'Levitus':
                XC = f.variables['lon'][99:290].data # lon=100~290 deg
                YC = f.variables['lat'][68:150].data # lat=-15~66 deg
                ZC = f.variables['depth_bnds'][:][:,0].data
            elif gridname == 'SIO':
                XC = f.variables['LONGITUDE'][79:271].data # lon=99.5~290.5 deg
                YC = f.variables['LATITUDE'][49:132].data # lat=-15.5~66.5 deg
                ZC = f.variables['PRESSURE'][:].data
            elif gridname=='glosea5':
                XC = np.concatenate([f.variables['x'][1120:-1].data, 360+f.variables['x'][:441].data])
                YC = f.variables['y'][75*4:156*4+1].data
                ZC = f.variables['deptht'][:].data            
            
    return grd

def rmSeason(data,nt,ny,nx):
    """
    Remove seasonal cycle from given data.
    data.shape = (nt,ny,nx)
    """
    monMean = np.ma.empty([12,ny,nx])
    anmData = np.ma.empty([nt,ny,nx]) # SST anomaly
    for i in range(12):
        monMean[i,:,:] = np.average(data[i::12,:,:], axis=0)
        anmData[i::12,:,:] = data[i::12,:,:] - monMean[i,:,:]
    return anmData

def seasonalCycle(data,nt,ny,nx):
    """
    Get seasonal cycle from given data.
    data.shape = (nt,ny,nx)
    """
    monMean = np.ma.empty([12,ny,nx])
    for i in range(12):
        monMean[i,:,:] = np.average(data[i::12,:,:], axis=0)
    return monMean

def plot2D(data,XC,YC,cmap='jet',vmin=None,vmax=None,title=None):
    plt.figure(figsize=(10,8), facecolor='w')
    m = Basemap(projection='cea', llcrnrlat=-15,urcrnrlat=66,llcrnrlon=100,urcrnrlon=290, resolution='c')
    m.drawcoastlines()
    m.drawparallels(np.arange(0,61,30),labels=[1,0,0,0],fontsize=14) #present label in the left side
    m.drawmeridians(np.arange(120,281,40),labels=[0,0,0,1],fontsize=14) # present label in the bottom
    
    if YC.shape[0]%2==1:
        i=-1
    else :i=None
    x, y = np.meshgrid(XC, YC[:i])
    m.pcolormesh(x,y,data[:i,:],latlon=True, cmap=cmap, vmin=vmin, vmax=vmax)
    m.colorbar()
    plt.title(title, fontsize=16)
    
def plot2Dccrs(data,XC,YC,projection=ccrs.PlateCarree,central_longitude=195,pacific=True,
               cmap='jet',vmin=None,vmax=None,title=None, barLoc='right',shrink=.45,label=True, unit='', colorbar=True):
    plt.figure(figsize=(10,8), facecolor='w')
    m = plt.gca(projection=projection(central_longitude=central_longitude))
    if pacific == True:
        #m.set_extent([120,270,-15,60])
        m.set_extent([100,290,-15,66], crs=ccrs.PlateCarree())
    else:
        pass
    m.coastlines()
    cm = m.pcolormesh(XC, YC, data, transform=ccrs.PlateCarree(), cmap=cmap, vmin=vmin, vmax=vmax)
    m.add_feature(cart.feature.LAND, color='grey', zorder=100, edgecolor='k')
    
    if label==True:
        m.gridlines(xlocs=[100,120,160,200,240,280,290],ylocs=[-15,0,30,60,66])
        m.set_xticks(np.arange(120,281, 40), crs=ccrs.PlateCarree())
        m.set_yticks(np.arange(0,61,30), crs=ccrs.PlateCarree())
        m.xaxis.set_major_formatter(LongitudeFormatter(zero_direction_label=True))
        m.yaxis.set_major_formatter(LatitudeFormatter())
        m.tick_params('both', labelsize=15)
        
    if colorbar==True:
        cax, kw = matplotlib.colorbar.make_axes(m, location=barLoc, pad=0.05, shrink=shrink)
        b = plt.colorbar(cm, cax=cax, **kw)
        b.set_label(unit, fontsize=16, labelpad=-10, y=1.13, rotation=0)
        b.ax.tick_params(labelsize=16)
        b.ax.yaxis.get_offset_text().set(size=16)
    m.set_title(title, fontsize=16, fontweight = 'bold')

