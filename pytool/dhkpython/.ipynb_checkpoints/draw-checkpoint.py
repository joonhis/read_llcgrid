from dhkpython import draw
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import datetime
from dateutil import relativedelta
import netCDF4 as nc; from netCDF4 import Dataset 
import matplotlib.patches as patches 
import seaborn as sns ; sns.set(style='white')
import glob
import cartopy as cart
import cartopy.crs as ccrs
# from mpl_toolkits.basemap import Basemap
#from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter

''' 1d plot '''
def scatter(X,Y,fsize=(6,6),title=None,xlabel=None,ylabel=None,axhline=None,axvline=None,\
            label=False,annotate=False,time=None,rng=None,dpi=100):
    fm = plt.figure(figsize=fsize,dpi=dpi, facecolor='w')
    m = plt.gca()
    cm = m.scatter(X,Y)
    

    if label==True:
        m.set_xlabel(xlabel,size=15)
        m.set_ylabel(ylabel,size=15)
        m.set_title(title,size=18)
        m.set_xlim(-rng,rng)
        m.set_ylim(-rng,rng)
        plt.axhline(axhline)
        plt.axvline(axvline)
    if annotate==True:
        for i in range(len(time)):
            m.annotate(time[i],(X[i],Y[i]))
    return m

def scatter_yearcircle(X,Y,fsize=(6,6),title=None,xlabel=None,ylabel=None,axhline=None,axvline=None,\
            label=False,annotate=False,time=None,rng=None,dpi=100):
    fm = plt.figure(figsize=fsize,dpi=dpi, facecolor='w')
    m = plt.gca()
    cm = m.scatter(X,Y,s=300,alpha=1,c='white',edgecolors='b',linewidth=1.5)

    if label==True:
        m.set_xlabel(xlabel,size=15)
        m.set_ylabel(ylabel,size=15)
        m.set_title(title,size=18)
        m.set_xlim(-rng,rng)
        m.set_ylim(-rng,rng)
        plt.axhline(axhline,c='gray')
        plt.axvline(axvline,c='gray')
    if annotate==True:
        for i in range(len(time)):
            m.annotate(time[i],(X[i],Y[i]),ha='center',va='center')
# 빨간 선으로 칠하고 싶을때는 직접 입력.
            #plt.scatter(sst_std[39],sic_std[39],s=300,c='white',edgecolors='r',linewidth=1.5)
    return m

def timeseries(X,data,fsize=(10,4),title=None,xlabel=None,ylabel=None,axhline=None,axvline=None,scatter=None,grid=True,dpi=100,c=None,alpha=1):
    plt.figure(figsize=fsize,dpi=dpi, facecolor='w')
    m = plt.gca()
    cm = m.plot(X,data,c=c,alpha=alpha)
  ##  cm.set_axhline(axhline)
   # cm.set_avhline(axvline)
    m.set_title(title,size=18)
    m.set_xlabel(xlabel,size=15)
    m.set_ylabel(ylabel,size=15)
    m.grid(grid)
    if scatter!=None:
        m.scatter(X,data,c='gray',marker='o',alpha=0.5)
    return m

        
''' 2d plot '''        
def imshow(data,cmap='viridis',fsize=(8,5),vmin=None,vmax=None,title=None,xlabel=None,ylabel=None,\
           colorbar=True,box=None,orig='lower',dpi=100): # origin
    fm = plt.figure(figsize=fsize,dpi=dpi, facecolor='w')
    m = plt.gca()
    cm = m.imshow(data,cmap=cmap, vmin=vmin, vmax=vmax,origin=orig)
    if colorbar==True:
        cax, kw = matplotlib.colorbar.make_axes(m)
        b = plt.colorbar(cm, cax=cax, **kw,shrink=0.5)
    m.set_title(title,size=18)
    m.set_xlabel(xlabel,size=15)
    m.set_ylabel(ylabel,size=15)
    if box!=None:
        m.add_patch(patches.Rectangle((box[0],box[2]),box[1]-box[0],box[3]-box[2],fill=False,edgecolor='black'))
    return m

# 그래프. rectangle을 imshow와 cartopy에 둘다 그려주고 싶은데 좌표축이 안맞는다.

# 2열 : figsize는 2 대 1!!
def pcolorGridspec(XC,YC,data,projection=ccrs.PlateCarree,central_longitude=180,cmap='RdBu_r',form=None\
,vmin=None,vmax=None,title=None, barLoc='bottom',shrink=.8,colorbar=True,box=None,resolution=0.25,grid=None,ij=[1,1]):
    m = plt.subplot(grid[ij[0],ij[1]],projection=projection(central_longitude=central_longitude))
    m.coastlines()
    cm = m.pcolormesh(XC, YC, data, transform=ccrs.PlateCarree(), cmap=cmap, vmin=vmin, vmax=vmax)
    m.add_feature(cart.feature.LAND, color='grey', zorder=100, edgecolor='k')
        
    if colorbar==True:
        cax, kw = matplotlib.colorbar.make_axes(m, location=barLoc, pad=0.02, shrink=shrink,aspect=40)
        b = plt.colorbar(cm, cax=cax, **kw,format=form)
#        b.set_label(unit, fontsize=16, labelpad=-10, y=1.13, rotation=0)
        b.ax.tick_params(labelsize=16)
        b.formatter.set_powerlimits((-4, 4)) ### 지수 단위 설정
        b.ax.yaxis.get_offset_text().set(size=16)
    m.set_title(title, fontsize=17, fontweight = 'bold')
    if box!=None:
        BOX1,BOX2,BOX3,BOX4 = box[0]*resolution-180,box[1]*resolution-180,-box[2]*resolution+90,-box[3]*resolution+90
        m.add_patch(patches.Rectangle((BOX1,BOX3),BOX2-BOX1,BOX4-BOX3,fill=False,edgecolor='black'))
    return cm,m

def contourfGridspec(data,projection=ccrs.PlateCarree,central_longitude=180,cmap='RdBu_r'\
,vmin=None,vmax=None,title=None, barLoc='bottom',shrink=.8,colorbar=True,box=None,resolution=0.25,grid=None,ij=[1,1]):
    m = plt.subplot(grid[ij[0],ij[1]],projection=projection(central_longitude=central_longitude))
    m.coastlines()
    cm = m.contourf(data, transform=ccrs.PlateCarree(), cmap=cmap, vmin=vmin, vmax=vmax)
    m.add_feature(cart.feature.LAND, color='grey', zorder=100, edgecolor='k')
        
    if colorbar==True:
        cax, kw = matplotlib.colorbar.make_axes(m, location=barLoc, pad=0.02, shrink=shrink,aspect=40)
        b = plt.colorbar(cm, cax=cax, **kw)
#        b.set_label(unit, fontsize=16, labelpad=-10, y=1.13, rotation=0)
        b.ax.tick_params(labelsize=16)
        b.ax.yaxis.get_offset_text().set(size=16)
    m.set_title(title, fontsize=17, fontweight = 'bold')
    if box!=None:
        BOX1,BOX2,BOX3,BOX4 = box[0]*resolution-180,box[1]*resolution-180,-box[2]*resolution+90,-box[3]*resolution+90
        m.add_patch(patches.Rectangle((BOX1,BOX3),BOX2-BOX1,BOX4-BOX3,fill=False,edgecolor='black'))
    return cm,m

def timeseriesGridspec(X,data,year_str,title=None,grid=None,ij=[1,1],c='b',label=False,gridline=True,scatter=None,tick=2,\
                      name=None,axhline=None,axvline=None,multi=False):
    m = plt.subplot(grid[ij[0],ij[1]])
    if multi == False:
        cm = m.plot(X,data,c=c,alpha=1,label=str(name))
    elif multi == True:
        cm = m.plot(X,data,c='r',alpha=0.7,label=str(name))
        
    #m.axhline(axhline) ; m.avhline(axvline)
    if label==True:
        m.set_title(title,size=17);m.set_xlabel(xlabel,size=15);m.set_ylabel(ylabel,size=15)
    m.legend(loc=2,fontsize=12)
    if axhline!=None:
        m.axhline(axhline,c='k',alpha=0.5)
    if axvline!=None:
        m.axvline(axvline[0],c='r',alpha=0.5)
        m.axvline(axvline[1],c='r',alpha=0.5)
           
    m.grid(gridline)
    plt.xticks(X[12::12*tick],year_str[1::1*tick],size=13);plt.yticks(size=13)
    
    if scatter!=None:
        m.scatter(X,data,c='gray',marker='o',alpha=0.5)
    return cm,m  

# 실용성 없다.
def timeseriesmultiY(X,data,year_str,title=None,grid=None,ij=[1,1],c='b',label=False,gridline=True,scatter=None,tick=2,\
                      name=None,axhline=None,axvline=None,multi=False):
    m = plt.subplot(grid[ij[0],ij[1]])
    if multi == False:
        cm = m.plot(X,data,c=c,alpha=1,label=str(name))
    elif multi == True:
        cm = m.plot(X,data[0],c='r',ls='--',alpha=0.7,label=str(name[0]))
        cm = m.plot(X,data[1],c='g',ls=':',alpha=0.7,label=str(name[1]))
        cm = m.plot(X,data[2],c='k',alpha=0.7,label=str(name[2]))
        
    #m.axhline(axhline) ; m.avhline(axvline)
    if label==True:
        m.set_title(title,size=17);m.set_xlabel(xlabel,size=15);m.set_ylabel(ylabel,size=15)
    m.legend(loc=2,fontsize=12)
    if axhline!=None:
        m.axhline(axhline,c='k',alpha=0.5)
    if axvline!=None:
        m.axvline(axvline[0],c='r',alpha=0.5)
        m.axvline(axvline[1],c='r',alpha=0.5)
           
    m.grid(gridline)
    plt.xticks(X[12::12*tick],year_str[1::1*tick],size=13);plt.yticks(size=13)
    
    if scatter!=None:
        m.scatter(X,data,c='gray',marker='o',alpha=0.5)
    return cm,m  



def pcolorccrs(XC,YC,data,fsize=(10,8),projection=ccrs.PlateCarree,central_longitude=180,cmap='coolwarm'\
,vmin=None,vmax=None,title=None, barLoc='right',shrink=.35,colorbar=True,box=None,dpi=100):
    plt.figure(figsize=fsize, dpi=dpi,facecolor='w')
    m = plt.gca(projection=projection(central_longitude=central_longitude))
    m.coastlines()
    cm = m.pcolormesh(XC, YC, data, transform=ccrs.PlateCarree(), cmap=cmap, vmin=vmin, vmax=vmax)
    m.add_feature(cart.feature.LAND, color='grey', zorder=100, edgecolor='k')
 
    if colorbar==True:
        cax, kw = matplotlib.colorbar.make_axes(m, location=barLoc, pad=0.05, shrink=shrink)
        b = plt.colorbar(cm, cax=cax, **kw)
#        b.set_label(unit, fontsize=16, labelpad=-10, y=1.13, rotation=0)
        b.ax.tick_params(labelsize=16)
        b.ax.yaxis.get_offset_text().set(size=16)
    m.set_title(title, fontsize=16, fontweight = 'bold')
    if box!=None:
        BOX1,BOX2,BOX3,BOX4 = box[0]*resolution-180,box[1]*resolution-180,-box[2]*resolution+90,-box[3]*resolution+90
        m.add_patch(patches.Rectangle((BOX1,BOX3),BOX2-BOX1,BOX4-BOX3,fill=False,edgecolor='black'))
    return cm,m


def corrmap(XC,YC,data,fsize=(10,8),projection=ccrs.PlateCarree,central_longitude=180,cmap='RdBu_r'\
,vvalue=1, barLoc='right',shrink=.45,label=True, unit='', colorbar=True,box=None,dpi=100):
    plt.figure(figsize=fsize, dpi=dpi,facecolor='w')
    m = plt.gca(projection=projection(central_longitude=central_longitude))
    m.coastlines()
    cm = m.pcolormesh(XC, YC, data, transform=ccrs.PlateCarree(), cmap=cmap, vmin=-vvalue, vmax=vvalue)
    m.add_feature(cart.feature.LAND, zorder=100, edgecolor='k', color='grey')
    
    if label==True:
        #m.gridlines(xlocs=[100,120,160,200,240,280,290],ylocs=[-15,0,30,60,66])
        #m.set_xticks(np.arange(0,179, 50), crs=ccrs.PlateCarree())
        #m.set_yticks(np.arange(0,61,30), crs=ccrs.PlateCarree())
        #m.xaxis.set_major_formatter(LongitudeFormatter(zero_direction_label=True))
        #m.yaxis.set_major_formatter(LatitudeFormatter())
        m.tick_params('both', labelsize=15)
        #m.set_extent([-50,30,-30,30])
        
    if colorbar==True:
        cax, kw = matplotlib.colorbar.make_axes(m, location=barLoc, pad=0.05, shrink=shrink)
        b = plt.colorbar(cm, cax=cax, **kw)
        b.set_label(unit, fontsize=16, labelpad=-10, y=1.13, rotation=0)
        b.ax.tick_params(labelsize=16)
        b.ax.yaxis.get_offset_text().set(size=16)
    m.set_title(title, fontsize=16, fontweight = 'bold')
    if box!=None:
        m.add_patch(patches.Rectangle((box[0],box[2]),box[1]-box[0],box[3]-box[2],fill=False,edgecolor='black'))
    return m

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
    
    
    
    