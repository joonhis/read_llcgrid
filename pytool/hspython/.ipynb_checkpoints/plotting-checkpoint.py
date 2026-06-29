from mpl_toolkits.basemap import Basemap
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt

def dimesmap(cax,grd,fld,crng,intv=1):
    fldm = ma.masked_where(np.isnan(fld),fld)
    land = fld*1;
    land[np.isfinite(land)]=0;
    land[np.isnan(land)]=1;
    coastline=land*1;
    land[np.nonzero(land==0)]=np.nan;
    landm = ma.masked_where(np.isnan(land),land);
    # setup lambert conformal basemap.
    # lat_1 is first standard parallel.
    # lat_2 is second standard parallel (defaults to lat_1).
    # lon_0,lat_0 is central point.
    # rsphere=(6378137.00,6356752.3142) specifies WGS4 ellipsoid
    # area_thresh=1000 means don't plot coastline features less
    # than 1000 km^2 in area.
    m = Basemap(width=11500000,height=5800000,
                rsphere=(6378137.00,6356752.3142),\
                resolution='c',area_thresh=1000.,projection='lcc',\
                lat_0=-60,lon_0=-90., ax=cax)
    ff=m.pcolormesh(grd.XC[::intv,::intv],grd.YC[::intv,::intv],fldm[::intv,::intv],\
                    shading='flat',vmin=crng[0],vmax=crng[1],\
                    latlon=True, rasterized=True)
    m.contourf(grd.XC[::intv,::intv],grd.YC[::intv,::intv],landm[::intv,::intv],1,\
               colors=[0.6,0.6,0.6],latlon=True, rasterized=True)
    m.contour(grd.XC,grd.YC,coastline,[0],\
              colors='black',linewidths=0.2,latlon=True)
    m.plot([-160,-160],[-75,-35],'k',linewidth=.5,latlon=True);
    m.plot([-20,-20],[-75,-35],'k',linewidth=.5,latlon=True);
    m.plot(grd.XC[0,::intv],grd.YC[0,::intv],'k',linewidth=.5,latlon=True);
    m.plot(grd.XC[-1,::intv],grd.YC[-1,::intv],'k',linewidth=.5,latlon=True);
    # draw parallels and meridians.
    m.plot([-140,-140],[-75,-35],'k',linewidth=0.1,linestyle=':',latlon=True);
    m.plot([-120,-120],[-75,-35],'k',linewidth=0.1,linestyle=':',latlon=True);
    m.plot([-100,-100],[-75,-35],'k',linewidth=0.1,linestyle=':',latlon=True);
    m.plot([-80,-80],[-75,-35],'k',linewidth=0.1,linestyle=':',latlon=True);
    m.plot([-60,-60],[-75,-35],'k',linewidth=0.1,linestyle=':',latlon=True);
    m.plot([-40,-40],[-75,-35],'k',linewidth=0.1,linestyle=':',latlon=True);
    m.plot([-20,-20],[-75,-35],'k',linewidth=0.1,linestyle=':',latlon=True);
    m.plot(grd.XC[199,::intv],grd.YC[199,::intv],'k',linewidth=0.1,linestyle=':',latlon=True);
    m.plot(grd.XC[399,::intv],grd.YC[399,::intv],'k',linewidth=0.1,linestyle=':',latlon=True);
    m.plot(grd.XC[599,::intv],grd.YC[599,::intv],'k',linewidth=0.1,linestyle=':',latlon=True);

    x, y = m(-17.5,-76)
    plt.text(x, y, r'75$^{\circ}$S',fontsize=5,ha='left',va='top',color='k')
    x, y = m(-18,-66)
    plt.text(x, y, r'65$^{\circ}$S',fontsize=5,ha='left',va='top',color='k')
    x, y = m(-18.5,-56)
    plt.text(x, y, r'55$^{\circ}$S',fontsize=5,ha='left',va='top',color='k')
    x, y = m(-19,-46)
    plt.text(x, y, r'45$^{\circ}$S',fontsize=5,ha='left',va='top',color='k')
    x, y = m(-19.5,-36)
    plt.text(x, y, r'35$^{\circ}$S',fontsize=5,ha='left',va='top',color='k')

    x, y = m(-160,-76)
    plt.text(x, y, r'160$^{\circ}$W',fontsize=5,ha='left',va='top',color='k')
    x, y = m(-140,-76)
    plt.text(x, y, r'140$^{\circ}$W',fontsize=5,ha='left',va='top',color='k')
    x, y = m(-120,-76)
    plt.text(x, y, r'120$^{\circ}$W',fontsize=5,ha='center',va='top',color='k')
    x, y = m(-100,-76)
    plt.text(x, y, r'100$^{\circ}$W',fontsize=5,ha='center',va='top',color='k')
    x, y = m(-80,-76)
    plt.text(x, y, r'80$^{\circ}$W',fontsize=5,ha='center',va='top',color='k')
    x, y = m(-60,-76)
    plt.text(x, y, r'60$^{\circ}$W',fontsize=5,ha='center',va='top',color='k')
    x, y = m(-40,-76)
    plt.text(x, y, r'40$^{\circ}$W',fontsize=5,ha='right',va='top',color='k')
    x, y = m(-20,-76)
    plt.text(x, y, r'20$^{\circ}$W',fontsize=5,ha='right',va='top',color='k')

    #plt.savefig(filename,transparent=True,bbox_inches='tight')
    return ff

def dimesmapsave(cax,grd,fld,filename,crng,intv=1,res=300):
    fldm = ma.masked_where(np.isnan(fld),fld)
    land = fld*1;
    land[np.isfinite(land)]=0;
    land[np.isnan(land)]=1;
    coastline=land*1;
    land[np.nonzero(land==0)]=np.nan;
    landm = ma.masked_where(np.isnan(land),land);
    # setup lambert conformal basemap.
    # lat_1 is first standard parallel.
    # lat_2 is second standard parallel (defaults to lat_1).
    # lon_0,lat_0 is central point.
    # rsphere=(6378137.00,6356752.3142) specifies WGS4 ellipsoid
    # area_thresh=1000 means don't plot coastline features less
    # than 1000 km^2 in area.
    m = Basemap(width=11500000,height=5800000,
                rsphere=(6378137.00,6356752.3142),\
                resolution='c',area_thresh=1000.,projection='lcc',\
                lat_0=-60,lon_0=-90., ax=cax)
    ff=m.pcolormesh(grd.XC[::intv,::intv],grd.YC[::intv,::intv],fldm[::intv,::intv],\
                    shading='flat',vmin=crng[0],vmax=crng[1],\
                    latlon=True, rasterized=True)
    m.contourf(grd.XC[::intv,::intv],grd.YC[::intv,::intv],landm[::intv,::intv],1,\
               colors=[0.6,0.6,0.6],latlon=True, rasterized=True)
    m.contour(grd.XC,grd.YC,coastline,[0],\
              colors='black',linewidths=0.2,latlon=True)
    m.plot([-160,-160],[-75,-35],'k',linewidth=.5,latlon=True);
    m.plot([-20,-20],[-75,-35],'k',linewidth=.5,latlon=True);
    m.plot(grd.XC[0,::intv],grd.YC[0,::intv],'k',linewidth=.5,latlon=True);
    m.plot(grd.XC[-1,::intv],grd.YC[-1,::intv],'k',linewidth=.5,latlon=True);
    # draw parallels and meridians.
    m.plot([-140,-140],[-75,-35],'k',linewidth=0.1,linestyle=':',latlon=True);
    m.plot([-120,-120],[-75,-35],'k',linewidth=0.1,linestyle=':',latlon=True);
    m.plot([-100,-100],[-75,-35],'k',linewidth=0.1,linestyle=':',latlon=True);
    m.plot([-80,-80],[-75,-35],'k',linewidth=0.1,linestyle=':',latlon=True);
    m.plot([-60,-60],[-75,-35],'k',linewidth=0.1,linestyle=':',latlon=True);
    m.plot([-40,-40],[-75,-35],'k',linewidth=0.1,linestyle=':',latlon=True);
    m.plot([-20,-20],[-75,-35],'k',linewidth=0.1,linestyle=':',latlon=True);
    m.plot(grd.XC[199,::intv],grd.YC[199,::intv],'k',linewidth=0.1,linestyle=':',latlon=True);
    m.plot(grd.XC[399,::intv],grd.YC[399,::intv],'k',linewidth=0.1,linestyle=':',latlon=True);
    m.plot(grd.XC[599,::intv],grd.YC[599,::intv],'k',linewidth=0.1,linestyle=':',latlon=True);

    x, y = m(-17.5,-76)
    plt.text(x, y, r'75$^{\circ}$S',fontsize=5,ha='left',va='top',color='k')
    x, y = m(-18,-66)
    plt.text(x, y, r'65$^{\circ}$S',fontsize=5,ha='left',va='top',color='k')
    x, y = m(-18.5,-56)
    plt.text(x, y, r'55$^{\circ}$S',fontsize=5,ha='left',va='top',color='k')
    x, y = m(-19,-46)
    plt.text(x, y, r'45$^{\circ}$S',fontsize=5,ha='left',va='top',color='k')
    x, y = m(-19.5,-36)
    plt.text(x, y, r'35$^{\circ}$S',fontsize=5,ha='left',va='top',color='k')

    x, y = m(-160,-76)
    plt.text(x, y, r'160$^{\circ}$W',fontsize=5,ha='left',va='top',color='k')
    x, y = m(-140,-76)
    plt.text(x, y, r'140$^{\circ}$W',fontsize=5,ha='left',va='top',color='k')
    x, y = m(-120,-76)
    plt.text(x, y, r'120$^{\circ}$W',fontsize=5,ha='center',va='top',color='k')
    x, y = m(-100,-76)
    plt.text(x, y, r'100$^{\circ}$W',fontsize=5,ha='center',va='top',color='k')
    x, y = m(-80,-76)
    plt.text(x, y, r'80$^{\circ}$W',fontsize=5,ha='center',va='top',color='k')
    x, y = m(-60,-76)
    plt.text(x, y, r'60$^{\circ}$W',fontsize=5,ha='center',va='top',color='k')
    x, y = m(-40,-76)
    plt.text(x, y, r'40$^{\circ}$W',fontsize=5,ha='right',va='top',color='k')
    x, y = m(-20,-76)
    plt.text(x, y, r'20$^{\circ}$W',fontsize=5,ha='right',va='top',color='k')

    cax.axis('off');
    plt.savefig(filename,transparent=True,dpi=res,bbox_inches='tight')

def ecplot(cax,ecx,ecy,ssh,fld,crng):
    ang=np.arange(0,2*np.pi,0.01);
    xp=np.cos(ang);
    yp=np.sin(ang);
    #
    p=cax.pcolormesh(ecx,ecy,fld,cmap='RdBu_r',vmin=crng[0],vmax=crng[-1],\
                     shading='flat',rasterized=True);
    cax.contour(ecx,ecy,ssh,np.arange(-20,0,5),colors='b',linestyles='solid')
    cax.contour(ecx,ecy,ssh,np.arange(5,21,5),colors='r',linestyles='solid')
    cax.contour(ecx,ecy,ssh,[0],colors='k',linestyles='solid')
    cax.plot(xp,yp,'k',linewidth=1,linestyle=':')
    cax.plot([0,0],[-2,2],'silver',linestyle=':');
    cax.plot([-2,2],[0,0],'silver',linestyle=':');

    #p.axes.set_xticklabels([]);p.axes.set_yticklabels([])
    p.axes.set_aspect('equal')
    p.axes.set_xlim([ecx.min(),ecx.max()]);
    p.axes.set_ylim([ecy.min(),ecy.max()]);

    return p
