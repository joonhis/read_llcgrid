import numpy as np
import numpy.ma as ma

def smooth_yz(grd,fld,xdist,ydist):
    """
    object:    spatial smoothing with a given scale in both x and y direction
    inputs:    'grd' is the grid structure variable.
               'fld' is the field to be smoothed.
               'xdist' and 'ydist' is the distance (km) from the center 
                   that will be averaged.
    output:    'sfld' is the smoothed field.
    """
    
    ndim=len(fld.shape);
    if ndim==2:
        [ly,lx]=fld.shape;
        sfld=np.zeros((ly,lx))
    elif ndim==3:
        [lz,ly,lx]=fld.shape;
        sfld=np.zeros((lz,ly,lx))

    for iy in range(0,ly):
        dely = int(np.ceil(ydist/(grd.DYC[iy,0]/1000)));
        delx = int(np.ceil(xdist/(grd.DXC[iy,0]/1000)));
        yrng=range(iy-dely,iy+dely+1)
        #
        #  extract data to average in y direction
        #
        dy=min(iy,dely,ly-iy-1)
        avgrng=np.arange(iy-dy,iy+dy+1)
        if ndim==2:
            mdata=fld[avgrng,:]
            #
            #  Now, loop in x-direction, fill nan and compute the mean
            #
            tmp=[np.nanmean(ma.array(mdata[:,np.arange(\
                 ix-min(ix,delx,lx-ix-1),ix+min(ix,delx,lx-ix-1)+1)],\
                 mask=np.isnan(mdata[:,np.arange(\
                 ix-min(ix,delx,lx-ix-1),ix+min(ix,delx,lx-ix-1)+1)])).\
                 filled(fld[iy,ix])) for ix in range(0,lx)]
            sfld[iy,:]=np.asarray(tmp)
        elif ndim==3:
            for ik in range(0,lz):
                mdata=fld[ik,avgrng,:]
                #
                #  Now, loop in x-direction, fill nan and compute the mean
                #
                tmp=[np.nanmean(ma.array(mdata[:,np.arange(\
                     ix-min(ix,delx,lx-ix-1),ix+min(ix,delx,lx-ix-1)+1)],\
                     mask=np.isnan(mdata[:,np.arange(\
                     ix-min(ix,delx,lx-ix-1),ix+min(ix,delx,lx-ix-1)+1)])).\
                     filled(fld[ik,iy,ix])) for ix in range(0,lx)]
                sfld[ik,iy,:]=np.asarray(tmp)

    return sfld

def smooth_yz_new(grd,fld,xdist,ydist):
    """
    object:    spatial smoothing with a given scale in both x and y direction
    inputs:    'grd' is the grid structure variable.
               'fld' is the field to be smoothed.
               'xdist' and 'ydist' is the distance (km) from the center 
                   that will be averaged.
    output:    'sfld' is the smoothed field.
    """
    
    ndim=len(fld.shape);
    if ndim==2:
        [ly,lx]=fld.shape;
        sfld=np.zeros((ly,lx))
    elif ndim==3:
        [lz,ly,lx]=fld.shape;
        sfld=np.zeros((lz,ly,lx))

    for iy in range(0,ly):
        dely = int(np.ceil(ydist/(grd.DYC[iy,0]/1000)));
        delx = int(np.ceil(xdist/(grd.DXC[iy,0]/1000)));
        yrng=range(iy-dely,iy+dely+1)
        #
        #  extract data to average in y direction
        #
        dy=min(iy,dely,ly-iy-1)
        avgrng=np.arange(iy-dy,iy+dy+1)
        if ndim==2:
            mdata=fld[avgrng,:]
            #
            #  Now, loop in x-direction, fill nan and compute the mean
            #
            tmp=[np.nanmean(ma.array(mdata[:,np.arange(\
                 ix-min(ix,delx,lx-ix-1),ix+min(ix,delx,lx-ix-1)+1)],\
                 mask=np.isnan(mdata[:,np.arange(\
                 ix-min(ix,delx,lx-ix-1),ix+min(ix,delx,lx-ix-1)+1)])).\
                 filled(fld[iy,ix]),axis=1) for ix in range(0,lx)]
            sfld[iy,:]=np.mean(np.asarray(tmp).T,axis=0)
        elif ndim==3:
            tmp3d=np.zeros((lz,len(avgrng),lx));
            for ik in range(0,lz):
                mdata=fld[ik,avgrng,:]
                #
                #  Now, loop in x-direction, fill nan and compute the mean
                #
                tmp=[np.nanmean(ma.array(mdata[:,np.arange(\
                     ix-min(ix,delx,lx-ix-1),ix+min(ix,delx,lx-ix-1)+1)],\
                     mask=np.isnan(mdata[:,np.arange(\
                     ix-min(ix,delx,lx-ix-1),ix+min(ix,delx,lx-ix-1)+1)])).\
                     filled(fld[ik,iy,ix]),axis=1) for ix in range(0,lx)]
                tmp3d[ik,...]=np.asarray(tmp).T
            sfld=np.mean(tmp3d,axis=1)

    return sfld

