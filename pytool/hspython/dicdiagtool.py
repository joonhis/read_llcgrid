import numpy as np

def advection(grd,uflux,vflux,wflux,ik=[0],orient='h'):
    lz=len(ik)
    [ly,lx]=grd.XC.shape
    grdvol=np.tile(grd.RAC,[lz,1,1])*np.tile(grd.DRF[ik],[1,ly,lx])

    if orient=='h':
        tmp=np.diff(uflux[ik,:,:],axis=2)
        dufluxdx=np.append(tmp,np.zeros((lz,ly,1)),axis=2)
        tmp=np.diff(vflux[ik,...],axis=1)
        dvfluxdy=np.append(tmp,np.zeros((lz,1,lx)),axis=1)
        div=-(dufluxdx+dvfluxdy)/grdvol
    elif orient=='v':
        dwfluxdz=np.squeeze(np.diff(wflux,axis=0))
        # negative sign is cancelled out because it is positive upward
        div=(dwfluxdz[ik,:,:])/grdvol

    div=div*np.tile(grd.DRF[ik],[1,ly,lx])
    return div

def hadvection(grd,uflux,vflux,ik=[0]):
    lz=len(ik)
    [ly,lx]=grd.XC.shape
    if lz==1:
        grdvol=grd.RAC*grd.DRF[ik]
        tmp=np.diff(uflux,axis=1)
        dufluxdx=np.append(tmp,np.zeros((ly,1)),axis=1)
        tmp=np.diff(vflux,axis=0)
        dvfluxdy=np.append(tmp,np.zeros((1,lx)),axis=0)
        div=-(dufluxdx+dvfluxdy)/grdvol
        div=div*grd.DRF[ik]
    else:
        [ly,lx]=grd.XC.shape
        grdvol=np.tile(grd.RAC,[lz,1,1])*np.tile(grd.DRF[ik],[1,ly,lx])
        tmp=np.diff(uflux[ik,:,:],axis=2)
        dufluxdx=np.append(tmp,np.zeros((lz,ly,1)),axis=2)
        tmp=np.diff(vflux[ik,...],axis=1)
        dvfluxdy=np.append(tmp,np.zeros((lz,1,lx)),axis=1)
        div=-(dufluxdx+dvfluxdy)/grdvol
        div=div*np.tile(grd.DRF[ik],[1,ly,lx])
    return div


def vadvection(grd,wflux,ik=[0]):
    [ly,lx]=grd.XC.shape
    dwfluxdz=np.diff(wflux,axis=0)
    # negative sign is cancelled out because it is positive upward
    div=np.squeeze(dwfluxdz[ik,:,:])/grd.RAC
    return div
