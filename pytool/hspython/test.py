from MITgcmutils import rdmds
import numpy as np

def loadgrid(gridname,xrng=None,yrng=None,varname=None,flag=1):
        
    if gridname=="dimes50":
        dirGrid="/Users/hajsong/School/DIMES/DIMES_OCCA/Data/"
    elif gridname=="arctic4":
        dirGrid="/Users/hajsong/School/Arctic/4km/Data/"
    elif gridname=="arctic36":
        dirGrid="/Users/hajsong/School/Arctic/36km/Data/"
    elif gridname=="dimes1deg":
        dirGrid="/Volumes/BigLacie/DIMES_OCCA/dimes_occa_dic_1deg/"

    tmp=rdmds(dirGrid+"XC")
    [Ly,Lx]=tmp.shape
    if xrng is None:
        xrng=np.arange(Lx)
    if yrng is None:
        yrng=np.arange(Ly)
    if varname is None:
        varname=['XC','YC','RAC','DXC','DYC','hFacC','hFacW','hFacS','Depth',\
                 'RC','RF','DRC','DRF'];
        if flag==2
            varname=np.append(varname,['XG','YG','RAZ','DXG','DYG'])

    class grd(object):
        for iv in xrange(len(varname)):
            exec('tmpvar=rdmds(dirGrid+"'+varname[iv]+'")');
            tmpvar=tmpvar.squeeze();
            dim=len(tmpvar.shape)
            if dim==3:
                exec(varname[iv]+'=tmpvar[:,yrng][:,:,xrng]')
            elif dim==2:
                exec(varname[iv]+'=tmpvar[yrng][:,xrng]')
            if varname[iv]=='hFacC':
                mskC=hFacC.copy()
                mskC[mskC==0]=np.nan
                mskC[np.isfinite(mskC)]=1.
                grd.mskC=mskC;
            if varname[iv]=='hFacS':
                mskW=hFacW.copy()
                mskW[mskW==0]=np.nan
                mskW[np.isfinite(mskW)]=1. 
                grd.mskS=mskS;
            if varname[iv]=='hFacW':
                mskS=hFacS.copy()
                mskS[mskS==0]=np.nan
                mskS[np.isfinite(mskS)]=1.
                grd.mskW=mskW;

    del grd.iv
    del grd.mro
    return grd
