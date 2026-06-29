from MITgcmutils import rdmds
import numpy as np

def loadgrid_all(gridname):
        
    if gridname=="dimes50":
        dirGrid="/net/fs09/d0/hajsong/DIMES_OCCA_DIC/grid/"
    elif gridname=="arctic4":
        dirGrid="/net/fs09/d0/hajsong/Arctic/4km/Data/"
    elif gridname=="arctic36":
        #dirGrid="/Users/hajsong/School/Arctic/36km/Data/"
        dirGrid="/d1/jscott/ARCTIC/run_79-13_combo6_newdiags_layersT/"
    elif gridname=="dimes1deg":
        dirGrid="/Volumes/BigLacie/DIMES_OCCA/dimes_occa_dic_1deg/"
    elif gridname=="so_box":
        dirGrid="/home/hajsong/so_box_biogeo/run/"
    elif gridname=="speedy":
        dirGrid="/home/hajsong/Research/CPL_MITgcm_Speedy/speedy_box/run/"

    class grd(object):
        XC=rdmds(dirGrid+"XC")
        YC=rdmds(dirGrid+"YC")
        RAC=rdmds(dirGrid+"RAC")
        DXC=rdmds(dirGrid+"DXC")
        DYC=rdmds(dirGrid+"DYC")
        hFacC=rdmds(dirGrid+"hFacC")
        hFacW=rdmds(dirGrid+"hFacW")
        hFacS=rdmds(dirGrid+"hFacS")
        Depth=rdmds(dirGrid+"Depth")
        RC=rdmds(dirGrid+"RC")
        RF=rdmds(dirGrid+"RF")
        DRC=rdmds(dirGrid+"DRC")
        DRF=rdmds(dirGrid+"DRF")
        XG=rdmds(dirGrid+"XG")
        YG=rdmds(dirGrid+"YG")
        RAZ=rdmds(dirGrid+"RAZ")
        DXG=rdmds(dirGrid+"DXG")
        DYG=rdmds(dirGrid+"DYG")

    mskC=grd.hFacC.copy()
    mskC[mskC==0]=np.nan
    mskC[np.isfinite(mskC)]=1.
    mskW=grd.hFacW.copy()
    mskW[mskW==0]=np.nan
    mskW[np.isfinite(mskW)]=1.
    mskS=grd.hFacS.copy()
    mskS[mskS==0]=np.nan
    mskS[np.isfinite(mskS)]=1.
    grd.mskC=mskC;
    grd.mskW=mskW;
    grd.mskS=mskS;

    return grd

def loadgrid(gridname,region=None,varname=None,flag=1):

    if gridname=="dimes50":
        dirGrid="/home/hajsong/Data/dimes_occa_dic/"
    elif gridname=="arctic4":
        dirGrid="/net/fs09/d0/hajsong/Arctic/4km/Data/"
    elif gridname=="arctic36":
        #dirGrid="/Users/hajsong/School/Arctic/36km/Data/"
        dirGrid="/d1/jscott/ARCTIC/run_79-13_combo6_newdiags_layersT/"
    elif gridname=="dimes1deg":
        dirGrid="/Volumes/BigLacie/DIMES_OCCA/dimes_occa_dic_1deg/"
    elif gridname=="so_box":
        dirGrid="/home/hajsong/so_box_biogeo/run/"
    elif gridname=="sochannel":
        dirGrid="/home/hajsong/Research/SOchannel/Exp2d/run/"
    elif gridname=="sochannel3":
        dirGrid="/net/fs09/d0/siciak/so_chan_3d/run_1/res_t07/grid/"
    elif gridname=='soch1km':
        dirGrid="/net/fs09/d1/hajsong/SOCH/Exp1km/run/"
    elif gridname=='goarabia':
        dirGrid="/net/fs09/d1/jm_c/exp/ArG_setup/short_run/res_001/"
    elif gridname=="speedy":
        dirGrid="/home/hajsong/Research/CPL_MITgcm_Speedy/speedy_box/run/"
    elif gridname=="speedyLatLon":
        dirGrid="/home/silvan/research/speedyBox_LatLon/run/"
    elif gridname=="speedycs":
        dirGrid="/home/silvan/research/speedyBox_cs/run/"
    elif gridname=="ctrans3d":
        dirGrid="/data/SOchannel/Exp3d/run_ctrans/" ### 이부분 추가해봄
    elif gridname=="rw3d":
        dirGrid="/data/SOchannel/Exp3d/run_rw/" ### 이부분 추가해봄
        
    tmp=rdmds(dirGrid+"XC")
    [Ly,Lx]=tmp.shape
    if varname is None:
        varname=['XC','YC','RAC','DXC','DYC','hFacC','hFacW','hFacS','Depth',\
                 'RC','RF','DRC','DRF'];
        if flag==2:
            varname=np.append(varname,['XG','YG','RAZ','DXG','DYG'])

    class grd(object):
        for iv,vname in enumerate(varname):
            if region is None:
                exec('tmpvar=rdmds("'+dirGrid+varname[iv]+'")');
                tmpvar=tmpvar.squeeze();
                exec(varname[iv]+'=tmpvar')
            else:
                if vname is 'RC' or vname is 'RF' or vname is 'DRC' or vname is 'DRF':
                    exec('tmpvar=rdmds("'+dirGrid+varname[iv]+'")');
                else:
                    exec('tmpvar=rdmds("'+dirGrid+varname[iv]+'",region='+str(region)+')');
                tmpvar=tmpvar.squeeze();
                exec(varname[iv]+'=tmpvar')
            if vname=='hFacC':
                mskC=hFacC.copy()
                mskC[mskC==0]=np.nan
                mskC[np.isfinite(mskC)]=1.
            if vname=='hFacW':
                mskW=hFacW.copy()
                mskW[mskW==0]=np.nan
                mskW[np.isfinite(mskW)]=1.
            if vname=='hFacS':
                mskS=hFacS.copy()
                mskS[mskS==0]=np.nan
                mskS[np.isfinite(mskS)]=1.
            del tmpvar
    del grd.iv,grd.vname

    return grd
