import numpy as np

def calmld(grd,T,S,delrho=0.03,LK=50):
    """
    Computing MLD
    NEED TO GENERALIZE THE CODE!!
    """
    import numpy as np
    from MITgcmutils import densjmd95
    [ly,lx]=grd.XC.shape

    # compute density
    rho=0*grd.hFacC
    for k in xrange(LK):
      rho[k,:,:]=densjmd95(S[k,:,:],T[k,:,:],0)

    rho=rho*grd.mskC - 1000.0

    # Compute MLD
    mld=np.zeros((ly, lx))
    mld_dens = rho[0,:,:]+delrho
    for ix in range(0,lx):
        for iy in range(0,ly):
            mdens = mld_dens[iy,ix]
            if np.isnan(mdens):
                mld[iy,ix]=float('nan')
            else:
                rhocol = rho[:,iy,ix]
                id=np.where(rhocol>=mdens)[0]
                if len(id)==0:
                    mld[iy,ix]=grd.Depth[iy,ix]
                elif len(id)==len(rhocol):
                    mld[iy,ix]=-grd.RC[0]
                else:
                    alpha=(mdens-rhocol[id[0]-1])/(rhocol[id[0]]-rhocol[id[0]-1])
                    mld[iy,ix]=-(alpha*grd.RC[id[0]]+(1-alpha)*grd.RC[id[0]-1])

    return mld


def invert_map(fldZ, zC, zF, fldVec):
    """
    [ Function ] 
    [kFld,zFld,k0Fld]=invert_map(fldZ,zC,zF,fldVec)

    [Description]
    Find vertical index (kFld) and depth (zFld)
    for each element of the vector fldVec using fldZ, zC and zF

    [Input]
    fldZ   : Potential density on the model's vertical grid, 2D array (nr x nx) 
    zC     : Depth at the grid center, vector (nr)
    zF     : Depth at the interface, vector (nr+1)
    fldVec : potential density initially defined in Layers package, vector(ny)

    [Output]
    kFld  : vertical index at the layers interfaces, 2D array (ny, nx)
    zFld  : depth at the layers interfaces 2D array (ny, nx)
    k0Fld : vertical index before accounting fraction, 2D array (nx,ny)
    """
    ny=len(fldVec)
    [nr,nx] = fldZ.shape
    kFld=np.zeros([ny,nx]) 
    zFld=np.zeros([ny,nx]) 
    k0Fld=np.zeros([ny,nx])
    #- extended vert. res:
    zExt = np.append(zC,zF[-1])

    #%- assume field=0 @ land-pts and field > 0 elsewhere:
    MxV=fldZ.max() 
    botV=fldVec.max()
    botV=np.maximum(botV,MxV)+1
    var=fldZ.copy()
    var[np.nonzero(fldZ==0)] = botV

    fldV=np.ones([nr+1,nx])*botV 
    fldV[:-1,:] = var
    mnV=np.min(fldV)
    var1d = fldV.flatten()

    for j,ff in enumerate(fldVec):
        kk=np.zeros(nx, dtype='int32') 
        zz=np.zeros(nx)
        var=(fldV[1:nr+1,:]-ff)*(fldV[:-1,:]-ff)
        [Ks, Is] = np.nonzero(var < 0)
        kk[Is] = Ks

        [Ke, Ie] = np.nonzero(fldV[:-1,:] == ff)
        nEx=len(Ie) 
        nu=0
        if nEx > 0:
            Iu = np.unique(Ie) 
            nu = len(Iu)
            if nu == nEx:
                kk[Ie]=Ke
            else:
              print("warning 1")
              Ku = Iu.copy()
              for l in xrange(nu):
                  L = np.nonzero(Ie == Iu[l])
                  Ku[l] = Ke[L[0]]
              kk[Iu] = Ku

        k1 = np.maximum(kk,0)
        k2 = np.minimum(kk+1,nr)
        ik1 = k1*nx 
        ik1 = ik1 + np.arange(nx)
        ik2 = k2*nx 
        ik2 = ik2 + np.arange(nx)
        k1 = k1.astype('int32') 
        k2 = k2.astype('int32')
        ik1 = ik1.astype('int32') 
        ik2 = ik2.astype('int32')

        dfld = var1d[ik2] - var1d[ik1] 
        dfld[np.nonzero(kk==0)] = 1.
        [J]=np.nonzero(dfld <= 0 )

        dfld = 1./dfld
        frac = ff - var1d[ik1] 
        frac = frac * dfld
        zz = zExt[k1] + frac*(zExt[k2]-zExt[k1])
        dk = frac * (k2 - k1)
        zz[np.nonzero(kk==0)] = 0.
        dk[np.nonzero(kk==0)] = 0.

        [I1]=np.nonzero(fldV[0,:] > ff)

        if len(I1) > 0:
            zz[I1]=zF[0] 
            kk[I1]=0

        if (len(I1)+len(Is)+nu != nx):
            #fprintf(' error for j= %i, ff= %f : nSt= %i , nEx= %i , nL1= %i\n', ...
            #j,ff,length(Is),nu,length(I1))
            print("warning 2")
        #- save in output array:
        kFld[j,:]=kk+dk
        zFld[j,:]=zz
        k0Fld[j,:]=kk

    return kFld, zFld, k0Fld

def map2z(psiRHO,kzLay,zkLay,hFacS):
    """
    Map Psi from (lat-rho) to (lat-dep)
    """
    [nr,ny,nx] = hFacS.shape
    msk = np.ceil(hFacS)
    mskV = msk[:,:,0]
    
    [nlp,ny] = psiRHO.shape
    nLay = nlp - 1
    
    kBot = np.round(np.sum(mskV,axis=0))
    
    psiRes = np.zeros([nr+1,ny])
    for i in range(1,ny):
        [J] = np.nonzero(kzLay[:,i]>0)
        if len(J)*kBot[i] > 0:
            j1 = J.min()
            j1 = np.maximum(1,j1-1)
            z1 = kzLay[j1:nlp,i]
            v1 = psiRHO[j1:nlp,i]
            zg = np.arange(kBot[i])
            vg = np.interp(zg,z1,v1)
            psiRes[:kBot[i],i] = vg
            
    return psiRes
