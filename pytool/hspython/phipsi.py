import numpy as np
def calc_conv2D(grd,fldU,fldV,xrng=None,yrng=None,ik=0,list_factors=None):
    """
    object:    compute flow field convergent part (i.e. minus the divergence)
    inputs:    fldU and fldV are transport or velocity fields
    list_factors the list of factors that need to
                 be applied to fldU,fldV. By default it is empty (i.e. {}).
                 The most complete list wouldrbe ['dh','dz','hfac'].
    output:    fldDIV is the convergence
               (integrated, not averaged, over grid cell area)

    notes:     fldU,fldV  that may be either
               [A] a 3D vector field or
               [B] a 2D vector field

    in case [A], layer thicknesses = mygrid.DRF; in case [B] layer thickness = 1
    in any case, the global variable mygrid is supposed to be available
    """

    # Check the x/y range
    [Ly,Lx]=grd.XC.shape;
    if xrng==None:
        lx=Lx;ly=Ly;
        xs=0;xe=lx;
        ys=0;ye=ly;
        xend=lx-1;
        yend=ly-1;
    else:
        ly=len(yrng)
        lx=len(xrng)
        xs=xrng[0];xe=xs+lx;
        ys=yrng[0];ye=ys+ly;
        xend=xrng[-1]+1;
        yend=yrng[-1]+1;
        if xend==Lx: xend=Lx-1;
        if yend==Ly: yend=Ly-1;

    # prepare fldU/fldV:
    fldU[np.isnan(fldU)]=0.; fldV[np.isnan(fldV)]=0.;

    if list_factors!=None:
        facW=np.ones([Ly,Lx])
        facS=np.ones([Ly,Lx])
        for tmpstr in list_factors:
            if tmpstr=='dh': facW=facW*grd.DYG; facS=facS*grd.DXG;
            elif tmpstr=='dz': facW=facW*grd.DRF[ik]; facS=facS*grd.DRF[ik];
            elif tmpstr=='hfac':
                facW=facW*grd.hFacW[ik,...];
                facS=facS*grd.hFacS[ik,...];
        fldU=fldU*facW
        fldV=fldV*facS

    FLDU=np.empty([ly,lx+1]);
    FLDV=np.empty([ly+1,lx]);
    fldDIV=np.empty([ly,lx]);
    #FLDU=fldU[yrng][:,np.append(xrng,xend)]
    #FLDV=fldV[np.append(yrng,yend)][:,xrng]
    FLDU=np.hstack((fldU[ys:ye,xs:xe],fldU[ys:ye,xend:xend+1]))
    FLDV=np.vstack((fldV[ys:ye,xs:xe],fldV[yend:yend+1,xs:xe]))
    fldDIV=FLDU[:,0:-1]-FLDU[:,1:]+FLDV[0:-1,:]-FLDV[1:,:];
    if xend==Lx-1:
        fldDIV[:,-1]=fldDIV[:,-2];
    if yend==Ly-1:
        fldDIV[-1,:]=fldDIV[-2,:]
    return fldDIV

def calc_UV_grad(grd,fld,xrng=None,yrng=None):
    """
    object:    compute horizontal gradients
    inputs:    fld is the 'tracer' field of interest
    output:    [dFLDdx,dFLDdy] are the gradient fields
    """

    # Check the x/y range
    [Ly,Lx]=grd.XC.shape;
    if xrng==None:
        xs=0;xsn=xs;
        ys=0;ysn=ys;
        xe=Lx;xep=Lx-1;
        ye=Ly;yep=Ly-1;
        lx=Lx;ly=Ly;
    else:
        ly=len(yrng)
        lx=len(xrng)
        xs=xrng[0];
        ys=yrng[0];
        xe=xs+lx;
        ye=ys+ly;
        xsn=max(0,xs-1);
        ysn=max(0,ys-1);
        xep=min(Lx-1,xe);
        yep=min(Ly-1,ye);

    FLD=np.empty([ly+2,lx+2]);
    FLD[1:-1,1:-1]=fld[ys:ye,xs:xe];
    FLD[1:-1,0]=fld[ys:ye,xsn];
    FLD[1:-1,-1]=fld[ys:ye,xep]
    FLD[0,1:-1]=fld[ysn,xs:xe];
    FLD[-1,1:-1]=fld[yep,xs:xe]
    FLD[0,0]=fld[ysn,xsn];FLD[-1,-1]=fld[yep,xep]
    FLD[0,-1]=fld[ysn,xep];FLD[-1,0]=fld[yep,xsn]

    dFLDdx=np.empty((ly,lx))
    dFLDdy=np.empty((ly,lx))
    tmpA=FLD[1:-1,1:-1];
    tmpB=FLD[0:-2,1:-1];
    dFLDdy=(tmpA-tmpB)/grd.DYC[ys:ye,xs:xe];
    tmpA=FLD[1:-1,1:-1];
    tmpB=FLD[1:-1,0:-2];
    dFLDdx=(tmpA-tmpB)/grd.DXC[ys:ye,xs:xe];

    return dFLDdx, dFLDdy

def calc_dFLDdt(grd,fld,msk,ik):
    """ compute dFLDdt """
    [ly,lx]=fld.shape
    #facx=np.ones((ly,lx))
    #facy=np.ones((ly,lx))
    facx=grd.DRF[ik]*grd.DYG*grd.hFacW[ik,...]/grd.DXC
    facy=grd.DRF[ik]*grd.DXG*grd.hFacS[ik,...]/grd.DYC
    #facx=grd.DRF[ik]*grd.DYG/grd.DXC
    #facy=grd.DRF[ik]*grd.DXG/grd.DYC

    fld=fld*msk;
    FLD=np.empty([ly+2,lx+2])
    FLD[1:-1,1:-1]=fld;
    FLD[:,0]=FLD[:,1];FLD[:,-1]=FLD[:,-2]
    FLD[0,:]=FLD[1,:];FLD[-1,:]=FLD[-2,:]
    FacX=np.empty([ly+2,lx+2])
    FacX[1:-1,1:-1]=facx;
    FacX[:,0]=FacX[:,1];FacX[:,-1]=FacX[:,-2]
    FacX[0,:]=FacX[1,:];FacX[-1,:]=FacX[-2,:]
    FacY=np.empty([ly+2,lx+2])
    FacY[1:-1,1:-1]=facy;
    FacY[:,0]=FacY[:,1];FacY[:,-1]=FacY[:,-2]
    FacY[0,:]=FacY[1,:];FacY[-1,:]=FacY[-2,:]

    dFLDdt=(fld-FLD[1:-1,0:-2])*FacX[1:-1,1:-1]
    dFLDdt[np.isnan(dFLDdt)]=0.
    atmp=(fld-FLD[1:-1,2:])*FacX[1:-1,2:];
    atmp[np.isnan(atmp)]=0.
    dFLDdt=dFLDdt+atmp
    atmp=(fld-FLD[2:,1:-1])*FacY[2:,1:-1];
    atmp[np.isnan(atmp)]=0.
    dFLDdt=dFLDdt+atmp
    atmp=(fld-FLD[0:-2,1:-1])*FacY[1:-1,1:-1];
    atmp[np.isnan(atmp)]=0.
    dFLDdt=dFLDdt+atmp

    return dFLDdt
    
def calc_velpot(grd,fldU,fldV,xrng=None,yrng=None,ik=0):
    """
    object: compute the divergent part of a tranport vector field
            by solving a Laplace equation with Neumann B.C.
 
    Some math: u=dPhi/dx and v=dPhi/dy for the nondivergent part
               dyg(i,j)*drf*hFacW(i,j)*u(i,j) 
                = dyg(i,j)*drf*hFacW(i,j)*(Phi(i,j)-Phi(i-1,j))/dxc(i,j)
               dxg(i,j)*drf*hFacS(i,j)*v(i,j) 
                = dxg(i,j)*drf*hFacS(i,j)*(Phi(i,j)-Phi(i,j-1))/dyc(i,j)

    input:  fldU,fldV       2D transport vectors (masked with NaNs)
            xrng,yrng       arrays for the subset of the data (optional)
            ik              vertical level
    output: fldUdiv,fldVdiv divergent part (m/s)
                            To convert to transport, fldUdiv*dyg*drf*hfacW
            fldDivPot       Potential for the divergent flow (optional)
                            (m^2/s)
    note 1: here we do everything in terms of transports, so we omit all
            grid factors (i.e. with dx=dy=1), and the resulting potential
            has usits of transports.
    note 2: this routine was derived from diffsmooth2D_extrap_inv.m but
        it differs in several aspects.
        It : 1) omits grid factors;
             2) uses of a neumann, rather than dirichlet, boundary condition;
             3) gets the mask directly from fld
    """
    from scipy.sparse import coo_matrix, spdiags
    from scipy.sparse.linalg import spsolve

    # Check the x/y range
    [Ly,Lx]=grd.XC.shape;
    if xrng==None:
        lx=Lx;ly=Ly;
        xs=0;xe=Lx;
        ys=0;ye=Ly;
        # no grid scaling factor in div. or transports
        fld=calc_conv2D(grd,fldU,fldV)
    else:
        ly=len(yrng)
        lx=len(xrng)
        xs=xrng[0];xe=xs+lx;
        ys=yrng[0];ye=ys+ly;
        # no grid scaling factor in div. or transports
        fld=calc_conv2D(grd,fldU,fldV,xrng=xrng,yrng=yrng)

    fld=fld*grd.mskC[ik,ys:ye,xs:xe];
    mskWet=np.empty((ly,lx));
    np.copyto(mskWet,grd.mskC[ik,ys:ye,xs:xe]);
    mskDry=1.*np.isnan(mskWet); mskDry[mskDry==0]=np.nan;

    # put 0 first guess if needed and switch land mask:
    fld[np.isnan(fld)]=0; fld=fld*mskWet;

    # define mapping from global array to (no nan points) global vector
    tmp1=mskWet.flatten()
    kk=np.where(np.isfinite(tmp1))[0]; nn=len(kk);
    KK=np.where(np.isfinite(tmp1),np.arange(0,ly*lx),np.nan);
    KK=KK.reshape(ly,lx);
    LL=np.zeros([ly*lx]); LL[kk]=np.arange(nn); LL=LL.reshape(ly,lx);

    A=coo_matrix(([],([],[])), shape=(nn,nn))
    for ii in xrange(0,3):
        for jj in xrange(0,3):
            #
            # 1) seed points (FLDones)
            #    and neighborhood of influence (FLDkkFROM)
            #
            FLDones=np.zeros([ly,lx])
            FLDones[jj::3,ii::3]=1;
            FLDones[np.isnan(KK)]=0;
            #
            FLDkkFROMtmp=np.zeros([ly,lx])
            FLDkkFROMtmp[jj::3,ii::3]=KK[jj::3,ii::3]*1
            #FLDkkFROMtmp[np.isnan(fld)]=0;
            #
            tmp1=np.zeros([ly+2,lx+2]);
            tmp1[1:-1,1:-1]=FLDkkFROMtmp*1
            FLDkkFROM=np.zeros([ly,lx]);
            for ii2 in xrange(0,3):
                for jj2 in xrange(0,3):
                    tmp2=tmp1[np.arange(0,ly)+jj2][:,np.arange(0,lx)+ii2]
                    FLDkkFROM=np.nansum(np.asarray([FLDkkFROM,tmp2]),axis=0)
            #
            # 2) compute effect of each point on neighboring target point:
            #
            dFLDdt=calc_dFLDdt(grd,FLDones,mskWet,ik)
            #
            # 3) include seed contributions in matrix:
            #   3.1) for wet points
            dFLDdtWet=dFLDdt*mskWet;
            #dFLDdtWet[np.isnan(dFLDdtWet)]=0.
            tmp1=np.where(np.isfinite(dFLDdtWet.flatten()))[0]
            dFLDdtWet=dFLDdtWet.flatten()[tmp1];
            FLDkkFROMtmp=FLDkkFROM.flatten()[tmp1];
            FLDkkTOtmp=KK.flatten()[tmp1];
            row=LL.flatten()[FLDkkTOtmp.astype(int)]
            col=LL.flatten()[FLDkkFROMtmp.astype(int)]
            A=A+coo_matrix((dFLDdtWet,(row,col)),shape=(nn,nn));
            #   3.2) for dry points
            #       (this part reflects the neumann boundary condition)
            dFLDdtDry=dFLDdt*mskDry;
            #dFLDdtDry[np.isnan(dFLDdtDry)]=0.
            tmp1=np.where(np.isfinite(dFLDdtDry.flatten()))[0]
            dFLDdtDry=dFLDdtDry.flatten()[tmp1];
            FLDkkFROMtmp=FLDkkFROM.flatten()[tmp1];
            col=LL.flatten()[FLDkkFROMtmp.astype(int)]
            A=A+coo_matrix((dFLDdtDry,(col,col)),shape=(nn,nn));

    # 4) solve for potential:
    # add noise to compute the inverse
    if nn>1e4:
        A=A+spdiags((np.random.rand(nn)-1)*1e-12,0,nn,nn)

    yy=fld.flatten()
    yy=yy[np.where(np.isfinite(KK.flatten()))]
    xx=spsolve(A,yy);

    # 5) Prepare output:
    fldDivPot=np.empty([ly,lx])
    fldDivPot[np.isfinite(KK)]=xx
    fldDivPot=fldDivPot.reshape(ly,lx)*grd.mskC[ik,ys:ye,xs:xe]
    FDP=np.empty([Ly,Lx])
    FDP[ys:ye,xs:xe]=fldDivPot;
    [fldUdiv,fldVdiv]=calc_UV_grad(grd,FDP,xrng=range(xs,xe),yrng=range(ys,ye));

    return fldUdiv,fldVdiv,fldDivPot

def calc_barostream(grd,fldU,fldV,ik=0,noDiv=1,list_factors=['dh','dz']):
    """
    <<< NOT TESTED THOROUGHLY, DO NOT USE YET! >>>
    object:    compute barotropic streamfunction
    inputs:    fldU and fldV are the fields of grid point transport
    optional:  noDiv (default is 1). If 1 then remove the divergent 
               part of the flow field first. If 0 then dont.
               list_factors (default is {'dh','dz'})
    output:    FLD is the streamfunction
    notes:     the result is converted to Sv
    """

    # 0) prepare fldU/fldV (transport fields):
    if len(fldU.shape)==3:
        [n3,ly,lx]=fldU.shape
        dxg=np.tile(grd.DXG,[n3,1,1]);
        dyg=np.tile(grd.DYG,[n3,1,1]);
    elif len(fldU.shape)==2:
        [ly,lx]=fldU.shape
        dxg=grd.DXG
        dyg=grd.DYG
        n3=1;

    fldU[np.isnan(fldU)]=0.
    fldV[np.isnan(fldV)]=0.

    if n3==len(grd.DRF):
        drf=np.tile(grd.DRF,[1,ly,lx])
        facW=np.ones([n3,ly,lx])
        facS=np.ones([n3,ly,lx])
    elif n3==1:
        drf=np.ones([ly,lx])
        facW=np.ones([ly,lx])
        facS=np.ones([ly,lx])

    for tmpstr in list_factors:
        if tmpstr=='dh': facW=facW*dyg; facS=facS*dxg;
        elif tmpstr=='dz': facW=facW*drf; facS=facS*drf;
        elif tmpstr=='hfac': facW=facW*mygrid.hFacW; facS=facS*mygrid.hFacS;

    # apply mask:
    if n3>1:
        fldU=fldU.sum(axis=0)*grd.mskW[0,...];
        fldV=fldV.sum(axis=0)*grd.mskS[0,...];

    # take out thedivergent part of the flow:
    if noDiv==1:
        [fldUdiv,fldVdiv,fldDivPot]=diffsmooth2D_div_inv(grd,fldU,fldV,0);
        fldU=fldU-fldUdiv; fldV=fldV-fldVdiv;

    # 1) compute streamfunction
    FLDV=np.zeros([ly,lx+1]); FLDU=np.zeros([ly+1,lx])
    FLDV[:,1:]=fldV; #FLDV[:,-1]=fldV[:,-1]
    FLDU[1:,:]=fldU; #FLDU[-1,:]=fldU[-1,:]
    FLDU[np.isnan(FLDU)]=0.; FLDV[np.isnan(FLDV)]=0.;
    tmp1=np.cumsum(FLDU,axis=0);
    tmp1=np.append(tmp1,tmp1[:,-1].reshape(len(tmp1),1),1)

    tmp2=np.diff(tmp1,axis=0)+FLDV;
    tmp3=np.cumsum(np.mean(tmp2,axis=1));

    # to check divergence implied errors:  
    # figure; for iF=1:fldU.nFaces; subplot(3,2,iF); plot(std(tmp2{iF},0,1));
    tmp1=tmp1-np.tile(np.hstack((tmp3,0)),[tmp1.shape[1],1]).transpose()
    bf_step1=tmp1;

    # 3) put streamfunction at cell center
    #tmp1=(bf_step1[:-1,:]+bf_step1[1:,:])/2.;
    #bf_step2=(tmp1[:,:-1]+tmp1[:,1:])/2.

    #return fldU,fldV,bf_step1,bf_step2

    # 4) set 0 on land on average:
    iland=np.where(np.isnan(grd.mskC[ik,...])&np.isfinite(bf_step2[:-1,:-1]))
    tmp2=np.nanmedian(bf_step2[iland])
    bf_step3=(bf_step2-tmp2)

    # 5) return the result:
    fldBAR=bf_step3;
    return fldU,fldV,fldDivPot,fldBAR

def calc_localpsi(grd,fldU,fldV,xrng=None,yrng=None):
    """
    Compute local streamfunction
    Input  
      + grd : grid structure variable
      + fldU, fldV : 2D nondivergent U and V
      + xrng, yrng : local x and y indices
    """
    #
    #  Select local U and V
    #
    [Ly,Lx]=grd.XC.shape
    dim=fldU.shape
    ndim=len(dim)
    #
    if xrng is None: xrng=np.array(Lx)
    if yrng is None: yrng=np.array(Ly)
    lx=len(xrng)
    ly=len(yrng)
    #
    #  Initialize Psi
    #
    if ndim==3:
        Psi=np.zeros((dim[0],ly+1,lx+1))
    elif ndim==2:
        Psi=np.zeros((ly+1,lx+1))
    if xrng[-1]>=Lx-1:
        if ndim==3:
            FLDU=fldU[:,yrng][:,:,xrng]
            FLDU=np.append(FLDU,np.zeros((dim[0],ly,1))\
                           ,axis=2)
        elif ndim==2:
            FLDU=fldU[yrng][:,xrng]
            FLDU=np.append(FLDU,np.zeros((ly,1)),axis=1)
    else:
        if ndim==3:
            FLDU=fldU[:,yrng][:,:,np.append(xrng,xrng[-1]+1)]
        elif ndim==2:
            FLDU=fldU[yrng][:,np.append(xrng,xrng[-1]+1)]

    if yrng[-1]>=Ly-1:
        if ndim==3:
            FLDV=fldV[:,yrng][:,:,xrng]
            FLDV=np.append(FLDV,np.zeros((dim[0],1,lx))\
                           ,axis=1)
        elif ndim==2:
            FLDV=fldV[yrng][:,xrng]
            FLDV=np.append(FLDV,np.zeros((1,FLDV.shape[1])),axis=0)
    else:
        if ndim==3:
            FLDV=fldV[:,np.append(yrng,yrng[-1]+1)][:,:,xrng]
        elif ndim==2:
            FLDV=fldV[np.append(yrng,yrng[-1]+1)][:,xrng]

    FLDU[np.isnan(FLDU)]=0.;
    FLDV[np.isnan(FLDV)]=0.;
    #
    #  start from x0 and y0
    #
    ix=lx/2
    iy=ly/2
    if ndim==3:
        psix=np.zeros((dim[0],1,lx+1))
        psix[:,0,1:]=np.cumsum(FLDV[:,iy,:],axis=1)
        psix[:,0,:]=psix[:,0,:]-psix[:,0,ix:ix+1]
        #Psi[:,iy,:]=Psi[:,iy,:]-np.tile(Psi[:,iy,ix:ix+1],[1,lx+1])
        psiy=np.zeros((dim[0],ly+1,1))
        psiy[:,1:,0]=np.cumsum(-FLDU[:,:,ix],axis=1)
        psiy[:,:,0]=psiy[:,:,0]-psiy[:,iy:iy+1,0]
    elif ndim==2:
        Psi[iy,1:]=np.cumsum(FLDV[iy,:],axis=0)
        Psi[iy,:]=Psi[iy,:]-Psi[iy,ix]
        Psi[1:,ix]=np.cumsum(-FLDU[:,ix],axis=0)
        Psi[:,ix]=Psi[:,ix]-Psi[iy,ix]
    #
    #  Fill Psi with the mean of two integrals (x and y direction)
    #
    if ndim==3:
        inty=np.zeros((dim[0],ly+1,lx+1))
        inty[:,1:,:]=np.cumsum(-FLDU,axis=1)
        inty=inty-inty[:,iy:iy+1,:]+psix
        intx=np.zeros((dim[0],ly+1,lx+1))
        intx[:,:,1:]=np.cumsum(FLDV,axis=2)
        intx=intx-intx[:,:,ix:ix+1]+psiy
    elif ndim==2:
        inty=np.cumsum(-FLDU[iy:,(ix+1):],axis=0)\
             +np.tile(Psi[iy,(ix+1):],[ly-iy,1])
        intx=np.cumsum(FLDV[(iy+1):,ix:],axis=1)\
             +np.tile(Psi[(iy+1):,ix].T,[lx-ix,1]).T
        Psi[(iy+1):,(ix+1):]=(intx+inty)/2

    Psi=(intx+inty)/2
    return Psi

