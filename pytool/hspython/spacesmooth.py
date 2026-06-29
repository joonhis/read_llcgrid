from numba import double,jit
import math

@jit("double[:,:](double[:,:],double[:,:],i4,i4,f4)",nopython=True)
def spsmooth(fld,sfld,xdist,ydist,thres):
    """
    object:    spatial smoothing with a given scale in both x and y direction
    inputs:    'fld' is the field to be smoothed.
               'sfld' is the array that will carry the smoothed field. 
                      As an input, it can be any array with the same size as the 'fld'
               'xdist' and 'ydist' is the distance (km) from the center 
                      that will be averaged.
    output:    'sfld' is the smoothed field.
    """
    [ly,lx]=fld.shape;
    for iy in xrange(0,ly):
        dely = int(ydist/4);    # assume that the grid spacing is roughly 4 km
        delx = int(xdist/4);
        dy=min(iy,dely,ly-iy-1)
        ys=iy-dy
        ye=iy+dy+1
        for ix in xrange(0,lx):
            #if math.isnan(fld[iy,ix])!=True:
                dx=min(ix,delx,lx-ix-1)
                xs=ix-dx;
                xe=ix+dx+1;
                #mdata=fld[ymin:ymax,xmin:xmax]
                z=0
                M=ye-ys
                N=xe-xs
                for j in range(M):
                    for i in range(N):
                        val=fld[ys+j,xs+i]
                        if thres>0:
                            if abs(val-fld[iy,ix])>thres:
                                val=fld[iy,ix]
                        if math.isnan(val):val=fld[iy,ix]
                        z += val
                sfld[iy,ix]=z/(M*N)
            #else:
            #    sfld[iy,ix]=0

    return sfld


@jit("double[:,:](double[:,:],double[:,:],i4,i4)",nopython=True)
def usmooth(fldU,umerge,delx,dely):
    """
    Compute the transport fields on a coarse grid
    after merging over grids defined by delx and dely

    input:
      fldU   - u transport
      umerge - empty array that will be filled in this function
      delx   - number of grid points from the center in x direction
      dely   - number of grid points from the center in y direction

    output:
      umerge - smoothed u transport
    """

    y0=int(math.ceil(float(dely)/2));
    [ly,lx]=umerge.shape;

    for ix in xrange(0,lx,delx):
        #
        # 1. smoothing u along the patch boundary in y direction
        #  (a) compute umerge by averaging the transport for the patch
        #
        ys=0;ye=y0;
        usum=0;
        for i in xrange(ys,ye):
            usum += fldU[i,ix]
        if math.isnan(usum) == False:
            for i in xrange(ys,ye):
                umerge[i,ix]=usum/(ye-ys);
        else:
            for i in xrange(ys,ye):
                umerge[i,ix]=fldU[i,ix]

        for iy in xrange(y0,ly,dely):
            ys=iy;
            ye=min(ys+dely,ly);
            usum=0.
            c=0
            for i in xrange(ys,ye):
                if math.isnan(fldU[i,ix])==False:
                    usum += fldU[i,ix]
                    c+=1
            for i in xrange(ys,ye):
                if math.isnan(fldU[i,ix])==False:
                    umerge[i,ix]=usum/c;
                else:
                    umerge[i,ix]=fldU[i,ix];
      
     
        #
        #  (b) umerge along iy is done. Smoothing vmerge
        #
#        ys=0;
#        ye=ys+dely;
#        for iy in xrange(0,ly):
#            if iy>=ye:
#                ys=iy;
#                ye=ys+dely;
#            usum=0.
#            for i in xrange(ys,ye):
#                usum += umerge[i,ix]
#            if math.isnan(usum) == False:
#                y3=ye-1;
#                a=umerge[ys,ix];
#                b=umerge[y3,ix];
#                n=y3-ys
#                umerge[iy,ix]=((y3-iy)*a+(iy-ys)*b)/n
#
#        #
#        # 2. smoothing u inside of the patch boundary in x direction
#        #
#        if ix==0:
#            xe=0;
#        elif ix>0:
#            xs=xe;xe=ix;
#            for iy in xrange(0,ly):
#                usum=0.;
#                for i in xrange(xs,ix):
#                    usum+=umerge[iy,i]
#                if math.isnan(usum) == False:
#                    n=xe-xs
#                    for i in xrange(xs+1,xe):
#                        umerge[iy,i]=\
#                            ((i-xs)*umerge[iy,xe]+(xe-i)*umerge[iy,xs])/n
#                else:
#                    for i in xrange(xs+1,xe):
#                        umerge[iy,i]=fldU[iy,i]
#    #
#    #  If the ix is not the last grid cell, then repeat the process for 
#    #  the grid points that were not covered before
#    # 
#    if ix<(lx-1):
#        xs=ix
#        ix=lx-1
#        ys=0;ye=y0;
#        usum=0;
#        for i in xrange(ys,ye):
#            usum += fldU[i,ix]
#        if math.isnan(usum) == False:
#            for i in xrange(ys,ye):
#                umerge[i,ix]=usum/(ye-ys);
#        else:
#            for i in xrange(ys,ye):
#                umerge[i,ix]=fldU[i,ix]
#
#        for iy in xrange(y0,ly,dely):
#            ys=iy;
#            ye=min(ys+dely,ly);
#            usum=0.
#            c=0
#            for i in xrange(ys,ye):
#                if math.isnan(fldU[i,ix])==False:
#                    usum += fldU[i,ix]
#                    c+=1
#            for i in xrange(ys,ye):
#                if math.isnan(fldU[i,ix])==False:
#                    umerge[i,ix]=usum/c;
#                else:
#                    umerge[i,ix]=fldU[i,ix];
#
#        ys=0;
#        ye=ys+dely;
#        for iy in xrange(0,ly):
#            if iy>=ye:
#                ys=iy;
#                ye=ys+dely;
#            usum=0.
#            for i in xrange(ys,ye):
#                usum += umerge[i,ix]
#            if math.isnan(usum) == False:
#                y3=ye-1;
#                a=umerge[ys,ix];
#                b=umerge[y3,ix];
#                n=y3-ys
#                umerge[iy,ix]=((y3-iy)*a+(iy-ys)*b)/n
#
#        for iy in xrange(0,ly):
#            usum=0.;
#            for i in xrange(xs,ix):
#                usum+=umerge[iy,i]
#            if math.isnan(usum) == False:
#                n=ix-xs
#                for i in xrange(xs+1,ix):
#                    umerge[iy,i]=\
#                        ((i-xs)*umerge[iy,ix]+(ix-i)*umerge[iy,xs])/n
#            else:
#                for i in xrange(xs+1,ix):
#                    umerge[iy,i]=fldU[iy,i]
    return umerge
    
    
@jit("double[:,:](double[:,:],double[:,:],i4,i4)",nopython=True)
def vsmooth(fldV,vmerge,delx,dely):
    """
    Compute the transport fields on a coarse grid
    after merging over grids defined by delx and dely

    input:
      fldV   - v transport
      vmerge - empty array that will be filled in this function
      delx   - number of grid points from the center in x direction
      dely   - number of grid points from the center in y direction

    output:
      vmerge - smoothed u transport
    """

    x0=int(math.ceil(float(delx)/2));
    [ly,lx]=vmerge.shape;

    for iy in xrange(0,ly,dely):
        #
        # 1. smoothing v along the patch boundary in x direction
        #  (a) compute vmerge by averaging the transport for the patch
        #
        xs=0;xe=x0;
        vsum=0;
        for i in xrange(xs,xe):
            vsum += fldV[iy,i]
        if math.isnan(vsum) == False:
            for i in xrange(xs,xe):
                vmerge[iy,i]=vsum/(xe-xs);

        for ix in xrange(x0,lx,delx):
            xs=ix;
            xe=min(xs+delx,lx);
            vsum=0.
            c=0
            for i in xrange(xs,xe):
                if math.isnan(fldV[iy,i])==False:
                    vsum += fldV[iy,i]
                    c+=1
            for i in xrange(xs,xe):
                if math.isnan(fldV[iy,i]) == False:
                    vmerge[iy,i]=vsum/c
                else:
                    vmerge[iy,i]=fldV[iy,i];
        #
        #  (b) vmerge along ix is done. Smoothing vmerge
        #
#        xs=0;
#        xe=xs+delx;
#        for ix in xrange(0,lx):
#            if ix>=xe:
#                xs=ix;
#                xe=xs+delx;
#            vsum=0.
#            for i in xrange(xs,xe):
#                vsum += vmerge[iy,i]
#            if math.isnan(vsum) == False:
#                x3=xe-1;
#                a=vmerge[iy,xs];
#                b=vmerge[iy,x3];
#                n=x3-xs
#                vmerge[iy,ix]=((x3-ix)*a+(ix-xs)*b)/n
#
#        #
#        # 2. smoothing u inside of the patch boundary in x direction
#        #
#        if iy==0:
#            ye=0;
#        elif iy>0:
#            ys=ye;ye=iy;
#            for ix in xrange(0,lx):
#                vsum=0.;
#                for i in xrange(ys,iy+1):
#                    vsum+=vmerge[i,ix]
#                if math.isnan(vsum) == False:
#                    n=ye-ys
#                    for i in xrange(ys+1,ye):
#                        vmerge[i,ix]=\
#                            ((i-ys)*vmerge[ye,ix]+(ye-i)*vmerge[ys,ix])/n
#               
#                else:
#                    for i in xrange(ys+1,ye):
#                        vmerge[i,ix]=fldV[i,ix]
#    #
#    #  If the iy is not the last grid cell, then repeat the process for
#    #  the grid points that were not covered before
#    #
#    if iy<(ly-1):
#        ys=iy
#        iy=ly-1
#        xs=0;xe=x0;
#        vsum=0;
#        for i in xrange(xs,xe):
#            vsum += fldV[iy,i]
#        if math.isnan(vsum) == False:
#            for i in xrange(xs,xe):
#                vmerge[iy,i]=vsum/(xe-xs);
#
#        for ix in xrange(x0,lx,delx):
#            xs=ix;
#            xe=min(xs+delx,lx);
#            vsum=0.
#            c=0
#            for i in xrange(xs,xe):
#                if math.isnan(fldV[iy,i])==False:
#                    vsum += fldV[iy,i]
#                    c+=1
#            for i in xrange(xs,xe):
#                if math.isnan(fldV[iy,i]) == False:
#                    vmerge[iy,i]=vsum/c
#                else:
#                    vmerge[iy,i]=fldV[iy,i];
#         
#        xs=0;
#        xe=xs+delx;
#        for ix in xrange(0,lx):
#            if ix>=xe:
#                xs=ix;
#                xe=xs+delx;
#            vsum=0.
#            for i in xrange(xs,xe):
#                vsum += vmerge[iy,i]
#            if math.isnan(vsum) == False:
#                x3=xe-1;
#                a=vmerge[iy,xs];
#                b=vmerge[iy,x3];
#                n=x3-xs
#                vmerge[iy,ix]=((x3-ix)*a+(ix-xs)*b)/n
#
#        for ix in xrange(0,lx):
#            vsum=0.;
#            for i in xrange(ys,iy+1):
#                vsum+=vmerge[i,ix]
#            if math.isnan(vsum) == False:
#                n=iy-ys
#                for i in xrange(ys+1,iy):
#                    vmerge[i,ix]=\
#                        ((i-ys)*vmerge[iy,ix]+(iy-i)*vmerge[ys,ix])/n
#
#            else:
#                for i in xrange(ys+1,iy):
#                    vmerge[i,ix]=fldV[i,ix]
    return vmerge
    
    
