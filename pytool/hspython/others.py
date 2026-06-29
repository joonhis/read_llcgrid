import numpy as np
import scipy.signal

def inpaint_nans(im):
    import scipy.signal
    ipn_kernel = np.array([[1,1,1],[1,0,1],[1,1,1]]) # kernel for inpaint_nans
    nans = np.isnan(im)
    while np.sum(nans)>0:
        im[nans] = 0
        vNeighbors = scipy.signal.convolve2d((nans==False),ipn_kernel,mode='same',boundary='symm')
        im2 = scipy.signal.convolve2d(im,ipn_kernel,mode='same',boundary='symm')
        im2[vNeighbors>0] = im2[vNeighbors>0]/vNeighbors[vNeighbors>0]
        im2[vNeighbors==0] = np.nan
        im2[(nans==False)] = im[(nans==False)]
        im = im2
        nans = np.isnan(im)
    return im

def inpaintnans(im, N):
    import scipy.signal
    ipn_kernel = np.array([[1,1,1],[1,0,1],[1,1,1]]) # kernel for inpaint_nans
    nans = np.isnan(im)
    for i in xrange(N):
        im[nans] = 0
        vNeighbors = scipy.signal.convolve2d((nans==False),ipn_kernel,mode='same',boundary='symm')
        im2 = scipy.signal.convolve2d(im,ipn_kernel,mode='same',boundary='symm')
        im2[vNeighbors>0] = im2[vNeighbors>0]/vNeighbors[vNeighbors>0]
        im2[vNeighbors==0] = np.nan
        im2[(nans==False)] = im[(nans==False)]
        im = im2
        nans = np.isnan(im)

    return im

def detrend2d(Z):
    N,M=Z.shape
    [X,Y] = np.meshgrid(range(M),range(N));
    #
    # Make the 2D data as 1D vector
    #
    Xcolv = X.reshape(M*N,1);
    Ycolv = Y.reshape(M*N,1);
    Zcolv = Z.reshape(M*N,1);
    Const = np.ones((len(Xcolv),1))
    #
    # find the coeffcients of the best plane fit
    #
    A=np.concatenate((Xcolv,Ycolv,Const),axis=1)
    Coeff=np.linalg.lstsq(A,Zcolv)

    Z_p = Coeff[0][0] * X + Coeff[0][1] * Y + Coeff[0][2];
    Z_f = Z - Z_p;
    return Z_f,Z_p
