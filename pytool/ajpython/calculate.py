import numpy as np

def geostrophicCurrent(ssh, X, Y, g=9.81, a=6.371e6, obs=False):
    """
    + object : get geostrophic current with sea surface height data
    + input
      - ssh : sea surface height
              ssh.shape = (nt,ny,nx)
      - X, Y : grid data with unit of degree which is divided uniformly
    + output
      - u, v : geostrophic current velocity of zonal and meridional direction 
              u.shape = v.shape = (nt,ny,nx)
    + cf. Central-difference method is used for derivative. 
          So the output data doesn't have value in the boundary.
          Equator has been masked since Coriolis parameter is zero. 
    """
    nt, ny, nx = ssh.shape
    omega = 7.29e-5
    
    f = 2*omega*np.sin(Y*(np.pi/180))
    if obs==False:
        eqIdx = np.where(Y==0.)[0][0]
        f[eqIdx] = None
    elif obs==True:
        eqIdx = np.where(Y==-0.5)[0][0]
        f[eqIdx] = None
        f[eqIdx+1] = None
    
    dx = a*np.cos(Y*(np.pi/180))*(X[2]-X[0])*(np.pi/180) # 1D array
    dy = a*(Y[2]-Y[0])*(np.pi/180) # constant
    
    u = np.ma.empty(ssh.shape)
    v = np.ma.empty(ssh.shape)
    for iy in range(ny-2):
        u[:,iy+1,:] = -(g/f[iy+1])*((ssh[:,iy+2,:]-ssh[:,iy,:])/dy)
        for ix in range(nx-2):
            v[:,iy+1,ix+1] = (g/f[iy+1])*((ssh[:,iy+1,ix+2]-ssh[:,iy+1,ix])/dx[iy+1])

    return u, v
