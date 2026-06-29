import numpy as np

def readbin(grd,filename):
    f=open(filename,'rb')
    data = np.fromfile(f,'>f4')
    f.close()
    lendata=len(data)
    [ly,lx]=grd.XC.shape
    lz=lendata/lx/ly
    if lz==1:
        var = np.reshape(data,[ly,lx])
    else:
        if lz==np.floor(lz):
            var = np.reshape(data,[lz,ly,lx])
        else:
            raise IntersectException("Check the dimension again.")
           
    return var

def read_data(filename,lx,ly):
    f=open(filename,'rb')
    data = np.fromfile(f,'>f4')
    lendata=len(data)
    lz=lendata/lx/ly
    if lz==1:
        var = np.reshape(data,[ly,lx])
    else:
        if lz==np.floor(lz):
            var = np.reshape(data,[lz,ly,lx])
        else:
            raise IntersectException("Check the dimension again.")

    return var

def tic():
    #Homemade version of matlab tic and toc functions
    import time
    global startTime_for_tictoc
    startTime_for_tictoc = time.time()

def toc():
    import time
    if 'startTime_for_tictoc' in globals():
        print("Elapsed time is "+str(time.time() - startTime_for_tictoc) + " seconds.")
    else:
        print("Toc: start time not set")
