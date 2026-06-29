# 코드를 복사해서 주피터에서 직접 사용하는 것을 추천: 상황이 다 달라서 그게 편함
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import datetime
from dateutil import relativedelta
import netCDF4 as nc; from netCDF4 import Dataset 
import matplotlib.patches as patches 
import seaborn as sns ; sns.set(style='white')
import glob
from dhkpython import draw

def monthly2onefile(files,variable,nt,ny,nx):
    mvimd = np.empty([nt,ny,nx])
    for i in range(nt):
        f = Dataset(files[i])
        mvimd[i] = f[str(variable)][0]
        if np.mod(i,30)==0:
            print(i)
    return mvimd
#np.save('/home/dhkim/data/original data/calculated npy/ERA5_mvimd.npy',mvimd)


def standardized_anomaly(data,nt,ny,nx): # 너무 CPU 사용량 많이 잡아먹는다. 다른 변수들 저장시 메모리 차지하므로 다 지우고 import이후 바로 돌리기.
    data_mean=np.empty([nt,ny,nx]);data_std=np.copy(data_mean)
    for i in range(12):
        data_mean[i] = np.mean(data[i::12,:,:],axis=0)
        data_std[i] = np.std(data[i::12],axis=0)
    rmdata = np.empty([nt,ny,nx],dtype=np.float32)
    for i in range(nt):
        rmdata[i,:,:] = (data[i,:,:] - data_mean[np.mod(i,12)])/data_std[np.mod(i,12)]
    return rmdata

def standardized_anomaly_1d(data):
    data_mean=np.empty([len(data)]);data_std=np.copy(data_mean)
    for i in range(12):
        data_mean[i] = np.mean(data[i::12])
        data_std[i] = np.std(data[i::12])
    rmdata = np.empty([len(data)])
    for i in range(len(data)):
        rmdata[i] = (data[i] - data_mean[np.mod(i,12)])/data_std[np.mod(i,12)]
    return rmdata

### 할 때마다 변수이름 바꿔주기
def daily_to_monthly(filelist,year,variable, mask=False,region=None): # np.array 로 나타내기!
    filelist.append(filelist[0]) # 마지막 연도 평균하기 위해서! 추가
    f = Dataset(filelist[0],'r')
    shape = f[variable][0].shape
    #sst = np.empty([len(year),region[1]-region[0],region[3]-region[2]]) # year는 그냥 time range
    sst = np.empty([len(year),shape[0],shape[1]])
    daily_sst = []
    count = 0
    
    f = Dataset(filelist[0],'r')
    date_test = nc.num2date(f.variables['time'][0],f.variables['time'].units) 
    
    for i in range(len(filelist)):
        f = Dataset(filelist[i],'r')
        file_date = nc.num2date(f.variables['time'][0],f.variables['time'].units) #읽어줘야 편하다.
        
        if file_date.month != date_test.month:## 한달치를 모아 월평균
            ############ 2월만 하면 month는 다르고 year이 같아야 한다.!!
            if mask == True:
                sst[count,:,:] = np.ma.masked_less(np.mean(daily_sst,axis=0),-10)
            else:
                sst[count,:,:] = np.mean(daily_sst,axis=0)
            daily_sst=[]; daily_sic=[]    
            count = count + 1
        if region != None:
            daily_sst.append(f[variable][0][region[0]:region[1],region[2]:region[3]])
        else:
            daily_sst.append(f[variable][0])
        date_test = file_date
        f.close()
    return sst


def rmSeason(data,nt,ny,nx):
    """
    Remove seasonal cycle from given data.
    data.shape = (nt,ny,nx)
    """
    monMean = np.ma.empty([12,ny,nx])
    anmData = np.ma.empty([nt,ny,nx]) # SST anomaly
    for i in range(12):
        monMean[i,:,:] = np.average(data[i::12,:,:], axis=0)
        anmData[i::12,:,:] = data[i::12,:,:] - monMean[i,:,:]
    return anmData

def seasonalCycle(data,nt,ny,nx): #1월평균,2월평균...
    """
    Get seasonal cycle from given data.
    data.shape = (nt,ny,nx)
    """
    monMean = np.ma.empty([12,ny,nx])
    for i in range(12):
        monMean[i,:,:] = np.average(data[i::12,:,:], axis=0)
    return monMean


           
def lagged_corr_map_rmSeason(X,field,lag): 
    '''계절변동성 제거한 데이터 전용. X라는 것으로 field(SST,SLP...)의 변화 살피기'''
    corr = np.full([field.shape[1],field.shape[2]],1e3);
    for i in range(field.shape[1]):
        for j in range(field.shape[2]):
            corr[i,j] = np.corrcoef(field[lag:,i,j],X[:-lag])[0,1] 
        if np.mod(i,100)==0:
            print(i)
    return corr

    #np.save('/home/dhkim/data/research/QBO/npy/ERA5_corr/lag/IOB_lag'+str(lag)+'_SST_12mon_corr.npy',SST_corr)