# from utility import month , dx, dy로 변수를 바로 사용하자

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import datetime
####################### 영문 월 이름 계산기 ################################################
m12 = np.arange(1,13) ; month = []
for i in range(12):
    dt = datetime.datetime.strptime(str(m12[i]),'%m')
    month.append(dt.strftime('%B')[:3])
    
####################### dx, dy 거리 구하기 ##############
'''1도당 거리를 해상도에 맞게 또 나눠줘야 한다'''
circum = 2*np.pi*6400000/360 # 1도 단위의 거리
dy = 2*np.pi*6400000/360
def dx(lat):
    dx = 2*np.pi*6400000*np.cos(np.deg2rad(lat))/360
    return dx

# 날짜 계산
from dhkpython import utility
from dateutil.relativedelta import relativedelta
months = utility.month
date = datetime.datetime(1979,1,1)

years = [] ; years_last2=[]
while date < datetime.datetime(2021,1,1):
    years.append(str(date.year)[2:]+'/'+str(date.month))
    if date.month==1:
        years_last2.append(str(date.year)[2:])
    date = date + relativedelta(months=1)
years_last2.append('21')


