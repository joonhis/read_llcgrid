#'''원하는 시간대만 골라내기'''
#files[153:181]
files_Feb=[]
for i in range(len(files)):
    DATE = datetime.date(int(files[i][34:38]),int(files[i][38:40]),int(files[i][40:42]))#.strftime('%Y%m%d')
    if DATE.month == 2:# or DATE.month == 1:
        files_Feb.append(files[i])
files_Feb.append('/home/dhkim/data/QBOdata/mur_data/20020901090000-JPL-L4_GHRSST-SSTfnd-MUR25-GLOB-v02.0-fv04.2.nc')
files_Feb ### 가짜 데이터 마지막줄에 추가



'''daliy mean: datetime으로 month넘어가는지 확인
- 대륙의 값이 해양 최솟값보다 커서 snapshot에서는 마스킹 불가=> sum해서 마스킹
- np.array 붙이면 마스킹이 깨져버린다!!!'''

def daily_to_monthly(filelist, ):
    sst=[]; sic=[]; # sst: monthly
    daily_sst=[]; daily_sic=[]

    f = Dataset(filelist[0],'r')
    date_test = nc.num2date(f.variables['time'][0],f.variables['time'].units) 
    
    for i in range(len(filelist)):
        f = Dataset(filelist[i],'r')
        file_date = nc.num2date(f.variables['time'][0],f.variables['time'].units) #읽어줘야 편하다.
        
        if file_date.year != date_test.year:## 한달치를 모아 월평균
            if mask == True:
                sst.append(np.ma.masked_less(np.mean(daily_sst,axis=0),-10))
            else:
                sst.append(np.mean(daily_sst,axis=0))
                daily_sst=[]; daily_sic=[]    
            
        
        daily_sst.append(f['analysed_sst'][0])
        daily_sic.append(f['sea_ice_fraction'][0])
        date_test = file_date
        f.close()
    return sst, sic


'''시간대 분리: 한달만 할때.'''
f = Dataset(fdir+'HadISST_sst.nc','r')
f2 = Dataset(fdir+'HadISST_ice.nc','r')
sst1=[]; sst=[];sic1=[]; sic=[]

index_mm = 5
for i in range(len(year)-2): #### len year = 40까지
    sst1.append(f['sst'][109*12+i*12+index_mm]) ## 79년 1월부터
    sic1.append(f2['sic'][109*12+i*12+index_mm])
    
if index_mm <=4: #### 5월 [4], len year = 41까지
    sst1.append(f['sst'][109*12+(len(year)-1)*12+index_mm])
    sic1.append(f2['sic'][109*12+(len(year)-1)*12+index_mm])
    
sst = np.ma.masked_less(sst1,-20)
sic = np.ma.masked_less(sic1,-0.1)

len(sst)
#''' 여태까지 20년 10월부터 어떻게 된거지? == 20년은 append로 붙이기'''

'''파일 없는부분 검출기!! .'''
DATEdelta = datetime.date(int(files[0][31:35]),int(files[0][35:37]),int(files[0][37:39]))
for i in range(len(files)):
    DATE = datetime.date(int(files[i][31:35]),int(files[i][35:37]),int(files[i][37:39]))
    if DATE != DATEdelta:
        print(DATEdelta)
        DATEdelta = DATE
    DATEdelta = DATEdelta + datetime.timedelta(days=1)
    
    
# 여러달 평균할때 (삭제 금지)
'''시간대 분리: index만 나눠주면된다.'''
sst1=[]; sst=[];sic1=[]; sic=[]
for i in range(len(year)):
    sst1.append(f['sst'][110*12+i*12-3])   # 10
    sst1.append(f['sst'][110*12+i*12-2])  # 11
#    sst1.append(f['sst'][110*12+i*12-1])  # 12
#    sst1.append(f['sst'][110*12+i*12])    # 1월
    sst.append(np.mean(sst1,axis=0))
    sst1=[]
    sic1.append(f2['sic'][110*12+i*12-3])
    sic1.append(f2['sic'][110*12+i*12-2])
#    sic1.append(f2['sic'][110*12+i*12-1])
#    sic1.append(f2['sic'][110*12+i*12])
    sic.append(np.mean(sic1,axis=0))
    sic1=[]
sst = np.ma.masked_less(sst,-20)
sic = np.ma.masked_less(sic,-0.1)

len(sst)

