# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 11:18:59 2023

@author: Glacier
"""

#%% Import different libraries..

"Importing Libraries..."

from zipfile import ZipFile
import os
import glob
from osgeo import gdal
import numpy as np
import datetime as dt
import numpy as np
import datetime as dt
import glob
from osgeo import gdal
import cv2
from matplotlib import pyplot as plt
from matplotlib import image as mpimg
from osgeo import gdalconst
from osgeo import gdalnumeric
import pandas as pd
import os
import math
import sys
sys.path
sys.path.append('C:\ProgramData\Anaconda3\Scripts')
# from osgeo import gdal_merge as gm
import subprocess
from netCDF4 import Dataset
import netCDF4
from netCDF4 import MFDataset
from osgeo import ogr
import json
import seaborn as sn
from mpl_toolkits.basemap import Basemap
from numpy import linspace
from numpy import meshgrid
import calendar
import shutil
import geopandas as gpd
from sklearn.metrics import mean_squared_error
import rasterio as rs
from scipy.interpolate import griddata
import earthpy.plot as ep
import earthpy.clip as ec
from rasterio.transform import from_origin
from haversine import haversine, Unit

#%% Merging of shapefiles of SARAL Data over Antarctica Region...

## MERGING OF CYCLES FOR EACH YEAR, Eg. 2013, 2014, 2015....

path=r"H:\SARAL Data Shapefile - Antarctica\Merged_Cycles"
lst=os.listdir(path)
lst.sort()

for file in lst[0:9]:
    print(file)
    directory=path+"//"+file
    files=glob.glob(directory+"//"+"*.shp")
    os.makedirs(r"H:\SARAL Data Shapefile - Antarctica\Merged_Cycles\2013/"+file)
    df1=pd.DataFrame(columns={'Latitude', 'Longitude', 'Altitude', 'Range', 'DTC', 'WTC', 'IC',
       'SET', 'PT', 'Ice2_Sigma', 'Sigma0', 'Elevation','geometry'})
    for i in files:
        # try:
        print(i)
        df=gpd.read_file(i)
        # df=df.drop(columns={'geometry'})
        df1 = pd.concat([df1,df],axis=0)
        # except:
            # print(i)
    df3 = gpd.GeoDataFrame(df1)
    # df3.set_geometry(gpd.points_from_xy(df3["Longitude"], df3["Latitude"]), inplace = True, crs = "EPSG:4326")
    df3.set_geometry(df3['geometry'], inplace = True, crs = "EPSG:4326")
    df3.to_file(r"H:\SARAL Data Shapefile - Antarctica\Merged_Cycles\2013/"+file+"//"+file+"2013_merged.shp", driver = "ESRI Shapefile")



#%% Direct Method

dem = rs.open(r"E:\New folder\Python\NSIDC_Ant500m_wgs84_elev_m.tif")
slp = rs.open(r"E:\New folder\Python\NSIDC_Ant500m_wgs84_slp_m.tif")
asp = rs.open(r"H:\NSIDC_Antarctica_DEM\NSIDC_Ant500m_wgs84_Aspect.tif")

dem_array = dem.read(1)

x0 = dem.read_transform()[0]
y0 = dem.read_transform()[3]

latitude = np.zeros(dem.shape)
for i in range(dem.shape[0]):
    latitude[i,:] = (y0 - 250) + i*-500

longitude = np.zeros(dem.shape)
for j in range(dem.shape[1]):
    longitude[:,j] = (x0 + 250) + j*500

#%%

def dm_method(array1,minx,maxx,miny,maxy):
    
    P_dem = array1[minx:maxx,miny:maxy]
    
    P_Lon = np.zeros(P_dem.shape)
    for j in range(int(P_dem.shape[0]/2)*-1,int(P_dem.shape[0]/2)+1):
        P_Lon[:,j+int(P_dem.shape[0]/2)] = j*500
    
    P_Lat = np.zeros(P_dem.shape)
    for j in range(int(P_dem.shape[0]/2)*-1,int(P_dem.shape[0]/2)+1):
        P_Lat[j+int(P_dem.shape[0]/2),:] = j*500
    
    X,Y,Z = P_Lon,P_Lat,P_dem
    
    A = np.ones((len(X.flatten()),6))
    A[:,1] = X.flatten()
    A[:,2] = Y.flatten()
    A[:,3] = A[:,1]**2
    A[:,4] = A[:,1]*A[:,2]
    A[:,5] = A[:,2]**2

    B =np.ones((len(Z.flatten()),1))
    B[:,0] = Z.flatten()

    A_T = np.matrix.getT(A)
    AT_A = np.matmul(A_T,A)
    A_INV = np.linalg.inv(AT_A)
    
    C = np.matmul(np.matmul(A_INV,A_T),B)
    theta = np.sqrt((C[1,0]**2)+(C[2,0]**2))
    h_bar = C[0,0]
    intercept = (2*(((C[1,0]**2)*C[3,0])+(C[1,0]*C[2,0]*C[4,0])+((C[2,0]**2)*C[5,0])))/(theta**2)
    
    return theta, h_bar, intercept

path = r"E:\New folder\cycle_008_shapefile"
lst = glob.glob(path+"/*.shp")
lst.sort()

name=[]
r1,r2,r3,r4 = [],[],[],[]
b1,b2,b3,b4 = [],[],[],[]

for p in lst:
    # p = r"E:\New folder\cycle_008_shapefile\cycle_008_merged\cycle_008.shp"
    print(p)
    df = gpd.read_file(p)
    # print(df.columns)
    # try:
    #     df = df.drop(columns={'DEM_1', 'Slope_1', 'Aspect_1', 'Lat', 'Lon', 'Slope_Angl', 'Curvature', 
    #                           'h_bar', 'hdm', 'Rc', 'ELE_DM1', 'ELE_DM2', 'ELE_DM3'})
    # except:
    #     pass
    
    # lake = gpd.read_file(r"E:\New folder\Vostok_Lake\vostok_lake.shp")
    
    # df = gpd.clip(df,lake)
    
    df1 = df.dropna()
    
    df1 = df1.set_index(np.arange(0,df1.shape[0],1))
    
    try:
        df1 = df1.to_crs(epsg=3976)
        df1['Lat'] = df1.geometry.y
        df1['Lon'] = df1.geometry.x
    except:
        df1 = df1.drop(columns={"geometry"})
        df4 = gpd.GeoDataFrame(df1)
        df4.set_geometry(gpd.points_from_xy(df4["Longitude"], df4["Latitude"]), inplace = True, crs = "EPSG:4326")
        df1 = df4
        df1 = df1.to_crs(epsg=3976)
        df1['Lat'] = df1.geometry.y
        df1['Lon'] = df1.geometry.x
        
    df1['Row']=dem.index(df1.Lon,df1.Lat)[0]
        
    df1['Col']=dem.index(df1.Lon,df1.Lat)[1]
            
    df1 = df1[(df1['Row']>=0)&(df1['Row']<dem_array.shape[0])&(df1['Col']>=0)&(df1['Col']<dem_array.shape[1])]
    
    df1['DEM_1'] = dem.read(1)[df1['Row'],df1['Col']]
    df1['Slope_1'] = slp.read(1)[df1['Row'],df1['Col']]
        
    df1['DEM_Lat'] = latitude[df1['Row'],df1['Col']]
    
    df1['DEM_Lon'] = longitude[df1['Row'],df1['Col']]
    
    df1['Distance'] = np.sqrt((df1['Lon']-df1['DEM_Lon'])**2+(df1['Lat']-df1['DEM_Lat'])**2)
    
    df1['Pixel'] = df1.apply(lambda row: 1 if round(row['Distance']/500)==0 else round(row['Distance']/500), axis=1)
    
    df1['X0'] = df1['Row']-df1['Pixel']
    
    df1['Xn'] = df1['Row']+df1['Pixel']+1
    
    df1['Y0'] = df1['Col']-df1['Pixel']
    
    df1['Yn'] = df1['Col']+df1['Pixel']+1
    
    df1['Xn'] = pd.to_numeric(df1['Xn'],downcast="integer")
    
    df1['X0'] = pd.to_numeric(df1['X0'],downcast="integer")
    
    df1['Y0'] = pd.to_numeric(df1['Y0'],downcast="integer")
    
    df1['Yn'] = pd.to_numeric(df1['Yn'],downcast="integer")
    
    df1 = df1[(df1['X0']>0)&(df1['Xn']<dem_array.shape[0])&(df1['Y0']>0)&(df1['Yn']<dem_array.shape[1])]

    df1['Constants'] = df1.apply(lambda row: dm_method(dem_array,row['X0'],row['Xn'],row['Y0'],row['Yn']),axis=1)
    
    df1['Slope_Angle'] = df1.apply(lambda row: row["Constants"][0],axis=1)
    
    df1['h_bar'] = df1.apply(lambda row: row["Constants"][1],axis=1)
    
    df1['Curvature'] = df1.apply(lambda row: row["Constants"][2]/100,axis=1)
    
    df1['hdm'] = (df1['Range']*(df1['Slope_Angle']**2))/(2*(1-(df1['Range']*(df1['Curvature']-(1/(6371000+df1['h_bar']))))))
    
    df1['Rc'] = df1['Range'] + df1['hdm']
    
    df1['ELE_DM1'] = df1['Altitude'] - df1['Rc'] - df1['DTC']-df1['WTC']-df1['IC']-df1['SET']-df1['PT']
    
    # df1['Difference'] = df1['ELE_DM1'] - df1['DEM_1']

    df2 = df1[(df1['Slope_1']>=0)&(df1['Slope_1']<=0.85)]
    df2 = df2.dropna()

    
    rmse1 = math.sqrt(mean_squared_error(df2['DEM_1'], df2['ELE_DM1']))
    
    bias1 = np.nanmean(np.abs(df2['DEM_1']- df2['ELE_DM1']))

    print(bias1)

    name.append(os.path.basename(p)[:-4])
    
    b2.append(bias1)
    
    r2.append(rmse1)
    
    df3 = gpd.GeoDataFrame(df1)
    df3.set_geometry(gpd.points_from_xy(df3["Longitude"], df3["Latitude"]), inplace = True, crs = "EPSG:4326")
    df3 = df3.drop(columns={'Ice2_Sigma', 'Sigma0','Constants'})
    df3.to_file(r"E:\New folder\cycle_008_shapefile\Corrected_Elevation/"+os.path.basename(p), driver = "ESRI Shapefile")
    
df_stats = pd.DataFrame({"Name":name,"DM_Bias":b1,"DM1_Bias":b2,"DM2_Bias":b3,"DM3_Bias":b4,"DM_RMSE":r1,"DM1_RMSE":r2,"DM2_RMSE":r3,"DM3_RMSE":r4})
df_stats.to_excel(r"E:\New folder\Stats_1.xlsx")

#%%

# ln0,lnx,lt0,lty = df2.DEM_Lon.min()-250,df2.DEM_Lon.max()+250,df2.DEM_Lat.max()+250,df2.DEM_Lat.min()-250
# m,n = int((lnx-ln0)/500),int((lt0-lty)/500)

# yy = np.zeros((n,m))
# for i in range(n):
#     yy[i,:] = (lt0 - 250) + i*-500

# xx = np.zeros((n,m))
# for j in range(m):
#     xx[:,j] = (ln0 + 250) + j*500

ln0,lt0 = x0,y0
xx,yy = longitude, latitude
n,m = dem.shape

df4 = df3[df3['ELE_DM1']>=0][['Lon','Lat','ELE_DM1']]
knew_xy_coord = df4[['Lon', 'Lat']].values
knew_values1 = df4['ELE_DM1'].values
result1 = griddata(points=knew_xy_coord, values=knew_values1, xi=(xx, yy), method='nearest')
with rs.open(r'E:\New folder\Images\AIS_DM1.tif', 'w', driver='GTiff', dtype='float64',width=m, height=n, count=1, 
             crs='EPSG:3976',transform=from_origin(ln0, lt0, 500, 500), nodata=0) as dataset:
    dataset.write(result1,1)
    
# df4 = df3[df3['ELE_DM2']>=0][['Lon','Lat','ELE_DM2']]
# knew_xy_coord = df4[['Lon', 'Lat']].values
# knew_values1 = df4['ELE_DM2'].values
# result2 = griddata(points=knew_xy_coord, values=knew_values1, xi=(xx, yy), method='nearest')
# with rs.open(r'E:\New folder\Images\AIS_DM2.tif', 'w', driver='GTiff', dtype='float64',width=m, height=n, count=1, 
#              crs='EPSG:3976',transform=from_origin(ln0, lt0, 500, 500), nodata=0) as dataset:
#     dataset.write(result2,1)
    
img1 = rs.open(r"E:\New folder\Images\AIS_DM1.tif").read(1)
# img2 = rs.open(r"E:\New folder\Images\AIS_DM2.tif")


img2 = rs.open(r"E:\New folder\Python\NSIDC_Ant500m_wgs84_elev_m.tif").read(1)

difference = img1 - img2


a1 = gdalnumeric.SaveArray(difference,r"E:\New folder\Images\difference.tif",format="GTiff",prototype=r"E:\New folder\Images\AIS_DM1.tif")
a1 = None



# ais = gpd.read_file(r"E:\New folder\Shapefile\GroundingLine_Antarctica_v02.shp")
# ais = ais.to_crs(4326)

# atm_all = gpd.read_file(r"E:\New folder\ATM\ILATM2_20131118_20131128_smooth_nadir3seg_50pt.shp",mask=ais)
# print(atm_all.columns)

# atm = atm_all

# # atm = gpd.clip(atm_all,lake)

# atm = atm.to_crs(3976)

# atm['Lat'] = atm.geometry.y
# atm['Lon'] = atm.geometry.x

# atm['DEM_1'] = dem.read(1)[dem.index(atm.Lon,atm.Lat)[0],dem.index(atm.Lon,atm.Lat)[1]]
# atm['Slope_1'] = slp.read(1)[slp.index(atm.Lon,atm.Lat)[0],slp.index(atm.Lon,atm.Lat)[1]]
# atm['ELE_DM1'] = img1.read(1)[img1.index(atm.Lon,atm.Lat)[0],img1.index(atm.Lon,atm.Lat)[1]]
# # atm['ELE_DM2'] = img2.read(1)[img2.index(atm.Lon,atm.Lat)[0],img2.index(atm.Lon,atm.Lat)[1]]

# atm = atm[(atm['Slope_1']>=0)&(atm['Slope_1']<=1)]
# atm = atm.dropna()

# # atm = atm[((atm['DEM_1']-atm['ELE_DM1'])>=-1)&((atm['DEM_1']-atm['ELE_DM1'])<=1)&
# #           ((atm['DEM_1']-atm['ELE_DM2'])>=-1)&((atm['DEM_1']-atm['ELE_DM2'])<=1)&
# #           ((atm['DEM_1']-atm['ELE_DM3'])>=-1)&((atm['DEM_1']-atm['ELE_DM3'])<=1)&
# #           ((atm['DEM_1']-atm['Corr_Ele'])>=-1)&((atm['DEM_1']-atm['Corr_Ele'])<=1)]

# rm = math.sqrt(mean_squared_error(atm[' WGS84_Ell'], atm['DEM_1']))
# rm1 = math.sqrt(mean_squared_error(atm[' WGS84_Ell'], atm['ELE_DM1']))
# # rm2 = math.sqrt(mean_squared_error(atm[' WGS84_Ell'], atm['ELE_DM2']))
# bs = np.nanmean(np.abs(atm[' WGS84_Ell']- atm['DEM_1']))
# bs1 = np.nanmean(np.abs(atm[' WGS84_Ell']- atm['ELE_DM1']))
# # bs2 = np.nanmean(np.abs(atm[' WGS84_Ell']- atm['ELE_DM2']))
# print(bs)
# print(bs1)
# # print(bs2)


