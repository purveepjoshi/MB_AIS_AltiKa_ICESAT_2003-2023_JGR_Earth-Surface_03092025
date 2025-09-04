# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 11:12:33 2023

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

#%%

"""
ANTARCTICA
"""

#%% Conversion of NetCDF Files of SARAL Data to Shapefile over Antarctica...

path = r"H:\SARAL Data"
lst = os.listdir(path)
lst.sort()

for file in lst[104:105]:
    print(file)
    directory=path+"//"+file
    files = glob.glob(directory+"//*.nc")
    files.sort()
    # sdt = os.path.basename(files[0]).split("_")[4]
    # edt = os.path.basename(files[len(files)-1]).split("_")[4]
    # if os.path.exists(path+"/shapefile_winter_2023_SRL/"+file)==False:
    #     os.makedirs(path+"/shapefile_winter_2023_SRL/"+file)
    df5 = pd.DataFrame()
    for i in files:
        print(i)
        try:
            ft = Dataset(i)
            lat = ft.variables['lat_40hz'][:]
            lon = ft.variables['lon_40hz'][:]
            alt = ft.variables["alt_40hz"][:]
            rnge = ft.variables["ice2_range_40hz"][:]
            dtc = ft.variables['model_dry_tropo_corr'][:]
            dtc = cv2.resize(dtc,(1,rnge.shape[0]*rnge.shape[1]),interpolation=cv2.INTER_LINEAR)
            wtc = ft.variables['model_wet_tropo_corr'][:]
            wtc = cv2.resize(wtc,(1,rnge.shape[0]*rnge.shape[1]),interpolation=cv2.INTER_LINEAR)
            ion_corr = ft.variables['iono_corr_gim'][:]
            ion_corr = cv2.resize(ion_corr,(1,rnge.shape[0]*rnge.shape[1]),interpolation=cv2.INTER_LINEAR)
            sol_ear_tide = ft.variables['solid_earth_tide'][:]
            sol_ear_tide = cv2.resize(sol_ear_tide,(1,rnge.shape[0]*rnge.shape[1]),interpolation=cv2.INTER_LINEAR)
            pole_tide = ft.variables['pole_tide'][:]
            pole_tide = cv2.resize(pole_tide,(1,rnge.shape[0]*rnge.shape[1]),interpolation=cv2.INTER_LINEAR)
            ice2_sigma = ft.variables['ice2_sig0_40hz'][:]
            sigma = ft.variables['sig0_40hz'][:]
            
            df = pd.DataFrame({'Latitude':lat.flatten(),
                                'Longitude':lon.flatten(),
                                'Altitude':alt.flatten(),
                                'Range':rnge.flatten(),
                                'DTC':dtc.flatten(),
                                'WTC':wtc.flatten(),
                                'IC':ion_corr.flatten(),
                                'SET':sol_ear_tide.flatten(),
                                'PT':pole_tide.flatten(),
                                'Ice2_Sigma0':ice2_sigma.flatten(),
                                'Sigma0':sigma.flatten()})
            
            df['Elevation']=df['Altitude']-df['Range']-df['DTC']-df['WTC']-df['IC']-df['SET']-df['PT']
            
            df1=df[(df['Longitude']>=180)]
            df2=df[(df['Longitude']<=180)]
            df1['Longitude']=df1['Longitude']-360
            
            df4 = pd.concat([df1,df2],axis=0)
            
            df4 = df4[df4['Latitude']<=-59]
            
            Boundary = gpd.read_file(r"I:\44th_Antarctic_Expedition\Validation\polgons\Maitri_polygon.shp")
            df_clip = gpd.GeoDataFrame(df4)
            df_clip.set_geometry(gpd.points_from_xy(df_clip["Longitude"], df_clip["Latitude"]), inplace = True, crs = "EPSG:4326")
            df_clip = df_clip.to_crs(3412)
            df_1 = gpd.clip(df_clip,Boundary)

            
            # df5 = pd.concat([df5,df4],axis=0)
            
            # BELOW code is for saving individual passes into shapefile
            df3 = gpd.GeoDataFrame(df_1)
            df3.set_geometry(gpd.points_from_xy(df3["Longitude"], df3["Latitude"]), inplace = True, crs = "EPSG:4326")
            df3 = df3.to_crs(3412)
            df3.to_file(r"H:\SARAL Data\SARAL_shapefile_Maitri/"+file+"//"+os.path.basename(i)[:-8] + ".shp", driver = "ESRI Shapefile")
        except:
            print(i)

    # df3 = gpd.GeoDataFrame(df5)
    # df3.set_geometry(gpd.points_from_xy(df3["Longitude"], df3["Latitude"]), inplace = True, crs = "EPSG:4326")
    # df3 = df3.to_crs(3412)
    # df3.to_file(path+"/shapefile_winter_2023_SRL/"+file+"//"+ "AIS_continental_ice_3412_%s_%s.shp"%(sdt,edt), driver = "ESRI Shapefile")

#%% Merging of all cycles for individual year 

path = r"I:\AltiKa\SARAL Data Shapefile - Antarctica\Merged_Cycles\2013"
text = "AIS_continental_ice_3412_20130314_20140123"

df1 = pd.DataFrame()

for root, dirs, files in sorted(os.walk(path)):
    for file in files:
        if file.endswith('.shp'):
            print(file)
            p = os.path.join(root,file)
            df = gpd.read_file(p)
            df1 = pd.concat([df1,df],axis=0)
            del df

df1 = df1.set_index(np.arange(0,df1.shape[0],1))

df3 = gpd.GeoDataFrame(df1)
df3.set_geometry(gpd.points_from_xy(df3["Longitude"], df3["Latitude"]), inplace = True, crs = "EPSG:4326")
df3 = df3.to_crs(3412)
df3.to_file(path + "//" + text + ".shp", driver = "ESRI Shapefile")

#%% Direct Method

dem = rs.open(r"H:\NSIDC_Antarctica_DEM\DEM\NSIDC_Ant500m_wgs84_elev_m.tif")
slp = rs.open(r"H:\NSIDC_Antarctica_DEM\Slope\NSIDC_Ant500m_wgs84_slope_m.tif")

dem_array = dem.read(1)

x0 = dem.read_transform()[0]
y0 = dem.read_transform()[3]

latitude = np.zeros(dem.shape)
for i in range(dem.shape[0]):
    latitude[i,:] = (y0 - 250) + i*-500

longitude = np.zeros(dem.shape)
for j in range(dem.shape[1]):
    longitude[:,j] = (x0 + 250) + j*500

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

p = r"I:\AltiKa\SARAL Data Shapefile - Antarctica\Merged_Cycles\2013_Merged.shp"

df = gpd.read_file(p)
print(df.columns)
try:
    df = df.drop(columns={'layer', 'path'})
except:
    pass

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

df2 = df1[(df1['Slope_1']>=0)&(df1['Slope_1']<=1)]
df2 = df2.dropna()

df2 = df2[(df2['ELE_DM1']>=0)&(df2['ELE_DM1']<=4093)]

rmse1 = math.sqrt(mean_squared_error(df2['DEM_1'], df2['ELE_DM1']))
bias1 = np.nanmean(np.abs(df2['DEM_1']- df2['ELE_DM1']))

print(bias1)
print("\n")
print(rmse1)

df3 = gpd.GeoDataFrame(df2)
df3.set_geometry(gpd.points_from_xy(df3["Longitude"], df3["Latitude"]), inplace = True, crs = "EPSG:4326")
df3 = df3.drop(columns={'Constants'})
df3 = df3.to_crs(3412)
df3.to_file(r"I:\H_Hard_drive\SARAL_DATA\Merged_Cycles\167_shapefile.shp", driver = "ESRI Shapefile")

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

df4 = df2[df2['ELE_DM1']>=0][['Lon','Lat','ELE_DM1']]
knew_xy_coord = df4[['Lon', 'Lat']].values
knew_values1 = df4['ELE_DM1'].values
result1 = griddata(points=knew_xy_coord, values=knew_values1, xi=(xx, yy), method='nearest')
with rs.open(r"I:\AltiKa\SARAL Data Shapefile - Antarctica\Merged_Cycles\2013\AIS_continental_ice_3412_20130314_20140128.tif", 
             'w', driver='GTiff', dtype='float64',width=m, height=n, count=1, 
              crs='EPSG:3976',transform=from_origin(ln0, lt0, 500, 500), nodata=0) as dataset:
    dataset.write(result1,1)

result1 = cv2.medianBlur(result1,5)
with rs.open(r"I:\AltiKa\SARAL Data Shapefile - Antarctica\Merged_Cycles\2013\AIS_continental_ice_3412_20130314_20140128_Median_Filter.tif",
             'w', driver='GTiff', dtype='float64',width=m, height=n, count=1, 
              crs='EPSG:3976',transform=from_origin(ln0, lt0, 500, 500), nodata=0) as dataset:
    dataset.write(result1,1)


"""
CLIPPING OF INTERPOLATED IMAGES
"""

shp_clip = r"H:\antarctica_boundary-shapefile\Antarctica_3412.shp"
for band in sorted(glob.glob(r"I:\AltiKa\SARAL Data Shapefile - Antarctica\Merged_Cycles\2013/"+"*.tif")):
    print(band)
    options = gdal.WarpOptions(dstSRS="EPSG:3412",xRes=500,yRes=500,cutlineDSName=shp_clip,cropToCutline=True,dstNodata=np.nan)
    outBand = gdal.Warp(srcDSOrSrcDSTab=band,
                        destNameOrDestDS=band[:-4]+"_CLIP.tif",
                        options=options)
    outBand = None

    os.remove(band)
    os.rename(band[:-4]+"_CLIP.tif",band)

#%% Clipping of a shapefile over a polygon
 

df = pd.DataFrame()
polygon = gpd.read_file(r"I:\Rignot_Basins_AIS-GrIS\Antarctica\ANT_Basins_IMBIE2_v1.6_Without_Iceland_AltiKa_Covered_Rignot_Basins_EA_3412\ANT_Basins_IMBIE2_v1.6_Without_Iceland_AltiKa_Covered_Rignot_Basins_EA_04_D-Dp_3412.shp")
for i in sorted(glob.glob(r"I:\DGPS_Insitu_Data_AIS\cycle_167\cycle_167"+"/*.shp")):
    print(i)
    saral = gpd.read_file(i)
    clipped = gpd.clip(saral,polygon)
    df = pd.concat([df,clipped],axis= 0)
    
df.to_file(r"I:\DGPS_Insitu_Data_AIS\clipped_saral_c-cp\clipped_C-cp_Saral_167.shp")
