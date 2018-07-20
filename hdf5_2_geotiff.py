###############################################################################
# imports
###############################################################################
from __future__ import division
import numpy as np
from glob import glob
import h5py
import gdal
from osgeo import osr
import os
from pyproj import Proj, transform

gdal.AllRegister()
###############################################################################
# read data
###############################################################################
#dir_name = r"f:\data\smap\download\SMAP_Twente\\hdf5\\"
#fNames = glob(dir_name+"*.h5")

#fNames = [r'SMAP_L3_SM_P_20180703_R16010_001.h5']
#fNames = [r'SMAP_L3_SM_P_E_20180703_R16010_001.h5']

fNames = glob('*.h5')
for i in fNames:
    with h5py.File(i, mode="r") as f:
        name = '/Soil_Moisture_Retrieval_Data_AM/soil_moisture'
        data = f[name][:]
        units = f[name].attrs['units']
        longname = f[name].attrs['long_name']
        _FillValue = f[name].attrs['_FillValue']
        valid_max = f[name].attrs['valid_max']
        valid_min = f[name].attrs['valid_min']
        invalid = np.logical_or(data > valid_max,
                                data < valid_min)
        invalid = np.logical_or(invalid, data == _FillValue)
        data[invalid] = np.nan
        data =  np.ma.masked_where(np.isnan(data), data)
        
        # get geodata
        latitude = f['/Soil_Moisture_Retrieval_Data_AM/latitude'][:]   #y
        longitude = f['/Soil_Moisture_Retrieval_Data_AM/longitude'][:] #x

        latitude_masked = np.ma.masked_where(latitude==-9999, latitude)
        longitude_masked = np.ma.masked_where(longitude==-9999, longitude)

#%%
        # reproject coordinates


        inProj = Proj(init='epsg:4326')
#        outProj = Proj(init='epgs:6933')
        outProj = Proj('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')
        
        
        
        longitude,latitude = transform(inProj, outProj, longitude_masked, latitude_masked)
        
        
        
        # setup geotransform
        num_cols = data.shape[1]
        num_rows = data.shape[0]        
        
        xmin = longitude[longitude!=-9999].min()
        xmax = longitude[longitude!=-9999].max()
        ymin = latitude[latitude!=-9999].min()
        ymax = latitude[latitude!=-9999].max()
    
        xres = (xmax-xmin)/num_cols
        yres = (ymax-ymin)/num_rows
        
        geotransform = (xmin, xres, 0, ymax, 0, -yres)
        
        #        # create geotiff
        directory, filename = os.path.split(i)
        filename_ext = os.path.splitext(filename)[0]
        driver = gdal.GetDriverByName("Gtiff")
        
#        src_ds = gdal.Open("SMAP_L3_SM_P_E_20150401_R14010_001.tif")
        
#        raster = driver.CreateCopy(r"geotiffs\\"+filename_ext+".tif", src_ds, strict=0)
#                src_ds = None

        raster = driver.Create(os.path.splitext(i)[0] + ".tif",
                               int(num_cols), int(num_rows),
                               1, gdal.GDT_Float32)
                
        # set geotransform and projection
        raster.SetGeoTransform(geotransform)
        srs = osr.SpatialReference()
#        srs.ImportFromEPSG(4326)
        wkt_text = r'PROJCS["unnamed",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9108"]],AUTHORITY["EPSG","4326"]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["standard_parallel_1",30],PARAMETER["central_meridian",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1],AUTHORITY["epsg","6933"]]'
        
        srs.ImportFromWkt(wkt_text)
        
         
        raster.SetProjection(srs.ExportToWkt())
        
        # write data array to raster
        raster.GetRasterBand(1).WriteArray(data.data)
        raster.GetRasterBand(1).SetNoDataValue(-9999)
        raster.FlushCache()
        raster = None        