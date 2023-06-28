import os
from pathlib import Path, WindowsPath
import xarray as xr
import rioxarray as rio
import pandas as pd
from datetime import datetime
import glob
import pymannkendall as mk

# Calculate slope and p-value of Mann-Kendall test
# 1) Load all yearly LINT files
# 2) 


""" Calculate slope and p-value of Mann-Kendall test

aoi_name: area of interest (must be defined in aoi_coords). Necessary in order to avoid memory error
aoi_coords: bounding box coordinates of AOI
base_path: path object (path to folder that contains LINT folder)
base1, base2: start and end of baseline period
"""
start_year = 2014
end_year = 2021 # last year for which we have data

base_path = WindowsPath('//cwsfileserver/projects/lulucf/f02_data')
base_path = os.path.join(base_path, 'GDMP')
base1 = 2014
base2 = 2021

aoi_bbox = dict(
        west=(1500000.0000, 900000.0000, 1500000.0000+(7400000.0000-1500000.0000)/2,5500000.0000 ),
        east=(1500000.0000+(7400000.0000-1500000.0000)/2,  900000.0000, 7400000.0000,5500000.0000 ),
        # luxembourg = (4014674.4725, 2933830.0285, 4071341.1143999994, 3015531.6271) ## used for testing
        )

# helper function to compute slope and p-value of MK test
def mk_trend_slope(vec, verbose=False):
    trend=mk.original_test(vec)
    slope=trend.slope
    p = trend.p
    if verbose:
        print(slope,p)
    return (slope,p)

for aoi_name in list(aoi_bbox):
    print("loading data for aoi area: "+aoi_name)
    aoi_coords = aoi_bbox[aoi_name]
    # AOI: subset of Europe, e.g. Iberian peninsula (in order to get everything into memory)
    # (adding 1x spacing to make sure SMA covers entire AOI)
    spacing = 1000
    xmin = aoi_coords[0] - spacing
    ymin = aoi_coords[1] - spacing
    xmax = aoi_coords[2] + spacing
    ymax = aoi_coords[3] + spacing

    baseline_period = (base1, base2)
    
    #################################################################
    #################################################################
    # select tiffs
    resolution = "1km" ### change this to point to the right folder
    in_path = os.path.join(base_path, f'GDMP{resolution}')

    out_path = os.path.join(base_path, f'GDMP{resolution}_mk_{aoi_name}')
    if not os.path.exists(out_path):
        os.mkdir(out_path)
    #################################################################
    #################################################################

    # load GDMP into xarray
    # tif_list = [f for f in os.listdir(in_path) if 'GDMP' in f]
    tif_list = [os.path.basename(f) for f in glob.glob(in_path + '/clip_GDMP_*.tif')]
    # Create variable used for time axis
    time_var = xr.Variable('time', pd.to_datetime([f'{fname[5+5:9+5]}-01-01' for fname in tif_list]))
    # Load in and concatenate all individual GeoTIFFs
    gdmp_ds= xr.concat([xr.open_dataset(os.path.join(in_path, i), engine='rasterio').sel(y=slice(ymax, ymin), x=slice(xmin, xmax)) for i in tif_list], dim=time_var)
    # gdmp_ds= xr.concat([xr.open_dataset(os.path.join(in_path, i), engine='rasterio').sel(y=slice(ymin+10000, ymin), x=slice(xmin, xmin+10000)) for i in tif_list], dim=time_var)
    # Rename the variable to a more useful name
    gdmp_ds = gdmp_ds.rename({'band_data': 'GDMP'})
    gdmp_ds = gdmp_ds*0.02
    print(gdmp_ds)
    print("start computation")
    
    # compute slope and p-value of Mann-Kendall test
    slope,pvalue = xr.apply_ufunc(mk_trend_slope,gdmp_ds,input_core_dims=[['time']], output_core_dims=[[],[]], vectorize=True)

    # compute standard deviation to create mask
    std = gdmp_ds.sel(time=slice(datetime(baseline_period[0], 1, 1), datetime(baseline_period[1], 12, 31))).std(dim='time')
    print("saving results")
    ## slope
    slope = slope.fillna(-999)
    slope = xr.where(std == 0, -999, slope)    # mask pixels where stdev = 0, i.e. anom is inf
    # write tiff file
    slope['GDMP'].rio.write_nodata(-999, inplace=True)
    slope.rio.write_crs("epsg:3035", inplace=True)
    slope['GDMP'].rio.to_raster(os.path.join(out_path, f'GDMP_{resolution}_trend_slope.tif'), compress='LZW')
    
    ## pvalue
    pvalue = pvalue.fillna(-999)
    pvalue = xr.where(std == 0, -999, pvalue)    # mask pixels where stdev = 0, i.e. anom is inf

    # write tiff file
    pvalue['GDMP'].rio.write_nodata(-999, inplace=True)
    pvalue.rio.write_crs("epsg:3035", inplace=True)
    pvalue['GDMP'].rio.to_raster(os.path.join(out_path, f'GDMP_{resolution}_trend_pvalue.tif'), compress='LZW')
    print("end")



