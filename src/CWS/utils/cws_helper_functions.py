import logging
from loguru import logger
import os
import rasterio
import rasterio.mask
from rasterio.windows import from_bounds
from rasterio.windows import Window
from rasterio import windows
from pathlib import Path
from osgeo import gdal
import numpy as np
from osgeo import gdal
import subprocess
from osgeo import osr
import pandas as pd
import glob
import geopandas as gpd
from rasterstats import zonal_stats
import datetime
import math

import pyodbc 
import sqlalchemy as sa
from sqlalchemy import create_engine, event
from sqlalchemy.engine.url import URL
import json




#https://stackoverflow.com/questions/21257899/writing-a-csv-file-into-sql-server-database-using-python
def greenmonkey(in_df, out_db, schema_name): 
       
        engine_greenmonkey = sa.create_engine('mssql+pyodbc://' + "GREENMONKEY" + '/' + "Carbon_Mapping" + '?trusted_connection=yes&driver=ODBC+Driver+17+for+SQL+Server')
        connection = engine_greenmonkey.raw_connection()
        cursor = connection.cursor()
        connection.commit()
        #export
        #########################################
        input_dataframe=in_df
        name_of_out_table =out_db
        schema_name = schema_name
        #########################################
        input_dataframe.to_sql(name_of_out_table, engine_greenmonkey, if_exists='replace', index = False, schema=schema_name)
        cursor.close()
        connection.commit()
    

        print ("----------------------------------------------------------------")































