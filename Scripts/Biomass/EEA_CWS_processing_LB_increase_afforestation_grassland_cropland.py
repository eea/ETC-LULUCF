

# Reading libaries and start a test-connection to different MS-SQL Server:
import os
import pyodbc 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sqlalchemy as sa
from sqlalchemy import create_engine, event , text
from sqlalchemy.engine.url import URL

import json
#check processing time:##
import time  
import datetime
from configparser import ConfigParser
start_time = time.time()


##########################

### SET conection to MS-sql server:
################################################## SET postgre-sql connection:

################################################## read database keys:
def config(filename, section='EEA-MS-SERVER_GREENMONKEY'):
    # create a parser
    parser = ConfigParser()
    # read config file
    parser.read(filename)

    # get section, default to postgresql
    db = {}
    if parser.has_section(section):
        params = parser.items(section)
        for param in params:
            db[param[0]] = param[1]
    else:
        raise Exception(
            'Section {0} not found in the {1} file'.format(section, filename))

    return db

keys = config(filename='database.ini')
SERVER=                keys['server']
Database_name =        keys['database']
USER =                 keys['user']
PSW =                  keys['password']


DRIVER = "ODBC Driver 17 for SQL Server"
SCHEMA = 'szenario_afforestation'

################################################## read database keys:

print ("connect to engine GREENMONKEY")  
engine_GREENMONKEY = sa.create_engine('mssql+pyodbc://' + SERVER+ '/' + Database_name + '?trusted_connection=yes&driver=ODBC+Driver+17+for+SQL+Server')
print ("connect to greenmonkey engine.....GREENMONKEY and store sql table")  




################################################################################## sql to dataframe:
query_check=('''
             SELECT * 
             FROM [Carbon_Mapping].[szenario_afforestation].[LUT_trees]
             ''')  
with engine_GREENMONKEY.begin() as conn:
    query =   text(query_check) 
    readable_database = pd.read_sql_query(query, conn)
#print(readable_database)
##################################################################################

print ("----------------------------------------------------------------")




#############part1 reading the LB results by NUTS3 region from the scenario:

####set outputfolder of BIOMASS - aff - scenario:
LB_afforestation_results_folder = r'L:\f02_data\carbon_model_data\output\LB_NUTS_stat'



df_list =[]
folder_path = LB_afforestation_results_folder
with os.scandir(folder_path) as entries:
    for entry in entries:
        if entry.is_file():
            file_for_df = folder_path+"\\" + entry.name
            #print(file_for_df)  # Print file names

            ## csv to pandas:
            df = pd.read_csv(file_for_df)    ## csv file to dataframe
            df_list.append(df)   ## append single df to df_list
            

combined_df = pd.concat(df_list, ignore_index=True)    # combine all dataframes
print (combined_df)


#############part2 storing results on MS-SQL:#############################


###################################################
#name_of_table = name_of_table
export_df_to_sql = combined_df  # dataframe to be exported
schmema_name = SCHEMA
name_of_table = 'LB_increase_afforestation_grassland_cropland_nuts3'
###################################################
export_df_to_sql.to_sql(name_of_table, engine_GREENMONKEY,  schema=schmema_name,if_exists='replace')
print ("end storing on SQL")






### QC:#### to check the
'''
Select * from 
  FROM [Carbon_Mapping].[szenario_afforestation].[LB_increase_afforestation_grassland_cropland_nuts3_draft

  where 
  Nuts_ID = 'AT111' and  
  [Year_potential] = 2035 and 
  [Tree_species_factor] ='Fraxinus_excelsior'and 
  [land_use_selection] = 'grassland'and 
  [perc_reforest] = 10
'''

####  update attributes
'''
 CONCAT(  'NUTS_ID_', [NUTS_ID],
	'_Tree_species_',[Tree_species_factor],
	'_perc_reforest_',[perc_reforest],
	'_RCP_',[RCP],
	'_Year_potential_',[Year_potential],
	'_land_use_selection_',[land_use_selection]) as model_parameter
'''


## add:join forest zone to output:
'''
,[FOREST_ZONE]
,[NUTS_LEVEL0_ID]
	into [Carbon_Mapping].[szenario_afforestation].[LB_increase_afforestation_grassland_cropland_nuts3_2]
  FROM [Carbon_Mapping].[szenario_afforestation].[LB_increase_afforestation_grassland_cropland_nuts3]

  left join 
[szenario_afforestation].[LUT_FOREST_ZONE] on[LUT_FOREST_ZONE].[NUTS_LEVEL3_ID] = [LB_increase_afforestation_grassland_cropland_nuts3].[NUTS_ID]


'''


## add: group by tree gourps and FGS:
'''
/****** Script for SelectTopNRows command from SSMS  ******/
drop table if exists [Carbon_Mapping].[szenario_afforestation].[LB_increase_afforestation_grassland_cropland_nuts3_forest_group]
go


SELECT 
     
       [LB_increase_afforestation_grassland_cropland_nuts3].[NUTS_ID]
      ,[LB_increase_afforestation_grassland_cropland_nuts3].[NUTS_LEVEL]
      ,[LB_increase_afforestation_grassland_cropland_nuts3].[Slope_factor]
      ,[LB_increase_afforestation_grassland_cropland_nuts3].[Slope_src]

      ,[LB_increase_afforestation_grassland_cropland_nuts3].[perc_reforest]
      ,[LB_increase_afforestation_grassland_cropland_nuts3].[perc_reforest_src]
      ,[LB_increase_afforestation_grassland_cropland_nuts3].[RCP]
      ,[LB_increase_afforestation_grassland_cropland_nuts3].[Year_potential]
      ,[LB_increase_afforestation_grassland_cropland_nuts3].[FT]
      ,[LB_increase_afforestation_grassland_cropland_nuts3].[FGS]
      ,[LB_increase_afforestation_grassland_cropland_nuts3].[land_use_selection]
      ,[LB_increase_afforestation_grassland_cropland_nuts3].[FOREST_ZONE]
      ,[LB_increase_afforestation_grassland_cropland_nuts3].[NUTS_LEVEL0_ID]

	  ,[LUT_C_SEQ_AFFOR_JRC_V4_GROUPED_FGS].[yrl_C_seq_age_0_20_avg]
      ,[LUT_C_SEQ_AFFOR_JRC_V4_GROUPED_FGS].[yrl_C_seq_age_21_30_avg]
      ,[LUT_C_SEQ_AFFOR_JRC_V4_GROUPED_FGS].[yrl_C_seq_age_0_20_min]
      ,[LUT_C_SEQ_AFFOR_JRC_V4_GROUPED_FGS].[yrl_C_seq_age_21_30_min]
      ,[LUT_C_SEQ_AFFOR_JRC_V4_GROUPED_FGS].[yrl_C_seq_age_0_20_max]
      ,[LUT_C_SEQ_AFFOR_JRC_V4_GROUPED_FGS].[yrl_C_seq_age_21_30_max]


	  ,sum([nr_pixels]) as      [nr_pixels]
  into [Carbon_Mapping].[szenario_afforestation].[LB_increase_afforestation_grassland_cropland_nuts3_forest_group]

  FROM [Carbon_Mapping].[szenario_afforestation].[LB_increase_afforestation_grassland_cropland_nuts3]

  left join [szenario_afforestation].[LUT_C_SEQ_AFFOR_JRC_V4_GROUPED_FGS]
  on [LUT_C_SEQ_AFFOR_JRC_V4_GROUPED_FGS].[FGS]  =[LB_increase_afforestation_grassland_cropland_nuts3].[FGS] AND
  [LUT_C_SEQ_AFFOR_JRC_V4_GROUPED_FGS].[FOREST_ZONE]  =[LB_increase_afforestation_grassland_cropland_nuts3].[FOREST_ZONE] 

  group by 
    
       [LB_increase_afforestation_grassland_cropland_nuts3].[NUTS_ID]
      ,[LB_increase_afforestation_grassland_cropland_nuts3].[NUTS_LEVEL]
      ,[LB_increase_afforestation_grassland_cropland_nuts3].[Slope_factor]
      ,[LB_increase_afforestation_grassland_cropland_nuts3].[Slope_src]

      ,[LB_increase_afforestation_grassland_cropland_nuts3].[perc_reforest]
      ,[LB_increase_afforestation_grassland_cropland_nuts3].[perc_reforest_src]
      ,[LB_increase_afforestation_grassland_cropland_nuts3].[RCP]
      ,[LB_increase_afforestation_grassland_cropland_nuts3].[Year_potential]
      ,[LB_increase_afforestation_grassland_cropland_nuts3].[FT]
      ,[LB_increase_afforestation_grassland_cropland_nuts3].[FGS]
      ,[LB_increase_afforestation_grassland_cropland_nuts3].[land_use_selection]
  
      ,[LB_increase_afforestation_grassland_cropland_nuts3].[FOREST_ZONE]
      ,[LB_increase_afforestation_grassland_cropland_nuts3].[NUTS_LEVEL0_ID]
	  
	  ,[LUT_C_SEQ_AFFOR_JRC_V4_GROUPED_FGS].[yrl_C_seq_age_0_20_avg]
      ,[LUT_C_SEQ_AFFOR_JRC_V4_GROUPED_FGS].[yrl_C_seq_age_21_30_avg]
      ,[LUT_C_SEQ_AFFOR_JRC_V4_GROUPED_FGS].[yrl_C_seq_age_0_20_min]
      ,[LUT_C_SEQ_AFFOR_JRC_V4_GROUPED_FGS].[yrl_C_seq_age_21_30_min]
      ,[LUT_C_SEQ_AFFOR_JRC_V4_GROUPED_FGS].[yrl_C_seq_age_0_20_max]
      ,[LUT_C_SEQ_AFFOR_JRC_V4_GROUPED_FGS].[yrl_C_seq_age_21_30_max]


'''



### update : change group with max area:

'''
SELECT
       [NUTS_ID]
      ,[Tree_species_factor]
      ,[Year_potential]
      ,[FT]
      ,[FGS]
      ,[land_use_selection]
      ,[FOREST_ZONE]
 
	  ,  [nr_pixels]
  FROM [Carbon_Mapping].[szenario_afforestation].[LB_increase_afforestation_grassland_cropland_nuts3]

  where [NUTS_ID] = 'AT124'  order by     [NUTS_ID]
   
      ,[Year_potential]
      ,[FT]
      ,[FGS]
      ,[land_use_selection]
      ,[FOREST_ZONE], [nr_pixels] desc

'''