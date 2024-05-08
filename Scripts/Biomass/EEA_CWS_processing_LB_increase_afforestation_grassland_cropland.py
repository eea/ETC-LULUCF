

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
LB_afforestation_results_folder = r'L:\f02_data\carbon_model_data\output\LB_NUTS_stats'





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
name_of_table = 'LB_increase_afforestation_grassland_cropland_nuts3_EL_updated'
###################################################
export_df_to_sql.to_sql(name_of_table, engine_GREENMONKEY,  schema=schmema_name,if_exists='replace')
print ("end storing on SQL")









### update tables:



################################################################################## sql to dataframe:
engine_GREENMONKEY = sa.create_engine('mssql+pyodbc://' + SERVER+ '/' + Database_name + '?trusted_connection=yes&driver=ODBC+Driver+17+for+SQL+Server')
query_1=('''
             SELECT * 
         into [Carbon_Mapping].[szenario_afforestation].[LUT_trees_test]
             FROM [Carbon_Mapping].[szenario_afforestation].[LUT_trees]
             ''')  
with engine_GREENMONKEY.begin() as conn:
    query_1 =   text(query_check) 
    conn.execute(query_1)
  

#print(readable_database)
##################################################################################
    


# (1) update attributes:

## new col:
column_name = 'model_parameter'
query_drop_column=('''      
        ALTER TABLE '''+ schmema_name+'.'+name_of_table +''' 
        DROP COLUMN if exists '''+column_name )  
print (query_drop_column)
with engine_GREENMONKEY.begin() as conn:
    query_txt =   text(query_drop_column) 
    conn.execute(query_txt)


print ("--------------------------")
query_1=('''      
        ALTER TABLE '''+ schmema_name+'.'+name_of_table +''' 
        ADD '''+column_name+''' varchar(2550);

             ''')  
print (query_1)
with engine_GREENMONKEY.begin() as conn:
    query_txt =   text(query_1) 
    conn.execute(query_txt)



## update col:

query_1=('''         UPDATE ''' + schmema_name+'.'+name_of_table +''' 
        set model_parameter =
            CONCAT(  'NUTS_ID_', [NUTS_ID],
                '_Tree_species_',[Tree_species_factor],
                '_perc_reforest_',[perc_reforest],
                '_RCP_',[RCP],
                '_Year_potential_',[Year_potential],
                '_land_use_selection_',[land_use_selection]) 
             ''')  
with engine_GREENMONKEY.begin() as conn:
    query_txt =   text(query_1) 
    conn.execute(query_txt)







## add:join forest zone to output:
name_of_table_new  =  name_of_table+'_with_forest_zone'
query_1=('''        
      Drop table if exists ''' + schmema_name+'.'+name_of_table_new +''' 
             ''')  
with engine_GREENMONKEY.begin() as conn:
    query_txt =   text(query_1) 
    conn.execute(query_txt)
# into  ''' + schmema_name+'.'+name_of_table_new +''' 
## add:join forest zone to output:
query_1=('''        
      Select 
             '''+name_of_table+'''.[index]
            ,'''+name_of_table+'''.[LB_mean_yrly_age_0_20]
            ,'''+name_of_table+'''.[LB_mean_yrly_age_21_30]
            ,'''+name_of_table+'''.[LB_min_yrly_age_0_20]
            ,'''+name_of_table+'''.[LB_min_yrly_age_21_30]
            ,'''+name_of_table+'''.[LB_max_yrly_age_0_20]
            ,'''+name_of_table+'''.[LB_max_yrly_age_21_30]
            ,'''+name_of_table+'''.[nr_pixels]
            ,'''+name_of_table+'''.[LB_total_2035_avg]
            ,'''+name_of_table+'''.[LB_total_2035_min]
            ,'''+name_of_table+'''.[LB_total_2035_max]
            ,'''+name_of_table+'''.[NUTS_ID]
            ,'''+name_of_table+'''.[NUTS_LEVEL]
            ,'''+name_of_table+'''.[Slope_factor]
            ,'''+name_of_table+'''.[Slope_src]
            ,'''+name_of_table+'''.[Tree_prob]
            ,'''+name_of_table+'''.[Tree_prob_src]
            ,'''+name_of_table+'''.[Tree_species_factor]
            ,'''+name_of_table+'''.[Tree_species_src]
            ,'''+name_of_table+'''.[perc_reforest]
            ,'''+name_of_table+'''.[perc_reforest_src]
            ,'''+name_of_table+'''.[RCP]
            ,'''+name_of_table+'''.[Year_potential]
            ,'''+name_of_table+'''.[FT]
            ,'''+name_of_table+'''.[FGS]
            ,'''+name_of_table+'''.[land_use_selection]
            ,'''+name_of_table+'''.[scenario_name]
            ,'''+name_of_table+'''.[LB_total_2050_avg]
            ,'''+name_of_table+'''.[LB_total_2050_min]
            ,'''+name_of_table+'''.[LB_total_2050_max]
            ,'''+name_of_table+'''.[model_parameter]
                
         ,LUT_FOREST_ZONE_v4.FOREST_ZONE 
         into  ''' + schmema_name+'.'+name_of_table_new +''' 
         from ''' +  schmema_name+'.'+name_of_table      +'''
      left join  [szenario_afforestation].[LUT_FOREST_ZONE_V4] 
         on [LUT_FOREST_ZONE_v4].[NUTS_LEVEL3_ID] = [LB_increase_afforestation_grassland_cropland_nuts3_EL_updated].[NUTS_ID] 
             ''')  
with engine_GREENMONKEY.begin() as conn:
    query_txt =   text(query_1) 
    conn.execute(query_txt)



#     update u
#   set u.[FGS] = s.[FGS]
#   from [Carbon_Mapping].[szenario_afforestation].[LB_increase_afforestation_grassland_cropland_nuts3_EL_updated_with_forest_zone] u
#       left join [szenario_afforestation].[LUT_C_SEQ_AFFOR_JRC_V6]  s on
#           u.[Tree_species_factor] = s.[SPECIES_NAME]







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

### get max area of afforestation for group:

'''

	 

	


	  SELECT
       [LB_increase_afforestation_grassland_cropland_nuts3].[NUTS_ID]
      ,[LB_increase_afforestation_grassland_cropland_nuts3].[Year_potential]

      ,[LB_increase_afforestation_grassland_cropland_nuts3].[FGS]
      ,[LB_increase_afforestation_grassland_cropland_nuts3].[land_use_selection]
      ,[LB_increase_afforestation_grassland_cropland_nuts3].[FOREST_ZONE]
	  ,max([LB_increase_afforestation_grassland_cropland_nuts3].[nr_pixels]) as max_nr_pixels

      ,[LUT_C_SEQ_AFFOR_JRC_V4_GROUPED_FGS].[yrl_C_seq_age_0_20_avg]
      ,[LUT_C_SEQ_AFFOR_JRC_V4_GROUPED_FGS].[yrl_C_seq_age_21_30_avg]
      ,[LUT_C_SEQ_AFFOR_JRC_V4_GROUPED_FGS].[yrl_C_seq_age_0_20_min]
      ,[LUT_C_SEQ_AFFOR_JRC_V4_GROUPED_FGS].[yrl_C_seq_age_21_30_min]
      ,[LUT_C_SEQ_AFFOR_JRC_V4_GROUPED_FGS].[yrl_C_seq_age_0_20_max]
      ,[LUT_C_SEQ_AFFOR_JRC_V4_GROUPED_FGS].[yrl_C_seq_age_21_30_max]


  FROM [Carbon_Mapping].[szenario_afforestation].[LB_increase_afforestation_grassland_cropland_nuts3]


  left join [szenario_afforestation].[LUT_C_SEQ_AFFOR_JRC_V4_GROUPED_FGS] on 

[LUT_C_SEQ_AFFOR_JRC_V4_GROUPED_FGS].[FGS]         =[LB_increase_afforestation_grassland_cropland_nuts3].[FGS] AND
[LUT_C_SEQ_AFFOR_JRC_V4_GROUPED_FGS].[FOREST_ZONE]= [LB_increase_afforestation_grassland_cropland_nuts3].[FOREST_ZONE]

  ----where [NUTS_ID] = 'AT124' ----and [LB_increase_afforestation_grassland_cropland_nuts3].[FGS] = 'F'

group by 
	   [LB_increase_afforestation_grassland_cropland_nuts3].[NUTS_ID]
      ,[LB_increase_afforestation_grassland_cropland_nuts3].[Year_potential]
      ,[LB_increase_afforestation_grassland_cropland_nuts3].[FGS]
      ,[LB_increase_afforestation_grassland_cropland_nuts3].[land_use_selection]
      ,[LB_increase_afforestation_grassland_cropland_nuts3].[FOREST_ZONE]
	  ,[LUT_C_SEQ_AFFOR_JRC_V4_GROUPED_FGS].[yrl_C_seq_age_0_20_avg]
      ,[LUT_C_SEQ_AFFOR_JRC_V4_GROUPED_FGS].[yrl_C_seq_age_21_30_avg]
      ,[LUT_C_SEQ_AFFOR_JRC_V4_GROUPED_FGS].[yrl_C_seq_age_0_20_min]
      ,[LUT_C_SEQ_AFFOR_JRC_V4_GROUPED_FGS].[yrl_C_seq_age_21_30_min]
      ,[LUT_C_SEQ_AFFOR_JRC_V4_GROUPED_FGS].[yrl_C_seq_age_0_20_max]
      ,[LUT_C_SEQ_AFFOR_JRC_V4_GROUPED_FGS].[yrl_C_seq_age_21_30_max]
	   order by     [NUTS_ID]
'''



print(
'''
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡠⣄⡀⠀⠀⡠⠞⠛⢦⣠⢤⡀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡀⢠⠏⠀⠀⢱⡀⣸⠁⠀⡴⠋⠀⠀⣹⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡴⠋⠉⢿⢀⡤⠶⣴⠇⣯⠀⣼⠁⠀⢀⡴⠷⣄
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⠞⠁⠀⣀⡾⠋⠀⠀⢹⣼⠁⢠⡇⠀⡴⠋⠀⠀⡼
⠀⠀⠀⠀⢠⠊⠑⢦⠀⡴⠋⢀⣠⠞⠉⠀⠀⠀⣠⣿⠧⣄⡾⠁⡼⠁⣀⣤⠾⡁
⠀⠀⠀⠀⢸⠀⠀⣨⠟⠁⢠⡞⠁⠀⠀⠀⣠⡾⠛⠁⠀⣿⠃⣰⠃⣴⠋⠀⠀⣷
⠀⠀⠀⠀⣸⢠⠞⠁⠀⢠⠏⠀⠀⢀⡴⠋⠁⠀⢀⣠⡴⠿⣶⡇⢰⠇⠀⠀⢠⠇
⠀⠀⠀⢠⢿⠏⠀⠀⠀⠉⠀⠀⣠⠞⠁⠀⡴⠚⠉⠁⠀⢀⡟⠀⣼⠀⠀⠀⢸⠀
⠀⠀⠀⡾⣼⢀⠀⠀⠀⠀⠀⠈⠉⠀⣠⠞⠁⠀⠀⢀⡴⠋⠙⢼⠃⠀⠀⠀⣸⠀
⠀⠀⠀⡇⠉⡎⠀⣰⠃⠀⠀⠀⠀⠀⠁⠀⠀⠀⡼⠉⠀⠀⠀⠘⠂⠀⠀⣠⠇⠀
⠀⠀⠀⡇⢸⠀⣰⠃⠀⡴⠀⠀⠀⠀⠀⠀⣠⠞⠁⠀⠀⠀⠀⠀⠀⣠⠖⠁⠀⠀
⠀⠀⢸⠁⡏⢠⠃⢀⠞⠀⠀⠀⠀⠀⠀⢸⠁⠀⠀⠀⠀⢀⣠⠖⠋⠁⠀⠀⠀⠀
⠀⠀⡞⠀⠃⡎⢀⠏⠀⠀⠀⠀⠀⠀⢀⡏⠀⣀⡤⠴⠚⠉⠀⠀⠀⠀⠀⠀⠀⠀
⡴⢺⠇⠀⠀⠀⠞⠀⠀⠀⠀⠀⠀⢀⡾⠒⠋⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⡇⠘⣆⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⠞⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⢳⡀⠘⢦⡀⠀⠀⠀⠀⠀⠀⡰⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀

⠀⠳⣄⠀⠙⠲⣤⣀⣠⠴⠊⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠈⠓⠦⣄⣀⡠⠎⠀
''')

print ("job done")