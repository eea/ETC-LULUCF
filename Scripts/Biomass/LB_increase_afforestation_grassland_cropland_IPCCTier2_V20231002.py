
# # IPCC Tier 2 method LB (Living Biomass) increase simulation for afferestation on cropland and grassland

# #### Load and (install) used packages

import os
import glob
from pathlib import Path
import geopandas as gpd
from loguru import logger
from pathlib import Path, WindowsPath
#import geopandas as gpd
import sys
import matplotlib.pyplot as plt
import logging
import rasterio
import rasterio.mask
import numpy as np
from osgeo import gdal
import subprocess
from osgeo import osr
import pandas as pd
import pprint
import datetime
import pyodbc 
import sqlalchemy as sa
from sqlalchemy import create_engine, event
from sqlalchemy.engine.url import URL
import json
import dask_geopandas
from dask.distributed import LocalCluster, Client, performance_report

### location of the code to run the scenarios
sys.path.append(Path(os.getcwd()).joinpath('src').as_posix())
from SOC_scenarios.utils.soc_helper_functions import *
from constants import (
    type_method,
    dict_C_pool_mng_options,
    dir_signature
)

from Biomass.utils.biom_helper_functions import *
sys.path.append(Path(os.getcwd()).joinpath('Scripts').as_posix())
from SOC.SOC_stratification import SOC_strat_IPCC_block_proc
from Biomass.run_afforestation_scenario import afforestation_LUT_block_proc

@logger.catch()
def get_settings():
    logger.info('Get settings')
    # Year_potential = 2035
    Year_baseline = 2024

    # define which management action is applied
    mng_option = 'afforestation'

    # load the carbon pool that is targeted here
    carbon_pool = dict_C_pool_mng_options.get(mng_option)

    # List of CLC classes (LEVEL-3) that can be used for afforestation
    lst_CLC_afforestation = [211, 212, 213, 231] #TODO land use selection cropland or  grassland
    # lst_CLC_afforestation_cropland = [211, 212, 213] #TODO land use selection cropland or  grassland
    # lst_CLC_afforestation_grassland = [231]

    dict_default_afforestation_factors = {
            'Slope': 0,
            'RCP': 'rcp45',
            'Tree_prob': 70,
            'Tree_species': 'Betula_pendula', #TODO 36 tree species
            'Perc_reforest': 10, # TODO afforestation intensity
            # 'Year_potential': Year_potential,
            'input_source': 'EEA39'}

    SCENARIO_SPECS = {
            'mng_option': mng_option, 
            # 'Year_potential': Year_potential,
            'Year_baseline': Year_baseline,
            'carbon_pool': carbon_pool, # Living biomass pool
            'lst_CLC_affor': lst_CLC_afforestation,
            'afforestation_config': dict_default_afforestation_factors
        }


    # ## PART2: CONFIGURE PROCESSING SETTINGS

    # Define if the results should be created by block-based
    # processing, meaning only the window for the AOI (i.e., NUTS region) will
    # be loaded

    # Only this type of processing is supported for the moment!!!
    block_based_processing = True

    # Country_running
    Country = None  # set to None if want to run entire EEA39 extent 

    # suffix of the output raster and ID that is
    # used to distinguish different scenario runs
    # TODO fix
    scenario_name = f'Scenario_JRCV3_{str(2035)}_fix'

    # if want to run with some NUTS specific factors
    # set the below parameter to true
    # if a country is defined only the NUTS regions
    # in that country will be considered
    run_NUTS_specific_scenario = False

    # Indicate if stats at NUTS LEVEL need to be provided
    # this will be witten in shapefile format and saved
    add_stats_NUTS_level = True


    # define if the result might be overwritten
    overwrite = True

    # the location of the basefolder in which all the input data needed for the simulation is provided
    Basefolder_input_data = os.path.join(dir_signature, 'input')

    ### the basefolder in which all the output results will be stored
    Basefolder_output =  os.path.join(dir_signature, 'output')

    # the name of the drive on which all the data can be accessed:
    # dir_signature = os.getcwd()

    # location where the NUTS specific configuration for afforesation is stored
    NUTS_LUT_factors_folder = os.path.join(Basefolder_output,'NUTS_LUT_afforestation_scenario')

    # define the yield table LUT and forest zone LUT
    # location that should be used for processing
    name_yield_table_LUT = 'LUT_C_SEQ_AFFOR_JRC_V3.csv'
    name_LUT_forest_zones = 'LUT_FOREST_ZONE.csv'
    folder_JRC_table = os.path.join(Basefolder_output, 'NUTS_LUT_afforestation_scenario',
                                    'JRC_yield_table')
    yield_table_LUT_dir = os.path.join(folder_JRC_table, name_yield_table_LUT)
    forest_zone_dir = os.path.join(folder_JRC_table, name_LUT_forest_zones)


    CONFIGURATION_SPECS = {
        'block_based_processing': block_based_processing,
        'Country': Country,
        'scenario_name': scenario_name,
        'yield_table_dir': yield_table_LUT_dir,
        'forest_zone_dir': forest_zone_dir,
        'run_NUTS_SPECIFIC': run_NUTS_specific_scenario,
        'add_stats_NUTS': add_stats_NUTS_level,
        'scaling': 100,
        'NUTS_LUT_folder': NUTS_LUT_factors_folder,
        'overwrite': overwrite,
        'dir_signature': dir_signature,
        'type_method': type_method,
        'Basefolder_output': Basefolder_output,
        'Basefolder_input': Basefolder_input_data
    }


    # ## PART3: DEFINE REQUIRED DATASETS

    # Define processing extent
    ## location with the EEA39 extent
    # EEA_extent_layer = gpd.read_file(os.path.join(Basefolder_input_data, 'AOI', 'EEA39_extent_noDOM.shp'))
    VECTORS = {'NUTS':os.path.join(Basefolder_input_data, 'NUTS', 'NUTS_RG_20M_2021_3035.shp')}
            # 'EEA_extent': EEA_extent_layer}

    # Below define the locations of the rasters needed for the afforestation mask

    DEM_location = os.path.join(Basefolder_input_data, 'DEM', 'dt1_slope_INCA_100m_EPSG3035.tif')
    N2K_location = os.path.join(Basefolder_input_data, 'N2K', 'N2K_2022_100m.tif')      ###s4e: we added a new N2K dataset with 100m resolution
    # if a region specific external mask is available,
    # please insert it below otherwise set it to 'None'
    # !! Please ensure that all the rasters are on the same grid
    external_mask = None


    ########## TODO: ETC/DI #############

    """
    Following additional dtaatsets should be implemented to identify those areas 
    that should be protected for afforestation:
    * Nationally designated areas (CDDA): Protected areas should be EXCLUDED for afforestation
    * High nature value farmland 2012: High nature value farmland should be EXCLUDED for afforestation
    * Extended wetland layer (2018): Any wetland should be EXCLUDED for afforestation
    * Peatland map (2017): Tanneberger peatland map --> peatland should be EXCLUDED for afforestation


    The abandoned cropland locations should be PRIORITIZED for afforestation if is has been found suitable
    for afforestation based on the afforestation mask
    """

    CLC_dir = os.path.join(Basefolder_input_data, 'CLC_ACC', 'CLC2018ACC_V2018_20.tif')

    DATASETS_afforestation = {'DEM': DEM_location,
                            'N2K': N2K_location,
                            'external_mask': external_mask,
                            'CLC': CLC_dir}

    """Load the datasets needed for defining afforestation suitability of tree species"""

    # Below define the raster(s) for defining the suitability
    #  of planting a certain tree species

    ## Load the basefolder of the EU4Trees dataset
    EU4Trees_base = os.path.join(Basefolder_input_data, 'EU-Trees4F', 'ens_sdms','prob','resampled')

    DATASETS_AFFORESTATION_SUITABILTY = {
        'EU4Trees': EU4Trees_base
    }

    """Load IPCC datasets needed to assign strata specific increment values"""

    # Put all the IPCC related datasets used for strata based assigned using the LUT

    ## the directory of the IPCC climate  
    IPCC_climate_raster = os.path.join(Basefolder_input_data, 'IPCC_layers', 'climate', 'ipcc_climate_zone_100m_EPSG3035_EEA39.tif') 

    DATASETS_IPCC = {
        'CLIMATE_ZONES': IPCC_climate_raster,
    }

    DATASETS_SPECS = {
        'VECTORS': VECTORS,
        'AFFORESTATION_MASK': DATASETS_afforestation,
        'AFFORESTATION_SUITABILITY': DATASETS_AFFORESTATION_SUITABILTY,
        'IPCC_STRATIFICATION': DATASETS_IPCC
    }

    # ### TODO Specify area of interest
    # (running on the entire EEA39 extent is not recommended in a Notebook environment due to memory limitations!) 
    # NUTS based processing of entire extent is possible, but in that case the **'run_NUTS_specific_scenario'** parameter should be set to True.
    # 
    # A NUTS3 region should be selected for processing (e.g. 'PL432').The original script located here: https://github.com/VITObelgium/ETC-CCA-LULUCF/blob/master/Scripts/Biomass/run_afforestation_scenario.py  will allow to run all different NUTS regions at country or EU level at once. This notebook examplifies just the workflow for one NUTS region. 

    #select now the country of interest; if set to None the entire EEA39 extent will be processed
    NUTS_region = None #'CZ020' #LU000' # 'CZ020'

    # ASSING ALL THE ABOVE DEFINED CONDIGURATION INTO THE SETTINGS DICTIONARY WHICH WILL BE USED IN THE AFFORESTATION FUNCTIONS

    settings = {'SCENARIO_SPECS': SCENARIO_SPECS,
                'CONFIG_SPECS': CONFIGURATION_SPECS,
                'DATASETS': DATASETS_SPECS,
                'commit_id': 'cb72cc9159d438f3e235e3e9c10e60dd2853f95c'}

                # 'NUTS_region': NUTS_region,
                # 'NUTS3_info': NUTS3_region}
    return settings
    

def nuts_wrapper(row, settings):
    logger.info(f'{row.NUTS_ID} Define settings for NUTS' )
    # define settings for this nuts
    settings_nuts = settings.copy()
    # NUTS_layer = gpd.read_file(settings.get('DATASETS').get('VECTORS').get('NUTS'))
    # NUTS3_region = NUTS_layer.loc[((NUTS_layer.LEVL_CODE == 3) & (NUTS_layer.NUTS_ID == row.NUTS_ID))].iloc[0,:]
    settings_nuts['NUTS_region'] = row.NUTS_ID
    settings_nuts['NUTS3_info'] = row

    outfolder = Path(settings.get('CONFIG_SPECS').get('Basefolder_output')).joinpath(
                    f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}_NUTS_stats')
    outname = f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}_stats_NUTS_EEA39_{settings.get("CONFIG_SPECS").get("scenario_name")}_{row.NUTS_ID}.csv'
    if(os.path.exists(Path(outfolder).joinpath(outname))):
        logger.info(f'{row.NUTS_ID} Statistics already exist')
        return pd.DataFrame(data=[])

    # from https://github.com/VITObelgium/ETC-CCA-LULUCF/blob/master/notebooks/output/NUTS_LUT_afforestation_scenario/JRC_yield_table/LUT_C_SEQ_AFFOR_JRC_V3.csv
    tree_species = ['Abies_alba', 'Acer_campestre',
    'Alnus_glutinosa', 'Alnus_incana',
       'Acer_opalus', 'Acer_platanoides', 'Acer_pseudoplatanus',
       'Betula_pendula', 'Carpinus_betulus', 'Carpinus_orientalis',
       'Castanea_sativa', 'Fraxinus_angustifolia', 'Fraxinus_excelsior',
       'Fraxinus_ornus', 'Fagus_sylvatica', 'Juglans_regia',
       'Larix_decidua', 'Picea_abies', 'Pinus_cembra', 'Pinus_pinaster',
       'Pinus_nigra', 'Populus_nigra', 'Pinus_sylvestris', 'Pinus_pinea',
       'Prunus_avium', 'Quercus_cerris', 'Quercus_faginea',
       'Quercus_frainetto', 'Quercus_ilex', 'Quercus_robur',
       'Quercus_suber', 'Quercus_petraea', 'Quercus_pubescens',
       'Quercus_coccifera', 'Salix_alba', 'Ulmus_minor'] 
    # # TODO define settings variants
    lst_NUTS_stats = []
    for classes, land_use_selection in [[[211, 212, 213],'cropland'], [[231], 'grassland']]:  # cropland or grassland
        settings_nuts['SCENARIO_SPECS']['lst_CLC_affor'] = classes
        settings_nuts['SCENARIO_SPECS']['lst_CLC_affor_name'] = land_use_selection
        # afforestation mask depends only on land use
        logger.info(f'{row.NUTS_ID} Compute afforestation mask for land use selection: {land_use_selection}')
        affor_mask_layer = define_affor_areas(settings_nuts)

        # for intensity in [10, 20, 50]:
        # settings_nuts['SCENARIO_SPECS']['afforestation_config']['Perc_reforest'] = intensity
        # logger.info(f'{row.NUTS_ID} Compute scenario with % reforestation = {intensity}')
        for tree in tree_species:
            logger.info(f'{row.NUTS_ID}  Compute scenario with tree species = {tree}')
            settings_nuts['SCENARIO_SPECS']['afforestation_config']['Tree_species'] = tree
            
            logger.info(f'{row.NUTS_ID} Compute afforestation potential layer')
            affor_potential_layer, outname_affor_pot = create_affor_potential(settings_nuts, affor_mask_layer)

            if(affor_potential_layer is None):
                logger.warning(f'{row.NUTS_ID} tree species: {tree} land use: {land_use_selection} skipped')
                
            elif(outname_affor_pot == -1): 
                logger.warning(f'Skipping region {row.NUTS_ID}')
                pd.DataFrame(data=[]).to_csv(Path(outfolder).joinpath(outname))
                return pd.DataFrame(data=[])
            else:
                logger.info(f'{row.NUTS_ID} Compute stats')
                for intensity in [10, 20, 50]:
                    settings_nuts['SCENARIO_SPECS']['afforestation_config']['Perc_reforest'] = intensity
                    logger.info(f'{row.NUTS_ID} Compute scenario with % reforestation = {intensity}')
                    for year_potential in [2035, 2050]:
                        logger.info(f'{row.NUTS_ID} Compute scenario with year_potential = {year_potential}')
                        settings_nuts["SCENARIO_SPECS"]["Year_potential"] = year_potential
                        df_stats_NUTS = calc_stats_biomass_NUTS(os.path.join(
                            settings_nuts.get('CONFIG_SPECS').get('Basefolder_output'),
                            'LB_scenario', 
                            settings_nuts.get('NUTS3_info')['CNTR_CODE'], 
                            outname_affor_pot),
                            settings_nuts.get('NUTS3_info'),
                            settings_nuts)
                        df_stats_NUTS['land_use_selection'] = land_use_selection
                        df_stats_NUTS['scenario_name'] = f"{intensity}_{tree}_{land_use_selection}"
                        lst_NUTS_stats.append(df_stats_NUTS)

    # merge all the statistics of each NUTS-3 region
    # in a single summary dataframe
    if(len(lst_NUTS_stats)>0):
        df_stats_all_NUTS3 = pd.concat(lst_NUTS_stats)

        # # now that the NUTS3 values are calculated also the
        # # weighted average at NUTS LEVEL0 will be derived
        # df_stats_NUTS_final = calc_weighted_average_NUTS(
        #     df_stats_all_NUTS3, gpd.read_file(settings.get('DATASETS').get('VECTORS').get('NUTS')))
        outfolder = Path(settings.get('CONFIG_SPECS').get('Basefolder_output')).joinpath(
                    f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}_NUTS_stats')
        outname = f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}_stats_NUTS_EEA39_{settings.get("CONFIG_SPECS").get("scenario_name")}_{row.NUTS_ID}.csv'
        outfolder.mkdir(parents=True, exist_ok=True)
        df_stats_all_NUTS3.to_csv(Path(outfolder).joinpath(outname))
        logger.info(f'Write results in {Path(outfolder).joinpath(outname)}')
        # return df_stats_NUTS
        return df_stats_all_NUTS3
    else:
        pd.DataFrame(data=[]).to_csv(Path(outfolder).joinpath(outname))
        return pd.DataFrame(data=[])

if __name__ == '__main__':
    # logger.add('LB_increase_afforestation_grassland_cropland_IPCCTier2_V20231002_dask.log')
    logger.info("Get general settings")
    settings = get_settings()


    # # get NUTS list
    logger.info("Get NUTS")
    shp_NUTS = gpd.read_file(settings.get(
        'DATASETS').get('VECTORS').get('NUTS'))

    # filter only on the NUTS levels of interest
    shp_NUTS = shp_NUTS.loc[shp_NUTS.LEVL_CODE == 3] #only select NUTS3
    # TODO filter NUTS with FOREST ZONE LUT 
    df_LUT_forest_zones = pd.read_csv(os.path.join(settings.get('CONFIG_SPECS').get('forest_zone_dir')))
    shp_NUTS_filtered = shp_NUTS.loc[shp_NUTS.NUTS_ID.isin(df_LUT_forest_zones.NUTS_LEVEL3_ID) ]
    ddf = dask_geopandas.from_geopandas(shp_NUTS_filtered, npartitions=8)
    # ddf = ddf.spatial_shuffle()
 
    client=Client()
    print(client.dashboard_link)
    with performance_report(filename="afforestation_dask-report.html"):
        stats = ddf.apply(lambda row: nuts_wrapper(row, settings), axis=1, meta=('x', 'object')).compute()
    # stats = shp_NUTS_filtered[:5].apply(lambda row: nuts_wrapper(row, settings), axis=1)
    client.shutdown()
