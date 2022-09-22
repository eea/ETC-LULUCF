"""In this code a stratification on the EEA39 extent will be applied
in such a way that the IPCC default values can be used for making some projections
of the SOC change for different scenarios. """
import os
import glob
from pathlib import Path
import geopandas as gpd
import pandas as pd
from loguru import logger
import numpy as np

from SOC_scenarios.utils.soc_helper_functions import mask_raster_extent\
    , create_SOC_REF_raster, create_FLU_layer, create_factor_layer\
    ,create_SOC_scenario_layer, resample_raster_to_ref, mosaic_raster\
    , calc_stats_SOC_NUTS, calc_weighted_average_NUTS, \
    create_metadata_description_SOC, add_metadata_raster


def SOC_strat_IPCC_block_proc(settings: dict):
    """
    Function that will calculate the SOC scenario by block processing the process per NUTS LEVEL. This function should
    be enabled if want to run the SOC with scenarios at NUTS level.
    :param settings: dictionary with the settings needed to run the SOC scenario
    :return: Output map with the requested SOC scenario
    """
    shp_NUTS = settings.get('NUTS_extent_map')
    ### filter only on the NUTS levels of interest
    shp_NUTS = shp_NUTS.loc[shp_NUTS.LEVL_CODE.isin([0,3])]

    ### reproject to LAEA of not yet done
    crs = shp_NUTS.geometry.crs
    if crs is None:
        shp_NUTS.crs = shp_NUTS.crs(epsg=3035)
    elif crs.srs[-4:] != '3035':
        shp_NUTS = shp_NUTS.to_crs(epsg=3035)

    if settings.get('Country') is not None:
        shp_NUTS = shp_NUTS.loc[((shp_NUTS.CNTR_CODE == settings.get('Country')[0:2])& (shp_NUTS.LEVL_CODE == 3))]
    else:
        shp_NUTS = shp_NUTS.loc[shp_NUTS.LEVL_CODE == 3]

    ### LOOP NOW over the available NUTS3 regions in the area that is specified

    #TODO add maybe check if the EEA39 data on FLU, SOCREF, IPCC SOIL & CLIMATE is available
    shp_NUTS = shp_NUTS.reset_index(drop = True)
    lst_raster_files = []
    lst_NUTS_stats = []
    for index, NUTS3_region in shp_NUTS.iterrows():
        logger.info(f'CALCULATING SOC FOR NUTS3 REGION: {NUTS3_region.NUTS_ID} \n'
                    f'{np.round((index/shp_NUTS.shape[0]) *100,2)}% COMPLETED')
        settings.update({'NUTS3_info': NUTS3_region})


        ### define the output name of the SOC scenarios per NUTS LEVEL
        outfolder = Path(settings.get('Basefolder_output')).joinpath('SOC_scenario')\
            .joinpath(settings.get('NUTS3_info')['CNTR_CODE']).as_posix()
        outname_SOC_scenario = f'SOC_{settings.get("Scenario_name")}' \
                               f'_{settings.get("NUTS3_info")["NUTS_ID"]}.tif'
        settings.update({'outdir_SOC_scenario': os.path.join(outfolder, outname_SOC_scenario)})

        lst_raster_files.append(os.path.join(outfolder, outname_SOC_scenario))

        if not os.path.exists(os.path.join(outfolder, outname_SOC_scenario)) or settings.get('overwrite'):

            if settings.get('Fixed_factor_FLU'):
                FLU_layer, meta_FLU_layer = create_FLU_layer(settings, fixed_factor_creation=True)
                settings.update({'FLU_layer': FLU_layer,
                                 'FLU_meta_layer': meta_FLU_layer})
                if FLU_layer is None and meta_FLU_layer is None:
                    continue

            ### create stock change factors specified for NUTS region
            FMG_layer, meta_FMG_layer = create_factor_layer(settings, type_factor='FMG')
            if FMG_layer is None and meta_FMG_layer is None:
                continue
            FI_layer, meta_FI_layer = create_factor_layer(settings, type_factor='FI',
                                           fixed_factor_creation=True)
            if FI_layer is None and meta_FI_layer is None:
                continue

            ### add this loaded layer to the settings_file
            settings.update({'FMG_layer': FMG_layer,
                             'FMG_meta_layer': meta_FMG_layer,
                             'FI_layer': FI_layer,
                             'FI_meta_layer': meta_FI_layer})
            ### NOW THE FINAL SOC BASED ON THE DEFINED SCENARIO CAN BE CALCULATED

            create_SOC_scenario_layer(settings)


            #### ADD SOME METADATA FOR SOC SCENARIO
            dict_general, dict_band = create_metadata_description_SOC(settings, extent='NUTS')
            add_metadata_raster(dict_general, dict_band, os.path.join(outfolder, outname_SOC_scenario))




        ### check if some stats of the map need to be provided
        if settings.get('add_stats_NUTS_level') is not None:

            df_stats_NUTS = calc_stats_SOC_NUTS(os.path.join(outfolder, outname_SOC_scenario)
                                                , NUTS3_region,settings)
            lst_NUTS_stats.append(df_stats_NUTS)


    if settings.get('add_stats_NUTS_level') is not None:

        ### define output folder for writing the result
        if settings.get('Country') is not None:
            outfolder = Path(settings.get('Basefolder_output')).joinpath('SOC_NUTS_stats').joinpath(settings.get('Country'))
            outname = f'SOC_stats_NUTS_{settings.get("Country")}_{settings.get("Scenario_name")}.shp'

        else:
            outfolder = Path(settings.get('Basefolder_output')).joinpath('SOC_NUTS_stats')
            outname = f'SOC_stats_NUTS_EEA39_{settings.get("Scenario_name")}.shp'

        outfolder.mkdir(parents=True, exist_ok=True)

        if not os.path.exists(Path(outfolder).joinpath(outname)) or settings.get('overwrite'):
            df_stats_all_NUTS3 = pd.concat(lst_NUTS_stats)

            ## now that the NUTS3 is calculated also the weighted average at NUTS LEVEL0 can be derived
            df_stats_NUTS_final = calc_weighted_average_NUTS(df_stats_all_NUTS3, settings.get('NUTS_extent_map'))


            ## create now geodataframe out of it and write out the result
            gpd_NUTS_stats = gpd.GeoDataFrame(df_stats_NUTS_final, geometry=df_stats_NUTS_final.geometry)

            ## write out the result
            gpd_NUTS_stats.crs = shp_NUTS.crs
            gpd_NUTS_stats.to_file(Path(outfolder).joinpath(outname))



    ### Now mosaic all the pieces together
    dir_mosaic_raster = mosaic_raster(lst_raster_files, settings, return_outdir=True)


    dict_general, dict_band = create_metadata_description_SOC(settings)
    add_metadata_raster(dict_general, dict_band, dir_mosaic_raster)



def SOC_strat_IPCC_full_extent(settings):
    """"
    Function that will calculate the SOC for a specific scenario for the full EEA39 extent
    So even if country scale processing is enabled first the full extent will be opened to calculate
    country specific SOC.
    :param settings: dictionary with the settings needed to run the SOC scenario
    :return: Output map with the requested SOC scenario
    """

    ### Load the needed variables from the settings file
    Basefolder_input_data = settings.get('Basefolder_input_data')
    dir_signature = settings.get('dir_signature')
    overwrite = settings.get('overwrite')
    shp_extent_map = settings.get('EEA_extent_map')



    #### SOIL, CLIMATE layer should be uploaded

    IPCC_folder = os.path.join(Basefolder_input_data,'IPCC', 'IPCC_data_converted')


    ### search for soil data
    IPCC_files = glob.glob(os.path.join(IPCC_folder, '*.tif'))
    soil_file = [item for item in IPCC_files if 'soil' in Path(item).stem]

    ### search for climate data
    climate_file = [item for item in IPCC_files if 'climate' in Path(item).stem]

    ### apply a resampling and clipping to the data

    ### first do the clipping to limit data redundancy
    if settings.get('Country') is None:
        mask_raster_extent(soil_file,shp_extent_map, Path(IPCC_folder).joinpath('soil_final_map')
                           , settings,overwrite=False)
        mask_raster_extent(climate_file, shp_extent_map, Path(IPCC_folder).joinpath('climate_final_map'),settings,
                           overwrite=False)
        Country = 'EEA39'
    else:
        mask_raster_extent(soil_file,settings.get('NUTS_extent_map'), Path(IPCC_folder).joinpath('soil_final_map')
                           , settings,overwrite=False, Country_clipping=True)
        mask_raster_extent(climate_file, settings.get('NUTS_extent_map'), Path(IPCC_folder).joinpath('climate_final_map'),settings,
                           overwrite=False, Country_clipping=True)
        Country = settings.get('Country')

    ## now apply the resampling of the files to 100m
    if Country != 'EEA39':
        soil_map_to_resample = glob.glob(os.path.join(IPCC_folder, 'soil_final_map',Country, '*clipped.tif'))
        climate_map_to_resample = glob.glob(os.path.join(IPCC_folder, 'climate_final_map',Country, '*clipped.tif'))
    else:
        soil_map_to_resample = glob.glob(os.path.join(IPCC_folder, 'soil_final_map', '*clipped.tif'))
        climate_map_to_resample = glob.glob(os.path.join(IPCC_folder, 'climate_final_map', '*clipped.tif'))



    ### IPCC soil resampling

    ### need a reference file with the proper resolution and extent
    CLC_ref_file = os.path.join(dir_signature,'input_data','general','CLCACC', 'CLC2018ACC_V2018_20.tif')
    if Country != 'EEA39': ### convert it to the proper extent
        mask_raster_extent([CLC_ref_file], settings.get('NUTS_extent_map')
                           , Path(settings.get('Basefolder_input_data')).joinpath('CLC_ACC'),settings,
                           overwrite=False, Country_clipping=True)
        CLC_ref_file = Path(settings.get('Basefolder_input_data')).joinpath('CLC_ACC')\
                        .joinpath(Country).joinpath('CLC2018ACC_V2018_20_clipped.tif').as_posix()



    outfile_name_soil_map_resampled = settings.get('path_IPCC_soil_resampled')
    resample_raster_to_ref(soil_map_to_resample,CLC_ref_file, Country,outfile_name_soil_map_resampled.parent.as_posix(),
                           resample_factor=1, overwrite=overwrite,
                           resampling=True,
                           outname=outfile_name_soil_map_resampled.name)

    ### IPCC climate resampling
    outfile_name_climate_resampled = settings.get('path_IPCC_climate_resampled')
    resample_raster_to_ref(climate_map_to_resample,CLC_ref_file, Country,Path(outfile_name_climate_resampled).parent.as_posix(),
                           resample_factor=1, overwrite=overwrite,
                           resampling=True,
                           outname=outfile_name_climate_resampled.name)

    ### Now the SOC REF will be created based on IPCC soil and climate data
    create_SOC_REF_raster(settings)


    ### CREATING OF STOCK CHANGE FACTORS MAPS RELATED TO
    # MANAGEMENT AND LAND USE (FLU, FI, FMG)

    ### FLU generation
    create_FLU_layer(settings)

    ###### FMG LAYER GENERATION
    create_factor_layer(settings, type_factor='FMG'
                        ,fixed_factor_creation=settings.get('Fixed_factor_FMG'))

    #### FI LAYER GENERATION
    create_factor_layer(settings, type_factor='FI'
                        ,fixed_factor_creation=settings.get('Fixed_factor_FI'))

    ### NOW THE FINAL SOC BASED ON THE DEFINED SCENARIO CAN BE CALCULATED
    create_SOC_scenario_layer(settings)







def main_stratification(settings):

    #### In the first place the needed datasets will be
    # uploaded and will be checked if they are all in the same format

    SOC_method = settings.get('SOC_method')
    block_based_processing = settings.get('block_based_processing')

    if SOC_method == 'IPCC' and not block_based_processing:
        SOC_strat_IPCC_full_extent(settings)
    elif SOC_method == 'IPCC' and block_based_processing:
        SOC_strat_IPCC_block_proc(settings)


if __name__ == '__main__':

    logger.info('*' * 50)

    ### Some default settings

    dir_signature = 'L:'
    overwrite = False
    Basefolder_strata = os.path.join(dir_signature, 'etc', 'lulucf', 'strata')
    SOC_LUT_folder = os.path.join(Basefolder_strata, 'SOC_LUT')


    Basefolder_input_data = os.path.join(dir_signature,'etc','lulucf','input_data')
    CLC_ACC_folder = os.path.join(dir_signature, 'input_data', 'general', 'CLCACC')
    SOC_method = 'IPCC' ### if the scenarios should be made based on IPCC defaul values

    shp_extent_map = gpd.read_file(os.path.join(dir_signature, 'etc','lulucf','AOI','EEA39_extent_noDOM.shp'))


    path_IPCC_climate_resampled = os.path.join(Basefolder_input_data,'IPCC', 'IPCC_data_converted', 'climate_final_map',
                                                                  'ipcc_climate_zone_100m_EPSG3035_EEA39.tif')
    path_IPCC_soil_resampled = os.path.join(Basefolder_input_data,'IPCC', 'IPCC_data_converted', 'soil_final_map',
                                                              'ipcc_soil_type_100m_EPSG3035_EEA39.tif')


    #LUT CLC IPCC MAPPING

    CLC_IPCC_mapping = {
        'Forest': {'CLC_classes': [311,312,313,324,334],
                   'ID_class': 1},
        'Cropland': {'CLC_classes': [211,212,213,221,222,223,
                                    241,242,243,244],
                     'ID_class': 2},
        'Grassland': {'CLC_classes': [231,321,322,323],
                      'ID_class': 3},
        'Wetlands': {'CLC_classes': [411,412,421,422,423
            ,511,512,521,522,523],
                     'ID_class':4},
        'Settlements': {'CLC_classes': [111,112,121,
                                        122,123,124,
                                        131,132,133,
                                        141,142],
                        'ID_class': 5},
        'Other_land': {'CLC_classes': [331,332,333,335],
                       'ID_class': 6}}

    ###### SPECIFY SOME SCENARIO (FLU & FMG & FI) WITH DIFFERENT STOCK CHANGE
    # FACTORS FROM IPCC TO MAP THE SOC RESULT

    ##### FLU (Land Use) Cropland options:
    """
    1: Longterm cultivated
    2: Paddy rice
    3: Perrennial/Tree crop
    4: Set aside (<20 yrs)
    """

    ##### FMG (Management) Cropland options:
    """
    1: Full tillage
    2: Reduced
    3: No-till
    """
    ##### FI (Input) Cropland options:
    """
    1: Low
    2: Medium
    3: High without manure
    4: High with manure
    """

    ##### FMG Grassland options:
    """
    1: Nominally managed 
    2: Moderately degraded
    3: Severely degraded
    4: Improved grassland
    """

    ##### FI (Input) Grassland options:
    """
    1: Nominal
    2: High
    """

    #### Now define which factors you want to apply for SOC scenario

    ## One of the LU categories can be removed from the dictionary
    # if only want to do the estimation for one of the classes

    #### use some default factors for EU LEVEL if no NUTS specific factors are provided
    dict_default_stock_change_factors = {
        'Cropland': {'FLU': 1, 'FMG': 1, 'FI': 1}
        }
    #'Grassland': {'FMG': 1, 'FI': 1}


    ### Define if the results should be created by block-based
    ### processing, meaning only the window for the AOI will
    ### be loaded
    block_based_processing = True

    ### extension that is used to distinguish different scenario runs
    scenario_name = 'Scenario_1_baseline'


    ### Country_running
    Country = None #set to None if want to run entire EEA39 extent
    shp_NUTS_borders = gpd.read_file(os.path.join(dir_signature, 'etc','lulucf','AOI','NUTS_RG_20M_2021_3035.shp'))



    ### if want to run with some NUTS specific factors
    ### set the below parameter to true
    ### if a country is defined only the NUTS regions in that country will be considered
    run_NUTS_specific_scenario = False
    SOC_NUTS_factors_folder = os.path.join(Basefolder_strata,'NUTS_LUT_SOC_scenario')


    #### Indicate of stats at NUTS LEVEL need to be provided
    add_stats_NUTS_level = True

    #### if you have an input layer for the FMG or FI please set the below parameter to False:
    Fixed_factor_FMG = True
    Fixed_factor_FI = True
    Fixed_factor_FLU = True ## in casea baseline is used instead of the current FLU based on CLC

    scaling = 100 # the scaling that is used on the factors to avoid working in float

    settings = {'year_focus': 2018,
                'dir_signature': dir_signature,
                'overwrite': overwrite,
                'Basefolder_input_data': Basefolder_input_data,
                'Basefolder_output': Basefolder_strata,
                'SOC_method': SOC_method,
                'EEA_extent_map': shp_extent_map,
                'NUTS_extent_map': shp_NUTS_borders,
                'CLC_ACC_folder': CLC_ACC_folder,
                'SOC_NUTS_factors_folder': SOC_NUTS_factors_folder,
                'run_NUTS_specific_scenario': run_NUTS_specific_scenario,
                'Stock_change_scenario': dict_default_stock_change_factors,
                'SOC_LUT_folder': SOC_LUT_folder,
                'Fixed_factor_FMG': Fixed_factor_FMG,
                'Fixed_factor_FI': Fixed_factor_FI,
                'Fixed_factor_FLU': Fixed_factor_FLU,
                'Scaling': scaling,
                'Scenario_name': scenario_name,
                'Country': Country,
                'block_based_processing': block_based_processing,
                'path_IPCC_climate_resampled': path_IPCC_climate_resampled,
                'path_IPCC_soil_resampled': path_IPCC_soil_resampled,
                'add_stats_NUTS_level': add_stats_NUTS_level,
                'commit_id': 'bcfbc8fe0eb6318c51f34ebc7a57e5f37e1d7230'}

    main_stratification(settings)




