"""In this code a stratification on the EEA39 extent will be applied
in such a way that the IPCC default values can be used for making some projections
of the SOC change for different scenarios. """
import os
import glob
from pathlib import Path
import geopandas as gpd
from loguru import logger

from SOC_scenarios.utils.soc_helper_functions import mask_raster_extent\
    , create_SOC_REF_raster, create_FLU_layer, create_factor_layer\
    ,create_SOC_scenario_layer, resample_raster_to_ref

def SOC_strat_IPCC(settings):

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
        mask_raster_extent(soil_file,settings.get('Country_extent_map'), Path(IPCC_folder).joinpath('soil_final_map')
                           , settings,overwrite=False, Country_clipping=True)
        mask_raster_extent(climate_file, settings.get('Country_extent_map'), Path(IPCC_folder).joinpath('climate_final_map'),settings,
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
        mask_raster_extent([CLC_ref_file], settings.get('Country_extent_map')
                           , Path(settings.get('Basefolder_input_data')).joinpath('CLC_ACC'),settings,
                           overwrite=False, Country_clipping=True)
        CLC_ref_file = Path(settings.get('Basefolder_input_data')).joinpath('CLC_ACC')\
                        .joinpath(Country).joinpath('CLC2018ACC_V2018_20_clipped.tif').as_posix()



    outfile_name_soil_map_resampled = Path(soil_map_to_resample[0]).parent.joinpath(Path(soil_map_to_resample[0]).stem.replace('clipped', Country).replace('1000','100')+'.tif')
    resample_raster_to_ref(soil_map_to_resample,CLC_ref_file, Country,outfile_name_soil_map_resampled.parent.as_posix(),
                           resample_factor=1, overwrite=overwrite,
                           resampling=True,
                           outname = outfile_name_soil_map_resampled.name)

    ### IPCC climate resampling
    outfile_name_climate_resampled = Path(climate_map_to_resample[0]).parent.joinpath(Path(climate_map_to_resample[0]).stem.replace('clipped', Country).replace('1000','100')+'.tif')
    resample_raster_to_ref(climate_map_to_resample,CLC_ref_file, Country,Path(outfile_name_climate_resampled).parent.as_posix(),
                           resample_factor=1, overwrite=overwrite,
                           resampling=True,
                           outname=outfile_name_climate_resampled.name)


    ## store the location of the resampled soil and climate in the settings dictionary
    settings.update({'path_IPCC_climate_resampled': outfile_name_climate_resampled,
                     'path_IPCC_soil_resampled': outfile_name_soil_map_resampled})

    ### Now the SOC REF will be created based on IPCC soil and climate data
    create_SOC_REF_raster(settings)


    ### CREATING OF STOCK CHANGE FACTORS MAPS RELATED TO
    # MANAGEMENT AND LAND USE (FLU, FI, FMG)

    ### FLU generation
    create_FLU_layer(SOC_LUT_folder,settings)

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
    if SOC_method == 'IPCC':
        SOC_strat_IPCC(settings)


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

    ###### SPECIFY SOME SCENARIO (FMG & FI) WITH DIFFERENT STOCK CHANGE
    # FACTORS FROM IPCC TO MAP THE SOC RESULT

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
    dict_stock_change_factors = {
        'Cropland': {'FMG': 1, 'FI': 2},
        'Grassland': {'FMG': 1, 'FI': 1}}

    scenario_name = 'Scenario1'


    ### Country_running
    Country = 'FR' #set to None if want to run entire EEA39 extent
    shp_Country_borders = gpd.read_file(os.path.join(dir_signature, 'etc','lulucf','AOI','NUTS_RG_20M_2021_3035.shp'))

    #### if you have an input layer for the FMG or FI please set the below parameter to False:

    Fixed_factor_FMG = True
    Fixed_factor_FI = True

    scaling = 100 # the scaling that is used on the factors to avoid working in float

    settings = {'year_focus': 2018,
                'dir_signature': dir_signature,
                'overwrite': overwrite,
                'Basefolder_input_data': Basefolder_input_data,
                'Basefolder_output': Basefolder_strata,
                'SOC_method': SOC_method,
                'EEA_extent_map': shp_extent_map,
                'Country_extent_map': shp_Country_borders,
                'CLC_ACC_folder': CLC_ACC_folder,
                'Stock_change_scenario': dict_stock_change_factors,
                'SOC_LUT_folder': SOC_LUT_folder,
                'Fixed_factor_FMG': Fixed_factor_FMG,
                'Fixed_factor_FI': Fixed_factor_FI,
                'Scaling': scaling,
                'Scenario_name': scenario_name,
                'Country': Country}


    main_stratification(settings)




