"""
Script that will allow to run an afforestation scenario on grassland/cropland based on EO4Trees
"""



#### Import needed packages
import os
import glob
from pathlib import Path
import geopandas as gpd
import pandas as pd
from loguru import logger
import numpy as np

### load the needed functions
from Biomass.utils.biom_helper_functions import \
    (define_affor_areas,
     create_affor_potential,
    align_raster_data_reference,
     create_metadata_description_afforestation,
     calc_stats_biomass_NUTS)

from SOC_scenarios.utils.soc_helper_functions import \
    (add_metadata_raster,
    calc_weighted_average_NUTS,
     add_atrributes_stats)


def afforestation_LUT_block_proc(settings: dict):
    """
    Function that will predict the potential biomass increase if afforestation is applied
    on cropland or grassland. This function will be computed in a block-chained way per NUTS level. Furthermore,
    the function will allow to compute different afforestation scenario's per NUTS region.

    :param settings: dictionary with the settings needed to run the afforestation scenario
    :return: Output map with the requested afforestation scenario
    """

    mng_option = 'afforestation'
    ## load the carbon pool that should be adressed
    carbon_pool = dict_C_pool_mng_options.get(mng_option)
    settings.update({'carbon_pool': carbon_pool})

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

    """RESAMPLING THE EU4TREES DATA TO THE REFERENCE RESOLUTION AND FORMAT"""
    ###### Ensure that the needed EU4Trees raster is in the proper resolution and extent
    # If the EU-Trees4F data is not yet in the proper projection system and resolution --> first resample
    Basefolder_input_data = settings.get('Basefolder_input_data')
    Folder_EUtrees4F_prob = os.path.join(Basefolder_input_data, 'EU-Trees4F', 'ens_sdms', 'prob')
    Prob_files = glob.glob(os.path.join(Folder_EUtrees4F_prob, '*.tif'))

    ## filter only on the tree species for which we have LUT
    LUT_tree_species = list(pd.read_csv(os.path.join(settings.get('Basefolder_output'),
                                                       'NUTS_LUT_afforestation_scenario',
                                                       f'EU4_trees_LUT_biom_increment.csv'), sep=';').Tree_species)
    Prob_files_filtered = [item for item in Prob_files if Path(item).stem.split('_ens')[0] in LUT_tree_species]

    ### need a reference file with the proper resolution and extent
    CLC_ref_file = os.path.join(settings.get('dir_signature'),'input_data','general','CLCACC', 'CLC2018ACC_V2018_20.tif')
    for prob_file in Prob_files_filtered:
        align_raster_data_reference(prob_file, CLC_ref_file, settings.get('EEA_extent_map'), settings)

    # TODO ADD SOME CODE THAT AUTOMATICALLY GENERATES THE CLC ACC IPCC LAYER IF NOT PRESENT
    ## check if the IPCC LU CATEGORIES MAPPING IS AVAILABLE
    file_IPCC_LU_CAT = glob.glob(os.path.join(settings.get('Basefolder_output'), 'CLC_ACC_IPCC',
                                              'CLC2018ACC_V2018_20_IPCC_LU_Categories*.tif'))
    if len(file_IPCC_LU_CAT) == 0:
        logger.error('IPCC LU CAT FILE is not available --> please generate it')
        raise
    if len(file_IPCC_LU_CAT) > 1:
        logger.error('MULTIPLE TIFF FILES IN FOLDER OF CLC IPCC LU CAT --> please have a look')
        raise
    settings.update({'file_IPCC_LU': file_IPCC_LU_CAT[0]})

    ### LOOP NOW over the available NUTS3 regions in the area that is specified

    shp_NUTS = shp_NUTS.reset_index(drop = True)
    lst_raster_files = []
    lst_NUTS_stats = []
    for index, NUTS3_region in shp_NUTS.iterrows():
        logger.info(f'CALCULATING {carbon_pool} FOR NUTS3 REGION: {NUTS3_region.NUTS_ID} \n'
                    f'{np.round((index/shp_NUTS.shape[0]) *100,2)}% COMPLETED')
        settings.update({'NUTS3_info': NUTS3_region})

        ### define the output name of the scenarios per NUTS LEVEL
        outfolder = Path(settings.get('Basefolder_output')).joinpath(f'{carbon_pool}_scenario') \
            .joinpath(settings.get('NUTS3_info')['CNTR_CODE']).as_posix()
        outname_scenario = f'{carbon_pool}_{settings.get("Scenario_name")}' \
                               f'_{settings.get("NUTS3_info")["NUTS_ID"]}.tif'
        settings.update({f'outdir_{carbon_pool}_scenario': os.path.join(outfolder, outname_scenario)})

        lst_raster_files.append(os.path.join(outfolder, outname_scenario))

        if not os.path.exists(os.path.join(outfolder, outname_scenario)) or settings.get('overwrite'):


            #TODO do not write mask output but just calculate the outcome because the slope is a flexible parameter
            ### STEP 1: Create the afforestation mask
            affor_mask_layer = define_affor_areas(settings)

            if affor_mask_layer is None:
                logger.info(f'AFFORESTATION MASK COULD NOT BE GENERATED FOR NUTS REGION {NUTS3_region} --> skipped')
                continue


            ### STEP 2: Create the afforesation carbon potential layer
            affor_potential_layer, outname_affor_pot = create_affor_potential(settings, affor_mask_layer)

            #### ADD SOME METADATA FOR SOC SCENARIO
            dict_general, dict_band = create_metadata_description_afforestation(settings, extent='NUTS')
            add_metadata_raster(dict_general, dict_band, os.path.join(outfolder, outname_affor_pot))

        ### check if some stats of the map need to be provided
        if settings.get('add_stats_NUTS_level') is not None:

            df_stats_NUTS = calc_stats_biomass_NUTS(os.path.join(outfolder, outname_scenario)
                                                , NUTS3_region,settings)
            lst_NUTS_stats.append(df_stats_NUTS)

    if settings.get('add_stats_NUTS_level') is not None:

        ### define output folder for writing the result
        if settings.get('Country') is not None:
            outfolder = Path(settings.get('Basefolder_output')).joinpath(f'{settings.get("carbon_pool")}NUTS_stats').joinpath(settings.get('Country'))
            outname = f'SOC_stats_NUTS_{settings.get("Country")}_{settings.get("Scenario_name")}.geojson'

        else:
            outfolder = Path(settings.get('Basefolder_output')).joinpath(f'{settings.get("carbon_pool")}_NUTS_stats')
            outname = f'{settings.get("carbon_pool")}_stats_NUTS_EEA39_{settings.get("Scenario_name")}.geojson'

        outfolder.mkdir(parents=True, exist_ok=True)

        if not os.path.exists(Path(outfolder).joinpath(outname)) or settings.get('overwrite'):
            df_stats_all_NUTS3 = pd.concat(lst_NUTS_stats)

            ## now that the NUTS3 is calculated also the weighted average at NUTS LEVEL0 can be derived
            df_stats_NUTS_final = calc_weighted_average_NUTS(df_stats_all_NUTS3, settings.get('NUTS_extent_map'))


            ## create now geodataframe out of it and write out the result
            gpd_NUTS_stats = gpd.GeoDataFrame(df_stats_NUTS_final, geometry=df_stats_NUTS_final.geometry)

            ## add additional attribute tables that enable to interpret the results
            gpd_NUTS_stats = add_atrributes_stats(gpd_NUTS_stats, level_NUTS_focus=3)


            ## write out the result
            gpd_NUTS_stats.crs = shp_NUTS.crs
            if outname.endswith('geojson'):
                gpd_NUTS_stats.to_file(Path(outfolder).joinpath(outname), driver='GeoJSON')
            else:
                gpd_NUTS_stats.to_file(Path(outfolder).joinpath(outname))









def main_afforestation(settings):
    mng_option = 'afforestation'

    #### In the first place the needed datasets will be
    # uploaded and will be checked if they are all in the same format

    type_method = settings.get('type_method')
    block_based_processing = settings.get('block_based_processing')

    ## load the carbon pool that should be adressed
    carbon_pool = dict_C_pool_mng_options.get(mng_option)


    if (type_method == 'LUT' and block_based_processing) and carbon_pool == 'ABGbiomass':
        afforestation_LUT_block_proc(settings)




if __name__ == '__main__':

    logger.info('*' * 50)

    ### Import the default parameters
    from constants import *
    Biom_LUT_folder = os.path.join(Basefolder_strata, 'biomass_LUT')
    os.makedirs(Biom_LUT_folder, exist_ok=True)

    ### Define if the results should be created by block-based
    ### processing, meaning only the window for the AOI will
    ### be loaded
    block_based_processing = True

    ### Country_running
    Country = None #set to None if want to run entire EEA39 extent
    shp_NUTS_borders = gpd.read_file(os.path.join(dir_signature, 'etc','lulucf','AOI','NUTS_RG_20M_2021_3035.shp'))


    ## define which type of management you want to apply for the future potential carbon seq scenarion
    ## if list contains multiple options, the percentage per option should be provided
    ## the management should be defined per LU type

    dict_scenario_management = {
        'Cropland': ['afforestation'],
        'Grassland': ['afforestation']
    }

    #### use some default factors for EU LEVEL if no NUTS specific factors are provided for afforestation

    ###### Description flexibe parameters for afforestation scenario

    """
    1: Slope: threshold on slope for afforestation. Steeper slope than the thresholds will be afforested
    2: Tree_species: The tree species you want to use for the afforestation
    3: Tree_prob: The probability that a certain tree species should have before it can be afforested 
    4: Perc_reforest: The percentage of grassland/cropland that should be afforested
    5: Year_potential: The year for which the carbon potential calculation should be done
    """
    Year_potential = 2035
    dict_default_afforestation_factors = {
        'Cropland': {'Slope': 5,
                     'RCP': 'rcp45',
                     'Tree_prob': 70,
                     'Tree_species': 'Quercus_robur',
                     'Perc_reforest': 50,
                     'Year_potential': Year_potential,
                     'input_source': 'EEA39'},
        'Grassland': {'Slope': 5,
                      'RCP': 'rcp45',
                      'Tree_prob': 70,
                      'Tree_species': 'Quercus_robur',
                      'Perc_reforest': 50,
                      'Year_potential': Year_potential,
                      'input_source': 'EEA39'}
    }

    ### extension that is used to distinguish different scenario runs
    scenario_name = f'Scenario_1_{str(Year_potential)}'

    ## Below define the locations of the different rasters needed for the afforestation scenario
    DEM_location = os.path.join(dir_signature, 'input_data', 'wood_provision',
                                'input_intermediate', 'DEM', 'dt1_slope_INCA_100m_EPSG3035.tif')
    N2K_location = os.path.join(dir_signature, 'input_data', 'wood_provision',
                                'input_intermediate', 'N2K', 'N2K_2018_100m.tif')
    external_mask = None
    dict_afforestation_masking_layers = {'DEM': DEM_location,
                                         'N2K': N2K_location,
                                         'external_mask': external_mask}


    ### if want to run with some NUTS specific factors
    ### set the below parameter to true
    ### if a country is defined only the NUTS regions in that country will be considered
    run_NUTS_specific_scenario = False
    Afforestation_NUTS_factors_folder = os.path.join(Basefolder_strata,'NUTS_LUT_afforestation_scenario')

    #### Indicate of stats at NUTS LEVEL need to be provided
    add_stats_NUTS_level = True


    scaling = 100 # the scaling that is used on the factors to avoid working in float

    settings = {'year_baseline': 2018,
                'year_potential': Year_potential,
                'dir_signature': dir_signature,
                'overwrite': overwrite,
                'Basefolder_input_data': Basefolder_input_data,
                'Basefolder_output': Basefolder_strata,
                'type_method': type_method,
                'EEA_extent_map': shp_extent_map,
                'NUTS_extent_map': shp_NUTS_borders,
                'CLC_ACC_folder': CLC_ACC_folder,
                'afforestation_NUTS_factors_folder': Afforestation_NUTS_factors_folder,
                'run_NUTS_specific_scenario': run_NUTS_specific_scenario,
                'Scaling': scaling,
                'Scenario_name': scenario_name,
                'Country': Country,
                'block_based_processing': block_based_processing,
                'path_IPCC_climate_resampled': path_IPCC_climate_resampled,
                'path_IPCC_soil_resampled': path_IPCC_soil_resampled,
                'add_stats_NUTS_level': add_stats_NUTS_level,
                'mng_options': dict_scenario_management,
                'afforestation_scenario': dict_default_afforestation_factors,
                'dict_afforestation_masking_layers': dict_afforestation_masking_layers,
                'LUT_folder': Biom_LUT_folder,
                'commit_id': 'f73a0dcd4f963d2a06e4574c5e72b8b2c12d9a6b'}

    main_afforestation(settings)




