"""
Script that will allow to run an afforestation scenario on grassland/cropland based on EO4Trees
"""


# Import needed packages
import os
import glob
from pathlib import Path
import geopandas as gpd
import pandas as pd
from loguru import logger
import numpy as np

# load the needed functions
from Biomass.utils.biom_helper_functions import \
    (define_affor_areas,
     create_affor_potential,
     align_raster_data_reference,
     create_metadata_description_afforestation,
     calc_stats_biomass_NUTS,
     add_atrributes_stats_afforestation)

from SOC_scenarios.utils.soc_helper_functions import \
    (add_metadata_raster,
     calc_weighted_average_NUTS)


def afforestation_LUT_block_proc(settings: dict):
    """
    Function that will predict the potential biomass increase if afforestation is applied
    on cropland or grassland. This function will be computed in a block-chained way per NUTS level. Furthermore,
    the function will allow to compute different afforestation scenario's per NUTS region.

    :param settings: dictionary with the settings needed to run the afforestation scenario
    :return: Output map with the requested afforestation scenario
    """

    shp_NUTS = gpd.read_file(settings.get(
        'DATASETS').get('VECTORS').get('NUTS'))
    # filter only on the NUTS levels of interest
    shp_NUTS = shp_NUTS.loc[shp_NUTS.LEVL_CODE.isin([0, 3])]

    # reproject to LAEA of not yet done
    crs = shp_NUTS.geometry.crs
    if crs is None:
        shp_NUTS.crs = shp_NUTS.crs(epsg=3035)
    elif crs.srs[-4:] != '3035':
        shp_NUTS = shp_NUTS.to_crs(epsg=3035)

    if settings.get('CONFIG_SPECS').get('Country') is not None:
        shp_NUTS = shp_NUTS.loc[((shp_NUTS.CNTR_CODE == settings.get('CONFIG_SPECS')
                                  .get('Country')[0:2]) & (shp_NUTS.LEVL_CODE == 3))]
    else:
        shp_NUTS = shp_NUTS.loc[shp_NUTS.LEVL_CODE == 3]

    # If one single NUTS region is defined for processing
    # it will be filtered on that one
    if settings.get('NUTS_region'):
        shp_NUTS = shp_NUTS.loc[shp_NUTS.NUTS_ID == settings.get('NUTS_region')]

    """RESAMPLING THE EU4TREES DATA TO THE REFERENCE RESOLUTION AND FORMAT"""
    # Ensure that the needed EU4Trees raster is
    # in the proper resolution and extent
    # If the EU-Trees4F data is not yet in the
    # proper projection system and resolution --> first resample
    Basefolder_input_data = settings.get(
        'CONFIG_SPECS').get('Basefolder_input')
    Folder_EUtrees4F_prob = os.path.join(
        Basefolder_input_data, 'EU-Trees4F', 'ens_sdms', 'prob')
    Prob_files = glob.glob(os.path.join(Folder_EUtrees4F_prob, '*.tif'))

    # filter only on the tree species for which we have an IPCC LUT
    LUT_tree_species = list(set(list(pd.read_csv(os.path.join(settings.get('CONFIG_SPECS').get('Basefolder_output'),
                                                     'NUTS_LUT_afforestation_scenario', 'JRC_yield_table',
                                                     f'LUT_C_SEQ_AFFOR_JRC_V2.csv')).SPECIES_NAME)))
    Prob_files_filtered = [item for item in Prob_files if Path(
        item).stem.split('_ens')[0] in LUT_tree_species]

    # need a reference file with the proper resolution and extent
    # this is actually hard coded now, but the CLC layer of 2018 should be eventually used in the processing chain
    # the code expects that the CLC layer is already in the proper extent.
    CLC_ref_file = settings.get('DATASETS').get('AFFORESTATION_MASK').get('CLC')
    # for prob_file in Prob_files_filtered:
    #     align_raster_data_reference(
    #         prob_file, CLC_ref_file, shp_NUTS, settings)

    # LOOP NOW over the available NUTS3 regions
    # in the area that is specified
    # Per NUTS region the afforestation
    # mask will be generated and as well the annual as the total
    # biomas increment over the prediction period

    shp_NUTS = shp_NUTS.reset_index(drop=True)
    lst_raster_files = []
    lst_NUTS_stats = []
    for index, NUTS3_region in shp_NUTS.iterrows():
        logger.info(f'CALCULATING {carbon_pool} FOR NUTS3 REGION: {NUTS3_region.NUTS_ID} \n'
                    f'{np.round((index/shp_NUTS.shape[0]) *100,2)}% COMPLETED')
        settings.update({'NUTS3_info': NUTS3_region})

        # define the output name of the scenarios per NUTS LEVEL
        outfolder = Path(settings.get('CONFIG_SPECS').get('Basefolder_output')).joinpath(f'{carbon_pool}_scenario') \
            .joinpath(settings.get('NUTS3_info')['CNTR_CODE']).as_posix()
        outname_scenario = f'{carbon_pool}_{settings.get("CONFIG_SPECS").get("scenario_name")}' \
            f'_{settings.get("NUTS3_info")["NUTS_ID"]}.tif'
        settings.update(
            {f'outdir_{carbon_pool}_scenario': os.path.join(outfolder, outname_scenario)})

        lst_raster_files.append(os.path.join(outfolder, outname_scenario))

        if not os.path.exists(os.path.join(outfolder, outname_scenario)) or settings.get('CONFIG_SPECS').get('overwrite'):

            # TODO do not write mask output but just calculate
            # the outcome because the slope is a flexible parameter
            # STEP 1: Create the afforestation mask based on
            # the defined criteria in the settings

            affor_mask_layer = define_affor_areas(settings)

            if affor_mask_layer is None:
                logger.info(
                    f'AFFORESTATION MASK COULD NOT BE GENERATED FOR NUTS REGION {NUTS3_region} --> skipped')
                continue

            # STEP : Create the afforestation carbon potential layer,
            # which includes the annual increment
            # of living biomass per suitable afforestation pixel per year
            affor_potential_layer, outname_affor_pot = create_affor_potential(
                settings, affor_mask_layer)

            if affor_potential_layer is None:
                continue

            # ADD SOME METADATA to this afforestation potential
            # layer such that it is clear
            # which configuration was used for running.
            dict_general, dict_band = create_metadata_description_afforestation(
                settings, extent='NUTS')
            add_metadata_raster(dict_general, dict_band,
                                os.path.join(outfolder, outname_affor_pot))

        # check if some stats of the map need to be provided
        if settings.get('CONFIG_SPECS').get('add_stats_NUTS') is not None:
            # based on the raster with the annual increment of LB,
            # derive now the total increment
            # of LB based on the amount of pixels that will be
            # actually afforested and the period between
            # the baseline year and the year of the prediction
            df_stats_NUTS = calc_stats_biomass_NUTS(os.path.join(
                outfolder, outname_scenario), NUTS3_region, settings)
            lst_NUTS_stats.append(df_stats_NUTS)

    if settings.get('CONFIG_SPECS').get('add_stats_NUTS') is not None:

        # define output folder for writing the result
        if settings.get('CONFIG_SPECS').get('Country') is not None:
            outfolder = Path(settings.get('CONFIG_SPECS').get('Basefolder_output')).joinpath(
                f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}_NUTS_stats').joinpath(settings.get("CONFIG_SPECS").get('Country'))
            outname = f'SOC_stats_NUTS_{settings.get("CONFIG_SPECS").get("Country")}_{settings.get("CONFIG_SPECS").get("scenario_name")}.geojson'

        else:
            outfolder = Path(settings.get('CONFIG_SPECS').get('Basefolder_output')).joinpath(
                f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}_NUTS_stats')
            outname = f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}_stats_NUTS_EEA39_{settings.get("CONFIG_SPECS").get("scenario_name")}.geojson'

        outfolder.mkdir(parents=True, exist_ok=True)

        if not os.path.exists(Path(outfolder).joinpath(outname)) or settings.get("CONFIG_SPECS").get('overwrite'):
            # merge all the statistics of each NUTS-3 region
            # in a single summary dataframe
            df_stats_all_NUTS3 = pd.concat(lst_NUTS_stats)

            # now that the NUTS3 values are calculated also the
            # weighted average at NUTS LEVEL0 will be derived
            df_stats_NUTS_final = calc_weighted_average_NUTS(
                df_stats_all_NUTS3, gpd.read_file(settings.get('DATASETS').get('VECTORS').get('NUTS')))

            # create now geodataframe out of it and write out the result
            gpd_NUTS_stats = gpd.GeoDataFrame(
                df_stats_NUTS_final, geometry=df_stats_NUTS_final.geometry)

            # add additional attribute tables (total area NUTS, spaital
            # coverage crop and grassland in NUTS region)
            # that enables to interpret the results
            gpd_NUTS_stats = add_atrributes_stats_afforestation(
                gpd_NUTS_stats, level_NUTS_focus=None)

            # write out the result
            gpd_NUTS_stats.crs = shp_NUTS.crs
            if outname.endswith('geojson'):
                gpd_NUTS_stats.to_file(
                    Path(outfolder).joinpath(outname), driver='GeoJSON')
            else:
                gpd_NUTS_stats.to_file(Path(outfolder).joinpath(outname))


def __main__(settings):
    from constants import dict_C_pool_mng_options

    type_method = settings.get('CONFIG_SPECS').get('type_method')
    block_based_processing = settings.get(
        'CONFIG_SPECS').get('block_based_processing')

    # load the carbon pool that should be adressed
    carbon_pool = settings.get('SCENARIO_SPECS').get('carbon_pool')

    if (type_method == 'LUT' and block_based_processing) and carbon_pool == 'LB':
        # if the afforestation scenario should be run,
        # this function will be called
        afforestation_LUT_block_proc(settings)


if __name__ == '__main__':

    logger.info('*' * 50)

    # Import the default parameters
    from constants import (
        Basefolder_strata,
        dir_signature,
        Basefolder_input_data,
        shp_extent_map,
        path_IPCC_climate_resampled,
        path_IPCC_soil_resampled,
        type_method,
        dict_C_pool_mng_options

    )

    ##########################################################
    ##### PART1: CONFIGURE AFFORESTATION SCENARIO ############
    ##########################################################

    # define which type of management you want to apply
    # for the future potential carbon seq scenarion
    mng_option = 'afforestation'

    carbon_pool = dict_C_pool_mng_options.get(mng_option)

    # Description of the configured parameters for afforestation scenario

    """
    Define the configuration for the flexible parameters 
    that will determine the afforestation scenario
    
    1: Slope: threshold on slope for afforestation. 
        Steeper slope than the thresholds will be afforested


    2: RCP: the emission scenario of the EUTrees4F that 
       will be used for the future occurrence probability
    3: Tree_prob: The probability of occurrence that a 
       certain tree species should have in a certain target year
       before it can be used for afforestation. 
    4: Tree_species: The tree species you want 
       to use for the afforestation
    5: Perc_reforest: The percentage of suitable 
                     (determined by the afforestation mask) 
                     grassland/cropland that should be afforested
    5: Year_potential: The year for which the carbon potential 
                        calculation should be done
    
    Note that this are the default EU factors, if a specific LUT 
    is provided with these factors at NUTS level this will
    have priority. To enable to possibilty to use NUTS specific factors, 
    please set the run_NUTS_specific_scenario 
    (see below) to True'
    """
    Year_potential = 2050
    # starting year for afforestation
    Year_baseline = 2023

   # List of CLC classes that can be used for afforestation

    lst_CLC_afforestation = [211, 212, 213, 231]

    dict_default_afforestation_factors = {
        'Slope': 0,
        'RCP': 'rcp45',
        'Tree_prob': 70,
        'Tree_species': 'Betula_pendula',
        'Perc_reforest': 10,
        'Year_potential': Year_potential,
        'input_source': 'EEA39'}

    SCENARIO_SPECS = {
        'mng_option': mng_option,
        'Year_potential': Year_potential,
        'Year_baseline': Year_baseline,
        'carbon_pool': carbon_pool,
        'lst_CLC_affor': lst_CLC_afforestation,
        'afforestation_config': dict_default_afforestation_factors
    }

    ##########################################################
    ##### PART2: CONFIGURE PROCESSING SETTINGS ###############
    ##########################################################

    # Define if the results should be created by block-based
    # processing, meaning only the window for the AOI will
    # be loaded

    # Only this type of processing is supported for the moment!!!
    block_based_processing = True

    # Country_running
    Country = None  # set to None if want to run entire EEA39 extent

    # extension of the output raster and stats that is
    # used to distinguish different scenario runs
    scenario_name = f'Scenario_JRCV3_{str(Year_potential)}_fix'

    # define the yield table LUT and forest zone LUT
    # location that should be used for processing
    name_yield_table_LUT = 'LUT_C_SEQ_AFFOR_JRC_V3.csv'
    name_LUT_forest_zones = 'LUT_FOREST_ZONE.csv'
    folder_JRC_table = os.path.join(Basefolder_strata, 
                                    'NUTS_LUT_afforestation_scenario',
                                    'JRC_yield_table')
    yield_table_LUT_dir = os.path.join(folder_JRC_table, name_yield_table_LUT)
    forest_zone_dir = os.path.join(folder_JRC_table, name_LUT_forest_zones)

    # if want to run with some NUTS specific factors
    # set the below parameter to true
    # if a country is defined only the NUTS regions
    # in that country will be considered
    run_NUTS_specific_scenario = False

    # Indicate of stats at NUTS LEVEL need to be provided
    add_stats_NUTS_level = True

    # the scaling that is used on the factors to avoid working in float
    scaling = 100

    # define the folder with the needed LUT for NUTS specific parameter loading
    NUTS_LUT_factors_folder = os.path.join(
        Basefolder_strata, 'NUTS_LUT_afforestation_scenario')

    # define if the result might be overwritten
    overwrite = False

    CONFIGURATION_SPECS = {
        'block_based_processing': block_based_processing,
        'Country': Country,
        'scenario_name': scenario_name,
        'yield_table_dir': yield_table_LUT_dir,
        'forest_zone_dir': forest_zone_dir,
        'run_NUTS_SPECIFIC': run_NUTS_specific_scenario,
        'add_stats_NUTS': add_stats_NUTS_level,
        'scaling': scaling,
        'NUTS_LUT_folder': NUTS_LUT_factors_folder,
        'overwrite': overwrite,
        'dir_signature': dir_signature,
        'type_method': type_method,
        'Basefolder_output': Basefolder_strata,
        'Basefolder_input': Basefolder_input_data
    }

    ##########################################################
    ##### PART3: DEFINE REQUIRED DATASETS  ###################
    ##########################################################

    """Load the datasets needed for defining processing extent"""

    # Define processing extent
    # The shapefile containing the NUTS borders
    shp_extent_map = gpd.read_file(shp_extent_map)
    VECTORS = {'NUTS': os.path.join(dir_signature, 'etc',
                                    'lulucf', 'AOI',
                                    'NUTS_RG_20M_2021_3035.shp'),
               'EEA_extent': shp_extent_map}

    """Load the datasets needed for afforestation mask"""

    # Below define the locations of the
    # some rasters needed for the afforestation mask

    DEM_location = os.path.join(dir_signature, 'etc', 'lulucf', 'refs', 'dem',
                                'DEM_slope_3035_100m_warped.tif')
    N2K_location = os.path.join(dir_signature, 'input_data', 'wood_provision',
                                'input_intermediate', 'N2K', 'N2K_2018_100m.tif')
    # if a region specific external mask is available,
    # please insert it below otherwise set it to 'None'
    # !! Please ensure that all the rasters are on the same grid
    external_mask = None

    ########## TODO: ETC/DI #############

    """
    Following additional datasets should be implemented to identify those areas 
    that should be protected for afforestation:
    * Nationally designated areas (CDDA): Protected areas should be EXCLUDED for afforestation
    * High nature value farmland 2012: High nature value farmland should be EXCLUDED for afforestation
    * Extended wetland layer (2018): Any wetland should be EXCLUDED for afforestation
    * Peatland map (2017): Tanneberger peatland map --> peatland should be EXCLUDED for afforestation

    The abandoned cropland locations should be PRIORITIZED for afforestation if is has been found suitable
    for afforestation based on the afforestation mask
    """

    DATASETS_afforestation = {'DEM': DEM_location,
                              'N2K': N2K_location,
                              'external_mask': external_mask,
                              'CLC': os.path.join(dir_signature, 'input_data', 'general',
                                                  'CLCACC', 'CLC2018ACC_V2018_20.tif')}

    """Load the datasets needed for defining afforestation suitability of tree species"""

    # Below define the raster(s) for defining the suitability
    #  of planting a certain tree species

    # Load the basefolder of the EU4Trees dataset
    EU4Trees_base = os.path.join(
        Basefolder_input_data, 'EU-Trees4F', 'ens_sdms', 'prob', 'resampled')
    DATASETS_AFFORESTATION_SUITABILTY = {
        'EU4Trees': EU4Trees_base
    }

    """Load IPCC datasets needed to assign strata specific increment values"""

    # Put all the IPCC related datasets used for
    # strata based assigned using the LUT
    DATASETS_IPCC = {
        'CLIMATE_ZONES': path_IPCC_climate_resampled,
        'SOIL_ZONES':  path_IPCC_soil_resampled
    }

    DATASETS_SPECS = {
        'VECTORS': VECTORS,
        'AFFORESTATION_MASK': DATASETS_afforestation,
        'AFFORESTATION_SUITABILITY': DATASETS_AFFORESTATION_SUITABILTY,
        'IPCC_STRATIFICATION': DATASETS_IPCC
    }

    ################################################################
    ## PART4: COMPILE THE ENTIRE CONFIGURATION IN A SINGLE DICT ####
    ################################################################

    # Put all the defined settings in a single dictionary such that
    # it can be easily loaded in the afforestation
    # scenario workflow

    settings = {'SCENARIO_SPECS': SCENARIO_SPECS,
                'CONFIG_SPECS': CONFIGURATION_SPECS,
                'DATASETS': DATASETS_SPECS,
                'commit_id': 'cb72cc9159d438f3e235e3e9c10e60dd2853f95c'}

    __main__(settings)
