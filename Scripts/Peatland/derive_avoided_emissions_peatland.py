""""
This script will help to identify degraded peatland based
one the procedure suggested by Daniel and correspondignly estimate
the avoided emissions:

Degraded peatland mapping:
    Overlay between Tanneberger peatland map and extended wetland map (2018)
    should allow to identify the degraded peatland


Determine which is the source IPCC LUC (6 options):
    Overlay the CLC with the degraded peatland map

Create pixel-based map with the avoided emissions in case of rewetting

    Use the IPCC CLIMATE zones to determine in which climatic zone the
    degraded peatland is located

    Use the IPCC tables to estimate thea avoided emissions on an annual basis
"""

# import required packages
import os
from loguru import logger as log
from pathlib import Path
import pandas as pd
from SOC_scenarios.utils.soc_helper_functions import (
    define_processing_grid,
    add_metadata_raster)
import rasterio
import numpy as np
import datetime
from SOC_scenarios.utils.spark import get_spark_sql


def _get_degraded_peatland_window(window, settings, peat_id):

    peat_occur = rasterio.open(settings.get(
        'DATASETS').get('Peatland_map')).read(1, window=window)

    ext_wetl = rasterio.open(settings.get('DATASETS').
                             get('EXT_WETL')).read(1,
                                                   window=window)
    # Hard coding needed as not correctly retrieved
    nodata_wet = 65535

    # degraded peat window
    degr_peat = np.zeros((settings.get('Kernel'),
                          settings.get('Kernel'))).astype(np.uint8)

    if np.all(ext_wetl == nodata_wet):
        no_data_window = np.ones((settings.get(
            'Kernel'), settings.get(
            'Kernel'))).astype(np.uint8) * 255
        return window, (no_data_window)

    # If peatland and not
    # in extended wetland --> degraded
    loc_degraded = np.where((peat_occur == peat_id) &
                            (ext_wetl == 0))

    # If no degraded pixels --> set all zero in window
    if loc_degraded[0].size == 0:
        return window, (degr_peat)
    else:
        degr_peat[loc_degraded] = 1
        return window, (degr_peat)


def get_degraded_peatland(settings, df_grid,
                          outfolder, outname, sql=None):
    log.info('START MAPPING DEGRADED PEATLAND')

    # open metadata from reference raster
    meta_ref = rasterio.open(settings.get(
        'DATASETS').get('Peatland_map')).meta
    meta_ref.update({'dtype': 'uint8',
                     'nodata': 255,
                     'compress': 'LZW'})

    peat_id = settings.get('peat_id')

    # Now create a byte like raster which
    # in which the degraded pixels will be stored

    with rasterio.open(os.path.join(outfolder,
                                    outname), 'w',
                       **meta_ref) as dst_out:

        # Open the needed layers
        # for identifying degraded peatland
        if sql is not None:
            df_spark = sql.createDataFrame(df_grid).persist()

            log.info('Start processing patch extraction on the executors ...')
            sc_output = df_spark.repartition(len(df_grid)) \
                .rdd.map(lambda row: (row.window_id,
                                      _get_degraded_peatland_window(((row.x0, row.x1),
                                                                     (row.y0, row.y1)),
                                                                    settings, peat_id)))\
                .collectAsMap()

            df_spark.unpersist()

            log.success('All kernel rerieval done!')

            for window_id in sc_output.keys():
                window, degr_peat = sc_output[window_id]
                dst_out.write(degr_peat, window=window,
                              indexes=1)

        else:
            for i, row in df_grid.iterrows():
                log.info(
                    f'Processing window {str(i)} out of {str(df_grid.shape[0])}')

                window, degr_peat = _get_degraded_peatland_window(((row.x0, row.x1),
                                                                   (row.y0, row.y1)),
                                                                  settings, peat_id)

                dst_out.write(degr_peat, window=window,
                              indexes=1)
    dst_out.close()


def _get_avoided_emiss_window(window, settings):
    # Open needed LUT
    LUT_drained = pd.read_csv(settings.get('TABLES').get('DRAINED_EMISS'),
                              sep=';')
    LUT_rewet = pd.read_csv(settings.get('TABLES').get('REWETTED_EMISS'),
                            sep=';')

    # Open needed DATASETS

    # Classification of CLC categories
    CLC = rasterio.open(settings.get(
        'DATASETS').get('LUC')).read(1, window=window)
    nodata_CLC = rasterio.open(settings.get(
        'DATASETS').get('LUC')).profile.get('nodata')
    IPCC_CLIMATE = rasterio.open(settings.get(
        'DATASETS').get('CLIMATE_ZONES')).read(1, window=window)
    nodata_climate = rasterio.open(settings.get(
        'DATASETS').get('CLIMATE_ZONES')).profile.get('nodata')
    DEGR_PEAT = rasterio.open(settings.get(
        'DATASETS').get('DEGRADED_PEAT')).read(1, window=window)
    nodata_degr = rasterio.open(settings.get(
        'DATASETS').get('DEGRADED_PEAT')).profile.get('nodata')

    # if all nodata just skip this window
    if np.all(CLC == nodata_CLC) or np.all(IPCC_CLIMATE == nodata_climate) \
            or np.all(DEGR_PEAT == nodata_degr):
        no_data_window = np.ones((settings.get(
            'Kernel'), settings.get(
            'Kernel'))).astype(np.int16) * 32767
        return window, (no_data_window)

    emission_window = np.zeros((settings.get(
        'Kernel'), settings.get(
        'Kernel'))).astype(np.int16)

    CLC_LULUCF_mapping = settings.get('CLC_mapping_LULUCF')

    # Iterate over the LUT to determine location
    # specific avoided emission

    column_rewet_emis = [item for item in LUT_rewet.columns if 'CO2' in item]

    for i, cat in LUT_drained.iterrows():
        column_drained_emis = [item for item in cat.index if 'CO2' in item]

        # Check which LU class and climate zone the emissions refer to
        CLIMATE_CLASS = cat['IPCC_climate_id']
        LULUCF_NAME = cat['IPCC_LU']

        # Check if for the rewetting emission factors are
        # available for the corresponding climate id
        rewet_emission = LUT_rewet.loc[LUT_rewet.IPCC_climate_id ==
                                       CLIMATE_CLASS][column_rewet_emis].values[0][0]
        if 'None' in rewet_emission:
            continue

        rewet_emission = float(rewet_emission)

        # Check also that for the drained LULUCF CAT and
        # climate ID emission factors are available
        drained_emission = cat[column_drained_emis].values[0]

        if 'None' in drained_emission:
            continue

        drained_emission = float(drained_emission)

        # Calculate now the difference in emissions
        # Recale to int
        emission_diff = drained_emission - rewet_emission
        emission_diff = int(emission_diff * 1/settings.get('scale_factor'))

        # Now based on LULUCF class use LUT table to know
        # to which CLC classes it corresponds
        CLC_CLASSES = CLC_LULUCF_mapping.get('LEVEL1').get(LULUCF_NAME)

        # find the location to attribute the emission factor
        # based on the different datasets
        lst_loc = []

        loc_CLC_classes = np.where(np.isin(CLC, CLC_CLASSES))
        lst_loc.append(loc_CLC_classes)
        loc_CLIMATE = np.where(IPCC_CLIMATE == CLIMATE_CLASS)
        lst_loc.append(loc_CLIMATE)
        loc_degr = np.where(DEGR_PEAT == 1)
        lst_loc.append(loc_degr)

        # Create now an nd array with the same shape as the window
        # to check where all conditions are met

        counter_match_strata = np.zeros(emission_window.shape)
        for loc_match in lst_loc:
            counter_match_strata[loc_match] += 1

        # find locations where all conditions apply
        loc_merged = np.where(counter_match_strata == len(lst_loc))

        if loc_merged[0].size == 0:
            continue

        emission_window[loc_merged] = emission_diff
    return window, (emission_window)


def estimate_avoided_emissions(settings, df_grid,
                               outfolder, outname,
                               sql=None):

    # first prepare a raster for writing the output

    # open metadata from reference raster
    meta_ref = rasterio.open(settings.get(
        'DATASETS').get('DEGRADED_PEAT')).meta
    shape_ref = rasterio.open(settings.get(
        'DATASETS').get('DEGRADED_PEAT')).shape
    meta_ref.update({'dtype': 'int16',
                     'nodata': 32767,
                     'compress': 'LZW'})

    # Defines which factor should be applied
    # when reading the raster

    # Some metadata that will be assigned to the raster
    # such the user has some background information
    dict_band = {
        'nodata': 32767,
        'long_name': 'Emissions (tCO2e/ha/yr)',
        'scale': settings.get('scale_factor')
    }
    dict_general = {'copyright': 'DO NOT DISTRIBUTE',
                    'doi': 'Prototype dataset, no registered DOI',
                    'institution': 'EEA',
                    'name': 'Pixel-based avoided emissions',
                    'processing_date': datetime.datetime.now().date(),
                    'source': 'Derived from IPCC tables and degraded peatland map'
                    }

    # Now create a raster
    # in which the avoided emissions will be stored

    with rasterio.open(os.path.join(outfolder,
                                    outname), 'w',
                       **meta_ref) as dst_out:

        # Open the needed layers
        # for estimating yearly avoided emissions
        if sql is not None:
            df_spark = sql.createDataFrame(df_grid).persist()

            log.info('Start processing patch extraction on the executors ...')
            sc_output = df_spark.repartition(len(df_grid)) \
                .rdd.map(lambda row: (row.window_id,
                                      _get_avoided_emiss_window(((row.x0, row.x1),
                                                                 (row.y0, row.y1)),
                                                                settings)))\
                .collectAsMap()

            df_spark.unpersist()

            log.success('All kernel rerieval done!')

            for window_id in sc_output.keys():
                window, avoid_emiss = sc_output[window_id]
                dst_out.write(avoid_emiss, window=window,
                              indexes=1)

        else:
            for i, row in df_grid.iterrows():
                log.info(
                    f'Processing window {str(i)} out of {str(df_grid.shape[0])}')

                # if ((row.x0, row.x1), (row.y0, row.y1)) != ((512, 1024), (38912, 39424)):
                #     continue

                window, avoid_emiss = _get_avoided_emiss_window(((row.x0, row.x1),
                                                                 (row.y0, row.y1)),
                                                                settings)
                dst_out.write(avoid_emiss, window=window,
                              indexes=1)

    dst_out.close()

    # Now assign some metadata description to the raster
    add_metadata_raster(dict_general,
                        dict_band,
                        os.path.join(outfolder, outname))


def main_avoided_emissions_peatland(settings, sql):
    #######################################
    # PART1: Determine degraded peatland
    #######################################

    # first obtain the degraded peatland map

    # Define the processing grid to allow
    # parallellisation of the process
    keys_dataset = list(settings.get('DATASETS').keys())
    dir_ref = settings.get('DATASETS').get(keys_dataset[0])
    xdim, ydim = rasterio.open(dir_ref).shape
    df_grid_proc = define_processing_grid(settings,
                                          xdim, ydim)
    log.info((f'A total of {df_grid_proc.shape[0]} windows have'
              ' been defined in this patch ...'))

    # Check if can access the file in the first place,
    # otherwise regeneration will be needed

    outfolder_degr_peat = os.path.join(settings.get('BASEFOLDER'),
                                       'raster')
    os.makedirs(outfolder_degr_peat, exist_ok=True)

    outname_degr_peat = 'Degraded_peatland.tif'

    # if not os.access(os.path.join(outfolder_degr_peat, outname_degr_peat), os.R_OK) \
    #         and os.path.exists(os.path.join(outfolder_degr_peat, outname_degr_peat)):  # NOQA
    #     log.warning('Cannot open degradation file --> so do processing again!')
    #     os.unlink(os.path.join(outfolder_degr_peat, outname_degr_peat))

    # if not os.path.exists(os.path.join(outfolder_degr_peat, outname_degr_peat)) or settings.get('overwrite'):

    #     if os.path.exists(os.path.join(outfolder_degr_peat, outname_degr_peat)):
    #         os.unlink(os.path.join(outfolder_degr_peat, outname_degr_peat))
    #     get_degraded_peatland(settings, df_grid_proc, outfolder_degr_peat,
    #                           outname_degr_peat, sql=sql)

    #     log.success('Degraded peatland mapping finished')

    # add patch of degraded peatland to settings
    DATASETS = settings.get('DATASETS')
    DATASETS.update({'DEGRADED_PEAT': os.path.join(outfolder_degr_peat,
                                                   outname_degr_peat)})
    settings.update({'DATASETS': DATASETS})

    ######################################################
    # PART2: Determine avoided emissions degraded peatland
    ######################################################

    # IPCC EMISSION FACTORS ARE USED
    # IN CASE OF MULTIPLE EMISSIONS FACTORS
    # DEPENDING ON DRAINAGE OR NUTRIENT STATUS
    # THE MEAN VALUE WAS TAKEN

    # FOR THE REWETTING EMISSIONS ALSO
    # THE MEAN VALULE BETWEEN NUTRIENT
    # RICH AND POOR WAS TAKEN

    # Pixel-based calculation for avoided emissions
    # on degraded peatland
    log.info('Start mapping pixel-based avoided emissions')

    outfolder_emis = os.path.join(settings.get('BASEFOLDER'),
                                  'raster')
    outname_emis = 'Avoided_emissions_yrl.tiff'

    if not os.access(os.path.join(outfolder_emis, outname_emis), os.R_OK) \
            and os.path.exists(os.path.join(outfolder_emis, outname_emis)):  # NOQA
        log.warning('Cannot open degradation file --> so do processing again!')
        os.unlink(os.path.join(outfolder_emis, outname_emis))

    if not os.path.exists(os.path.join(outfolder_emis, outname_emis)) or settings.get('overwrite'):

        if os.path.exists(os.path.join(outfolder_emis, outname_emis)):
            os.unlink(os.path.join(outfolder_emis, outname_emis))

        estimate_avoided_emissions(settings, df_grid_proc,
                                   outfolder_emis, outname_emis,
                                   sql)


if __name__ == '__main__':
    # for setting permission correctly
    oldmask = os.umask(0o002)

    log.info('*' * 50)

    from constants import (
        dir_signature,
        CLC_IPCC_mapping_refinement,
        path_IPCC_climate_resampled)

    # below the id for peatland
    # occurance is defined based
    #  on the peatland map
    peatland_id = 104

    # Define locations of raster based datasets

    DATASETS = {
        'LUC': os.path.join(dir_signature, 'input_data', 'general',
                            'CLCACC', 'CLC2018ACC_V2018_20.tif'),
        'EXT_WETL': os.path.join(dir_signature, 'etc', 'lulucf',
                                 'input_data', 'Ext_wetland',
                                 'ext_wetland_2018_v2021_nw.tif'),
        'CLIMATE_ZONES': path_IPCC_climate_resampled,
        'Peatland_map': os.path.join(dir_signature, 'etc', 'lulucf',
                                     'input_data', 'Peatland',
                                     'PeatlandMap_EEA39_compress.tif'),

    }

    # Define locations of required LUT

    TABLES = {
        'DRAINED_EMISS': os.path.join(dir_signature, 'etc', 'lulucf',
                                      'strata', 'Peatland_LUT',
                                      'IPCC_emissions_drained_peat.csv'),
        'REWETTED_EMISS': os.path.join(dir_signature, 'etc', 'lulucf',
                                       'strata', 'Peatland_LUT',
                                       'IPCC_emissions_rewetted_peat.csv')
    }

    overwrite = True
    Basefolder_output = os.path.join(dir_signature, 'etc', 'lulucf',
                                     'strata', 'Peatland')
    os.makedirs(Basefolder_output, exist_ok=True)

    # Define if the script should be parallellized or processed locally
    # First deal with the spark stuff
    run_local = False
    if not run_local:
        sql = get_spark_sql(local=run_local)
    else:
        sql = None

    # store all processing settings in dictionary
    settings = {'DATASETS': DATASETS,
                'TABLES': TABLES,
                'BASEFOLDER': Basefolder_output,
                'peat_id': peatland_id,
                'scale_factor': 0.01,
                'overwrite': overwrite,
                'CLC_mapping_LULUCF': CLC_IPCC_mapping_refinement,
                'Kernel': 1024}

    main_avoided_emissions_peatland(settings, sql=sql)
