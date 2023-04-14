"""
In this script an analysis will be conducted to define distinct classed for generating
a LUT that defines the SOC for a certain stratum

Considered datasets for stratification:

* DEM (100 m)
* CLC 2018 (100 m)
* Environmental zones (100 m)
"""


# Import the necessary packages
import os
import numpy as np
from pathlib import Path
import rasterio
import pandas as pd
from SOC_scenarios.utils.soc_helper_functions import resample_raster_to_ref
from SOC_scenarios.utils.spark import get_spark_sql
from loguru import logger as log
from itertools import product
from ast import literal_eval


def window_generation(xdim: int, ydim: int, windowsize: int, stride: int,
                      force_match_grid: bool = True):
    """
    Function that will generate a list of (unique) windows that could be processed along the defined
    dimensions
    :param xdim: the size of the x-dimension
    :param ydim: the size of the y-dimension
    :param windowsize: the size that each window should have (windowsize x windowsize).
    :param stride: the overlap that the windows may have
    :param force_match_grid: define if the windows should cover the full xdim, ydim extent even this causes an overlap
    of some of the windows. If this is set to False together with a stide of zero, there will be no overlap between
    the windows. Hence, it might happen that the windowlist does not fully cover the dimensions.
    :return: list of windows that could be processed
    """
    #
    # force_match_grid: determines that wen a window would fall outside the xdim, ydim, a reshufelling of that
    # window should be forced such that is nicely ends at the bounds of the image or not. If set to false this last
    # window will not be created and as such the last pixels in the grid are not covered
    #
    # Get the windows
    windowlist = []
    for xStart in range(0, xdim, windowsize - 2 * stride):
        for yStart in range(0, ydim, windowsize - 2 * stride):
            # We need to check if we're at the end of the master image
            # We have to make sure we have a full subtile
            # so we need to expand such tile and the resulting overlap
            # with previous subtile is not an issue
            if xStart + windowsize > xdim:
                if force_match_grid or stride > 0:
                    xStart = xdim - windowsize
                    xEnd = xdim
                else:
                    continue

            else:
                xEnd = xStart + windowsize
            if yStart + windowsize > ydim:
                if force_match_grid or stride > 0:
                    yStart = ydim - windowsize
                    yEnd = ydim
                else:
                    continue
            else:
                yEnd = yStart + windowsize

            windowlist.append(((xStart, xEnd), (yStart, yEnd)))

    return windowlist


def _get_stats_window(window, settings, df_stratification):
    # Iterate through the datasets and store the corresponding stats in a NPZ format
    dict_stats = {}

    window_id = f'{window[0][0]}_{window[0][1]}_' \
                f'{window[1][0]}_{window[1][1]}'

    for i, strata in df_stratification.iterrows():
        lst_loc_strata = []
        for dataset in settings.get('DATASETS'):
            idx_dataset = list(settings.get('DATASETS')).index(dataset)
            if dataset == 'SOC':
                continue
            nodatavalue = rasterio.open(settings.get(
                'DATASETS').get(dataset)).profile.get('nodata')
            if nodatavalue is None:
                nodatavalue = 0
            values_dataset = rasterio.open(settings.get(
                'DATASETS').get(dataset)).read(1, window=window)
            if np.all([values_dataset == nodatavalue]) and idx_dataset == 0:
                log.warning(
                    f'NO STRATIFICATION COULD BE DONE FOR WINDOW {window} DUE TO NO VALID VALUES')
                return
            # Check if the window is not partially located at a border for where no data is available
            if idx_dataset == 0:
                loc_invalid = np.where(values_dataset == nodatavalue)

            # Retrieve information on how the specific opened dataset should be stratified
            if type(strata[f'{dataset}_RANGE']) != int:
                lst_filtering = literal_eval(strata[f'{dataset}_RANGE'])
            else:
                lst_filtering = [strata[f'{dataset}_RANGE']]
            if dataset == 'SLOPE':
                loc_strata = np.where(np.logical_and(np.greater_equal(values_dataset, lst_filtering[0]),
                                                     np.less(values_dataset, lst_filtering[-1])))
            else:
                loc_strata = np.where(np.isin(values_dataset, lst_filtering))

            # no match with one of the datasets --> skip to next strata
            if loc_strata[0].size == 0:
                lst_loc_strata = []
                break

            lst_loc_strata.append(loc_strata)
        if lst_loc_strata:
            nodata_SOC = rasterio.open(settings.get(
                'DATASETS').get('SOC')).profile.get('nodata')
            values_SOC = rasterio.open(settings.get(
                'DATASETS').get('SOC')).read(1, window=window)
            if np.all([values_SOC == nodata_SOC]):
                continue

            # create an empty array to define the locations where
            # all the conditions match based
            # on the stratification datasets
            values_SOC_match = np.zeros(values_SOC.shape)
            for loc_filtering in lst_loc_strata:
                values_SOC_match[loc_filtering] += 1

            values_SOC_match[loc_invalid] = 0
            values_SOC_match[values_SOC == nodatavalue] = 0

            sample_SOC = values_SOC[values_SOC_match == len(lst_loc_strata)]
            if sample_SOC.size == 0:
                continue
            dict_stats[f'STRATA_{str(strata["STRATA_ID"])}'] = sample_SOC

    # Write to a NPZ file in the end per window
    if dict_stats:
        outfolder = os.path.join(settings.get("outfolder"),
                                 "stratification", "inputs")
        os.makedirs(outfolder, exist_ok=True)
        log.info(f'Writing sample for: {window_id}')
        np.savez(os.path.join(outfolder, window_id + '.npz'), **dict_stats)


def create_df_strata(lst_comb, cols, type_column=None):
    df_strata = pd.DataFrame(lst_comb)
    df_strata.columns = cols
    df_strata.index = [item + 1 for item in df_strata.index]
    df_strata.index.name = 'STRATA_ID'
    df_strata.reset_index(inplace=True)

    return df_strata


def main_SOC_analysis(settings, sql=None):

    # First check if stratification classes are already generated
    Level_LUC = settings.get('Level_LUC_classes')
    outname = f'Stratification_SOC_LUC_classes_LEVEL{str(Level_LUC)}.csv'
    if not os.path.exists(os.path.join(settings.get('outfolder'), 'stratification', outname))  \
            or settings.get('overwrite'):
        # Now for each dataset layer the different classes will be defined
        log.info(
            f'{os.path.join(settings.get("outfolder"), "stratification", outname)} does not exist')  # NOQA

        # 1 DEM SLOPE
        SLOPE_CAT = list(settings.get('Slope_classes').keys())
        SLOPE_RANGE = [settings.get('Slope_classes').get(item)
                       for item in SLOPE_CAT]

        # 2 Environmental zones:
        ENV_CAT = list(settings.get('Env_zones_mapping').keys())
        ENV_RANGE = [settings.get('Env_zones_mapping').get(item)
                     for item in ENV_CAT]
        # Translation table env classes

        # 3 IPCC LULUCF classes based on CLC
        LUC_CAT = list(settings.get('CLC_cross_walk').get(
            f'LEVEL{str(Level_LUC)}').keys())
        LUC_RANGE = [settings.get('CLC_cross_walk').get(
            f'LEVEL{str(Level_LUC)}').get(item) for item in LUC_CAT]

        # Get now a list of all possible combinations for stratification
        All_combinations_CAT = list(product(SLOPE_CAT, ENV_CAT, LUC_CAT))
        All_combinations_RANGE = list(
            product(SLOPE_RANGE, ENV_RANGE, LUC_RANGE))

        # Now turn both list to a pandas dataframe where the ID of the strata will be defined
        # This ID will be used to join both lists together to obtain a full dataframe expressing
        # the entire defined stratification
        df_all_range = create_df_strata(All_combinations_RANGE, [
            'SLOPE_RANGE', 'ENV_RANGE', 'LUC_RANGE'])
        df_all_cat = create_df_strata(
            All_combinations_CAT, ['SLOPE_CAT', 'ENV_CAT', 'LUC_CAT'])

        df_full_stratification = pd.merge(
            df_all_cat, df_all_range, on='STRATA_ID')
        df_full_stratification['LEVEL_LUC'] = [
            Level_LUC] * df_full_stratification.shape[0]
        df_full_stratification.to_csv(os.path.join(settings.get(
            'outfolder'), 'stratification', outname), index=False)

    # Open the stratification file needed for extracting the values per window
    df_full_stratification = pd.read_csv(os.path.join(
        settings.get('outfolder'), 'stratification', outname))

    # check if the dimensions of all the datasets used for stratication are the same
    for dataset_name in settings.get('DATASETS'):
        ix = list(settings.get('DATASETS').keys()).index(dataset_name)
        if ix == 0:
            # Open the grid of one of the raster layers
            xdim, ydim = rasterio.open(settings.get(
                'DATASETS').get(dataset_name)).shape
        else:
            if rasterio.open(settings.get('DATASETS').get(dataset_name)).shape != (xdim, ydim):
                raise ValueError(
                    f'Dimensions of datasets do not match {dataset_name}')

    # Get list of windows that should be processed
    windowlist = window_generation(xdim, ydim,
                                   settings.get('Kernel'), stride=0,
                                   force_match_grid=False)

    log.info((f'A total of {len(windowlist)} windows have'
              ' been defined in this patch ...'))

    # We create unique window IDs, based on the window position itself,
    # and the EEA grid
    window_ids = ['_'.join([str(window[0][0]),
                            str(window[0][1]), str(window[1][0]),
                            str(window[1][1])]) for window in windowlist]

    # Put it together in a dataframe
    df = pd.DataFrame(window_ids).rename(columns={0: 'window_id'})
    df['x0'] = [window[0][0] for window in windowlist]
    df['x1'] = [window[0][1] for window in windowlist]
    df['y0'] = [window[1][0] for window in windowlist]
    df['y1'] = [window[1][1] for window in windowlist]

    # Below the script will be parallellized
    # if not locally processed

    if sql is None:
        for i, row in df.iterrows():
            log.info(f'Processing window {str(i)} out of {str(df.shape[0])}')
            # if ((row.x0, row.x1), (row.y0, row.y1)) == ((768, 896), (40448, 40576)):
            _get_stats_window(
                ((row.x0, row.x1),
                 (row.y0, row.y1)),
                settings,
                df_full_stratification
            )
    else:
        df_spark = sql.createDataFrame(df).persist()

        log.info('Start processing patch extraction on the executors ...')
        df_spark.repartition(len(df)) \
            .rdd.map(lambda row: (row.window_id,
                                  _get_stats_window(((row.x0, row.x1),
                                                     (row.y0, row.y1)),
                                                    settings,
                                                    df_full_stratification)))\

        df_spark.unpersist()

    log.success('All done!')


if __name__ == '__main__':
    # for setting permission correctly
    oldmask = os.umask(0o002)

    log.info('*' * 50)

    from constants import (
        dir_signature,
        CLC_IPCC_mapping_refinement,
        Env_zones_mapping)

    # The DEM may not be set as the first dataset!!!!!
    DATASETS = {
        'ENV': os.path.join(dir_signature, 'etc', 'lulucf', 'input_data',
                            'EnvZones', 'eea_r_3035_100_m_EnvZ-Metzger_2020_v1_r00.tif'),
        'SLOPE':  os.path.join(dir_signature, 'etc', 'lulucf', 'refs', 'dem', 'DEM_slope_3035_100m_warped.tif'),
        'LUC': os.path.join(dir_signature, 'input_data', 'general',
                            'CLCACC', 'CLC2018ACC_V2018_20.tif'),
        'SOC': os.path.join(dir_signature, 'etc', 'lulucf', 'refs', 'isric', 'ocs_0-30cm_mean_3035_100m.tif')
    }

    # Below the slope categories are defined
    # For stratification

    SLOPE_CAT = {
        'FLAT': [0, 5],
        'MODERATE': [5, 20],
        'STEEP': [20, 1000]
    }

    # Below define level of IPCC CLC crosswalk table for SOC LUT
    Level_crosswalk = 1

    # Define the output folder where the statistics will be stored
    outfolder_SOC_LUC = os.path.join(dir_signature, 'etc', 'lulucf',
                                     'strata', 'LUC')

    overwrite = False
    # If set to False, will run on the cluster
    run_local = True

    # CLC IPCC mapping refinement contains the cross-walk between CLC and IPCC
    # at two defined levels

    settings = {'DATASETS': DATASETS,
                'CLC_cross_walk': CLC_IPCC_mapping_refinement,
                'Kernel': 128,
                'Level_LUC_classes': Level_crosswalk,
                'Env_zones_mapping': Env_zones_mapping,
                'Slope_classes': SLOPE_CAT,
                'outfolder': outfolder_SOC_LUC,
                'overwrite': overwrite}

    # resample_raster_to_ref([DATASETS.get('DEM')], DATASETS.get('CLC'), None,
    #                        r'L:\etc\lulucf\refs\dem',
    #                        resample_method='bilinear',
    #                        resample_factor=1, overwrite=False,
    #                        resampling=True,
    #                        dtype='Byte',
    #                        outname='DEM_slope_3035_100m_warped.tif')

    ###########################################
    # First deal with the spark stuff
    if not run_local:
        sql = get_spark_sql(local=run_local)
    else:
        sql = None

    main_SOC_analysis(settings, sql=sql)
