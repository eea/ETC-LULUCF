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
from SOC_scenarios.utils.spark import get_spark_sql
from loguru import logger as log
from itertools import product
from ast import literal_eval
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')


def create_hist(lst_input: list, dict_printing_stats: dict, outname: str,
                outfolder: str, label: str, label_x: str = 'DOY of prediction', title=False, title_str=None):

    fig, (ax1) = plt.subplots(1, figsize=(15, 10))
    ax1.hist(x=lst_input, bins=15, label=label, color='black', align='mid')
    ax1.legend(loc='upper right', fontsize=27)
    ax1.set_xlabel(label_x, fontsize=27)
    ax1.set_ylabel('# Pixels', fontsize=27)

    # add also the text to the fig
    lst_text = []
    for stat_type, stat_out in dict_printing_stats.items():
        lst_text.append(f'{stat_type}: {stat_out}')
    textstr = '\n'.join(lst_text)

    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

    # place a text box in upper left in axes coords
    ax1.text(0.05, 0.95, textstr, transform=ax1.transAxes, fontsize=20,
             verticalalignment='top', bbox=props)
    ax1.tick_params(axis="x", labelsize=25)
    ax1.tick_params(axis="y", labelsize=25)
    ax1.set_xlim([0, 150])
    if title:
        ax1.set_title(title_str,  fontdict={'fontsize': 25})

    plt.tight_layout()
    fig.savefig(os.path.join(outfolder, outname))
    plt.close()


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


def _get_LUT_strata(strata, lst_data,
                    df_strata_info, settings,
                    LUC_stats='SOC'):

    # check for which files info
    # on this strata is available

    LEVEL_LUC = f'LEVEL_{str(settings.get("Level_LUC_classes"))}'

    log.info(f'SEARCHING STATS FOR STRATA ID: {strata}')
    files_strata = [item for item in lst_data if f'STRATA_{str(strata)}'
                    in list(np.load(item, mmap_mode='r').keys())]

    log.info(f'FOUND {str(len(files_strata))} FILES FOR STRATA ID: {strata}')

    if not files_strata:
        return strata, (None)

    # Now concatenate the corresponding data
    data_strata = np.concatenate([np.load(item)[f'STRATA_{str(strata)}']
                                  for item in files_strata])
    dict_strata = {'data': data_strata}

    np.savez_compressed(os.path.join(settings.get('outfolder_compiled'),
                                     f'{LEVEL_LUC}_STRATA_{strata}.npz'),
                        **dict_strata)
    log.info('CONCATENATED ALL FILES IN SINGLE ARRAY')

    # get the corresponding stats of the strata
    mean_strata = np.round(np.nanmean(data_strata), 2)
    median_strata = np.round(np.nanmedian(data_strata), 2)
    std_strata = np.round(np.nanstd(data_strata), 2)
    perc_25, perc_75 = np.round(np.nanpercentile(data_strata, q=[25, 75]), 2)

    df_strata_stats = pd.DataFrame([mean_strata, median_strata,
                                    std_strata, perc_25, perc_75,
                                    data_strata.size]).T
    df_strata_stats.columns = [f'mean_{LUC_stats}', f'median_{LUC_stats}',
                               f'stdv_{LUC_stats}', f'perc25_{LUC_stats}',
                               f'perc75_{LUC_stats}', f'nr_px']
    df_strata_stats = df_strata_stats.reset_index(drop=True)

    # create now a histogram from the
    # distribution and the stats
    dict_printing_hist = {'median': median_strata,
                          'stdv': std_strata,
                          'nr_px': data_strata.size}
    outname_hist = f'HIST_{LEVEL_LUC}_STRATA_{strata}.png'
    strata_meta = df_strata_info.loc[df_strata_info.STRATA_ID == strata]
    ENV_ZONE = strata_meta['ENV_CAT'].values[0]
    SLOPE_ZONE = strata_meta['SLOPE_RANGE'].values[0]
    LUC_CAT = strata_meta['LUC_CAT'].values[0]
    title_str = f'ENV_{ENV_ZONE}_SLOPE_{SLOPE_ZONE}_LUC_{LUC_CAT}'

    create_hist(data_strata, dict_printing_hist,
                outname_hist, settings.get('outfolder_hist'),
                label='SOC [ton/ha]', label_x='SOC',
                title=True, title_str=title_str)

    return strata, (df_strata_stats)


def _get_stats_window(window, settings, df_stratification):
    # Iterate through the datasets and store the corresponding stats in a NPZ format
    dict_stats = {}

    window_id = f'{window[0][0]}_{window[0][1]}_' \
        f'{window[1][0]}_{window[1][1]}'

    Level_LUC = f'LEVEL_{str(settings.get("Level_LUC_classes"))}'

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
                                 "stratification", "inputs",
                                 Level_LUC)
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

    if settings.get('Retrieve_SOC_kernel'):

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
                log.info(
                    f'Processing window {str(i)} out of {str(df.shape[0])}')
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
            sc_output = df_spark.repartition(len(df)) \
                .rdd.map(lambda row: (row.window_id,
                                      _get_stats_window(((row.x0, row.x1),
                                                         (row.y0, row.y1)),
                                                        settings,
                                                        df_full_stratification)))\
                .collectAsMap()

            df_spark.unpersist()

        log.success('All kernel rerieval done!')

    if settings.get('Compile_SOC_LUT'):
        log.info('START COMPILING ALL DATA PER STRATA')

        infolder = os.path.join(settings.get("outfolder"),
                                "stratification", "inputs",
                                f'LEVEL_{str(settings.get("Level_LUC_classes"))}')

        # in this folder the compiled data per strata from
        # the input folder will be stored
        outfolder_compiled_data = Path(infolder).parent.joinpath(
            f'LEVEL_{str(settings.get("Level_LUC_classes"))}_compiled')
        outfolder_LUT = os.path.join(settings.get("outfolder"),
                                     "stratification", "LUT", 'csv')
        outfolder_hist = os.path.join(settings.get("outfolder"),
                                      "stratification", "LUT", 'hist')
        outname_LUT = f'LUT_SOC_LEVEL_{str(Level_LUC)}_V1.csv'

        # create all folders
        os.makedirs(outfolder_LUT, exist_ok=True)
        os.makedirs(outfolder_hist, exist_ok=True)
        os.makedirs(outfolder_compiled_data, exist_ok=True)

        # update the settings dict with the output folders
        settings.update({'outfolder_compiled': outfolder_compiled_data,
                         'outfolder_hist': outfolder_hist})

        # First obtain list of all data
        lst_SOC_data = [os.path.join(infolder, item) for item
                        in os.listdir(infolder)]

        if not os.path.exists(os.path.join(outfolder_LUT, outname_LUT)) or settings.get('overwrite'):

            # Below the script will be parallellized
            # if not locally processed
            lst_df_stats = []
            if sql is None:
                for i, row in df_full_stratification.iterrows():
                    log.info(
                        f'Processing window {str(i)} out of {str(df_full_stratification.shape[0])}')
                    strata_id, df_strata = _get_LUT_strata(
                        row.STRATA_ID,
                        lst_SOC_data,
                        df_full_stratification,
                        settings)

                    if df_strata is None:
                        continue

                    # Get summary stats for LUT
                    df_strata_meta = df_full_stratification.iloc[i]
                    df_strata_meta = pd.DataFrame(df_strata_meta.values).T
                    df_strata_meta.columns = row.index
                    for j in df_strata_meta.columns:
                        df_strata[j] = df_strata_meta[j]

                    lst_df_stats.append(df_strata)

            else:
                df_spark = sql.createDataFrame(df_full_stratification).persist()

                log.info('Start processing patch extraction on the executors ...')
                sc_output = df_spark.repartition(len(df_full_stratification)) \
                    .rdd.map(lambda row: (row.STRATA_ID,
                                          _get_LUT_strata(
                                              row.STRATA_ID,
                                              lst_SOC_data,
                                              df_full_stratification,
                                              settings
                                          )))\
                    .collectAsMap()

                df_spark.unpersist()

                # Now merge all the output in a single
                # dataframe
                for strata in sc_output.keys():
                    strata_id, df_strata = sc_output[strata]
                    if not df_strata is None:
                        df_strata_meta = df_full_stratification.loc[
                            df_full_stratification.STRATA_ID == strata_id]
                        df_strata_meta = pd.DataFrame(df_strata_meta.values).T
                        df_strata_meta.index = df_full_stratification.loc[
                            df_full_stratification.STRATA_ID == strata_id].columns
                        for j in df_strata_meta.index:
                            df_strata[j] = df_strata_meta.loc[j]

                        lst_df_stats.append(df_strata)

                df_all_stats = pd.concat(lst_df_stats)
                df_all_stats = df_all_stats.reset_index(drop=True)
                df_all_stats.to_csv(os.path.join(
                    outfolder_LUT, outname_LUT), index=False)


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

    overwrite = True
    # If set to False, will run on the cluster
    run_local = False

    # If the  SOC sample data should not be retrieved
    #  anymore set this to false
    retrieve_SOC_strata_kernel = False

    # Set to True if want to compile all the
    # SOC per strate to create a LUT
    # Please ensure that first the retrieval
    # of the data per kernel is finished

    compile_SOC_LUT = True

    # CLC IPCC mapping refinement contains the cross-walk between CLC and IPCC
    # at two defined levels

    settings = {'DATASETS': DATASETS,
                'CLC_cross_walk': CLC_IPCC_mapping_refinement,
                'Kernel': 128,
                'Level_LUC_classes': Level_crosswalk,
                'Env_zones_mapping': Env_zones_mapping,
                'Slope_classes': SLOPE_CAT,
                'Retrieve_SOC_kernel': retrieve_SOC_strata_kernel,
                'Compile_SOC_LUT': compile_SOC_LUT,
                'outfolder': outfolder_SOC_LUC,
                'overwrite': overwrite}

    ###########################################
    # First deal with the spark stuff
    if not run_local:
        sql = get_spark_sql(local=run_local)
    else:
        sql = None

    main_SOC_analysis(settings, sql=sql)
