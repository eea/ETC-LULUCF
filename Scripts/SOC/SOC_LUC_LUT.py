"""
In this script an analysis will be conducted to define distinct classed for generating
a LUT that defines the SOC for a certain stratum

Considered datasets for stratification:

* DEM (100 m)
* CLC 2018 (100 m)
* Environmental zones (100 m)
"""


# Import the necessary packages
import random
import os
import numpy as np
from pathlib import Path
import rasterio
import pandas as pd
from SOC_scenarios.utils.spark import get_spark_sql
from SOC_scenarios.utils.soc_helper_functions import (
    window_generation, define_processing_grid)
from loguru import logger as log
from itertools import product
from ast import literal_eval
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as mcolors
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
        return strata, [(None), (None)]

    # Now concatenate the corresponding data
    data_strata = np.concatenate([np.load(item)[f'STRATA_{str(strata)}']
                                  for item in files_strata])
    dict_strata = {'data': data_strata}

    np.savez_compressed(os.path.join(settings.get('outfolder_compiled'),
                                     f'{LEVEL_LUC}_STRATA_{strata}.npz'),
                        **dict_strata)
    log.info('CONCATENATED ALL FILES IN SINGLE ARRAY')

    # get cdf information for storing in dataframe
    bins = np.arange(0, 150, 1)
    cdf, bins_image = get_cdf(data_strata, bins)

    df_cdf = pd.DataFrame(cdf)
    df_cdf.columns = ['cdf']
    df_cdf = df_cdf.T
    df_cdf.index = [f'STRATA_{strata}']

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

    return strata, [(df_strata_stats), (df_cdf)]


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


def plot_cdf(dict_CDF, outfolder_hist, outname, param_asses):
    fig, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(15, 10))

    options_plotting = list(dict_CDF.keys())
    options_plotting = [item for item in options_plotting if 'cdf' in item]
    colors = mcolors.CSS4_COLORS
    random.seed(4)
    keys_colors = list(colors.keys())
    keys_colors = [item for item in keys_colors if not 'white' in item]
    random.shuffle(keys_colors)

    for option in options_plotting:
        ix = options_plotting.index(option)
        color = keys_colors[ix]
        data = dict_CDF.get(option)
        label = option.split('_cdf')[0]

        ax1.plot(dict_CDF.get(f'{label}_bins')[1:], data,
                 label=label, color=color,
                 linewidth=5)

    ax1.set_xlabel('SOC [ton/ha]', fontsize=27)
    ax1.set_xlim([0, 150])
    ax1.set_ylabel('CDF [0-1]', fontsize=27)
    ax1.legend(loc='center left', fontsize=27, bbox_to_anchor=(1, 0.5))
    ax1.tick_params(axis="x", labelsize=25)
    ax1.tick_params(axis="y", labelsize=25)
    plt.suptitle(f'IMPACT {param_asses} on SOC', fontsize=35)
    plt.tight_layout()
    fig.savefig(os.path.join(outfolder_hist, outname))
    plt.close()


def get_cdf(data, bins):
    # get now the CDF of the distribution
    count, bins_image = np.histogram(data, bins=bins)
    # finding the PDF of the histogram using count values
    pdf = count / sum(count)

    # using numpy np.cumsum to calculate the CDF
    # We can also find using the PDF values by looping and adding
    cdf = np.cumsum(pdf)

    return cdf, bins_image


def get_conversion_LUT_strata(df, to_class, stats_conv='median'):
    # create dataframe that provides an
    # overview of the LUC impact

    # check if the LUC to convert to is
    # available in the specific checked strata
    if not to_class in list(df.LUC_CAT.unique()):
        return pd.DataFrame()

    to_pool_tot = df.loc[df.LUC_CAT == to_class][f'{stats_conv}_SOC'].values[0]
    to_pool_unc = df.loc[df.LUC_CAT == to_class]['stdv_SOC'].values[0]
    to_nr_px = df.loc[df.LUC_CAT == to_class]['nr_px'].values[0]

    # only consider applicabale LUC conversions

    df_filter = df.loc[df.LUC_CAT != to_class]

    df_filter = df_filter.rename(columns={'LUC_CAT': 'from_LULUCF_cat',
                                          f'{stats_conv}_SOC': 'from_SOC',
                                          'stdv_SOC': 'from_SOC_stdv',
                                          'nr_px': 'nr_px_from'})
    df_filter = df_filter[['from_LULUCF_cat', 'from_SOC',
                           'from_SOC_stdv', 'STRATA_ID',
                           'SLOPE_CAT', 'ENV_CAT', 'SLOPE_RANGE',
                           'ENV_RANGE', 'nr_px_from']]
    df_filter['to_from_LULUCF_cat'] = to_class
    df_filter['to_SOC'] = to_pool_tot
    df_filter['to_SOC_stdv'] = to_pool_unc
    df_filter['SOC_seq'] = df_filter['to_SOC'] - df_filter['from_SOC']
    df_filter['nr_px_to'] = to_nr_px
    return df_filter


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

        # define the window on which
        # the processing should be done
        df = define_processing_grid(settings, xdim, ydim)

        log.info((f'A total of {df.shape[0]} windows have'
                  ' been defined in this patch ...'))

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
                                      "stratification", "LUT", 'hist',
                                      f'LEVEL_{str(settings.get("Level_LUC_classes"))}')
        outname_LUT = f'LUT_SOC_LEVEL_{str(Level_LUC)}_V1.csv'
        outname_CDF = f'CDF_SOC_LEVEL_{str(Level_LUC)}_V1.csv'

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

        if (not os.path.exists(os.path.join(outfolder_LUT, outname_LUT)) or not os.path.exists(os.path.join(outfolder_LUT, outname_CDF))) or settings.get('overwrite'):

            # Below the script will be parallellized
            # if not locally processed
            lst_df_stats = []
            lst_df_cdf = []
            if sql is None:
                for i, row in df_full_stratification.iterrows():
                    log.info(
                        f'Processing window {str(i)} out of {str(df_full_stratification.shape[0])}')
                    strata_id, lst_df_out = _get_LUT_strata(
                        row.STRATA_ID,
                        lst_SOC_data,
                        df_full_stratification,
                        settings)

                    df_strata = lst_df_out[0]
                    df_cdf = lst_df_out[1]

                    if df_strata is None:
                        continue

                    # Get summary stats for LUT
                    df_strata_meta = df_full_stratification.iloc[i]
                    df_strata_meta = pd.DataFrame(df_strata_meta.values).T
                    df_strata_meta.columns = row.index
                    for j in df_strata_meta.columns:
                        df_strata[j] = df_strata_meta[j]
                        df_cdf[j] = df_strata_meta[j].values[0]

                    lst_df_stats.append(df_strata)
                    lst_df_cdf.append(df_cdf)

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
                    strata_id, lst_df_out = sc_output[strata]
                    df_strata = lst_df_out[0]
                    df_cdf = lst_df_out[1]
                    if not df_strata is None:
                        df_strata_meta = df_full_stratification.loc[
                            df_full_stratification.STRATA_ID == strata_id]
                        df_strata_meta = pd.DataFrame(df_strata_meta.values)
                        df_strata_meta.columns = df_full_stratification.loc[
                            df_full_stratification.STRATA_ID == strata_id].columns

                        log.info(f'STRATA META IS {df_strata_meta}')
                        log.info(f'STRATA META SHAPE {df_strata_meta.shape}')
                        for j in df_strata_meta.columns:
                            log.info(f'STRATA to write {df_strata_meta[j]}')
                            try:
                                df_strata[j] = df_strata_meta[j].values[0]
                                df_cdf[j] = df_strata_meta[j].values[0]
                            except:
                                log.info(f'STRATA META IS {df_strata_meta}')

                        lst_df_stats.append(df_strata)
                        lst_df_cdf.append(df_cdf)

                df_all_stats = pd.concat(lst_df_stats)
                df_all_stats = df_all_stats.reset_index(drop=True)
                df_all_stats.to_csv(os.path.join(
                    outfolder_LUT, outname_LUT), index=False)

                df_all_cdf = pd.concat(lst_df_cdf)
                df_all_cdf = df_all_cdf.reset_index(drop=True)
                df_all_cdf.to_csv(os.path.join(
                    outfolder_LUT, outname_CDF), index=False)

    if settings.get('Assessment_SOC'):

        # for each dataset it will be
        # determined what is the corresponding
        # impact on SOC

        infolder = os.path.join(settings.get("outfolder"),
                                "stratification", "inputs",
                                f'LEVEL_{str(settings.get("Level_LUC_classes"))}')

        outfolder_compiled_data = Path(infolder).parent.joinpath(
            f'LEVEL_{str(settings.get("Level_LUC_classes"))}_compiled')

        outfolder_cdf = os.path.join(settings.get("outfolder"),
                                     "stratification", "LUT", 'CDF')
        os.makedirs(outfolder_cdf, exist_ok=True)

        data_folder = outfolder_compiled_data

        # load all the files with the data

        lst_data = os.listdir(data_folder)

        LEVEL_LUC = f'LEVEL_{str(settings.get("Level_LUC_classes"))}'

        bins = np.arange(0, 150, 1)

        for dataset in DATASETS:

            log.info(f'ASSESSING IMPACT FOR DATASET {dataset}')

            if dataset == 'SOC':
                continue

            # determine which strata should
            # be grouped to assess the impact

            # define unique classes for the dataset

            if not 'SLOPE' in dataset:
                dataset_options = list(
                    df_full_stratification[f'{dataset}_CAT'].unique())
            else:
                dataset_options = list(
                    df_full_stratification[f'{dataset}_RANGE'].unique())

            dict_data_option = {}

            for option in dataset_options:
                log.info(
                    f'LOADING ALL DATA FOR {option} {dataset_options.index(option)}/{len(dataset_options)}')

                if not 'SLOPE' in dataset:
                    strata_IDs_option = list(
                        df_full_stratification.loc[df_full_stratification[f'{dataset}_CAT'] == option]
                        .STRATA_ID.values)
                else:
                    strata_IDs_option = list(
                        df_full_stratification.loc[df_full_stratification[f'{dataset}_RANGE'] == option]
                        .STRATA_ID.values)

                strata_IDs_option = [
                    f'{LEVEL_LUC}_STRATA_{str(item)}.npz' for item in strata_IDs_option]

                # filter now on the files that should be loaded
                lst_data_option = [os.path.join(
                    data_folder, item) for item in lst_data if item in strata_IDs_option]

                if not lst_data_option:
                    continue

                data_option = np.concatenate([np.load(item)['data']
                                              for item in lst_data_option])

                cdf, bins_image = get_cdf(data_option, bins)

                dict_data_option.update({f'{option}_or': data_option,
                                         f'{option}_cdf': cdf,
                                         f'{option}_bins': bins_image})

            # plot the CDF of the specific assessed options
            outname = f'{LEVEL_LUC}_{dataset}_IMPACT_SOC.png'
            plot_cdf(dict_data_option, outfolder_cdf, outname, dataset)

    if settings.get('Conversion_table'):
        log.info('Start creating conversion table')

        LEVEL_LUC = f'LEVEL_{str(settings.get("Level_LUC_classes"))}'

        # for each LUC a conversion could be established
        outname_conv_table = f'LUT_CONVERSION_LUC_{LEVEL_LUC}.csv'
        # name of the file containing the SOC content per strata
        outname_LUT = f'LUT_SOC_LEVEL_{str(Level_LUC)}_V1.csv'

        outfolder_LUT = os.path.join(settings.get("outfolder"),
                                     "stratification", "LUT", 'csv')

        if not os.path.exists(os.path.join(outfolder_LUT, outname_conv_table)) or settings.get('overwrite'):
            df_LUT = pd.read_csv(os.path.join(outfolder_LUT, outname_LUT))
            conversion_options = list(df_LUT.LUC_CAT.unique())

            lst_df_conv = []

            for to_LUC in conversion_options:
                # it is no important to group the LUC in the same strata
                # to determine which will
                # be the end SOC when applying a LUC
                # in that strata
                df_conv_table = df_LUT.groupby(
                    ['ENV_CAT', 'SLOPE_RANGE']).apply(get_conversion_LUT_strata, to_LUC)
                df_conv_table = df_conv_table.reset_index(drop=True)
                lst_df_conv.append(df_conv_table)

            df_full_conv_table = pd.concat(lst_df_conv)
            df_full_conv_table = df_full_conv_table.reset_index(drop=True)
            df_full_conv_table.to_csv(os.path.join(outfolder_LUT, outname_conv_table),
                                      index=False)


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
    Level_crosswalk = 2

    # Define the output folder where the statistics will be stored
    outfolder_SOC_LUC = os.path.join(dir_signature, 'etc', 'lulucf',
                                     'strata', 'LUC')

    overwrite = False
    # If set to False, will run on the cluster
    run_local = True

    # If the  SOC sample data should not be retrieved
    #  anymore set this to false
    retrieve_SOC_strata_kernel = False

    # Set to True if want to compile all the
    # SOC per strate to create a LUT
    # Please ensure that first the retrieval
    # of the data per kernel is finished

    compile_SOC_LUT = False

    # Set to True if want to run assessment
    # on the impact of the three selected
    # stratification layers on SOC

    assess_impact_lyrs_SOC = False

    # The final step is to create a from to
    # conversion table that expresses the
    # impact of LUC on SOC (over a certain period)

    create_conversion_table = True

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
                'Assessment_SOC': assess_impact_lyrs_SOC,
                'Conversion_table': create_conversion_table,
                'outfolder': outfolder_SOC_LUC,
                'overwrite': overwrite}

    ###########################################
    # First deal with the spark stuff
    if not run_local:
        sql = get_spark_sql(local=run_local)
    else:
        sql = None

    main_SOC_analysis(settings, sql=sql)
