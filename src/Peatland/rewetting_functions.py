import os
from loguru import logger as log
from pathlib import Path
import pandas as pd
import geopandas as gpd
import rasterio
import numpy as np
import datetime
import shapely
from shapely.geometry.multipolygon import MultiPolygon


from SOC_scenarios.utils.soc_helper_functions import (
    add_metadata_raster)


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
        CLC_CLASSES = CLC_LULUCF_mapping.get('LEVEL_LULUCF').get(LULUCF_NAME)

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


def _get_stats_region(ID, region, dir_raster_avd_emis,
                      dir_raster_degraded, settings):

    if type(region) == gpd.GeoDataFrame:
        region = pd.DataFrame(region.loc[region.NUTS_ID == ID]).squeeze()
        region.geometry = shapely.wkt.loads(region.geometry)

    if region.geometry is None:
        return ID, pd.DataFrame()

    # Load the LULUCF mapping information
    CLC_mapping_LULUCF = settings.get('CLC_mapping_LULUCF')
    LULUCF_categories = CLC_mapping_LULUCF.get('LEVEL_LULUCF')
    CLC_dir_raster = settings.get('DATASETS').get('LUC')

    # Open first all needed rasters
    if type(region.geometry) == shapely.geometry.multipolygon.MultiPolygon:
        geom = MultiPolygon(region.geometry).geoms
    else:
        geom = MultiPolygon([region.geometry]).geoms

    with rasterio.open(CLC_dir_raster, 'r') as CLC_object:
        try:
            CLC_region, CLC_transform = rasterio.mask.mask(CLC_object,
                                                           geom,
                                                           crop=True)

        except:
            log.warning(
                f'No intersection with region {region.NUTS_ID} possible')
            return ID, pd.DataFrame()

    with rasterio.open(dir_raster_degraded, 'r') as degraded_object:
        degraded_region, degraded_transform = rasterio.mask.mask(degraded_object,
                                                                 geom,
                                                                 crop=True)

    with rasterio.open(dir_raster_avd_emis, 'r') as avd_emis_object:
        avd_emis_region, avd_emis_transform = rasterio.mask.mask(avd_emis_object,
                                                                 geom,
                                                                 crop=True)

    dict_info_per_category = {}
    for LULUCF_cat in LULUCF_categories:
        loc_LULUCF = np.where(
            np.isin(CLC_region, LULUCF_categories.get(LULUCF_cat)))

        if loc_LULUCF[0].size == 0:
            # write info in dictionary
            dict_info_per_category.update({LULUCF_cat: {
                'ha_degraded': 0,
                f'ha_{LULUCF_cat}': 0,
                'yrly_avoided_emiss_ha': None
            }})

        # Get overview of LULUCF category distribution
        # per region
        area_LULUCF = np.round(loc_LULUCF[0].size)

        # Count the nr of degraded pixels per LULUCF category
        area_degraded = np.round(
            np.where(degraded_region[loc_LULUCF] == 1)[0].size)
        loc_degraded = np.where(degraded_region == 1)

        # Get now average avoided emission per LULUCF category
        avd_emis_region_copy = np.zeros(avd_emis_region.shape)
        avd_emis_region_copy[loc_degraded] += 1
        avd_emis_region_copy[loc_LULUCF] += 1

        avd_emis_mean = np.round(
            np.mean(avd_emis_region[np.where(avd_emis_region_copy == 2)]), 2)

        # write info in dictionary
        dict_info_per_category.update({LULUCF_cat: {
            'degraded_[ha]': area_degraded,
            'LULUCF_cat_[ha]': area_LULUCF,
            'avoided_emiss_[tCO2/ha/yr]': avd_emis_mean
        }})

    df_region = pd.DataFrame.from_dict(dict_info_per_category).T
    df_region['LULUCF_cat'] = df_region.index
    df_region = df_region.reset_index(drop=True)

    region_info = pd.concat([region.to_frame().T]
                            * df_region.shape[0]).reset_index(drop=True)
    region_info = region_info[['NUTS_ID', 'LEVL_CODE',
                               'CNTR_CODE', 'NUTS_NAME', 'geometry']]
    region_info['CNTR_CODE'] = region_info['CNTR_CODE'].astype(str)
    region_info['LEVL_CODE'] = region_info['LEVL_CODE'].astype(str)

    df_merged = pd.concat([df_region, region_info], axis=1)

    return ID, df_merged


def stats_region(dir_vector, dir_raster_avd_emis,
                 dir_raster_degraded, outdir,
                 settings, sql=None):

    # Load the corresponding vector file

    VECTORS_geoms = gpd.read_file(dir_vector)

    lst_df_stats = []
    if sql is None:
        for k, region in VECTORS_geoms.iterrows():
            log.info(
                f'Processing region {str(k+1)} out of {VECTORS_geoms.shape[0]}')

            ID, df_stats = _get_stats_region(k, region, dir_raster_avd_emis,
                                             dir_raster_degraded, settings)

            lst_df_stats.append(df_stats)

    else:
        VECTORS_geoms_spark = pd.DataFrame(VECTORS_geoms)
        VECTORS_geoms_spark['geometry'] = VECTORS_geoms_spark['geometry'] = [
            item.wkt if item is not None else item for item in VECTORS_geoms_spark.geometry]
        df_spark = sql.createDataFrame(VECTORS_geoms_spark).persist()

        log.info('Start processing patch extraction on the executors ...')
        sc_output = df_spark.repartition(len(VECTORS_geoms_spark)) \
            .rdd.map(lambda row: (row.NUTS_ID,
                                  _get_stats_region(row.NUTS_ID, VECTORS_geoms, dir_raster_avd_emis,
                                                    dir_raster_degraded, settings)))\
            .collectAsMap()

        df_spark.unpersist()

        log.success('All kernel rerieval done!')

        for id in sc_output.keys():
            id_region, df_stats = sc_output[id]
            lst_df_stats.append(df_stats)

    # Merge all the stats together
    df_final = pd.concat(lst_df_stats)
    df_final = df_final.reset_index(drop=True)
    gpd_final = gpd.GeoDataFrame(df_final,
                                 geometry=df_final['geometry'])

    gpd_final.crs = 3035

    gpd_final.to_file(outdir, driver='GeoJSON')
