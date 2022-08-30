import logging
from loguru import logger
import os
import rasterio
import rasterio.mask
import geopandas as gpd
from pathlib import Path
import numpy as np
import gdal
import subprocess
import osr
import pandas as pd
import glob

def compress_file(infile, outfile):
    cmd = f'gdal_translate -co compress=LZW {infile} {outfile}'
    try:
        subprocess.check_call(cmd, shell=True, stderr= subprocess.STDOUT)
        os.unlink(infile)

    except subprocess.CalledProcessError as e:
        raise RuntimeError ('command {} return with error code {}: {}'.format(e.cmd, e.returncode, e.output))
def mask_raster_extent(raster_files, shp_extent, outdir,
                       settings,
                       overwrite = False, scale = 1,
                       CLC_layer = False,
                       Country_clipping = False):
    if Country_clipping:
        Country = settings.get('Country')
        outdir = os.path.join(outdir, Country)
        shp_extent = shp_extent.loc[((shp_extent.LEVL_CODE == 0) & (shp_extent.NUTS_ID == Country[0:2].upper()))]

    os.makedirs(outdir, exist_ok=True)
    for raster_file in raster_files:
        outname = Path(raster_file).stem + '_clipped.tif'

        if not os.path.exists(os.path.join(outdir,outname)) or overwrite:
            if not CLC_layer:
                with rasterio.open(raster_file, 'r') as CLC_object:
                    out_image, out_transform = rasterio.mask.mask(CLC_object, shp_extent.geometry.values, crop=True)
                    out_meta = CLC_object.meta
                    out_image[out_image != CLC_object.nodata] = out_image[out_image != CLC_object.nodata] * scale
                    write_raster(out_image, out_meta, out_transform, outdir, outname)
            else:
                upper_left_corner_x = shp_extent.geometry.bounds.minx.values[0]
                lower_right_corner_y = shp_extent.geometry.bounds.miny.values[0]
                lower_right_corner_x = shp_extent.geometry.bounds.maxx.values[0]
                upper_left_corner_y = shp_extent.geometry.bounds.maxy.values[0]

                cmd = f'gdal_translate -projwin {str(upper_left_corner_x)} {str(upper_left_corner_y)} {str(lower_right_corner_x)} {str(lower_right_corner_y)} ' \
                      f'-co COMPRESS=LZW -of GTiff {raster_file} {os.path.join(outdir, outname)}'
                try:
                    subprocess.check_call(cmd, shell=True, stderr= subprocess.STDOUT)
                except subprocess.CalledProcessError as e:
                    raise RuntimeError ('command {} return with error code {}: {}'.format(e.cmd, e.returncode, e.output))


def create_SOC_REF_raster(settings):

    ## load the default settings
    SOC_LUT_folder = settings.get('SOC_LUT_folder')
    outfile_name_soil_map_resampled = settings.get('path_IPCC_soil_resampled')
    outfile_name_climate_resampled = settings.get('path_IPCC_climate_resampled')
    overwrite = settings.get('overwrite')

    ### now the metadata to define the SOCref will be defined
    SOC_ref_file = [item for item in glob.glob(os.path.join(SOC_LUT_folder, '*.csv')) if 'SOCref_LUT'in item][0]
    df_SOCref_LUT = pd.read_csv(SOC_ref_file, sep = ';')

    ### load up the resampled soil and climate
    # to allow calculation of the SOC-ref on the fly
    if settings.get('Country') is None:
        outdir_SOCREF = Path(SOC_LUT_folder).parent.joinpath('SOC_REF').as_posix()
    else:
        outdir_SOCREF = Path(SOC_LUT_folder).parent.joinpath('SOC_REF').joinpath(settings.get('Country')).as_posix()
        os.makedirs(outdir_SOCREF, exist_ok=True)

    outname_SOC_REF_raster = 'SOC_REF_IPCC_climate_soil_100m.tif'


    if not os.path.exists(os.path.join(outdir_SOCREF, outname_SOC_REF_raster)) and not overwrite:

        with rasterio.open(outfile_name_soil_map_resampled, 'r') as file:
            soil_raster = file.read(1)
            out_meta = file.meta
            out_transform = file.transform
        climate_raster = rasterio.open(outfile_name_climate_resampled).read(1)

        SOCref_raster = np.full(soil_raster.shape,255).astype(out_meta.get('dtype'))

        ## only keep the classes for which SOCref can be calculated
        df_SOCref_LUT_only_valid = df_SOCref_LUT.mask(df_SOCref_LUT.eq('None')).dropna()

        for i, strata_group in df_SOCref_LUT_only_valid.iterrows():
            logger.info(f'Converting to SOCREF '
                        f' for soil:{strata_group["IPCC_soil_id"]}&'
                        f' climate: {strata_group["IPCC_climate_id"]}')
            soil_value = strata_group['IPCC_soil_id']
            climate_value = strata_group['IPCC_climate_id']

            ### find the pixels that match the condition for this strata

            loc_strata = np.where((climate_raster == climate_value) & (soil_raster == soil_value))

            SOCref_raster[loc_strata] = strata_group['SOC_REF (ton*C/ha)']

        out_meta.update({'nodata': 255})
        write_raster(np.expand_dims(SOCref_raster,0),out_meta,out_transform,Path(outdir_SOCREF)
                     , Path(outname_SOC_REF_raster).name)


def create_FLU_layer(SOC_LUT_folder, settings):

    ### load default parameters
    year_focus = settings.get('year_focus')
    CLC_ACC_folder_original = settings.get('CLC_ACC_folder')
    CLC_ACC_folder_Country = Path(settings.get('Basefolder_input_data')).joinpath('CLC_ACC')
    overwrite = settings.get('overwrite')
    path_climate_raster = settings.get('path_IPCC_climate_resampled').as_posix()
    Basefolder_output_data = settings.get('Basefolder_strata')

    if settings.get('Country') is None:
        CLC_ACC_file = [item for item in glob.glob(os.path.join(CLC_ACC_folder_original, '*.tif')) if 'CLC{}'.format(str(year_focus)) in Path(item).stem][0]
    else:
        CLC_ACC_file = [item for item in glob.glob(os.path.join(CLC_ACC_folder_Country, settings.get('Country'), '*.tif'))
                        if 'CLC{}'.format(str(year_focus)) in Path(item).stem][0]

    ### open LUT FLU mapping
    df_CLC_FLU_conversion = pd.read_csv(os.path.join(SOC_LUT_folder, 'IPCC_FLU_CLC_mapping_LUT.csv'), sep=';')
    #only consider non empty rows for mapping
    df_CLC_FLU_conversion_dropna = df_CLC_FLU_conversion.mask(df_CLC_FLU_conversion.eq('None'))\
                                                        .dropna(subset=['CLC_id'])

    ### open LUT FLU factors
    df_FLU_factors_LUT = pd.read_csv(os.path.join(SOC_LUT_folder, 'IPCC_LUT_factors_FLU.csv'), sep=';')
    ## drop locations with no factor given
    df_FLU_factors_LUT_dropna = df_FLU_factors_LUT.mask(df_FLU_factors_LUT.eq('None'))\
                                                  .dropna(subset=['Factor'])


    # remap CLC according to IPCC
    # land use stock change factors
    outname_CLC_FLU_mapped = Path(CLC_ACC_file).stem + '_IPCC_FLU_mapped.tif'
    if settings.get('Country') is None:
        outdir_CLC_FLU_mapped = Path(Basefolder_output_data).joinpath('SOC_FLU')
    else:
        outdir_CLC_FLU_mapped = Path(Basefolder_output_data).joinpath('SOC_FLU').joinpath(settings.get('Country'))

    outdir_CLC_FLU_mapped.mkdir(parents=True, exist_ok=True)

    ### Also create a layer that defines the IPCC land use categories
    ### from CLC
    outname_CLC_IPCC_LU_category = Path(CLC_ACC_file).stem + '_IPCC_LU_Categories_Grassland_Cropland.tif'
    if settings.get('Country') is None:
        outdir_CLC_IPCC_LU_category = Path(Basefolder_output_data).joinpath('CLC_ACC_IPCC')
    else:
        outdir_CLC_IPCC_LU_category = Path(Basefolder_output_data).joinpath('CLC_ACC_IPCC').joinpath(settings.get('Country'))

    outdir_CLC_IPCC_LU_category.mkdir(parents=True, exist_ok=True)

    if ((not os.path.exists(outdir_CLC_FLU_mapped.joinpath(outname_CLC_FLU_mapped)) and not overwrite)\
            or (not os.path.exists(outdir_CLC_IPCC_LU_category.joinpath(outname_CLC_IPCC_LU_category)) and not overwrite)):

        with rasterio.open(CLC_ACC_file, 'r') as CLC_file:
            CLC_raster = CLC_file.read(1)
            CLC_meta = CLC_file.meta
            CLC_transform = CLC_file.transform

        if not os.path.exists(outdir_CLC_FLU_mapped.joinpath(outname_CLC_FLU_mapped)) and not overwrite:
            ## open the IPCC climate raster which is also needed to calculate the factors
            climate_raster = rasterio.open(path_climate_raster).read(1)


            ## create FLU factors layer
            CLC_remap_FLU = np.full(CLC_raster.shape,255).astype(np.uint8)
            CLC_meta.update({'nodata': 255, 'dtype': 'uint8'})
            scaling = settings.get('Scaling') ## just to avoid the saving in floats of the raster


            for i, conversion_row in df_CLC_FLU_conversion_dropna.iterrows() :
                IPCC_FLU_id = int(conversion_row['IPCC_FLU_id'])
                for p, conversion_FLU_factor in df_FLU_factors_LUT_dropna.loc[df_FLU_factors_LUT_dropna['IPCC_FLU_id'] == IPCC_FLU_id].iterrows():
                    IPCC_climate_ID = int(conversion_FLU_factor['IPCC_climate_id'])
                    logger.info(f'CALCULATING FLU FACTOR FOR CLIMATE CATEGORY:{IPCC_climate_ID}\n'
                                f' & LU CATEGORY:{IPCC_FLU_id}')
                    IPCC_factor_strata = float(conversion_FLU_factor['Factor'])
                    CLC_remap_FLU[((CLC_raster == int(conversion_row['CLC_id']))&(climate_raster == IPCC_climate_ID))] = int(IPCC_factor_strata*scaling)

            write_raster(np.expand_dims(CLC_remap_FLU,0),CLC_meta,CLC_transform,outdir_CLC_FLU_mapped, outname_CLC_FLU_mapped)

        if not os.path.exists(outdir_CLC_IPCC_LU_category.joinpath(outname_CLC_IPCC_LU_category)) and not overwrite:
            ## create IPCC LU category mapping layer:
            CLC_remap_IPCC_LUCAT = np.full(CLC_raster.shape,255).astype(np.uint8)
            CLC_meta.update({'nodata': 255, 'dtype': 'uint8'})
            for i, conversion_row in df_CLC_FLU_conversion_dropna.iterrows():
                LUCAT_ID = int(conversion_row['IPCC_landuse_id'])
                CLC_ID = int(conversion_row['CLC_id'])
                logger.info(f'CONVERTING CLC: {CLC_ID} TO LU CATEGORY:{LUCAT_ID}')
                CLC_remap_IPCC_LUCAT[CLC_raster == CLC_ID] = LUCAT_ID
            write_raster(np.expand_dims(CLC_remap_IPCC_LUCAT,0),CLC_meta,CLC_transform,outdir_CLC_IPCC_LU_category, outname_CLC_IPCC_LU_category)




def create_factor_layer(settings, type_factor = 'FMG',fixed_factor_creation = True):
    """
    Returns a layer containing the IPCC default stock change
    factors for the defined management  or input activities activities
    :param SOC_LUT_folder: The folder in which the LUT are stored for SOC estimation based on IPCC tables
    :param settings: Dictionary with different settings that define which type of management or input is applied
    :param fixed_factor_creation: If the management layer is already available, just ignore the defined settings and load up this file
                      to derive IPCC default factors. Set to true in this case
    :param type_factor: indicates for which stock change factor the calculation sgould be done. Options: FI, FMG
    :return: Raster layer with the IPCC stock change factor based on the type of factor addressed
    """

    #TODO add part in which the user can plug-in their own management layer, which then will be pre-processed to derive the stock change factors

    ### Load the settings
    Basefolder_strata_output = Path(settings.get('SOC_LUT_folder')).parent.as_posix()
    factor_scenario = settings.get('Stock_change_scenario')
    path_climate_raster = settings.get('path_IPCC_climate_resampled').as_posix()
    overwrite = settings.get('overwrite')


    ### First open the FLU layer to know which crops are where located
    ## the combination of the LUT and the CLC FLU mapped layer allows this
    if settings.get('Country') is None:
        outdir_IPCC_LUCAT = Path(Basefolder_strata_output).joinpath('CLC_ACC_IPCC')
    else:
        outdir_IPCC_LUCAT = Path(Basefolder_strata_output).joinpath('CLC_ACC_IPCC').joinpath(settings.get('Country'))
    CLC_IPCC_LUCAT_dir = glob.glob(os.path.join(outdir_IPCC_LUCAT, 'CLC{}ACC*Grassland_Cropland.tif'.format(settings.get("year_focus"))))
    df_CLC_FLU_conversion = pd.read_csv(os.path.join(Basefolder_strata_output, 'SOC_LUT', 'IPCC_FLU_CLC_mapping_LUT.csv'), sep=';')

    ### Define the output directory where the result will be stored
    if not fixed_factor_creation:
        outname_factor_mapped = f'IPCC_{type_factor}_mapped_input_management_layer.tif'
    else:
        outname_factor_mapped = f'IPCC_{type_factor}_mapped_scenario'
        for Crop in factor_scenario.keys():
            df_LUT_factor = pd.read_csv(os.path.join(Basefolder_strata_output, 'SOC_LUT', f'IPCC_LUT_factors_{type_factor}_{Crop}.csv'), sep=';')
            factor_name = df_LUT_factor.loc[df_LUT_factor[f'IPCC_{type_factor}_id']==factor_scenario.get(Crop).get(type_factor)][[f'IPCC_{type_factor}_name']].values[0][0]
            outname_factor_mapped = outname_factor_mapped + '_'+ Crop + '_' + factor_name
        outname_factor_mapped = outname_factor_mapped + '.tif'

    if settings.get('Country') is None:
        outdir_CLC_FLU_mapped = Path(Basefolder_strata_output).joinpath(f'SOC_{type_factor}')
    else:
        outdir_CLC_FLU_mapped = Path(Basefolder_strata_output).joinpath(f'SOC_{type_factor}').joinpath(settings.get('Country'))

    outdir_CLC_FLU_mapped.mkdir(parents=True, exist_ok=True)

    if not os.path.exists(outdir_CLC_FLU_mapped.joinpath(outname_factor_mapped)) and not overwrite:
        with rasterio.open(CLC_IPCC_LUCAT_dir[0], 'r') as CLC_IPCC_LUCAT_file:
            CLC_IPCC_LUCAT_raster = CLC_IPCC_LUCAT_file.read(1)

        ## now also open the climate raster to know in which climate regime we are
        with rasterio.open(path_climate_raster, 'r') as climate_file:
            Climate_raster = climate_file.read(1)
            meta_raster = climate_file.meta
            transform_raster = climate_file.transform

        meta_raster.update({'nodata': 255, 'dtype': 'uint8'})



        ### Create factor raster to write out factors in it
        factor_raster = np.full(CLC_IPCC_LUCAT_raster.shape,255).astype('uint8')
        scaling = settings.get('Scaling') ## apply some scaling to the factors to allow working with integers instead of floats


        for Crop in factor_scenario.keys():
            df_factor_IPCC_factors = pd.read_csv(os.path.join(Basefolder_strata_output, 'SOC_LUT', f'IPCC_LUT_factors_{type_factor}_{Crop}.csv'), sep=';')
            ## drop locations with no factor given
            df_factor_IPCC_factors_dropna = df_factor_IPCC_factors.mask(df_factor_IPCC_factors.eq('None')) \
                .dropna(subset=['Factor'])

            if fixed_factor_creation:
                factor_id = factor_scenario.get(Crop).get(type_factor)

            ## get the FLU classes that represents the Crop of interest
            IPCC_landuse_id = list(set(df_CLC_FLU_conversion.loc[df_CLC_FLU_conversion.IPCC_landuse_name == Crop]['IPCC_landuse_id'].values))

            for land_use_id in IPCC_landuse_id:
                logger.info('CALCULATING factor FOR LAND USE CATEGORY: {}'
                            .format(df_CLC_FLU_conversion.loc[df_CLC_FLU_conversion.IPCC_landuse_id == land_use_id]['IPCC_landuse_name'].values[0]))

                for i, conversion_row in df_factor_IPCC_factors_dropna.iterrows():
                    if fixed_factor_creation:
                        if int(conversion_row[f'IPCC_{type_factor}_id']) != int(factor_id):
                            continue
                    IPCC_climate_ID = int(conversion_row['IPCC_climate_id'])
                    IPCC_factor_strata = float(conversion_row['Factor'])
                    logger.info(f'CALCULATING {type_factor} FACTOR FOR FACTOR CATEGORY: {conversion_row[f"IPCC_{type_factor}_id"]} & CLIMATE ID: {conversion_row["IPCC_climate_id"]}')
                    factor_raster[((Climate_raster == IPCC_climate_ID)&(CLC_IPCC_LUCAT_raster == land_use_id))] = IPCC_factor_strata*scaling

        write_raster(np.expand_dims(factor_raster,0),meta_raster,transform_raster,outdir_CLC_FLU_mapped, outname_factor_mapped)


def create_SOC_scenario_layer(settings):
    ### Load the settings
    SOC_LUT_folder = settings.get('SOC_LUT_folder')
    overwrite = settings.get('overwrite')
    factor_scenario = settings.get('Stock_change_scenario')
    scaling = settings.get('Scaling')


    ### Define output name and output dir for scenario SOC
    if settings.get('Country') is None:
        outfolder = Path(SOC_LUT_folder).parent.joinpath('SOC_scenario')
        outname_SOC_scenario = f'SOC_{settings.get("Scenario_name")}_EEA39.tif'

    else:
        outfolder = Path(SOC_LUT_folder).parent.joinpath('SOC_scenario').joinpath(settings.get('Country'))
        outname_SOC_scenario = f'SOC_{settings.get("Scenario_name")}_{settings.get("Country")}.tif'


    outfolder.mkdir(parents = True,exist_ok=True)

    if os.path.exists(os.path.join(outfolder, outname_SOC_scenario)) and not overwrite:
        logging.info('SOC scenario file already created --> SKIP')
        return

    ### Load now the needed datasets for the calculation

    ## SOCREF
    if settings.get('Country') is None:
        path_SOCREF = Path(SOC_LUT_folder).parent.joinpath('SOC_REF').joinpath('SOC_REF_IPCC_climate_soil_100m.tif').as_posix()
    else:
        path_SOCREF = Path(SOC_LUT_folder).parent.joinpath('SOC_REF').joinpath(settings.get('Country'))\
            .joinpath('SOC_REF_IPCC_climate_soil_100m.tif').as_posix()


    ## FMG & FI
    ## use this list to find the correct scenario

    def find_fixed_factor_layer(settings,type_factor = 'FMG'):
        lst_factors_selected = []
        for Crop in factor_scenario.keys():
            df_LUT_factor = pd.read_csv(os.path.join(SOC_LUT_folder, f'IPCC_LUT_factors_{type_factor}_{Crop}.csv'), sep=';')
            factor_name = df_LUT_factor.loc[df_LUT_factor[f'IPCC_{type_factor}_id']==factor_scenario.get(Crop).get(type_factor)][[f'IPCC_{type_factor}_name']].values[0][0]
            lst_factors_selected.append(factor_name)

        if settings.get('Country') is None:
            folder_factor = Path(SOC_LUT_folder).parent.joinpath(f'SOC_{type_factor}')
        else:
            folder_factor = Path(SOC_LUT_folder).parent.joinpath(f'SOC_{type_factor}').joinpath(settings.get('Country'))


        # Find now the proper file
        if len(lst_factors_selected) == 2:
            path_factor = [item for item in glob.glob(os.path.join(folder_factor, f'*{lst_factors_selected[0]}_*_{lst_factors_selected[1]}.tif'))]
        elif len(lst_factors_selected) == 1:
            path_factor = [item for item in glob.glob(os.path.join(folder_factor, f'*{lst_factors_selected[0]}.tif'))]
        if len(path_factor) >1:
            raise logger.error('More then one factor scenario file is available which is not possible --> Check')

        return path_factor[0]

    if settings.get('Fixed_factor_FMG'):
        path_FMG = find_fixed_factor_layer(settings,type_factor='FMG')
    if settings.get('Fixed_factor_FI'):
        path_FI = find_fixed_factor_layer(settings,type_factor='FI')

    #FLU
    if settings.get('Country') is None:
        path_FLU = glob.glob(os.path.join(Path(SOC_LUT_folder).parent.as_posix(), 'SOC_FLU', '*_IPCC_FLU_mapped.tif'))
    else:
        path_FLU = glob.glob(os.path.join(Path(SOC_LUT_folder).parent.as_posix(), 'SOC_FLU', settings.get('Country'), '*_IPCC_FLU_mapped.tif'))

    if len(path_FLU) > 1:
        raise logger.error('MULTIPEL FLU FILES AVAILABLE --> Please check')
    path_FLU = path_FLU[0]


    ### Now that all locations for the rasters are defined it can be opened

    FMG_layer = rasterio.open(path_FMG).read(1)
    no_data_FMG = rasterio.open(path_FMG).meta.get('nodata')
    FI_layer = rasterio.open(path_FI).read(1)
    no_data_FI = rasterio.open(path_FI).meta.get('nodata')
    SOCref_layer = rasterio.open(path_SOCREF).read(1)
    no_data_SOCref = rasterio.open(path_SOCREF).meta.get('nodata')
    FLU_layer = rasterio.open(path_FLU).read(1)
    no_data_FLU = rasterio.open(path_FLU).meta.get('nodata')


    out_SOC_transform = rasterio.open(path_FMG).transform
    out_SOC_meta = rasterio.open(path_FMG).meta
    spatial_resolution = out_SOC_transform[0]

    ### Now apply the formula for SOC (2.25 IPCC) (ton C/ha)
    logger.info('CALCULATE SOCREF')

    SOC_raster = np.full(SOCref_layer.shape,65535).astype('uint16')

    ### the area of each pixel in hectares
    A_pixel_ha = int(((spatial_resolution**2) / 10000))


    SOC = SOCref_layer * (FLU_layer/scaling) * (FMG_layer/scaling)*(FI_layer/scaling) * A_pixel_ha

    ### ensure that the nodata pixels in one of the layers are set back to no data
    SOC[SOCref_layer == no_data_SOCref] = 65535
    SOC[FMG_layer == no_data_FMG] = 65535
    SOC[FI_layer == no_data_FI] = 65535
    SOC[FLU_layer == no_data_FLU] = 65535
    SOC = SOC.astype('uint16')

    out_SOC_meta.update({'nodata': 65535,
                        'dtype': 'uint16'})

    write_raster(np.expand_dims(SOC,0),out_SOC_meta,out_SOC_transform,Path(outfolder)
                 ,Path(outname_SOC_scenario).name)





def write_raster(out_image,out_meta, out_transform, outdir, outname):
    out_meta.update({"driver": "GTiff",
                     "height": out_image.shape[1],
                     "width": out_image.shape[2],
                     "transform": out_transform,
                     'compress': 'lzw'})
    if out_meta.get('count') == 1:
        with rasterio.open(os.path.join(outdir, outname), "w", **out_meta) as dest:
            dest.write(out_image)
    else:
        with rasterio.open(os.path.join(outdir, outname), "w", **out_meta) as dest:
            for band in range(out_image.shape[0]):
                band_result = out_image[band,:,:]
                dest.write_band(band+1, band_result)

def downsize_raster(input_raster_file, overwrite):
    outfile = Path(input_raster_file).parent.joinpath(Path(input_raster_file).stem + '_downsized.tif').as_posix()

    if not os.path.exists(outfile) and not overwrite:

        with rasterio.open(input_raster_file, 'r') as file:
            ref_raster = file.read(1)
            out_meta = file.meta
            out_transform = file.transform
            nodata_value = file.nodata

        downsized_raster = np.full(ref_raster.shape,255).astype(out_meta.get('dtype'))
        downsized_raster[ref_raster != nodata_value] = ref_raster[ref_raster != nodata_value]

        out_meta.update({'nodata': 255})
        write_raster(np.expand_dims(downsized_raster,0),out_meta,out_transform,Path(outfile).parent, Path(outfile).name)
        os.unlink(input_raster_file)
        ## give again the original name to the file
        os.rename(outfile, input_raster_file)


def resample_raster_to_ref(raster_files_to_clip, raster_file_ref, MS_focus, outdir, overwrite,
                        outname, resample_method = 'near', resample_factor = 1, scale_factor = 1, resampling = False,
                           ):

    for raster_file in raster_files_to_clip:


        if os.path.exists(os.path.join(outdir, outname)) and not overwrite:
            continue

        #remove if file exists and overwrite is true, otherwise gdal will crash
        if os.path.exists(os.path.join(outdir, outname)):
            os.unlink(os.path.join(outdir, outname))

        ds_reference = gdal.Open(raster_file_ref)
        ds_transform = ds_reference.GetGeoTransform()
        proj = osr.SpatialReference(wkt =ds_reference.GetProjection())
        epsg = proj.GetAttrValue('AUTHORITY',1)
        upper_left_corner_x = ds_transform[0]
        upper_left_corner_y = ds_transform[3]
        lower_right_corner_x = upper_left_corner_x + (ds_reference.RasterXSize * ds_transform[1])
        lower_right_corner_y = upper_left_corner_y + (ds_reference.RasterYSize * ds_transform[-1])
        if resampling:
            cmd_warp = f'gdalwarp --config GDAL_CACHEMAX 256 -co COMPRESS=LZW -s_srs epsg:{epsg} -t_srs epsg:{epsg} -te {str(upper_left_corner_x)} {str(lower_right_corner_y)} {str(lower_right_corner_x)} {str(upper_left_corner_y)}' \
                       f' -tr {ds_transform[1]*resample_factor} {abs(ds_transform[-1])*resample_factor} -of GTiff {raster_file} {os.path.join(outdir, outname)} -r {resample_method}'
        else:
            cmd_warp = f'gdalwarp -s_srs epsg:{epsg} -t_srs epsg:{epsg} -te {str(upper_left_corner_x)} {str(lower_right_corner_y)} {str(lower_right_corner_x)} {str(upper_left_corner_y)}' \
                       f'  -ts -of GTiff {raster_file} {os.path.join(outdir, outname)} -r {resample_method}'

        try:
            subprocess.check_call(cmd_warp, shell = True, stderr= subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            raise RuntimeError ('command {} return with error code {}: {}'.format(e.cmd, e.returncode, e.output))

        if scale_factor != 1:
            outname_tmp = rescale_raster(os.path.join(outdir, outname),outdir, rescale_factor= scale_factor)
            os.unlink(os.path.join(outdir, outname))
            os.rename(outname_tmp, os.path.join(outdir, outname))

        ### reduce the size of the raster because gdal is not that efficient
        downsize_raster(os.path.join(outdir, outname), overwrite)






def downsample_raster(raster_files_to_downsample, dir_ref_spatial,mixed_pixel_threshold, name_LC_info,
                      GDMP_300M = False,overwrite = False, outdir_downsampled = None):

    with rasterio.open(dir_ref_spatial, 'r') as src_spatial:
        nodata_ref_spatial = src_spatial.nodata
        values_ref_spatial = src_spatial.read(1)
        meta_ref_spatial = src_spatial.meta
        transform_ref_spatial = src_spatial.transform

    for raster_file in raster_files_to_downsample:
        if GDMP_300M:
            outname_downsampled = Path(raster_file).stem + '_downsampled_300m.tif'
        else:
            outname_downsampled = Path(raster_file).stem + '_downsampled_1000m.tif'

        if outdir_downsampled is None:
            outdir_downsampled = Path(raster_file).parent.parent.joinpath(f'{name_LC_info}_downsampled')
        outdir_downsampled.mkdir(parents=True, exist_ok=True)

        if os.path.exists(os.path.join(outdir_downsampled, outname_downsampled)) and not overwrite:
            continue

        with rasterio.open(raster_file, 'r') as src:
            nodata_raster_downsample = src.nodata
            raster_high_res_values = src.read(1)
            meta_raster_downsample = src.meta

        downsampling_factor = src_spatial.transform[0]/src.transform[0]

        raster_downsampled = np.full(values_ref_spatial.shape, nodata_raster_downsample)

        iteration = 0
        for (x,y), value in np.ndenumerate(values_ref_spatial):
            iteration += 1
            #check empty GDMP location (highly certain also no CLC available at this location
            if value == nodata_ref_spatial:
                continue

            values_downsampled_pixel =  raster_high_res_values[int(np.round(x*downsampling_factor)):int(np.round((x+1)*downsampling_factor))
            ,int(np.round(y*downsampling_factor)): int(np.round((y+1)*downsampling_factor))]

            values_cell, counts = np.unique(values_downsampled_pixel
                                            , return_counts= True)
            print(f'DOWNSAMPLED {str((iteration/values_ref_spatial.size)*100)}%')
            if values_cell[np.argmax(counts)] == nodata_raster_downsample:
                continue
            else:
                mode_values_downsampled_pixel = values_cell[np.argmax(counts)]
                count_mode_downsampled_pixel = counts[np.argmax(counts)]

                #check if not too much mixed pixels within cell
                if (count_mode_downsampled_pixel/values_downsampled_pixel.size) >= (1-mixed_pixel_threshold):
                    raster_downsampled[x,y] = mode_values_downsampled_pixel

        raster_downsampled = raster_downsampled.astype(meta_raster_downsample.get('dtype'))
        meta_ref_spatial.update({'dtype': meta_raster_downsample.get('dtype'),
                                 'nodata': nodata_raster_downsample})
        write_raster(np.expand_dims(raster_downsampled,0), meta_ref_spatial, transform_ref_spatial
                     , outdir_downsampled, outname_downsampled)

def rescale_raster(raster_file,outdir, rescale_factor):
    outname = Path(raster_file).stem + '_tmp.tif'

    with rasterio.open(raster_file, 'r') as raster_object:
        out_image = raster_object.read(1)
        out_meta = raster_object.meta
        out_transform = raster_object.transform
        out_image[out_image != raster_object.nodata] = out_image[out_image != raster_object.nodata] * rescale_factor
        write_raster(np.expand_dims(out_image,0), out_meta, out_transform, outdir, outname)

    return os.path.join(outdir, outname)
































