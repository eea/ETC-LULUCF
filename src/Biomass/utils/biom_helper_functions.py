"""
In this script all the functions needed for the estimation of the biomass carbon pool will be stored
"""

### Import needed modules
import os
from loguru import logger
from SOC_scenarios.utils.soc_helper_functions import (
    get_factors_from_NUTS,
    open_raster_from_window,
    resample_raster_to_ref,
    write_raster,
    mask_raster_extent)

from pathlib import Path
import numpy as np
import pandas as pd
import geopandas as gpd

def define_affor_areas(settings, slope_max = 87.5):
    """
    Function that will help in finding the potential afforestation areas based on CLC input, DEM (slope), N2K &
    external mask
    :param settings: different parameters used to create the afforestation mask
    :param slope_max: the maximum slope under which afforestation can be applied
    :return: Map with a pixel-based indication where afforestation could be applied
    """

    ## load all the required rasters and check if they are present

    ## CLC_ACC raster
    file_IPCC_LU = settings.get('file_IPCC_LU')

    ## DEM
    DEM_dir = settings.get('dict_afforestation_masking_layers').get('DEM')
    if not os.path.exists(DEM_dir):
        raise Exception

    ## N2K
    N2K_dir = settings.get('dict_afforestation_masking_layers').get('N2K')
    if not os.path.exists(N2K_dir):
        N2K_dir_coarse_res = Path(N2K_dir).parent.joinpath(Path(N2K_dir).name.replace('_100m', ''))
        if os.path.exists(N2K_dir_coarse_res):
            ## try generating the 100m spatial resolution layer
            resample_raster_to_ref([N2K_dir_coarse_res], file_IPCC_LU, 'EEA39',
                                   Path(N2K_dir_coarse_res).parent.as_posix(),
                                   resample_factor=1, overwrite=settings.get('overwrite'),
                                   resampling=True,
                                   outname=Path(N2K_dir).name)
        else:
            raise Exception

    ## External mask (if some areas or not suitable for afforestation). 1 --> not suitable. 0 --> suitable
    External_mask_dir = settings.get('dict_afforestation_masking_layers').get('external_mask')

    if External_mask_dir is not None:
        if not os.path.exists(External_mask_dir):
            raise Exception

    ## Start the masking
    ## load the scenario settings
    factor_scenario = settings.get('afforestation_scenario')


    overwrite = False
    ### get now a dict with the scenario for the NUTS ID and replace it with the default scenario

    if settings.get('run_NUTS_specific_scenario'):  #otherwise use default factors
        factor_scenario = get_factors_from_NUTS(settings, factor_scenario, 'Slope'
                                                , id_LUT_carbon_pool='afforestation')

    if not settings.get('NUTS3_info').geometry is None:
        bounds = settings.get('NUTS3_info').geometry.bounds

    else:
        ## ignore empty regions
        logger.info(f'NO GEOMETRY PROVIDED FOR AREA {settings.get("NUTS3_info").NUTS_ID}')
        return None

    Basefolder_mask = Path(settings.get('Basefolder_output')).joinpath('Afforestation_mask')
    if settings.get('Country') is None and not settings.get('block_based_processing'):
        outdir_mask = Basefolder_mask
        outname_mask = f'Afforestation_mask_{str(settings.get("year_baseline"))}_EEA39.tif'
    else:
        if not settings.get('block_based_processing'):
            outdir_mask = Path(Basefolder_mask).joinpath(settings.get('Country'))
            outname_mask = f'Afforestation_mask_{str(settings.get("year_baseline"))}_{settings.get("Country")}.tif'
        else:
            outdir_mask = Path(Basefolder_mask).joinpath(settings.get('NUTS3_info')['CNTR_CODE'])
            outname_mask = f'Afforestation_mask_{str(settings.get("year_baseline"))}_' \
                           f'{settings.get("NUTS3_info")["NUTS_ID"]}.tif'

    outdir_mask.mkdir(parents=True, exist_ok=True)

    if not os.path.exists(outdir_mask.joinpath(outname_mask)) or overwrite:
        # TODO need to add the option if not block-based processing
        if settings.get('block_based_processing'):
            ### block based opening of raster(s)
            CLC_IPCC_LUCAT_raster, src_CLC_raster = open_raster_from_window(file_IPCC_LU, bounds)
            if CLC_IPCC_LUCAT_raster.size == 0:
                logger.info(f'NO PIXEL DATA AVAILABLE FOR NUTS REGION: {settings.get("NUTS3_info").NUTS_ID}')
                return None

            no_data_IPCC_LUCAT = src_CLC_raster.get('nodata')

            DEM_raster, src_DEM_raster = open_raster_from_window(DEM_dir, bounds)
            if DEM_raster.size == 0:
                logger.info(f'NO PIXEL DATA AVAILABLE FOR NUTS REGION: {settings.get("NUTS3_info").NUTS_ID}')
                return None
            no_data_DEM_raster = src_DEM_raster.get('nodata')


            N2K_raster, src_N2K_raster = open_raster_from_window(N2K_dir, bounds)
            if N2K_raster.size == 0:
                logger.info(f'NO PIXEL DATA AVAILABLE FOR NUTS REGION: {settings.get("NUTS3_info").NUTS_ID}')
                return None
            if N2K_raster.size != DEM_raster.size:
                logger.warning(f'N2K RASTER IS NOT FULLY AVALABLE AT NUTS REGION: {settings.get("NUTS3_info").NUTS_ID}')
                return None

            if External_mask_dir is not None:
                External_mask_raster, src_external_mask = open_raster_from_window(External_mask_dir, bounds)
                if External_mask_raster.size == 0:
                    logger.info(f'NO PIXEL DATA AVAILABLE FOR NUTS REGION: {settings.get("NUTS3_info").NUTS_ID}')
                    return None
                no_data_external_mask = src_external_mask.get('nodata')

            else:
                External_mask_raster = np.zeros(DEM_raster.shape)
                no_data_external_mask = 255


        ## store the info needed for creating output mask layer
        meta_raster = src_CLC_raster
        transform_raster = src_CLC_raster.get('transform')
        meta_raster.update({'nodata': 255, 'dtype': 'uint8'})

        ### Create factor raster to write out factors in it
        mask_raster = np.full(CLC_IPCC_LUCAT_raster.shape,255).astype('uint8')

        for Crop in factor_scenario.keys():
            ## load the way the CLC layer is translated to IPCC LUCAT
            df_factor_IPCC_factors = pd.read_csv(os.path.join(settings.get('Basefolder_output'),
                                                 'SOC_LUT',
                                                 f'IPCC_FLU_CLC_mapping_LUT.csv'), sep=';')
            ## get the value assigned to the specific IPCC LUCAT
            value_LUCAT = df_factor_IPCC_factors.loc[df_factor_IPCC_factors['IPCC_landuse_name'] == Crop]\
                                                    ['IPCC_landuse_id'].values[0]
            loc_strata = np.where(CLC_IPCC_LUCAT_raster == value_LUCAT)

            # the minimum slope at which afforestation may start
            slope_min = factor_scenario.get(Crop).get('Slope')

            # Now filter on the locations with a suitable slope
            loc_affor_option = np.where((DEM_raster > slope_min) & (DEM_raster < slope_max)
                                 & (CLC_IPCC_LUCAT_raster == value_LUCAT) & (N2K_raster == 255) &
                                        (External_mask_raster == 0))

            ## set all the LU pixels to no masking initially
            mask_raster[loc_strata] = 0

            ## set the ones that could be afforested to one
            mask_raster[loc_affor_option] = 1


        ### ensure that the nodata pixels in one of the layers are set back to no data
        mask_raster[CLC_IPCC_LUCAT_raster == no_data_IPCC_LUCAT] = 255
        mask_raster[DEM_raster == no_data_DEM_raster] = 255
        mask_raster[External_mask_raster == no_data_external_mask] = 255

        write_raster(np.expand_dims(mask_raster,0),meta_raster,transform_raster,Path(outdir_mask)
                     , outname_mask)


        ### apply a cropping of the raster so that cells outside the extent are set to nodata
        if settings.get('block_based_processing'):
            NUTS_gpd = gpd.GeoDataFrame(geometry=[settings.get('NUTS3_info').geometry])
            outdir_raster = Path(outdir_mask).joinpath(outname_mask).as_posix()
            mask_raster_extent([Path(outdir_mask).joinpath(outname_mask).as_posix()],NUTS_gpd,
                               Path(outdir_mask), settings, overwrite= overwrite)
            ## remove the unclipped file and the rename the clipped file to the final output name
            os.remove(outdir_raster)
            os.rename(Path(outdir_mask).joinpath(outname_mask.replace('.tif', '_clipped.tif')).as_posix(),
                      outdir_raster)

            mask_raster_clipped, _ = open_raster_from_window(outdir_raster, bounds)

            return mask_raster_clipped
        else:
            return mask_raster
    else:
        mask_raster, _ = open_raster_from_window(outdir_mask.joinpath(outname_mask), bounds)
        return mask_raster


def create_affor_potential(settings):
    """
    Function that will calculate the afforestation carbon potential and rasterize the output
    :param settings: Different parameters used to create the afforestation potential map
    :return: raster with the carbon potential if a certain pixel is afforested
    """


    ### First load the needed layers
    Basefolder_mask = Path(settings.get('Basefolder_output')).joinpath('Afforestation_mask')
    Basefolder_affor_potential = Path(settings.get('Basefolder_output')).joinpath('Afforestation_pot')
    ## CLC_ACC raster
    file_IPCC_LU = settings.get('file_IPCC_LU')

    ### Define the location of the afforestation mask layer and the output of the afforestation potential
    if settings.get('Country') is None and not settings.get('block_based_processing'):
        outdir_mask = Basefolder_mask
        outname_mask = f'Afforestation_mask_{str(settings.get("year_baseline"))}_EEA39.tif'
        outdir_affor_pot = Basefolder_affor_potential
        outname_affor_pot = f'Afforestation_pot_{str(settings.get("year_potential"))}_EEA39.tif'

    else:
        if not settings.get('block_based_processing'):
            outdir_mask = Path(Basefolder_mask).joinpath(settings.get('Country'))
            outname_mask = f'Afforestation_mask_{str(settings.get("year_baseline"))}_{settings.get("Country")}.tif'

            outdir_affor_pot = Path(Basefolder_affor_potential).joinpath(settings.get('Country'))
            outname_affor_pot = f'Afforestation_pot_{str(settings.get("year_potential"))}_{settings.get("Country")}.tif'
        else:
            outdir_mask = Path(Basefolder_mask).joinpath(settings.get('NUTS3_info')['CNTR_CODE'])
            outname_mask = f'Afforestation_mask_{str(settings.get("year_baseline"))}_' \
                           f'{settings.get("NUTS3_info")["NUTS_ID"]}.tif'

            outdir_affor_pot = Path(Basefolder_affor_potential).joinpath(settings.get('NUTS3_info')['CNTR_CODE'])
            outname_affor_pot = f'Afforestation_pot_{str(settings.get("year_potential"))}_{settings.get("NUTS3_info")["NUTS_ID"]}.tif'

    afforestation_mask_dir = os.path.join(outdir_mask, outname_mask)
    if not os.path.exists(afforestation_mask_dir):
        raise Exception

    os.makedirs(outdir_affor_pot,exist_ok=True)
    affor_pot_dir = os.path.join(outdir_affor_pot, outname_affor_pot)

    ## Load the specific defined parameters (at EU or NUTS scale) used to calculate the potential
    factor_scenario = settings.get('afforestation_scenario')

    overwrite = False

    if not settings.get('NUTS3_info').geometry is None:
        bounds = settings.get('NUTS3_info').geometry.bounds
    else:
        ## ignore empty regions
        logger.info(f'NO GEOMETRY PROVIDED FOR AREA {settings.get("NUTS3_info").NUTS_ID}')
        return None



    if not os.path.exists(affor_pot_dir) or overwrite:
        if settings.get('block_based_processing'):
            ### block based opening of raster(s)
            affor_mask_raster, src_affor_mask_raster = open_raster_from_window(afforestation_mask_dir, bounds)
            if affor_mask_raster.size == 0:
                logger.info(f'NO PIXEL DATA AVAILABLE FOR NUTS REGION: {settings.get("NUTS3_info").NUTS_ID}')
                return None

            CLC_IPCC_LUCAT_raster, src_CLC_raster = open_raster_from_window(file_IPCC_LU, bounds)
            if CLC_IPCC_LUCAT_raster.size == 0:
                logger.info(f'NO PIXEL DATA AVAILABLE FOR NUTS REGION: {settings.get("NUTS3_info").NUTS_ID}')
                return None

            ## Load the LUT with info on the EU4trees and the IPCC increment
            df_trees_biom_increment = pd.read_csv(os.path.join(settings.get('Basefolder_output'),
                                                               'NUTS_LUT_afforestation_scenario',
                                                               f'EU4_trees_LUT_biom_increment.csv'), sep=';')

            # create raster to write output in
            affor_pot_raster = np.full(affor_mask_raster.shape,65535).astype('uint16')


            for Crop in factor_scenario.keys():
                ### load the specific (flexible) parameters per IPCC LUCAT for afforestation
                if settings.get('run_NUTS_specific_scenario'):  #otherwise use default factors
                    factor_scenario_tree_species = get_factors_from_NUTS(settings, factor_scenario, 'Tree_species'
                                                            , id_LUT_carbon_pool='afforestation')
                    tree_species = factor_scenario_tree_species.get(Crop).get('Tree_species')

                    factor_scenario_tree_prob = get_factors_from_NUTS(settings, factor_scenario, 'Tree_prob'
                                                                      , id_LUT_carbon_pool='afforestation')
                    tree_prob = factor_scenario_tree_prob.get(Crop).get('Tree_prob')

                else:
                    tree_species = factor_scenario.get(Crop).get('Tree_species')
                    tree_prob = factor_scenario.get(Crop).get('Tree_prob')


                RCP_scenario = factor_scenario.get(Crop).get('RCP')


                ## Now that we know the tree species and the probability we can start looking for the available data


                ## First check if we have the LUT tables to calculate the biomass increment for the requested tree species
                if tree_species in df_trees_biom_increment.Tree_species.to_list():
                    ## load the way the CLC layer is translated to IPCC LUCAT
                    df_factor_IPCC_factors = pd.read_csv(os.path.join(settings.get('Basefolder_output'),
                                                                      'SOC_LUT',
                                                                      f'IPCC_FLU_CLC_mapping_LUT.csv'), sep=';')
                    ## get the value assigned to the specific IPCC LUCAT
                    value_LUCAT = df_factor_IPCC_factors.loc[df_factor_IPCC_factors['IPCC_landuse_name'] == Crop] \
                        ['IPCC_landuse_id'].values[0]
                    loc_strata = np.where(CLC_IPCC_LUCAT_raster == value_LUCAT)































