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
import glob
import osr
from osgeo import gdal
import rasterio
import datetime
from rasterstats import zonal_stats

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

    ## External mask (if some areas or not suitable for afforestation). 1 --> suitable. 0 --> not suitable
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

        for LUCAT in factor_scenario.keys():
            ## load the way the CLC layer is translated to IPCC LUCAT
            df_factor_IPCC_factors = pd.read_csv(os.path.join(settings.get('Basefolder_output'),
                                                 'SOC_LUT',
                                                 f'IPCC_FLU_CLC_mapping_LUT.csv'), sep=';')
            ## get the value assigned to the specific IPCC LUCAT
            value_LUCAT = df_factor_IPCC_factors.loc[df_factor_IPCC_factors['IPCC_landuse_name'] == LUCAT]\
                                                    ['IPCC_landuse_id'].values[0]
            loc_strata = np.where(CLC_IPCC_LUCAT_raster == value_LUCAT)

            # the minimum slope at which afforestation may start
            slope_min = factor_scenario.get(LUCAT).get('Slope')

            # Now filter on the locations with a suitable slope
            loc_affor_option = np.where((DEM_raster > slope_min) & (DEM_raster < slope_max)
                                 & (CLC_IPCC_LUCAT_raster == value_LUCAT) & (N2K_raster == 255))

            # the external mask will add some additional areas even if not suitable
            mask_raster[External_mask_raster == 1] = 1  # 1 is suitable

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


def create_affor_potential(settings, affor_mask_array):
    """
    Function that will calculate the afforestation carbon potential and rasterize the output
    :param settings: Different parameters used to create the afforestation potential map
    :param affor_mask_array: the array containing the suitable locations for afforestation
    :return: raster with the carbon potential if a certain pixel is afforested
    """

    carbon_pool = settings.get('carbon_pool')

    ### First load the needed layers
    Basefolder_affor_potential = Path(settings.get('Basefolder_output')).joinpath(f'{carbon_pool}_scenario')
    Basefolder_affor_potential.mkdir(exist_ok=True)
    Basefolder_input_data = settings.get('Basefolder_input_data')
    Files_EUtrees4F_prob = glob.glob(os.path.join(Basefolder_input_data, 'EU-Trees4F', 'ens_sdms', 'prob', 'resampled', '*tif'))
    ## CLC_ACC raster
    file_IPCC_LU = settings.get('file_IPCC_LU')


    ### Define the location of the afforestation mask layer and the output of the afforestation potential
    if settings.get('Country') is None and not settings.get('block_based_processing'):
        outdir_affor_pot = Basefolder_affor_potential
        outname_affor_pot = f'Afforestation_pot_{str(settings.get("year_potential"))}_EEA39.tif'

    else:
        if not settings.get('block_based_processing'):
            outdir_affor_pot = Path(Basefolder_affor_potential).joinpath('EE3A9')
            outname_affor_pot = f'{carbon_pool}_{settings.get("Scenario_name")}_{"EEA39"}.tif'
        else:
            outdir_affor_pot = Basefolder_affor_potential \
                .joinpath(settings.get('NUTS3_info')['CNTR_CODE']).as_posix()
            outname_affor_pot = f'{carbon_pool}_{settings.get("Scenario_name")}' \
                                f'_{settings.get("NUTS3_info")["NUTS_ID"]}.tif'


    os.makedirs(outdir_affor_pot,exist_ok=True)
    affor_pot_dir = os.path.join(outdir_affor_pot, outname_affor_pot)


    ## the years for which a tree potential is available
    years_EU4trees_pred = [2035,2065,2095]

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
        #TODO write scaling factor in metadata of raster
        scaling = settings.get('Scaling') ## just to avoid the saving in floats of the raster

        if settings.get('block_based_processing'):
            ### block based opening of raster(s)
            if affor_mask_array.size == 0:
                logger.info(f'NO PIXEL DATA AVAILABLE FOR NUTS REGION: {settings.get("NUTS3_info").NUTS_ID}')
                return None

            CLC_IPCC_LUCAT_raster, src_CLC_raster = open_raster_from_window(file_IPCC_LU, bounds)
            if CLC_IPCC_LUCAT_raster.size == 0:
                logger.info(f'NO PIXEL DATA AVAILABLE FOR NUTS REGION: {settings.get("NUTS3_info").NUTS_ID}')
                return None
            no_data_IPCC_LUCAT = src_CLC_raster.get('nodata')

            ## Load the LUT with info on the EU4trees and the IPCC increment
            df_trees_biom_increment = pd.read_csv(os.path.join(settings.get('Basefolder_output'),
                                                               'NUTS_LUT_afforestation_scenario',
                                                               f'EU4_trees_LUT_biom_increment.csv'), sep=';')


            ## store the info needed for creating output mask layer
            meta_raster = src_CLC_raster
            transform_raster = src_CLC_raster.get('transform')
            meta_raster.update({'nodata': 65535, 'dtype': 'uint16'})
            spatial_resolution = transform_raster[0]
            ### the area of each pixel in hectares
            A_pixel_ha = int(((spatial_resolution**2) / 10000))

            # create raster to write output in
            affor_yrly_pot_raster = np.full(affor_mask_array.shape, 65535).astype('uint16')


            for LUCAT in factor_scenario.keys():
                ### load the specific (flexible) parameters per IPCC LUCAT for afforestation
                if settings.get('run_NUTS_specific_scenario'):  #otherwise use default factors
                    factor_scenario_tree_species = get_factors_from_NUTS(settings, factor_scenario, 'Tree_species'
                                                            , id_LUT_carbon_pool='afforestation')
                    tree_species = factor_scenario_tree_species.get(LUCAT).get('Tree_species')

                    factor_scenario_tree_prob = get_factors_from_NUTS(settings, factor_scenario, 'Tree_prob'
                                                                      , id_LUT_carbon_pool='afforestation')
                    tree_prob = factor_scenario_tree_prob.get(LUCAT).get('Tree_prob')

                else:
                    tree_species = factor_scenario.get(LUCAT).get('Tree_species')
                    tree_prob = factor_scenario.get(LUCAT).get('Tree_prob')


                ## These could not be NUTS specific factors
                RCP_scenario = factor_scenario.get(LUCAT).get('RCP')
                Year_potential = factor_scenario.get(LUCAT).get('Year_potential')

                ## now we know the year for the potential, we can select the most closest EU4trees prob year
                Diff_yrs_prob_predictions = [abs(item - Year_potential) for item in years_EU4trees_pred]
                check_duplicates = set([x for x in Diff_yrs_prob_predictions if Diff_yrs_prob_predictions.count(x) > 1])

                ## if multiple years are of the same distance from a probability layer take the probability layer in the future
                if len(check_duplicates) > 1:
                    idx_year_select = max([i for i, e in enumerate(Diff_yrs_prob_predictions) if e == list(check_duplicates)[0]])
                else:
                    idx_year_select = np.argmin(Diff_yrs_prob_predictions)

                ## now we can define the year for which the probability layer should be loaded
                Year_select_EU4trees = years_EU4trees_pred[idx_year_select]

                ## now select the suitable potential layer based on the potential year

                ## Open the occurence probability of the selected tree species for the selected RCP and year


                ## First check if we have the LUT tables to calculate the biomass increment for the requested tree species
                if tree_species in df_trees_biom_increment.Tree_species.to_list():
                    ## load the probability layer of the tree species
                    EU4Trees_prob_select = [item for item in Files_EUtrees4F_prob if
                                            Path(item).stem == f'{tree_species}_ens-sdms_{RCP_scenario}_fut' \
                                                               f'{str(Year_select_EU4trees)}_prob_pot_' \
                                                               f'EPSG3035']
                    if len(EU4Trees_prob_select) != 1:
                        raise Exception(f'Problem with loading tree species prob layer for '
                                        f'{settings.get("NUTS3_info").NUTS_ID}')

                    Tree_prob_raster, _ = open_raster_from_window(EU4Trees_prob_select[0], bounds)

                    ## load the way the CLC layer is translated to IPCC LUCAT
                    df_factor_IPCC_factors = pd.read_csv(os.path.join(settings.get('Basefolder_output'),
                                                                      'SOC_LUT',
                                                                      f'IPCC_FLU_CLC_mapping_LUT.csv'), sep=';')
                    ann_increment_tree_species = df_trees_biom_increment.loc[df_trees_biom_increment
                                                                                 .Tree_species == tree_species]['Ann_increment'].values[0]
                    ## get the value assigned to the specific IPCC LUCAT
                    value_LUCAT = df_factor_IPCC_factors.loc[df_factor_IPCC_factors['IPCC_landuse_name'] == LUCAT] \
                        ['IPCC_landuse_id'].values[0]

                    #now find the location on which can be afforested (based on the mask) and the LUCAT is present
                    loc_affor = np.where((affor_mask_array == 1) & (CLC_IPCC_LUCAT_raster == value_LUCAT) &
                                         (Tree_prob_raster > (tree_prob * 10)))

                    ## The actual potential calculation can now be computed based on a formula as follows:
                    """
                    annual CO2 sink = Volume increment x Density x BEFx Rootoshoot x CF x kg_to_ton
                    Where:
                    density = 500
                    BEF= 1.2
                    R= 1.25
                    CF= 0.47
                    kg_to_ton = 1000
                    """
                    #divide by 1000 to convert to tonC/hayr
                    affor_yrly_pot_raster[loc_affor] = int(ann_increment_tree_species * 500 * 1.2 * 0.25 *
                                                           0.47*scaling * A_pixel_ha * 0.001)


            ### ensure that the nodata pixels in one of the layers are set back to no data
            affor_yrly_pot_raster[CLC_IPCC_LUCAT_raster == no_data_IPCC_LUCAT] = 65535

            write_raster(np.expand_dims(affor_yrly_pot_raster,0),meta_raster,transform_raster,Path(outdir_affor_pot)
                         , outname_affor_pot)


            ### apply a cropping of the raster so that cells outside the extent are set to nodata
            if settings.get('block_based_processing'):
                NUTS_gpd = gpd.GeoDataFrame(geometry=[settings.get('NUTS3_info').geometry])
                outdir_raster = Path(outdir_affor_pot).joinpath(outname_affor_pot).as_posix()
                mask_raster_extent([Path(outdir_affor_pot).joinpath(outname_affor_pot).as_posix()],NUTS_gpd,
                                   Path(outdir_affor_pot), settings, overwrite= overwrite)
                ## remove the unclipped file and the rename the clipped file to the final output name
                os.remove(outdir_raster)
                os.rename(Path(outdir_affor_pot).joinpath(outname_affor_pot.replace('.tif', '_clipped.tif')).as_posix(),
                          outdir_raster)

                affor_yrly_pot_raster_clipped, _ = open_raster_from_window(outdir_raster, bounds)

                return affor_yrly_pot_raster_clipped, outname_affor_pot
            else:
                return affor_yrly_pot_raster, outname_affor_pot


    else:
        affor_yrly_pot_raster, _ = open_raster_from_window(outdir_affor_pot.joinpath(outname_affor_pot), bounds)
        return affor_yrly_pot_raster, outname_affor_pot



def align_raster_data_reference(input_raster: str, reference_raster: str, BBOX_mask: gpd.GeoDataFrame, settings: dict,
                                suffix_raster: str = 'EEA39'):
    """
    Function that will mask (based on BBOX) and resample a raster based on a reference raster.
    :param input_raster: the input directory of the raster that should be resampled
    :param reference_raster: the input directory of the reference raster that will be used to resample the input raster:
    :param settings: dictionary with settins that are helpfull for the resampling of the raster
    :param geometry that will be used to mask the raster:
    :param suffix_raster: the suffix that should be added to the output raster
    :return:
    """
    # define the outout directory of the resampled raster
    Path(input_raster).parent.joinpath('resampled').mkdir(exist_ok=True)

    resample_raster_to_ref([input_raster], reference_raster, suffix_raster,
                           Path(input_raster).parent.joinpath('resampled'),
                           resample_factor=1, overwrite=False,
                           resampling=True,
                           outname=Path(input_raster).stem + '_EPSG3035.tif',
                           resample_method='near', dtype='UInt16')


def reproject_raster(input_raster: str,output_raster: str, src_nodata: int, dst_nodata: int,
                     crs_to: int = 3035, dtype_out: str = 'uint16', resample_method: str = 'bilinear'):
    """
   Function that will reproject a raster to a certain projection system.
   :param input_raster: the input directory of the raster that should be resampled
   :param reference_raster: the output directory of the reprojected raster
   :param crs_to: the projection system to reproject to
   :param dtype_out: the output type that should be created. Please note that now only uint16 is supported
   :param resample_method: the method that should be used to resampling
   :return:
   """
    input_raster_obj = gdal.Open(input_raster)
    output_raster = output_raster
    warp = gdal.Warp(output_raster,input_raster_obj, dstSRS=f'EPSG:{str(crs_to)}',
                     outputType=gdal.GDT_UInt16, resampleAlg=resample_method,
                     srcNodata = src_nodata, dstNodata = dst_nodata)
    warp = None # Closes the files


### add metadata to created raster
def create_metadata_description_afforestation(settings: dict, extent: str = 'EEA39') -> dict:
    """
    Function that will create a description that will be saved in the final SOC scenario layer and also add band metadata
    :param settings:dictionary that contains some settings used for writing out the metadata
    :param extent: define if the description should be created for 'EEA39' scale or 'NUTS' scale
    :return:description that will be written into the raster
    """


    dict_general ={'copyright':'DO NOT DISTRIBUTE',
                   'commit_id': settings.get('commit_id'),
                   'doi':'Prototype dataset, no registered DOI',
                   'institution':'EEA',
                   'name':settings.get('Scenario_name') + '_afforestation_yrl_increment',
                   'processing_date':datetime.datetime.now().date(),
                   'source':' Derived from IPCC LUT and EU4TREES projection',
                   'RCP': settings.get('afforestation_scenario').get(list(settings.get('afforestation_scenario')
                                                                          .keys())[0]).get('RCP'),
                   'Year_potential': settings.get('afforestation_scenario').get(list(settings.get('afforestation_scenario')
                                                                                     .keys())[0]).get('Year_potential'),
                   'title':settings.get('Scenario_name')+' scenario at 100 m'}


    if settings.get('run_NUTS_specific_scenario') and extent == 'EEA39':
        dict_general.update({'Slope': f'Slope afforestation threshold value {settings.get("afforestation_scenario").get(list(settings.get("afforestation_scenario").keys())[0]).get("Slope")}% & FROM {"EEA39"}',
                             'Tree prob': f'Tree occurence probability value {settings.get("afforestation_scenario").get(list(settings.get("afforestation_scenario").keys())[0]).get("Tree_prob")}% & FROM {"EEA39"}',
                             'Tree species': f'Tree species selected {settings.get("afforestation_scenario").get(list(settings.get("afforestation_scenario").keys())[0]).get("Tree_species")} & FROM {"EEA39"}',})



    elif settings.get('run_NUTS_specific_scenario') and extent == 'NUTS':
        dict_Slope_factors_info = get_factors_from_NUTS(settings, settings.get('afforestation_scenario'), 'Slope',
                                                        id_LUT_carbon_pool= 'afforestation')
        dict_Tree_prob_factors_info = get_factors_from_NUTS(settings, settings.get('afforestation_scenario'), 'Tree_prob',
                                                            id_LUT_carbon_pool= 'afforestation')
        dict_Tree_species_selection_info = get_factors_from_NUTS(settings, settings.get('afforestation_scenario'), 'Tree_species',
                                                                 id_LUT_carbon_pool= 'afforestation')

        if 'Cropland' in dict_Slope_factors_info.keys():
            dict_general.update({'Slope Cropland ': 'Slope afforesation threshold value {}% & FROM {}'.format(dict_Slope_factors_info.get("Cropland").get("Slope"),
                                                                                    dict_Slope_factors_info.get("Cropland").get("input_source"))})
        if 'Cropland' in dict_Tree_prob_factors_info.keys():
            dict_general.update({'Tree prob Cropland': 'Tree occurence probability value {}% & FROM {}'.format(dict_Tree_prob_factors_info.get("Cropland").get("Tree_prob"),
                                                                                                     dict_Tree_prob_factors_info.get("Cropland").get("input_source"))})

        if 'Cropland' in dict_Tree_species_selection_info.keys():
            dict_general.update({'Tree species Cropland': 'Tree species selected {} & FROM {}'.format(dict_Tree_species_selection_info.get("Cropland").get("Tree_species"),
                                                                                                     dict_Tree_species_selection_info.get("Cropland").get("input_source"))})

        if 'Grassland' in dict_Slope_factors_info.keys():
            dict_general.update({'Slope Grassland ': 'Slope afforesation threshold value {}% & FROM {}'.format(dict_Slope_factors_info.get("Grassland").get("Slope"),
                                                                                                             dict_Slope_factors_info.get("Grassland").get("input_source"))})
        if 'Grassland' in dict_Tree_prob_factors_info.keys():
            dict_general.update({'Tree prob Grassland': 'Tree occurence probability value {}% & FROM {}'.format(dict_Tree_prob_factors_info.get("Grassland").get("Tree_prob"),
                                                                                                              dict_Tree_prob_factors_info.get("Grassland").get("input_source"))})

        if 'Grassland' in dict_Tree_species_selection_info.keys():
            dict_general.update({'Tree species Grassland': 'Tree species selected {} & FROM {}'.format(dict_Tree_species_selection_info.get("Grassland").get("Tree_species"),
                                                                                                      dict_Tree_species_selection_info.get("Grassland").get("input_source"))})


    else:

        if 'Cropland' in settings.get("afforestation_scenario").keys():
            dict_general.update({'Slope Cropland': 'Slope afforesation threshold value {}% FROM DEFAULT SCENARIO EEA39'
                                .format(settings.get("afforestation_scenario").get("Cropland").get("Slope"))})

            dict_general.update({'Tree prob Cropland': 'Tree occurence probability value {} FROM DEFAULT SCENARIO EEA39'
                                .format(settings.get("afforestation_scenario").get("Cropland").get("Tree_prob"))})
            dict_general.update({'Tree species Cropland': 'Tree species selected {} FROM DEFAULT SCENARIO EEA39'
                                .format(settings.get("afforestation_scenario").get("Cropland").get("Tree_species"))})

        if 'Grassland' in settings.get("afforestation_scenario").keys():
            dict_general.update({'Slope Grassland': 'Slope afforesation threshold value {}% FROM DEFAULT SCENARIO EEA39'
                                .format(settings.get("afforestation_scenario").get("Grassland").get("Slope"))})

            dict_general.update({'Tree prob Grassland': 'Tree occurence probability value {} FROM DEFAULT SCENARIO EEA39'
                                .format(settings.get("afforestation_scenario").get("Grassland").get("Tree_prob"))})
            dict_general.update({'Tree species Grassland': 'Tree species selected {} FROM DEFAULT SCENARIO EEA39'
                                .format(settings.get("afforestation_scenario").get("Grassland").get("Tree_species"))})

    dict_general.update({'MAP_EXTENT':extent})

    dict_band = {
        'nodata': 65535,
        'long_name': 'ABGB (tonnes/C/ha) '
    }

    return dict_general, dict_band



def calc_stats_biomass_NUTS(raster_dir: str, spatial_layer: gpd,
                        settings: dict, stats_type: list = ['mean', 'count'], calc_BGB:bool=True
                            , RTS_ratio:float = 0.25, correct_low_prod_areas:bool = False) -> pd.DataFrame:
    """
    :param raster_dir: the directory of the raster on which the stats should be calculated
    :param spatial_layer: the geospatial layer that will be used to determine the areas for aggregation
    :param settings: dictionary with some settings that are needed to do the aggregation
    :param stats_type: indicate the type of statistics should be derived from the raster
    :param correct_low_prod_areas: if true for areas with lower production (colder climates) a correction
    factor will be applied
    :return: dataframe with the calculated stats
    """

    raster_object = rasterio.open(raster_dir)
    raster_values = raster_object.read(1)
    affine = raster_object.meta.get('transform')
    no_data = raster_object.meta.get('nodata')
    scaling = settings.get('Scaling')

    ## load some settings
    year_baseline = settings.get('year_baseline')
    year_potential = settings.get('year_potential')

    ### Open the FLU layer on which a distinction is made
    ### per IPCC LU type
    if settings.get('Country') is None or settings.get('block_based_processing'):
        outdir_IPCC_LUCAT = Path(settings.get('Basefolder_output')).joinpath('CLC_ACC_IPCC')
    else:
        outdir_IPCC_LUCAT = Path(settings.get('Basefolder_output')).joinpath('CLC_ACC_IPCC').joinpath(settings.get('Country'))

    CLC_IPCC_LUCAT_dir = glob.glob(os.path.join(outdir_IPCC_LUCAT, 'CLC{}ACC*Grassland_Cropland.tif'.format(settings.get("year_baseline"))))[0]


    FLU_raster, meta = open_raster_from_window(CLC_IPCC_LUCAT_dir, spatial_layer.geometry.bounds)

    ## load the way the CLC layer is translated to IPCC LUCAT
    df_FLU_mapping = pd.read_csv(os.path.join(settings.get('Basefolder_output'),
                                                      'SOC_LUT',
                                                      f'IPCC_FLU_CLC_mapping_LUT.csv'), sep=';')
    lst_df_stats_NUTS = []

    ## also take into account the percentage that should be reforest for the total increase calculation
    if settings.get('run_NUTS_specific_scenario'):  #otherwise use default factors
        dict_perc_reforest_info = get_factors_from_NUTS(settings, settings.get('afforestation_scenario')
                                                    , 'Perc_reforest', id_LUT_carbon_pool='afforestation')
    else:
        dict_perc_reforest_info = settings.get('afforestation_scenario')

    # ## Load the specific defined parameters (at EU or NUTS scale) used to calculate the potential
    # factor_scenario = settings.get('afforestation_scenario')


    for FLU_class in df_FLU_mapping['IPCC_landuse_id'].unique():
        raster_values_FLU_filter = np.copy(raster_values)
        raster_values_FLU_filter[np.where(FLU_raster != FLU_class)] = no_data
        stats = zonal_stats(spatial_layer.geometry, raster_values_FLU_filter, affine=affine,
                            nodata=no_data, stats=stats_type)
        df_stats = pd.DataFrame.from_dict(stats)
        df_stats.columns = [f'{settings.get("carbon_pool")}_mean_yrly', 'nr_pixels']
        if correct_low_prod_areas:
            ## only apply this factor for Scandinavia
            if spatial_layer.NUTS_ID[0:2] in ['NO', 'FI', 'SE']:
                df_stats[f'{settings.get("carbon_pool")}_mean_yrly'] = (df_stats[f'{settings.get("carbon_pool")}_mean_yrly']/scaling) * 0.5

            else:
                df_stats[f'{settings.get("carbon_pool")}_mean_yrly'] = df_stats[f'{settings.get("carbon_pool")}_mean_yrly']/scaling

        else:
            ## add the outcome to the proper scale
            df_stats[f'{settings.get("carbon_pool")}_mean_yrly'] = df_stats[f'{settings.get("carbon_pool")}_mean_yrly']/scaling
        df_stats = df_stats.round(2)

        # derive the IPCC_cat for which the stats are calculated
        IPCC_cat = df_FLU_mapping.loc[df_FLU_mapping['IPCC_landuse_id'] == FLU_class] \
            ['IPCC_landuse_name'].values[0]

        ## add the below ground biomass value based on the ABGB as well
        value_ABG_NUTS = df_stats[f'{settings.get("carbon_pool")}_mean_yrly'].values[0]
        if calc_BGB:
            if value_ABG_NUTS == no_data:
                df_stats['BGB_biomass_mean_yrly'] = no_data
            else:
                df_stats['BGB_biomass_mean_yrly'] = df_stats[f'{settings.get("carbon_pool")}_mean_yrly'] * RTS_ratio

        ## calculate now the total potential of carbon increase at the year based on the difference
        # between the potential year and the baseline year
        nr_years_future = int(year_potential - year_baseline)

        ## now calculat the increase based on the percentage to reforest
        df_stats['nr_pixels'] = df_stats['nr_pixels'] * dict_perc_reforest_info.get(IPCC_cat).get('Perc_reforest')

        if value_ABG_NUTS != no_data:
            df_stats[f'{settings.get("carbon_pool")}_total_{str(year_potential)}'] = df_stats[f'{settings.get("carbon_pool")}' \
                                                                                        f'_mean_yrly']*nr_years_future * \
                                                                                    df_stats['nr_pixels'].values[0]
            df_stats[f'BGB_biomass_total_{str(year_potential)}'] = df_stats['BGB_biomass_mean_yrly'] * nr_years_future * \
                                                                    df_stats['nr_pixels'].values[0]
        else:
            df_stats[f'{settings.get("carbon_pool")}_total_{str(year_potential)}'] = no_data
            df_stats[f'BGB_biomass_total_{str(year_potential)}'] = no_data



        df_stats['IPCC_cat'] = IPCC_cat
        df_stats['NUTS_ID'] = spatial_layer.NUTS_ID
        df_stats['NUTS_LEVEL'] = spatial_layer.LEVL_CODE

        ### Assign also used factors for the calculation to shapefile
        if settings.get('run_NUTS_specific_scenario'):
            dict_Slope_factors_info = get_factors_from_NUTS(settings, settings.get('afforestation_scenario'), 'Slope',
                                                            id_LUT_carbon_pool='afforestation')
            dict_Tree_prob_factors_info = get_factors_from_NUTS(settings, settings.get('afforestation_scenario'), 'Tree_prob',
                                                            id_LUT_carbon_pool='afforestation')
            dict_Tree_species_factors_info = get_factors_from_NUTS(settings, settings.get('afforestation_scenario'), 'Tree_species',
                                                            id_LUT_carbon_pool='afforestation')

        else:
            dict_Slope_factors_info = settings.get('afforestation_scenario')
            dict_Tree_prob_factors_info = settings.get('afforestation_scenario')
            dict_Tree_species_factors_info = settings.get('afforestation_scenario')

        if not IPCC_cat in dict_Slope_factors_info.keys():
                ### Factors are not defined for this crop type
                continue
        df_stats['Slope_factor'] = dict_Slope_factors_info.get(IPCC_cat).get('Slope')
        df_stats['Slope_src'] = dict_Slope_factors_info.get(IPCC_cat).get('input_source')

        df_stats['Tree_prob'] = dict_Tree_prob_factors_info.get(IPCC_cat).get('input_source')
        df_stats['Tree_prob_src'] = dict_Tree_prob_factors_info.get(IPCC_cat).get('input_source')

        df_stats['Tree_species_factor'] = dict_Tree_species_factors_info.get(IPCC_cat).get('Tree_species')
        df_stats['Tree_species_src'] = dict_Tree_species_factors_info.get(IPCC_cat).get('input_source')

        df_stats['perc_reforest'] = dict_perc_reforest_info.get(IPCC_cat).get('Perc_reforest')
        df_stats['perc_reforest_src'] = dict_perc_reforest_info.get(IPCC_cat).get('input_source')


        RCP_scenario = settings.get('afforestation_scenario').get(IPCC_cat).get('RCP')
        Year_potential = settings.get('afforestation_scenario').get(IPCC_cat).get('Year_potential')
        df_stats['RCP'] = RCP_scenario
        df_stats['Year_potential'] = Year_potential


        df_stats['geometry'] = [spatial_layer.geometry]

        ## set nan to 0
        df_stats = df_stats.fillna(0)

        lst_df_stats_NUTS.append(df_stats)

    return pd.concat(lst_df_stats_NUTS)



































