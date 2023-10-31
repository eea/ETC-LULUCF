"""
In this script all the functions needed for the estimation of the biomass carbon pool will be stored
"""

# Import needed modules
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
from osgeo import gdal
import rasterio
import datetime
from rasterstats import zonal_stats





def define_affor_areas(settings, slope_max=87.5):
    #logger.add("loguru_biom_helper_functions_define_affor_areas.log") # s4e
    """
    Function that will help in finding the potential afforestation areas based on CLC input, DEM (slope), N2K &
    external mask
    :param settings: different parameters used to create the afforestation mask
    :param slope_max: the maximum slope under which afforestation can be applied
    :return: Map with a pixel-based indication where afforestation could be applied
    """
    print("Feedback from define_affor_areas:.....................................")
    # load all the required rasters and check if they are present
    AFFORESTATION_MASK_DATASETS = settings.get(
        'DATASETS').get('AFFORESTATION_MASK')

    # load the defined settings of the afforestation area
    AFFORESTATION_CONFIG = settings.get(
        'SCENARIO_SPECS').get('afforestation_config')

    # CLC_ACC raster
    CLC_dir = AFFORESTATION_MASK_DATASETS.get('CLC')

    # DEM
    print("DEM:")
    DEM_dir = AFFORESTATION_MASK_DATASETS.get('DEM')
    if not os.path.exists(DEM_dir):
        print(DEM_dir)
        raise Exception
    print(DEM_dir)
    print("--------------------------------------------------------")
    # N2K
    print("N2k:")
    N2K_dir = AFFORESTATION_MASK_DATASETS.get('N2K')
    if not os.path.exists(N2K_dir):
        N2K_dir_coarse_res = Path(N2K_dir).parent.joinpath(
            Path(N2K_dir).name.replace('_100m', ''))
        print("if 1")
        if os.path.exists(N2K_dir_coarse_res):
            # warp and resample the layer to the proper extent
            # and resolution based on a reference raster
            # (CLC in this case)
            resample_raster_to_ref([N2K_dir_coarse_res], CLC_dir, 'EEA39',
                                   Path(N2K_dir_coarse_res).parent.as_posix(),
                                   resample_factor=1, overwrite=settings.get("CONFIG_SPECS").get('overwrite'),
                                   resampling=True,
                                   outname=Path(N2K_dir).name)
            print("if 2")
        else:
            print("else")
            print(N2K_dir_coarse_res)
            raise Exception
    
    print(N2K_dir)
    print("--------------------------------------------------------")
    # Exter
    # nal mask (if some areas or not suitable for afforestation). 1 --> suitable. 0 --> not suitable
    External_mask_dir = AFFORESTATION_MASK_DATASETS.get('external_mask')

    if External_mask_dir is not None:
        if not os.path.exists(External_mask_dir):
            print(External_mask_dir)
            raise Exception

    # Start the masking
    # load the configuration for the afforestation mask
    factor_scenario = settings.get('SCENARIO_SPECS').get('afforestation_config')

    overwrite = settings.get('CONFIG_SPECS').get('overwrite')

    # get now a dict with the NUTS specific configuration
    # and replace it with the default scenario.
    # if not NUTS specific configuration use the default
    # configuration applicable at EU scale
    if settings.get('CONFIG_SPECS').get('run_NUTS_SPECIFIC'):
        factor_scenario = get_factors_from_NUTS(
            settings, factor_scenario, 'Slope', id_LUT_carbon_pool='afforestation')

    # load the geometry of the NUTS region which is needed for block based opening of the rasters
    if not settings.get('NUTS3_info').geometry is None:
        bounds = settings.get('NUTS3_info').geometry.bounds

    else:
        # ignore empty regions
        logger.info(
            f'NO GEOMETRY PROVIDED FOR AREA {settings.get("NUTS3_info").NUTS_ID}')
        return None

    # Define the output folder in which the afforesation mask should be written
    Basefolder_mask = Path(settings.get('CONFIG_SPECS').get('Basefolder_output')
                           ).joinpath('Afforestation_mask')
    if settings.get('CONFIG_SPECS').get('Country') is None and not settings.get('CONFIG_SPECS').get('block_based_processing'):
        outdir_mask = Basefolder_mask
        outname_mask = f'Afforestation_mask_{str(settings.get("SCENARIO_SPECS").get("Year_baseline"))}_EEA39.tif'
    else:
        if not settings.get('CONFIG_SPECS').get('block_based_processing'):
            outdir_mask = Path(Basefolder_mask).joinpath(
                settings.get('CONFIG_SPECS').get('Country'))
            outname_mask = f'Afforestation_mask_{str(settings.get("SCENARIO_SPECS").get("Year_baseline"))}_{settings.get("CONFIG_SPECS").get("Country")}.tif'
        else:
            outdir_mask = Path(Basefolder_mask).joinpath(
                settings.get('NUTS3_info')['CNTR_CODE'])
            outname_mask = f'Afforestation_mask_{str(settings.get("SCENARIO_SPECS").get("Year_baseline"))}_' \
                           f'{settings.get("NUTS3_info")["NUTS_ID"]}.tif'

    outdir_mask.mkdir(parents=True, exist_ok=True)

    if not os.path.exists(outdir_mask.joinpath(outname_mask)) or overwrite:
        # TODO need to add the option if not block-based processing
        if settings.get('CONFIG_SPECS').get('block_based_processing'):
            # block based opening of raster(s)
            CLC_LUCAT_raster, src_CLC_raster = open_raster_from_window(
                CLC_dir, bounds)
            if CLC_LUCAT_raster.size == 0:
                logger.info(
                    f'NO PIXEL DATA AVAILABLE FOR NUTS REGION: {settings.get("NUTS3_info").NUTS_ID}')
                return None

            no_data_IPCC_LUCAT = src_CLC_raster.get('nodata')

            # mask out the non-considered valid CLC
            # classes for afforestation

            CLC_classes_filter = settings.get(
                'SCENARIO_SPECS').get('lst_CLC_affor')

            loc_CLC = np.where(np.isin(CLC_LUCAT_raster,
                                       CLC_classes_filter))
            # identify locations of valid CLC classes for afforestation
            CLC_LUCAT_raster[loc_CLC] = 1

            # open the DEM raster only for the window that covers the NUTS region
            DEM_raster, src_DEM_raster = open_raster_from_window(
                DEM_dir, bounds)
            if DEM_raster.size == 0:
                logger.info(
                    f'NO PIXEL DATA AVAILABLE FOR NUTS REGION: {settings.get("NUTS3_info").NUTS_ID}')
                return None
            no_data_DEM_raster = src_DEM_raster.get('nodata')

            # open the N2K raster only for the window that covers the NUTS region
            N2K_raster, src_N2K_raster = open_raster_from_window(
                N2K_dir, bounds)
            if N2K_raster.size == 0:
                logger.info(
                    f'NO PIXEL DATA AVAILABLE FOR NUTS REGION: {settings.get("NUTS3_info").NUTS_ID}')
                return None
            if N2K_raster.size != DEM_raster.size:
                logger.warning(
                    f'N2K RASTER IS NOT FULLY AVALABLE AT NUTS REGION: {settings.get("NUTS3_info").NUTS_ID}')
                return None

            if External_mask_dir is not None:
                # open the external mask raster only for the window that covers the NUTS region
                External_mask_raster, src_external_mask = open_raster_from_window(
                    External_mask_dir, bounds)
                if External_mask_raster.size == 0:
                    logger.info(
                        f'NO PIXEL DATA AVAILABLE FOR NUTS REGION: {settings.get("NUTS3_info").NUTS_ID}')
                    return None
                no_data_external_mask = src_external_mask.get('nodata')

            else:
                External_mask_raster = np.zeros(DEM_raster.shape)
                no_data_external_mask = 255

            ########## TODO: ETC/DI #############

            """
            Load the following datasets (once they are on the same grid as the other layers!!)
            * Nationally designated areas (CDDA)
            * High nature value farmland 2012
            * Extended wetland layer (2018)
            * Peatland map (2017)
            * Abandoned cropland
            """

        # store the info needed for creating output mask layer
        meta_raster = src_CLC_raster
        transform_raster = src_CLC_raster.get('transform')

        meta_raster.update({'nodata': 255, 'dtype': 'uint8'})

        # Create factor raster to write out factors in it
        mask_raster = np.full(CLC_LUCAT_raster.shape, 255).astype('uint8')

        # find location of this LUCAT
        loc_strata = np.where(CLC_LUCAT_raster == 1)

        # the minimum slope at which afforestation may start
        slope_min = AFFORESTATION_CONFIG.get('Slope')

        ########## TODO: ETC/DI #############

        """
        Use the following layers to define suitable locations. 
        Only those pixels where the corresponding datasets do not apply
        could be considered as suitable 
        * Nationally designated areas (CDDA)
        * High nature value farmland 2012
        * Extended wetland layer (2018)
        * Peatland map (2017)

        For the abandoned cropland layer add a different category (class ID), 
        such it can be favored for afforestation
        """

        # Now filter on the locations for afforestation

        loc_affor_option = np.asarray((DEM_raster >= slope_min) & (DEM_raster < slope_max)
                                    & (CLC_LUCAT_raster == 1) & (N2K_raster == 255)).nonzero()
        
        # set all the LU pixels to no afforestation initially
        # this helps to trace how much of the suitable LUCAT can
        # be actually used for afforestation
        mask_raster[loc_strata] = 0



        # the external mask will add some additional areas even if not suitable
        mask_raster[External_mask_raster == 1] = 1  # 1 is suitable

        # set the ones that could be afforested to one
        mask_raster[loc_affor_option] = 1
        ########## TODO: ETC/DI #############

        """
        Append the conditions for the below datasets as the external mask might 
        not identify afforestation if the corresponding datasets apply 
        * Nationally designated areas (CDDA)
        * High nature value farmland 2012
        * Extended wetland layer (2018)
        * Peatland map (2017)
        """
        # if N2K set back to zero (PRIORITY!!!!):
        # Only for those where the CLC class is valid
        loc_mask_N2K = np.where((N2K_raster == 1) & (CLC_LUCAT_raster == 1))
        mask_raster[loc_mask_N2K] = 0 ### OK unti here..66  # original: 0

        ########## TODO: ETC/DI #############

        """
        Assign for the suitable afforestation locations under abandoned 
        cropland a different ID, such that it can be differentiated later on in the workflow
        """

        # ensure that the nodata pixels in one
        # of the layers are set back to no data
        mask_raster[CLC_LUCAT_raster == no_data_IPCC_LUCAT] = 255
        mask_raster[DEM_raster == no_data_DEM_raster] = 255
        mask_raster[External_mask_raster == no_data_external_mask] = 255

        write_raster(np.expand_dims(mask_raster, 0), meta_raster,
                     transform_raster, Path(outdir_mask), outname_mask)

        # apply a cropping of the raster so that cells
        # outside the extent are set to nodata
        if settings.get('CONFIG_SPECS').get('block_based_processing'):
            NUTS_gpd = gpd.GeoDataFrame(
                geometry=[settings.get('NUTS3_info').geometry])
            outdir_raster = Path(outdir_mask).joinpath(outname_mask).as_posix()
            mask_raster_extent([Path(outdir_mask).joinpath(outname_mask).as_posix()], NUTS_gpd,
                               Path(outdir_mask), settings, overwrite=overwrite)
            # remove the unclipped file and the rename the clipped file to the final output name
            os.remove(outdir_raster)
            os.rename(Path(outdir_mask).joinpath(outname_mask.replace('.tif', '_clipped.tif')).as_posix(),
                      outdir_raster)

            mask_raster_clipped, _ = open_raster_from_window(
                outdir_raster, bounds)

            return mask_raster_clipped
        else:
            return mask_raster
    else:
        mask_raster, _ = open_raster_from_window(
            outdir_mask.joinpath(outname_mask), bounds)
        return mask_raster






def create_affor_potential(settings, affor_mask_array):
    """
    Function that will calculate the afforestation carbon potential and rasterize the output
    :param settings: Different parameters used to create the afforestation potential map
    :param affor_mask_array: the array containing the suitable locations for afforestation
    :return: raster with the carbon potential if a certain pixel is afforested
    """

    # load the defined settings of the afforestation area
    AFFORESTATION_CONFIG = settings.get(
        'SCENARIO_SPECS').get('afforestation_config')

    carbon_pool = settings.get('SCENARIO_SPECS').get('carbon_pool')

    # First load the needed layers
    Basefolder_affor_potential = Path(settings.get('CONFIG_SPECS').get(
        'Basefolder_output')).joinpath(f'{carbon_pool}_scenario')
    Basefolder_affor_potential.mkdir(exist_ok=True)
    Files_EUtrees4F_prob = glob.glob(os.path.join(
        settings.get('DATASETS').get('AFFORESTATION_SUITABILITY').get('EU4Trees'), '*tif'))
    # CLC_ACC raster
    CLC_dir = settings.get('DATASETS').get('AFFORESTATION_MASK').get('CLC')

    # Define the location of the afforestation mask layer
    # and the output of the afforestation potential
    if settings.get('CONFIG_SPECS').get('Country') is None and not settings.get('CONFIG_SPECS').get('block_based_processing'):
        outdir_affor_pot = Basefolder_affor_potential
        outname_affor_pot = f'Afforestation_pot_{str(settings.get("SCENARIO_SPECS").get("Year_potential"))}_EEA39.tif'

    else:
        if not settings.get('CONFIG_SPECS').get('block_based_processing'):
            outdir_affor_pot = Path(
                Basefolder_affor_potential).joinpath('EE3A9')
            outname_affor_pot = f'{carbon_pool}_{settings.get("CONFIG_SPECS").get("scenario_name")}_{"EEA39"}.tif'
        else:
            outdir_affor_pot = Basefolder_affor_potential \
                .joinpath(settings.get('NUTS3_info')['CNTR_CODE']).as_posix()
            outname_affor_pot = f'{carbon_pool}_{settings.get("CONFIG_SPECS").get("scenario_name")}' \
                                f'_{settings.get("NUTS3_info")["NUTS_ID"]}.tif'

    os.makedirs(outdir_affor_pot, exist_ok=True)
    affor_pot_dir = os.path.join(outdir_affor_pot, outname_affor_pot)

    overwrite = settings.get('CONFIG_SPECS').get('overwrite')

    # load the geometry of the NUTS region which is needed for
    #  block based opening of the rasters
    if not settings.get('NUTS3_info').geometry is None:
        bounds = settings.get('NUTS3_info').geometry.bounds
    else:
        # ignore empty regions
        logger.info(
            f'NO GEOMETRY PROVIDED FOR AREA {settings.get("NUTS3_info").NUTS_ID}')
        return None

    if not os.path.exists(affor_pot_dir) or overwrite:
        # TODO write scaling factor in metadata of raster
        # just to avoid the saving in floats of the raster
        scaling = settings.get('CONFIG_SPECS').get('scaling')

        if settings.get('CONFIG_SPECS').get('block_based_processing'):
            # block based opening of raster(s)
            if affor_mask_array.size == 0:
                logger.info(
                    f'NO PIXEL DATA AVAILABLE FOR NUTS REGION: {settings.get("NUTS3_info").NUTS_ID}')
                return None

            CLC_LUCAT_raster, src_CLC_raster = open_raster_from_window(
                CLC_dir, bounds)
            if CLC_LUCAT_raster.size == 0:
                logger.info(
                    f'NO PIXEL DATA AVAILABLE FOR NUTS REGION: {settings.get("NUTS3_info").NUTS_ID}')
                return None

            no_data_CLC_LUCAT = src_CLC_raster.get('nodata')

            # Load the LUT with info on the EU4trees and the IPCC related volume increment
            df_trees_biom_increment = pd.read_csv(os.path.join(settings.get('CONFIG_SPECS').get('yield_table_dir')))
            
            # Load the LUT that related to each NUTS region the coresponding forest zone
            df_LUT_forest_zones = pd.read_csv(os.path.join(settings.get('CONFIG_SPECS').get('forest_zone_dir')))

            # store the info needed for creating output mask layer
            meta_raster = src_CLC_raster
            transform_raster = src_CLC_raster.get('transform')
            meta_raster.update({'nodata': 65535, 'dtype': 'uint16', 
                                'count': 6})
            spatial_resolution = transform_raster[0]
            # the area of each pixel in hectares
            A_pixel_ha = int(((spatial_resolution**2) / 10000))

            # create raster to write output in
            affor_yrly_pot_raster = np.full(
                (6,) + affor_mask_array.shape, 65535).astype('uint16')

            # load the specific (flexible) configuration per
            # NUTS region for afforestation
            # otherwise use default factors
            if settings.get('CONFIG_SPECS').get('run_NUTS_SPECIFIC'):

                # load the configuration for the type of
                # tree species that should be afforested
                factor_scenario_tree_species = get_factors_from_NUTS(
                    settings, AFFORESTATION_CONFIG, 'Tree_species',
                    id_LUT_carbon_pool='afforestation',
                    LU_specific=False,
                    folder=settings.get('CONFIG_SPECS').get('NUTS_LUT_folder'))
                tree_species = factor_scenario_tree_species.get('Tree_species')

                # load the configuration on the tree probability to determine if it will be suitable to grow
                factor_scenario_tree_prob = get_factors_from_NUTS(
                    settings, AFFORESTATION_CONFIG, 'Tree_prob',
                    id_LUT_carbon_pool='afforestation',
                    LU_specific=False,
                    folder=settings.get('CONFIG_SPECS').get('NUTS_LUT_folder'))
                tree_prob = factor_scenario_tree_prob.get('Tree_prob')

            else:
                tree_species = AFFORESTATION_CONFIG.get('Tree_species')
                tree_prob = AFFORESTATION_CONFIG.get('Tree_prob')

            # check first if there is any data provided to calculate
            #  the afforestation potental for this specific tree species 
            # and forest zone in the NUTS3 region
            NUTS3_ID = settings.get('NUTS3_info')['NUTS_ID']
            # search the forest zone
            if not NUTS3_ID in df_LUT_forest_zones.NUTS_LEVEL3_ID.to_list():
                logger.warning(f'No forest zone information for {NUTS3_ID}')
                return None, None

            Forest_zone = df_LUT_forest_zones.loc[df_LUT_forest_zones['NUTS_LEVEL3_ID'] == NUTS3_ID]['FOREST_ZONE'].values[0]
            
            # search for any increment data
            no_avail = df_trees_biom_increment.loc[((df_trees_biom_increment['FOREST_ZONE'] == Forest_zone) & (df_trees_biom_increment['SPECIES_NAME'] == tree_species))].empty
            if no_avail:
                logger.warning(f'No yield information available for {NUTS3_ID} in FOREST ZONE {Forest_zone} and species {tree_species}')
                return None, None
            
            # These could not be NUTS specific factors and hence are
            # directly loaded from the default configuration
            RCP_scenario = AFFORESTATION_CONFIG.get('RCP')
            # the year for which the total biomass
            # increase should be calculated
            Year_potential = AFFORESTATION_CONFIG.get('Year_potential')

            # First check if we have the LUT tables to calculate the
            # volume increment for the requested tree species
            if tree_species in df_trees_biom_increment.SPECIES_NAME.to_list():
                # now select the suitable EUTrees-4F layers
                # Open the occurence probability of the
                # selected tree species for the selected RCP
                # do not consider the status of 2005 anymore
                # as the 2035 situation is closer to the current
                # situation 
                EU4Trees_prob_select = [item for item in Files_EUtrees4F_prob if Path(item).stem.startswith(tree_species)]
                EU4Trees_prob_select = [item for item in EU4Trees_prob_select if  RCP_scenario in Path(item).stem]

                if len(EU4Trees_prob_select) == 0:
                    raise Exception(f'Problem with loading tree species prob layer for '
                                    f'{settings.get("NUTS3_info").NUTS_ID}')

                # Open the tree probability rasters
                lst_prob_rasters = []
                for raster in EU4Trees_prob_select:
                    Tree_prob_raster_tmp, meta_tree_prob = open_raster_from_window(
                        raster, bounds)
                    
                    lst_prob_rasters.append(Tree_prob_raster_tmp)

                # concatenate the different probability rasters 
                # and check if the minimum probability exceeds the thresholds
                Tree_prob_raster = np.amin(np.dstack(lst_prob_rasters), axis=2)
                no_data_tree_prob = meta_tree_prob.get('nodata')

                # based on the LUT, load the annual C seq increase 
                # including the min, average and maximum sequestration
                avg_ann_C_seq_tree_species_age_0_20 = df_trees_biom_increment.loc[((df_trees_biom_increment['FOREST_ZONE'] == Forest_zone) 
                & (df_trees_biom_increment['SPECIES_NAME'] == tree_species))]['yrl_C_seq_age_0_20_avg'].values[0]
                avg_ann_C_seq_tree_species_age_21_30 = df_trees_biom_increment.loc[((df_trees_biom_increment['FOREST_ZONE'] == Forest_zone) 
                & (df_trees_biom_increment['SPECIES_NAME'] == tree_species))]['yrl_C_seq_age_21_30_avg'].values[0]
                max_ann_C_seq_tree_species_age_0_20 = df_trees_biom_increment.loc[((df_trees_biom_increment['FOREST_ZONE'] == Forest_zone) 
                & (df_trees_biom_increment['SPECIES_NAME'] == tree_species))]['yrl_C_seq_age_0_20_max'].values[0]
                max_ann_C_seq_tree_species_age_21_30 = df_trees_biom_increment.loc[((df_trees_biom_increment['FOREST_ZONE'] == Forest_zone) 
                & (df_trees_biom_increment['SPECIES_NAME'] == tree_species))]['yrl_C_seq_age_21_30_max'].values[0]
                min_ann_C_seq_tree_species_age_0_20 = df_trees_biom_increment.loc[((df_trees_biom_increment['FOREST_ZONE'] == Forest_zone) 
                & (df_trees_biom_increment['SPECIES_NAME'] == tree_species))]['yrl_C_seq_age_0_20_min'].values[0]
                min_ann_C_seq_tree_species_age_21_30 = df_trees_biom_increment.loc[((df_trees_biom_increment['FOREST_ZONE'] == Forest_zone) 
                & (df_trees_biom_increment['SPECIES_NAME'] == tree_species))]['yrl_C_seq_age_21_30_min'].values[0]

                # now find the location on which can be
                # afforested (based on the afforestation  mask)
                # the presence of the LUCAT and the pixel
                # for which the tree probability is high enough
                loc_affor = np.where((affor_mask_array == 1) &
                                     ((Tree_prob_raster > (tree_prob * 10)) & (Tree_prob_raster < no_data_tree_prob)))

                # The yearly living biomass (above and below ground) increase is be computed based on a  IPCC formula as follows:
                """
                annual C sink = Volume increment x Density x BCEFx (1+Rootoshoot) x CF 
                """
                # Multiple by the pixel area to know the increase per pixel
                # first band for average age 0-20
                # second band for average age 21-30
                # third band for min age 0-20
                # fourth band for min age 21-30
                # fifth band for max age 0-20
                # sixt band for max age 21-30
                affor_yrly_pot_raster[0,loc_affor[0], loc_affor[1]] = int(avg_ann_C_seq_tree_species_age_0_20 *scaling * A_pixel_ha)
                affor_yrly_pot_raster[1,loc_affor[0], loc_affor[1]] = int(avg_ann_C_seq_tree_species_age_21_30 *scaling * A_pixel_ha)
                affor_yrly_pot_raster[2,loc_affor[0], loc_affor[1]] = int(min_ann_C_seq_tree_species_age_0_20 *scaling * A_pixel_ha)
                affor_yrly_pot_raster[3,loc_affor[0], loc_affor[1]] = int(min_ann_C_seq_tree_species_age_21_30 *scaling * A_pixel_ha)
                affor_yrly_pot_raster[4,loc_affor[0], loc_affor[1]] = int(max_ann_C_seq_tree_species_age_0_20 *scaling * A_pixel_ha)
                affor_yrly_pot_raster[5,loc_affor[0], loc_affor[1]] = int(max_ann_C_seq_tree_species_age_21_30 *scaling * A_pixel_ha)
            # ensure that the nodata pixels in one
            # of the layers are set back to no data
            affor_yrly_pot_raster[:, CLC_LUCAT_raster ==
                                  no_data_CLC_LUCAT] = 65535

            write_raster(affor_yrly_pot_raster, meta_raster,
                         transform_raster, Path(outdir_affor_pot), outname_affor_pot)

            # apply a cropping of the raster so that cells
            # outside the extent are set to nodata
            if settings.get('CONFIG_SPECS').get('block_based_processing'):
                NUTS_gpd = gpd.GeoDataFrame(
                    geometry=[settings.get('NUTS3_info').geometry])
                outdir_raster = Path(outdir_affor_pot).joinpath(
                    outname_affor_pot).as_posix()
                mask_raster_extent([Path(outdir_affor_pot).joinpath(outname_affor_pot).as_posix()], NUTS_gpd,
                                   Path(outdir_affor_pot), settings, overwrite=overwrite)
                # remove the unclipped file and the rename
                # the clipped file to the final output name
                os.remove(outdir_raster)
                os.rename(Path(outdir_affor_pot).joinpath(outname_affor_pot.replace('.tif', '_clipped.tif')).as_posix(),
                          outdir_raster)

                affor_yrly_pot_raster_clipped, _ = open_raster_from_window(
                    outdir_raster, bounds)

                return affor_yrly_pot_raster_clipped, outname_affor_pot
            else:
                return affor_yrly_pot_raster, outname_affor_pot

    else:
        affor_yrly_pot_raster, _ = open_raster_from_window(
            Path(outdir_affor_pot).joinpath(outname_affor_pot), bounds)
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


def reproject_raster(input_raster: str, output_raster: str, src_nodata: int, dst_nodata: int,
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
    warp = gdal.Warp(output_raster, input_raster_obj, dstSRS=f'EPSG:{str(crs_to)}',
                     outputType=gdal.GDT_UInt16, resampleAlg=resample_method,
                     srcNodata=src_nodata, dstNodata=dst_nodata)
    warp = None  # Closes the files


# add metadata to created raster
def create_metadata_description_afforestation(settings: dict, extent: str = 'EEA39') -> dict:
    """
    Function that will create a description that will be saved in the final afforestation
     scenario layer and also add band metadata
    :param settings:dictionary that contains some settings used for writing out the metadata
    :param extent: define if the description should be created for 'EEA39' scale or 'NUTS' scale
    :return:description that will be written into the raster
    """

    AFFORESTATION_CONFIG = settings.get(
        'SCENARIO_SPECS').get('afforestation_config')

    dict_general = {'copyright': 'DO NOT DISTRIBUTE',
                    'commit_id': settings.get('commit_id'),
                    'doi': 'Prototype dataset, no registered DOI',
                    'institution': 'EEA',
                    'name': settings.get('CONFIG_SPECS').get('scenario_name') + '_afforestation_yrl_increment',
                    'processing_date': datetime.datetime.now().date(),
                    'source': ' Derived from JRC YIELD TABLE, IPCC COEFFICIENT and EU4TREES projection',
                    'RCP': AFFORESTATION_CONFIG.get('RCP'),
                    'Year_potential': settings.get('SCENARIO_SPECS').get('Year_potential'),
                    'title': settings.get('CONFIG_SPECS').get('scenario_name')+' scenario at 100 m'}

    if settings.get('CONFIG_SPECS').get('run_NUTS_SPECIFIC') and extent == 'EEA39':
        dict_general.update({'Slope': f'Slope afforestation threshold value {AFFORESTATION_CONFIG.get("Slope")}% & FROM {"EEA39"}',
                             'Tree prob': f'Tree occurence probability value {AFFORESTATION_CONFIG.get("Tree_prob")}% & FROM {"EEA39"}',
                             'Tree species': f'Tree species selected {AFFORESTATION_CONFIG.get("Tree_species")} & FROM {"EEA39"}', })

    elif settings.get('CONFIG_SPECS').get('run_NUTS_SPECIFIC') and extent == 'NUTS':

        dict_Slope_factors_info = dict_Tree_prob_factors_info = get_factors_from_NUTS(
            settings, AFFORESTATION_CONFIG, 'Slope',
            id_LUT_carbon_pool='afforestation',
            LU_specific=False,
            folder=settings.get('CONFIG_SPECS').get('NUTS_LUT_folder'))

        dict_Tree_prob_factors_info = get_factors_from_NUTS(
            settings, AFFORESTATION_CONFIG, 'Tree_prob',
            id_LUT_carbon_pool='afforestation',
            LU_specific=False,
            folder=settings.get('CONFIG_SPECS').get('NUTS_LUT_folder'))

        dict_Tree_species_selection_info = get_factors_from_NUTS(
            settings, AFFORESTATION_CONFIG, 'Tree_species',
            id_LUT_carbon_pool='afforestation',
            LU_specific=False,
            folder=settings.get('CONFIG_SPECS').get('NUTS_LUT_folder'))

        dict_general.update({'Slope': 'Slope afforesation threshold value {}% & FROM {}'.format(dict_Slope_factors_info.get("Slope"),
                                                                                                dict_Slope_factors_info.get("input_source"))})

        dict_general.update({'Tree prob': 'Tree occurence probability value {}% & FROM {}'.format(dict_Tree_prob_factors_info.get("Tree_prob"),
                                                                                                  dict_Tree_prob_factors_info.get("input_source"))})

        dict_general.update({'Tree species': 'Tree species selected {} & FROM {}'.format(dict_Tree_species_selection_info.get("Tree_species"),
                                                                                         dict_Tree_species_selection_info.get("input_source"))})

    else:

        dict_general.update({'Slope': 'Slope afforesation threshold value {}% FROM DEFAULT SCENARIO EEA39'
                            .format(AFFORESTATION_CONFIG.get("Slope"))})

        dict_general.update({'Tree prob': 'Tree occurence probability value {} FROM DEFAULT SCENARIO EEA39'
                            .format(AFFORESTATION_CONFIG.get("Tree_prob"))})
        dict_general.update({'Tree species': 'Tree species selected {} FROM DEFAULT SCENARIO EEA39'
                            .format(AFFORESTATION_CONFIG.get("Tree_species"))})

    dict_general.update({'MAP_EXTENT': extent})

    dict_band = {
        'nodata': 65535,
        'long_name': 'LB (tonnes/C/ha) '
    }

    return dict_general, dict_band


def calc_stats_biomass_NUTS(raster_dir: str, spatial_layer: gpd,
                            settings: dict, stats_type: list = ['mean', 'count']) -> pd.DataFrame:
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
    # AVG
    avg_raster_values_age_0_20 = raster_object.read(1)
    avg_raster_values_age_21_30 = raster_object.read(2)
    # MIN
    min_raster_values_age_0_20 = raster_object.read(3)
    min_raster_values_age_21_30 = raster_object.read(4)
    # MAX
    max_raster_values_age_0_20 = raster_object.read(5)
    max_raster_values_age_21_30 = raster_object.read(6)
    affine = raster_object.meta.get('transform')
    no_data = raster_object.meta.get('nodata')
    scaling = settings.get('CONFIG_SPECS').get('scaling')

    AFFORESTATION_CONFIG = settings.get(
        'SCENARIO_SPECS').get('afforestation_config')

    # load some the baseline year and the target year for the calculation
    year_baseline = settings.get('SCENARIO_SPECS').get('Year_baseline')
    year_potential = settings.get('SCENARIO_SPECS').get('Year_potential')

    lst_df_stats_NUTS = []

    # also take into account the percentage that should be
    # reforest for the total increase calculation
    if settings.get('CONFIG_SPECS').get('run_NUTS_SPECIFIC'):  # otherwise use default factors

        dict_perc_reforest_info = get_factors_from_NUTS(
            settings, AFFORESTATION_CONFIG, 'Perc_reforest',
            id_LUT_carbon_pool='afforestation',
            LU_specific=False,
            folder=settings.get('CONFIG_SPECS').get('NUTS_LUT_folder'))
    else:
        dict_perc_reforest_info = AFFORESTATION_CONFIG

    lst_stats = []
    if 'mean' in stats_type:
        idx_mean = stats_type.index('mean')
        # AVG
        avg_stats_age_0_20 = zonal_stats(spatial_layer.geometry, avg_raster_values_age_0_20, 
                                          affine=affine, nodata=no_data, stats=stats_type[idx_mean])
        avg_stats_age_21_30 = zonal_stats(spatial_layer.geometry, avg_raster_values_age_21_30, 
                                          affine=affine, nodata=no_data, stats=stats_type[idx_mean])
        lst_stats.append(pd.DataFrame.from_dict(avg_stats_age_0_20))
        lst_stats.append(pd.DataFrame.from_dict(avg_stats_age_21_30))
        # MIN
        min_stats_age_0_20 = zonal_stats(spatial_layer.geometry, min_raster_values_age_0_20, 
                                          affine=affine, nodata=no_data, stats=stats_type[idx_mean])
        min_stats_age_21_30 = zonal_stats(spatial_layer.geometry, min_raster_values_age_21_30, 
                                          affine=affine, nodata=no_data, stats=stats_type[idx_mean])
        lst_stats.append(pd.DataFrame.from_dict(min_stats_age_0_20))
        lst_stats.append(pd.DataFrame.from_dict(min_stats_age_21_30))
        # MAX
        max_stats_age_0_20 = zonal_stats(spatial_layer.geometry, max_raster_values_age_0_20, 
                                          affine=affine, nodata=no_data, stats=stats_type[idx_mean])
        max_stats_age_21_30 = zonal_stats(spatial_layer.geometry, max_raster_values_age_21_30, 
                                          affine=affine, nodata=no_data, stats=stats_type[idx_mean])
       
        lst_stats.append(pd.DataFrame.from_dict(max_stats_age_0_20))
        lst_stats.append(pd.DataFrame.from_dict(max_stats_age_21_30))
    if 'count' in stats_type:
        idx_count = stats_type.index('count')
        count_occur = np.copy(avg_raster_values_age_0_20)
        # if afforestation possible set to 1
        count_occur[count_occur != no_data] = 1
        stats_count = zonal_stats(spatial_layer.geometry, count_occur, affine=affine,
                                  nodata=no_data, stats=stats_type[idx_count])

        lst_stats.append(pd.DataFrame.from_dict(stats_count))

    df_stats = pd.concat(lst_stats, axis=1)
    df_stats.columns = [
        f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}_mean_yrly_age_0_20',
        f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}_mean_yrly_age_21_30', 
        f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}_min_yrly_age_0_20',
        f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}_min_yrly_age_21_30', 
        f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}_max_yrly_age_0_20',
        f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}_max_yrly_age_21_30', 
        'nr_pixels']
   
    # add the outcome to the proper scale
    # AVG
    df_stats[f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}_mean_yrly_age_0_20'] = df_stats[
        f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}_mean_yrly_age_0_20']/scaling
    df_stats[f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}_mean_yrly_age_21_30'] = df_stats[
        f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}_mean_yrly_age_21_30']/scaling
    # MIN
    df_stats[f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}_min_yrly_age_0_20'] = df_stats[
        f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}_min_yrly_age_0_20']/scaling
    df_stats[f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}_min_yrly_age_21_30'] = df_stats[
        f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}_min_yrly_age_21_30']/scaling
    
    # MAX
    df_stats[f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}_max_yrly_age_0_20'] = df_stats[
        f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}_max_yrly_age_0_20']/scaling
    df_stats[f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}_max_yrly_age_21_30'] = df_stats[
        f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}_max_yrly_age_21_30']/scaling

    df_stats = df_stats.round(2)

    # the average LB increase for the NUTS region
    avg_value_LB_NUTS_age_0_20 = df_stats[f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}_mean_yrly_age_0_20'].values[0]

    # calculate now the total potential of carbon increase at the
    # year based on the difference
    # between the potential year and the baseline year
    nr_years_future = int(year_potential - year_baseline)

    # now calculate the increase based on the percentage to reforest
    df_stats['nr_pixels'] = df_stats['nr_pixels'] * \
        (dict_perc_reforest_info.get('Perc_reforest')/100)
    
    
    if avg_value_LB_NUTS_age_0_20 != no_data:
        if nr_years_future <= 20:
            # AVG
            avg_tot_seq_C_final = df_stats[f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}'
                                           f'_mean_yrly_age_0_20'] * nr_years_future * df_stats['nr_pixels'].values[0]
            # MIN
            min_tot_seq_C_final = df_stats[f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}'
                                           f'_min_yrly_age_0_20'] * nr_years_future * df_stats['nr_pixels'].values[0]
            # MAX
            max_tot_seq_C_final = df_stats[f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}'
                                           f'_max_yrly_age_0_20'] * nr_years_future * df_stats['nr_pixels'].values[0]
        else:
            # also take the part after 20 years into account
            # AVG
            tot_seq_C_part1 = df_stats[f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}'
                                           f'_mean_yrly_age_0_20'] * 20 * df_stats['nr_pixels'].values[0]
            tot_seq_C_part2 = df_stats[f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}'
                                           f'_mean_yrly_age_21_30'] * (nr_years_future-20) * df_stats['nr_pixels'].values[0]
            avg_tot_seq_C_final  = tot_seq_C_part1 + tot_seq_C_part2

            # MIN
            tot_seq_C_part1 = df_stats[f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}'
                                           f'_min_yrly_age_0_20'] * 20 * df_stats['nr_pixels'].values[0]
            tot_seq_C_part2 = df_stats[f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}'
                                           f'_min_yrly_age_21_30'] * (nr_years_future-20) * df_stats['nr_pixels'].values[0]
            min_tot_seq_C_final  = tot_seq_C_part1 + tot_seq_C_part2

            # MAX
            tot_seq_C_part1 = df_stats[f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}'
                                           f'_max_yrly_age_0_20'] * 20 * df_stats['nr_pixels'].values[0]
            tot_seq_C_part2 = df_stats[f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}'
                                           f'_max_yrly_age_21_30'] * (nr_years_future-20) * df_stats['nr_pixels'].values[0]
            max_tot_seq_C_final  = tot_seq_C_part1 + tot_seq_C_part2

        # AVG
        df_stats[f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}_total_{str(year_potential)}_avg'] = avg_tot_seq_C_final
        # MIN
        df_stats[f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}_total_{str(year_potential)}_min'] = min_tot_seq_C_final
        # MAX
        df_stats[f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}_total_{str(year_potential)}_max'] = max_tot_seq_C_final

    else:
        df_stats[f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}_total_{str(year_potential)}_avg'] = no_data
        df_stats[f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}_total_{str(year_potential)}_min'] = no_data
        df_stats[f'{settings.get("SCENARIO_SPECS").get("carbon_pool")}_total_{str(year_potential)}_max'] = no_data

    # write some NUTS specific information to the dataframe
    df_stats['NUTS_ID'] = spatial_layer.NUTS_ID
    df_stats['NUTS_LEVEL'] = spatial_layer.LEVL_CODE

    # Assign also used configuration factors for the calculation to shapefile
    if settings.get('CONFIG_SPECS').get('run_NUTS_SPECIFIC'):
        dict_Slope_factors_info = get_factors_from_NUTS(
            settings, AFFORESTATION_CONFIG, 'Slope',
            id_LUT_carbon_pool='afforestation',
            LU_specific=False,
            folder=settings.get('CONFIG_SPECS').get('NUTS_LUT_folder'))
        dict_Tree_prob_factors_info = get_factors_from_NUTS(
            settings, AFFORESTATION_CONFIG, 'Tree_prob',
            id_LUT_carbon_pool='afforestation',
            LU_specific=False,
            folder=settings.get('CONFIG_SPECS').get('NUTS_LUT_folder'))
        dict_Tree_species_factors_info = get_factors_from_NUTS(
            settings, AFFORESTATION_CONFIG, 'Tree_species',
            id_LUT_carbon_pool='afforestation',
            LU_specific=False,
            folder=settings.get('CONFIG_SPECS').get('NUTS_LUT_folder'))

    else:
        dict_Slope_factors_info = AFFORESTATION_CONFIG
        dict_Tree_prob_factors_info = AFFORESTATION_CONFIG
        dict_Tree_species_factors_info = AFFORESTATION_CONFIG

    df_stats['Slope_factor'] = dict_Slope_factors_info.get('Slope')
    df_stats['Slope_src'] = dict_Slope_factors_info.get('input_source')

    df_stats['Tree_prob'] = dict_Tree_prob_factors_info.get('Tree_prob')
    df_stats['Tree_prob_src'] = dict_Tree_prob_factors_info.get('input_source')

    df_stats['Tree_species_factor'] = dict_Tree_species_factors_info.get(
        'Tree_species')
    df_stats['Tree_species_src'] = dict_Tree_species_factors_info.get(
        'input_source')

    df_stats['perc_reforest'] = dict_perc_reforest_info.get('Perc_reforest')
    df_stats['perc_reforest_src'] = dict_perc_reforest_info.get('input_source')

    RCP_scenario = AFFORESTATION_CONFIG.get('RCP')
    Year_potential = AFFORESTATION_CONFIG.get('Year_potential')
    df_stats['RCP'] = RCP_scenario
    df_stats['Year_potential'] = Year_potential

    # add also the forest type (D/C) and growth speed type (S/F)
    df_F_info = pd.read_csv(os.path.join(settings.get('CONFIG_SPECS').get('yield_table_dir')))
    # GET forest type (FT)
    FT = df_F_info.loc[df_F_info['SPECIES_NAME'] == dict_Tree_species_factors_info.get(
        'Tree_species')]['FT'].values[0]
    df_stats['FT'] = FT
    # Get growth speed type
    FGS = df_F_info.loc[df_F_info['SPECIES_NAME'] == dict_Tree_species_factors_info.get(
        'Tree_species')]['FGS'].values[0]
    df_stats['FGS'] = FGS


    df_stats['geometry'] = [spatial_layer.geometry]

    # set nan to 0
    df_stats = df_stats.fillna(0)

    lst_df_stats_NUTS.append(df_stats)

    return pd.concat(lst_df_stats_NUTS)


def add_atrributes_stats_afforestation(spatial_layer: gpd, level_NUTS_focus=None,
                                       spatial_resolution: int = 100) -> gpd:
    """
    Function that will insert some additional columns that could be used for the interpretation or analysis
    of the SOC stats.
    :param spatial_layer: Geopandas layer in which the SOC stats are written
    :param level_NUTS_focus: if needed specify the NUTS level of focus for the stats
    :param spatial_resolution: the spatial resolution on which the aggregation was done
    :return: update of the geodataframe with some additional attributes
    """
    if level_NUTS_focus != None:
        spatial_layer = spatial_layer.loc[spatial_layer.NUTS_LEVEL ==
                                          level_NUTS_focus]

    spatial_layer['area_NUTS_[ha]'] = (
        spatial_layer.geometry.area)/(spatial_resolution**2)
    spatial_layer = spatial_layer.rename(
        columns={'nr_pixels': 'area_afforestation_[ha]'})
    spatial_layer['%_cover_afforestation'] = (
        spatial_layer['area_afforestation_[ha]'] / spatial_layer['area_NUTS_[ha]']) * 100

    return spatial_layer
