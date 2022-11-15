"""
Constant settings that can be recycled among most of the scripts for Carbon seq mapping
"""

import os
import geopandas as gpd

### Some default settings

dir_signature = 'L:'
overwrite = False
Basefolder_strata = os.path.join(dir_signature, 'etc', 'lulucf', 'strata')


Basefolder_input_data = os.path.join(dir_signature,'etc','lulucf','input_data')
CLC_ACC_folder = os.path.join(dir_signature, 'input_data', 'general', 'CLCACC')
type_method = 'LUT' ### define the kind of method that shoudl be used to define the carbon seq basline and potential

shp_extent_map = gpd.read_file(os.path.join(dir_signature, 'etc','lulucf','AOI','EEA39_extent_noDOM.shp'))


path_IPCC_climate_resampled = os.path.join(Basefolder_input_data,'IPCC', 'IPCC_data_converted', 'climate_final_map',
                                           'ipcc_climate_zone_100m_EPSG3035_EEA39.tif')
path_IPCC_soil_resampled = os.path.join(Basefolder_input_data,'IPCC', 'IPCC_data_converted', 'soil_final_map',
                                        'ipcc_soil_type_100m_EPSG3035_EEA39.tif')


#LUT CLC IPCC MAPPING

CLC_IPCC_mapping = {
    'Forest': {'CLC_classes': [311,312,313,324,334],
               'ID_class': 1},
    'Cropland': {'CLC_classes': [211,212,213,221,222,223,
                                 241,242,243,244],
                 'ID_class': 2},
    'Grassland': {'CLC_classes': [231,321,322,323],
                  'ID_class': 3},
    'Wetlands': {'CLC_classes': [411,412,421,422,423
        ,511,512,521,522,523],
                 'ID_class':4},
    'Settlements': {'CLC_classes': [111,112,121,
                                    122,123,124,
                                    131,132,133,
                                    141,142],
                    'ID_class': 5},
    'Other_land': {'CLC_classes': [331,332,333,335],
                   'ID_class': 6}}


# LUT management activaties and carbon pool

dict_C_pool_mng_options = {'afforestation': 'ABGbiomass',
                           'cropland_management': 'SOC',
                           'grassland_management': 'SOC'}