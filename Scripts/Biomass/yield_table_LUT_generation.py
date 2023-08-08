""""
In this code the yield table from JRC and harmonized by Lucia
will be used to automatize the generation of a LUT that can be 
used to estimate the carbon sequestration potential for planting
EUTrees-4F tree species. 
"""


# Required packages
from loguru import logger
import os
import pandas as pd
from pathlib import Path
import geopandas as gpd


def get_sheet(dir_excel, sheet_name):
    df_sheet = pd.read_excel(dir_excel, sheet_name=sheet_name)
    return df_sheet

def get_df_coefficient(df, column_name, coefficient_name):
    column_coef = df[[column_name]]
    # if is nan check if in the 'Not specified'
    # column some values are given
    if coefficient_name == 'Iv':
        column_not_specified = df[['Not specified']]
        # create final column filled by the 
        # not specified column where needed  
        column_final = [column_coef.iloc[item].values[0] if not pd.isnull(column_coef.iloc[item]).values[0] else column_not_specified.iloc[item].values[0] for item in column_not_specified.index]
        df[coefficient_name] = column_final

    else:
        df[coefficient_name] = column_coef
    
        

    return df

def get_coefficient(species_code, forest_zone, df_coeff, coeff_name):
    """
    Based on the species code and forest zone 
    find the corresponding coefficient
    """

    # if increment also filter specifically on forest zone
    if coeff_name == 'Iv':
        coeff_row = df_coeff.loc[((df_coeff['Code'].isin([item for item in df_coeff['Code'] if species_code in item])) & (df_coeff['Forest_zone'] == forest_zone))]
    else:
        coeff_row = df_coeff.loc[df_coeff['Code'].isin([item for item in df_coeff['Code'] if species_code in item])]

    if coeff_row.empty:
        logger.warning(f'No {coeff_name} FOUND FOR SPECIES: {species_code} AND FOREST ZONE: {forest_zone}')
        return None

    if coeff_row.shape[0] > 1:
        logger.error(f'MORE THAN TWO OCCURENCES OF {coeff_name} FOR SPECIES CODE {species_code}: AND FOREST ZONE: {forest_zone}')
        # INFO: TAKE THE MEAN OF THE VALUES (ONLY OCCUR FOR RTS OF ONE FOREST TYPE)
        coeff_row = pd.DataFrame(coeff_row.mean()).T
    coeff = coeff_row[coeff_name].values[0]
    
    return coeff


def create_unique_comb_df(DATAFRAMES):
    #first ensure that each dataframe 
    # has the species code and forest type linked
    df_code = DATAFRAMES.get('Code')

    DATAFRAMES_LINKED = {}
    for dataset_name in DATAFRAMES.keys():
        if dataset_name == 'Code':
            continue
        df = DATAFRAMES.get(dataset_name)

        # add the code to the dataframe
        if not 'Code' in df.columns:
            df['Code'] = [list(df_code.loc[df_code['Forest type IPCC'] == item]['Code']) for item in df['Forest type IPCC'].values]
        
        if not 'Forest type IPCC' in df.columns:
            df['Forest type IPCC'] = [list(df_code.loc[df_code['Code'] == item]['Forest type IPCC']) for item in df['Code'].values]
        
        DATAFRAMES_LINKED.update({dataset_name: df})
   
    # The linked dataframes will be used to create
    # a summary dataframe with all the coefficients 
    # per forest zone and type

    lst_C_seq = []

    for species_code in df_code['Code'].unique():
        species_name = df_code.loc[df_code['Code'] == species_code]['NFI Species'].values[0]
        for forest_zone in DATAFRAMES_LINKED.get('I')['Forest_zone'].unique():
            # get now all the coefficients for the species and forest zone

            # CF
            df_CF = DATAFRAMES_LINKED.get('CF')
            CF_coeff = get_coefficient(species_code, forest_zone, 
                                      df_CF, 'CF')
            if CF_coeff is None:
                continue

            # BCEF
            df_BCEF = DATAFRAMES_LINKED.get('BCEF')
            BCEF_coeff = get_coefficient(species_code, forest_zone, 
                                      df_BCEF, 'BCEF')
            if BCEF_coeff is None:
                continue
            
            # RTS
            df_RTS = DATAFRAMES_LINKED.get('RTS')
            RTS_coeff = get_coefficient(species_code, forest_zone, 
                                        df_RTS, 'RTS')
            if RTS_coeff is None:
                continue
            
            # Iv
            df_I = DATAFRAMES_LINKED.get('I')
            I_coeff = get_coefficient(species_code, forest_zone, 
                                      df_I, 'Iv')
            if I_coeff is None:
                continue
            
            # Apply formula for annual carbon SEQ (tonnes*C/yr)
            yrly_carbon_seq = I_coeff * BCEF_coeff * (1+RTS_coeff) * CF_coeff

            # create now dataframe from all these parameters 
            # for the specific species and forest zone

            lst_C_seq.append(
                (
                species_code,
                species_name,
                forest_zone,
                CF_coeff,
                BCEF_coeff,
                RTS_coeff,
                I_coeff,
                yrly_carbon_seq
                )
            )
    # now create the final LUT for 
    # C seq by applying afforestation
    df_C_seq_LUT = pd.DataFrame(
        lst_C_seq, 
        columns=[
            'SPECIES_CODE',
            'SPECIES_NAME',
            'FOREST_ZONE',
            'CF',
            'BCEF',
            'RTS', 
            'Iv',
            'yrl_C_seq'

        ],
    )

    return df_C_seq_LUT




    
def get_annual_C_species(settings):
    dir_excel = settings.get('CONFIG_SPECS').get('yield_table_dir')

    # carbon fraction of dry matter
    df_CF = get_sheet(dir_excel, 'CF')

    # Biomass conversion and expansion factor (BCEF)
    df_BCEF = get_sheet(dir_excel, 'BCEFremoval')
    df_BCEF = get_df_coefficient(df_BCEF, '<=20', 'BCEF')

    # Root-to-shoot ratio (RTS)
    df_RTS = get_sheet(dir_excel, 'Root-to-shoot')
    df_RTS = get_df_coefficient(df_RTS, '<=75', 'RTS')

    # Net annual increment (I)
    df_I = get_sheet(dir_excel, 'Iv_final')
    df_I = get_df_coefficient(df_I, '<30', 'Iv')

    # Get also some information on the 
    # different forest species codes 
    # (linked with forest types from IPCC) 
    df_F_code = get_sheet(dir_excel, 'Species_code')

    # combine all the dataframes to get 
    # factors for every unique combination
    DATAFRAMES = {
        'CF': df_CF,
        'BCEF': df_BCEF,
        'RTS': df_RTS,
        'I': df_I,
        'Code': df_F_code
    }
    
    df_C_seq_LUT = create_unique_comb_df(DATAFRAMES)

    return df_C_seq_LUT


def get_LUT_forest_zone(settings):
    # open first the shapefile of NUTS regions
    shp_NUTS = gpd.read_file(settings.get('VECTORS').get('NUTS'))

    # get sheet with forest zone information per country
    dir_excel = settings.get('CONFIG_SPECS').get('yield_table_dir')
    df_zones = get_sheet(dir_excel, 'Forest_zone')

    lst_df_LUT_zones = []
    for CNTR in df_zones.CNTR_CODE.unique():
        NUTS3_regions = shp_NUTS.loc[((shp_NUTS.CNTR_CODE == CNTR) & (shp_NUTS.LEVL_CODE == 3))]['NUTS_ID'].to_list()
        ZONE = df_zones.loc[df_zones['CNTR_CODE'] == CNTR]['Forest_zone'].values[0]
        df_NUTS_zone = pd.DataFrame(NUTS3_regions, columns=['NUTS_LEVEL3_ID'])
        df_NUTS_zone['NUTS_LEVEL0_ID'] = [CNTR] * df_NUTS_zone.shape[0]
        df_NUTS_zone['Forest_zone'] = [ZONE] * df_NUTS_zone.shape[0]
        lst_df_LUT_zones.append(df_NUTS_zone)

    df_NUTS_LUT_final = pd.concat(lst_df_LUT_zones)
    return df_NUTS_LUT_final





    


def _main_(settings):

    # first define outfolder and outname 
    # and check if file not yet exists

    outfolder = settings.get('CONFIG_SPECS').get('outfolder')
    outname_C_seq_LUT = 'LUT_C_SEQ_AFFOR_JRC_V1.csv'
    outname_forest_zone_LUT = 'LUT_FOREST_ZONE.csv'
    overwrite = settings.get('CONFIG_SPECS').get('overwrite')

    
    # C_seq LUT
    if not os.path.exists(os.path.join(outfolder, outname_C_seq_LUT)) or overwrite:
        # function that will generate 
        # the annual C increase for all species
        df_C_seq = get_annual_C_species(settings)

        # write the corresponding dataframe 
        # to the output folder
        df_C_seq.to_csv(os.path.join(outfolder, outname_C_seq_LUT), index=False)

    # Forest zone LUT
    if not os.path.exists(os.path.join(outfolder, outname_forest_zone_LUT)) or overwrite:

        # Link to each NUTS region the corresponding forest zone
        df_NUTS_LUT_final = get_LUT_forest_zone(settings)

        df_NUTS_LUT_final.to_csv(os.path.join(outfolder, 
                                              outname_forest_zone_LUT), 
                                              index=False)












if __name__ == '__main__':

    logger.info('*' * 50)

    from constants import (
        Basefolder_strata, 
        dir_signature)

    # define the directory of the yield table:
    name_file = 'inventory_JRC_Data_Collection_100523_nolink_ETC_CA_.xlsx'
    folder_table = os.path.join(Basefolder_strata, 
                                'NUTS_LUT_afforestation_scenario',
                                'JRC_yield_table')
    yield_table_dir = os.path.join(folder_table, name_file)

    CONFIGURATION_SPECS = {
        'yield_table_dir': yield_table_dir,
        "outfolder": Path(yield_table_dir).parent,
        'overwrite': False
    }

    VECTORS = {
        'NUTS': os.path.join(dir_signature, 'etc',
                             'lulucf', 'AOI',
                             'NUTS_RG_20M_2021_3035.shp')
    }

    
    settings = {
                'CONFIG_SPECS': CONFIGURATION_SPECS,
                'VECTORS': VECTORS
               }
    _main_(settings)