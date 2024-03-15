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


def get_sheet(dir_excel, sheet_name,
              skiprows=0):
    df_sheet = pd.read_excel(dir_excel, sheet_name=sheet_name, 
                             skiprows=skiprows)
    return df_sheet


def Iv_aggregate(df,  
                 mgmt_strategy=['E', 'S'],
                 agg_columns =['EEA Forest code',
                               'FOREST_ZONE',
                               'Forest type (IPCC)']):
    """
    Function that will calculate an aggregated value
    of the increment per forest zone and forest type
    """
    # filter first the dataframe on the relevant
    # management
    df_filter = df.loc[df['mgmt_strategy'].isin(mgmt_strategy)]
    df_agg_avg = df_filter.groupby(by=agg_columns).mean(numeric_only=True)
    df_agg_avg = df_agg_avg.reset_index()
    df_agg_min = df_filter.groupby(by=agg_columns).quantile(q=0.1)
    df_agg_min = df_agg_min.reset_index()
    df_agg_max = df_filter.groupby(by=agg_columns).quantile(q=0.9)
    df_agg_max = df_agg_max.reset_index()
    df_agg_avg['AgeCL_2_min'] = df_agg_min['AgeCL_2']
    df_agg_avg['AgeCL_3_min'] = df_agg_min['AgeCL_3']
    df_agg_avg['AgeCL_2_max'] = df_agg_max['AgeCL_2']
    df_agg_avg['AgeCL_3_max'] = df_agg_max['AgeCL_3']
    return df_agg_avg

def get_df_coefficient(df, column_name, coefficient_name):
    if not type(column_name) == list:
        column_name = [column_name]
    column_coef = df[column_name]
    # if is nan check if in the 'Not specified'
    # column some values are given
    if coefficient_name == 'Iv':
        df_agg = Iv_aggregate(df)
        # drop ID EEA column
        if 'ID EEA-JRC' in df_agg:
            df_agg = df_agg.drop(columns=['ID EEA-JRC'])
        for col in column_coef:
            if 'AgeCL_2' in col:
                age_info = '0-20'
            elif 'AgeCL_3' in col:
                age_info = '21-30'
            coef_col = f'{coefficient_name}_{age_info}'
            df_agg[coef_col] = df_agg[col]
            df_agg[coef_col + '_max'] = df_agg[col+'_max']
            df_agg[coef_col + '_min'] = df_agg[col+'_min']
        # drop all age classes:
        age_cols = [item for item in df_agg.columns if 'Age' in item]
        df_agg = df_agg.drop(columns=age_cols)
        return df_agg
    elif coefficient_name == 'Iv_add':
        for col in column_coef:
            if 'Avg' in col:
                df['Iv_0-20'] = df[col]
                df['Iv_21-30'] = df[col]
            else:
                df[f'Iv_0-20_{col.lower()}'] = df[col]
                df[f'Iv_21-30_{col.lower()}'] = df[col]
        df = df.drop(columns=column_coef)
        return df
    else:
        for col in column_coef:
            if '<=20' in col:
                age_info = '0-20'
            else:
                age_info = col
            if 'RTS' in coefficient_name:
                coef_col = f'{coefficient_name}'
            elif 'FG' in coefficient_name:
                coef_col = f'{coefficient_name}'
            elif coefficient_name in ['FT', 'FGS']:
                coef_col = f'{coefficient_name}'
            else:
                coef_col = f'{coefficient_name}_{age_info}'
            df[coef_col] = column_coef[col]
        return df

def get_coefficient(species_code, forest_zone, df_coeff, coeff_name,
                    ZONE_DEP=False, SPECIES_DEP = False):
    """
    Based on the species code and forest zone 
    find the corresponding coefficient
    """
    if not ZONE_DEP and not SPECIES_DEP:
        coeff_row = df_coeff.loc[df_coeff['Code'].isin([item for item in df_coeff['Code'] if species_code in item])]
    elif SPECIES_DEP:
         coeff_row = df_coeff.loc[df_coeff['EEA Forest code'] == species_code]
    else:
        coeff_row = df_coeff.loc[((df_coeff['EEA Forest code'] == species_code)
                                  &(df_coeff['FOREST_ZONE'] == forest_zone))]

    if coeff_row.empty:
        logger.warning(f'No {coeff_name} FOUND FOR SPECIES: {species_code} AND FOREST ZONE: {forest_zone}')
        return None

    if coeff_row.shape[0] > 1 and not coeff_name in ['FT', 'FGS']:
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
            df['Code'] = [list(df_code.loc[df_code['Forest type (IPCC)'] == item]['Code']) for item in df['Forest type (IPCC)'].values]
        
        if not 'Forest type (IPCC)' in df.columns:
            df['Forest type (IPCC)'] = [list(df_code.loc[df_code['Code'] == item]['Forest type (IPCC)']) for item in df['Code'].values]
        
        DATAFRAMES_LINKED.update({dataset_name: df})
   
    # The linked dataframes will be used to create
    # a summary dataframe with all the coefficients 
    # per forest zone and type

    lst_C_seq = []

    for species_code in df_code['Code'].unique():
        species_name = df_code.loc[df_code['Code'] == species_code]['EUTrees4F name'].values[0]
        species_name = species_name.replace(' ', '_')
        for forest_zone in DATAFRAMES_LINKED.get('I')['FOREST_ZONE'].unique():
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
                                      df_BCEF, ['BCEF_0-20', 'BCEF_21-40'])
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
                                      df_I, ['Iv_0-20', 'Iv_21-30',
                                             'Iv_0-20_min', 'Iv_21-30_min',
                                             'Iv_0-20_max', 'Iv_21-30_max'],
                                            ZONE_DEP=True)
            if I_coeff is None:
                continue

            # FT
            df_FG = DATAFRAMES_LINKED.get('FG')
            FT_coeff = get_coefficient(species_code, forest_zone,
                                       df_FG, 'FT', SPECIES_DEP=True)
            
            # FGS
            df_FGS = DATAFRAMES_LINKED.get('FGS')
            FGS_coeff = get_coefficient(species_code, forest_zone,
                                       df_FGS, 'FGS', SPECIES_DEP=True)

            # Apply formula for annual carbon SEQ (tonnes*C/yr)
            # seperate for period 0-20 and 21-30 and addd also min max range
            yrly_carbon_seq_1_avg = I_coeff[0] * BCEF_coeff[0] * (1+RTS_coeff) * CF_coeff
            yrly_carbon_seq_2_avg = I_coeff[1] * BCEF_coeff[1] * (1+RTS_coeff) * CF_coeff
            yrly_carbon_seq_1_min = I_coeff[2] * BCEF_coeff[0] * (1+RTS_coeff) * CF_coeff
            yrly_carbon_seq_2_min = I_coeff[3] * BCEF_coeff[1] * (1+RTS_coeff) * CF_coeff
            yrly_carbon_seq_1_max = I_coeff[4] * BCEF_coeff[0] * (1+RTS_coeff) * CF_coeff
            yrly_carbon_seq_2_max = I_coeff[5] * BCEF_coeff[1] * (1+RTS_coeff) * CF_coeff

            # create now dataframe from all these parameters 
            # for the specific species and forest zone

            lst_C_seq.append((
                    species_code,
                    species_name,
                    forest_zone,
                    FGS_coeff,
                    FT_coeff,
                    CF_coeff,
                    BCEF_coeff[0],
                    BCEF_coeff[1],
                    RTS_coeff,
                    I_coeff[0],
                    I_coeff[1],
                    I_coeff[2],
                    I_coeff[3],
                    I_coeff[4],
                    I_coeff[5],
                    yrly_carbon_seq_1_avg,
                    yrly_carbon_seq_2_avg,
                    yrly_carbon_seq_1_min,
                    yrly_carbon_seq_2_min,
                    yrly_carbon_seq_1_max,
                    yrly_carbon_seq_2_max,
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
            'FGS',
            'FT',
            'CF',
            'BCEF_0_20',
            'BCEF_21_30',
            'RTS',
            'Iv_0_20_avg',
            'Iv_21-30_avg',
            'Iv_0_20_min',
            'Iv_21-30_min',
            'Iv_0_20_max',
            'Iv_21-30_max',
            'yrl_C_seq_age_0_20_avg',
            'yrl_C_seq_age_21_30_avg',
            'yrl_C_seq_age_0_20_min',
            'yrl_C_seq_age_21_30_min',
            'yrl_C_seq_age_0_20_max',
            'yrl_C_seq_age_21_30_max',
        ],
    )
    return df_C_seq_LUT


def compile_species(df, idx_valid=6):
    # keep only the five first columns
    cols_filter = df.columns[0:idx_valid]
    df_filter = df[cols_filter]
    # only consider the once that are availbel in EUTrees-4F
    df_filter_av = df_filter.loc[df_filter['Presence in the EUTrees4F'] == 'y']
    df_filter_av = df_filter_av.rename(columns={'New Code': 'Code'})
    return df_filter_av


def get_annual_C_species(settings):
    dir_excel = settings.get('CONFIG_SPECS').get('yield_table_dir')

    # carbon fraction of dry matter
    df_CF = get_sheet(dir_excel, 'CF')

    # Biomass conversion and expansion factor (BCEF)
    df_BCEF = get_sheet(dir_excel, 'BCEFi')
    df_BCEF = get_df_coefficient(df_BCEF, ['<=20', '21-40'], 'BCEF')

    # Root-to-shoot ratio (RTS)
    df_RTS = get_sheet(dir_excel, 'Root-to-shoot')
    df_RTS = get_df_coefficient(df_RTS, '<=75', 'RTS')

    # Net annual increment (I) based on JRC
    df_I_JRC = get_sheet(dir_excel, 'EEA NAI_Ev clean')
    df_I_JRC = get_df_coefficient(df_I_JRC, ['AgeCL_2', 'AgeCL_3'], 'Iv')
    df_I_JRC = df_I_JRC.reset_index(drop=True)
    df_I = df_I_JRC

    # # Net annual increment (I) based on growth curve calculation
    # df_I_added = get_sheet(dir_excel, 'NAI_even_EEA_add_species')
    # df_I_added = get_df_coefficient(df_I_added, ['Min', 'Max', 'Avg'], 'Iv_add')
    # df_I_added = df_I_added.reset_index(drop=True)

    # # merge the two increment dataframes
    # df_I = pd.concat([df_I_JRC, df_I_added], axis=0)
    # df_I = df_I.reset_index(drop=True)

    # Get forest group info (FT)
    df_FG = get_sheet(dir_excel, 'Species_EEA_type')
    df_FG = get_df_coefficient(df_FG, 'Forest type', 'FT')
    
    # Get forest type growth speed info (FGS)
    df_FGS = get_sheet(dir_excel, 'Species_EEA_type')
    df_FGS = get_df_coefficient(df_FGS, 'Type of specie', 'FGS')

    # Get also some information on the 
    # different forest species codes 
    # (linked with forest types from IPCC) 
    df_F_code_JRC = get_sheet(dir_excel, 'Species_EEA', skiprows=1)
    df_F_code_JRC = compile_species(df_F_code_JRC)
    df_F_code_JRC = df_F_code_JRC.reset_index(drop=True)
    # df_F_code_added = get_sheet(dir_excel, 'Species_EEA_added')
    # df_F_code_added = df_F_code_added.reset_index(drop=True)
    # df_F_code = pd.concat([df_F_code_JRC, df_F_code_added], axis=0)
    df_F_code = df_F_code_JRC
    # df_F_code = df_F_code.reset_index(drop=True)
    df_F_code = df_F_code.sort_values(by=['Code'])

    # combine all the dataframes to get 
    # factors for every unique combination
    DATAFRAMES = {
        'CF': df_CF,
        'BCEF': df_BCEF,
        'RTS': df_RTS,
        'I': df_I,
        'FG': df_FG,
        'FGS': df_FGS,
        'Code': df_F_code
    } 

    df_C_seq_LUT = create_unique_comb_df(DATAFRAMES)

    return df_C_seq_LUT


def add_ms_species(settings):
    """
    Function that will add missing species from JRC 
    based on growth curves calculation 
    to the LUT.
    """
    # Open the JRC LUT based on yield information
    df_JRC_LUT = settings.get('JRC_LUT')
    
    # now open the excel file containing the yield
    # information for some additional species
    dir_excel = settings.get('CONFIG_SPECS').get('yield_table_ms_sp')
    # Net annual increment (I)
    df_I = get_sheet(dir_excel, 'NAI_add_species')
    # Get average, min and max
    df_I = df_I.rename(columns={'Avg': 'Iv_0-20',
                                'Min': 'Iv_0-20_min',
                                'Max': 'Iv_0-20_max'})
    # the increment is not specified per age
    df_I['Iv_21-30'] = df_I['Iv_0-20'] * df_I.shape[0]
    df_I['Iv_21-30_min'] = df_I['Iv_0-20_min'] * df_I.shape[0]
    df_I['Iv_21-30_max'] = df_I['Iv_0-20_max'] * df_I.shape[0]

    # Get per species and forest zone the yearly increment
    

def get_LUT_forest_zone(settings):
    # open first the shapefile of NUTS regions
    shp_NUTS = gpd.read_file(settings.get('VECTORS').get('NUTS'))

    # get sheet with forest zone information per country
    dir_excel = settings.get('CONFIG_SPECS').get('yield_table_dir')
    df_zones = get_sheet(dir_excel, 'FOREST_ZONE')

    lst_df_LUT_zones = []
    for CNTR in df_zones.CNTR_CODE.unique():
        NUTS3_regions = shp_NUTS.loc[((shp_NUTS.CNTR_CODE == CNTR) & (shp_NUTS.LEVL_CODE == 3))]['NUTS_ID'].to_list()
        ZONE = df_zones.loc[df_zones['CNTR_CODE'] == CNTR]['FOREST_ZONE'].values[0]
        df_NUTS_zone = pd.DataFrame(NUTS3_regions, columns=['NUTS_LEVEL3_ID'])
        df_NUTS_zone['NUTS_LEVEL0_ID'] = [CNTR] * df_NUTS_zone.shape[0]
        df_NUTS_zone['FOREST_ZONE'] = [ZONE] * df_NUTS_zone.shape[0]
        lst_df_LUT_zones.append(df_NUTS_zone)

    df_NUTS_LUT_final = pd.concat(lst_df_LUT_zones)

    return df_NUTS_LUT_final


def group_LUT_species(settings, dir_LUT_species,
                      grouplevel=['FGS', 'FOREST_ZONE']):
    # need to groupby forest zone and species type
    df_species_LUT = pd.read_csv(dir_LUT_species)
    # take the average per grouylevel
    df_avg_LUT_group = df_species_LUT.groupby(by=grouplevel).mean(numeric_only=True)
    df_avg_LUT_group = df_avg_LUT_group.reset_index()
    # keep only relevant columns
    cols_keep = [item for item in df_avg_LUT_group.columns.values
                 if 'yrl_C_seq' in item]
    df_avg_LUT_group = df_avg_LUT_group[grouplevel + cols_keep]
    return df_avg_LUT_group
    

def _main_(settings):

    # first define outfolder and outname 
    # and check if file not yet exists

    outfolder = settings.get('CONFIG_SPECS').get('outfolder')
    outname_C_seq_LUT = 'LUT_C_SEQ_AFFOR_JRC_V4.csv'
    outname_C_seq_LUT_FGS = 'LUT_C_SEQ_AFFOR_JRC_V4_GROUPED_FGS.csv'
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
    
    # C seq LUT for grouped species (F/S)
    if not os.path.exists(os.path.join(outfolder, outname_C_seq_LUT_FGS)) or overwrite:
        # function that will generate 
        # the annual C increase for the grouped species
        df_C_seq = group_LUT_species(settings,
                                     os.path.join(outfolder,
                                                  outname_C_seq_LUT))
        # write the corresponding dataframe 
        # to the output folder
        df_C_seq.to_csv(os.path.join(outfolder, outname_C_seq_LUT_FGS),
                        index=False)

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
    name_file = 'Growth_curves_JRC_Data_Collection_140324_EEA_V3.xlsx'
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