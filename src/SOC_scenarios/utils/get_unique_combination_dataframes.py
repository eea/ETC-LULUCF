"""CODE USED TO COMBINE TWO DATAFRAMES TO PROVIDE THE UNIQUE COMBINATION
OF THEM IN A NEW DATAFRAME"""

import pandas as pd
from itertools import product
import os

def combine_dataframes(df1, df2):
    """
    Combine two dataframes in such a way that a new dataframe is created with
    all possible combinations of the columns in two dataframes
    :param df1:
    :param df2:
    :return:
    """
    columns_df1 = list(df1.columns.values)
    columns_df2 = list(df2.columns.values)

    df1_drop_duplicated = df1.drop_duplicates(subset = [columns_df1[0]])
    df2_drop_duplicated = df2.drop_duplicates(subset = [columns_df2[0]])


    grouped_columns_1 = list(product(df1_drop_duplicated[columns_df1[0]], df2_drop_duplicated[columns_df2[0]]))
    grouped_columns_2 = list(product(df1_drop_duplicated[columns_df1[1]], df2_drop_duplicated[columns_df2[1]]))

    df_combined_1 = pd.DataFrame(grouped_columns_1, columns = [columns_df1[0], columns_df2[0]])
    df_combined_2 = pd.DataFrame(grouped_columns_2, columns = [columns_df1[1], columns_df2[1]])
    df_combined_final = pd.concat([df_combined_1, df_combined_2], axis = 1)

    return df_combined_final

LUT_SOCREF = False
LUT_FLU = False
LUT_FMG = False
LUT_FI = True
IPCC_land_use_class = 'Grassland'

if LUT_SOCREF:
    df_soil = pd.read_csv(r"L:\etc\lulucf\strata\SOC_LUT\IPCC_soil_type_zones_LUT.csv", sep = ';')
    df_climate = pd.read_csv(r"L:\etc\lulucf\strata\SOC_LUT\IPCC_climate_zones_LUT.csv", sep = ';')
    outdir = r'L:\etc\lulucf\strata\SOC_LUT'
    outname = 'SOCref_LUT_SOIL_CLIMATE_tmp.xlsx'

    df_combined = combine_dataframes(df_soil, df_climate)
    df_combined.to_excel(os.path.join(outdir,outname), index = False)

if LUT_FLU:
    df_climate = pd.read_csv(r"L:\etc\lulucf\strata\SOC_LUT\IPCC_climate_zones_LUT.csv", sep = ';')
    df_FLU = pd.read_csv(r"L:\etc\lulucf\strata\SOC_LUT\IPCC_FLU_CLC_mapping_LUT.csv", sep=';')
    outdir = r'L:\etc\lulucf\strata\SOC_LUT'
    outname = r'IPCC_LUT_factors_FLU_tmp.xlsx'

    df_combined = combine_dataframes(df_FLU, df_climate)
    df_combined.to_excel(os.path.join(outdir,outname),index=False)

if LUT_FMG:
    df_climate = pd.read_csv(r"L:\etc\lulucf\strata\SOC_LUT\IPCC_climate_zones_LUT.csv", sep = ';')
    df_FMG = pd.read_csv(rf"L:\etc\lulucf\strata\SOC_LUT\IPCC_FMG_mapping_{IPCC_land_use_class}.csv", sep=';')
    outdir = r'L:\etc\lulucf\strata\SOC_LUT'
    outname = rf'IPCC_LUT_factors_FMG_tmp_{IPCC_land_use_class}.xlsx'
    df_combined = combine_dataframes(df_FMG, df_climate)
    df_combined.to_excel(os.path.join(outdir,outname),index=False)

if LUT_FI:
    df_climate = pd.read_csv(r"L:\etc\lulucf\strata\SOC_LUT\IPCC_climate_zones_LUT.csv", sep = ';')
    df_FMG = pd.read_csv(rf"L:\etc\lulucf\strata\SOC_LUT\IPCC_FI_mapping_{IPCC_land_use_class}.csv", sep=';')
    outdir = r'L:\etc\lulucf\strata\SOC_LUT'
    outname = rf'IPCC_LUT_factors_FI_tmp_{IPCC_land_use_class}.xlsx'
    df_combined = combine_dataframes(df_FMG, df_climate)
    df_combined.to_excel(os.path.join(outdir,outname),index=False)
