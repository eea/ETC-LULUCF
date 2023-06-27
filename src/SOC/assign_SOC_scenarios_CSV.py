""""
Script that will be used to assign the output of the SOC stats to the CSV table of NUTS regions used to upload
on the tableau environment of EEA
"""

import pandas as pd
import geopandas as gpd
import numpy as np
import datetime
import os

df_EEA_NUTS = pd.read_csv(r"L:\etc\lulucf\AOI\EEA\admlayer_eea38plus_v_june_126_65.csv")
gpd_NUTS_STATS = gpd.read_file(r"L:\etc\lulucf\strata\SOC_NUTS_stats\SOC_baseline_VS_scenarios_stats_NUTS.geojson")
outfolder = r'L:\etc\lulucf\strata\SOC_NUTS_stats\EEA_format'
os.makedirs(outfolder, exist_ok=True)
outname = f'IPCC_SOC_baseline_VS_scenarios_V_{str(datetime.datetime.now().date())}.csv'
### drop redundant columns:
columns = list(gpd_NUTS_STATS.columns.values)
columns_valid = [item for item in columns if not 'factor' in item and not 'src' in item]
columns_valid = [item for item in columns_valid if not 'geometry' in item and not 'NUTS_LEVEL' in item]
gpd_NUTS_STATS = gpd_NUTS_STATS[columns_valid]

lst_stats_NUTS = []
for i, NUTS_stats in gpd_NUTS_STATS.iterrows():
    NUTS_ID = NUTS_stats.NUTS_ID
    ### check if NUTS ID IN EEA TABLE
    if not NUTS_ID in df_EEA_NUTS.NUTS_EU.to_list():
        continue
    NUTS_stats_filtered = pd.DataFrame(NUTS_stats).T
    lst_stats_NUTS.append(NUTS_stats_filtered)

df_NUTS_stats_SOC = pd.concat(lst_stats_NUTS)

df_EEA_NUTS = df_EEA_NUTS.rename(columns={'NUTS_EU': 'NUTS_ID'})
df_merged_stats_EEA_format = pd.merge(df_EEA_NUTS, df_NUTS_stats_SOC, on='NUTS_ID')
df_merged_stats_EEA_format = df_merged_stats_EEA_format.rename(columns={'NUTS_ID': 'NUTS_EU'})
df_merged_stats_EEA_format['area_NUTS_[ha]'] = [np.round(item,2) for item in df_merged_stats_EEA_format['area_NUTS_[ha]'].values if not np.isnan(item)]
df_merged_stats_EEA_format['%_cover_IPCC_cat'] = [np.round(item,2) for item in df_merged_stats_EEA_format['%_cover_IPCC_cat'].values if not np.isnan(item)]
df_merged_stats_EEA_format.to_csv(os.path.join(outfolder, outname), index=False, encoding="utf-8-sig")
