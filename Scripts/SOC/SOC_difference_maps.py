""""
Script that will calculate the difference maps between the baseline/current and potential scenario
"""
import glob
import os
from pathlib import Path
import geopandas as gpd

#### Specify the extent of the SOC dataset. Use EEA39 or country abbrevation used in NUTS shapefile to filter on this

AOI = 'EEA39'
Basefolder_SOC_stats = r'L:\etc\lulucf\strata\SOC_NUTS_stats'
suffix_file = '.geojson' ## if only want to select certain versions of the stats
scenario_baseline = '1'
years_SOC_stable = 20  # the amount of years before SOC reach equilibrium

if AOI == 'EEA39':
    Basefolder_AOI = os.path.join(Basefolder_SOC_stats)
else:
    Basefolder_AOI = os.path.join(Basefolder_SOC_stats, AOI)

SOC_stats_files = glob.glob(os.path.join(Basefolder_AOI, f'*{suffix_file}'))
Baseline_file = [item for item in SOC_stats_files if ('baseline' in Path(item).stem) & ('Scenario_1' in Path(item).stem)]
if len(Baseline_file) >1:
    print("Multiple baseline files available --> don't know which to choose")
    assert Exception

Scenario_files = [item for item in SOC_stats_files if 'future' in Path(item).stem]


### now derive some stats of the difference between the baseline and the future scenario and assert it to the baseline
## file

gpd_baseline = gpd.read_file(Baseline_file[0])
Baseline_scenario_name = 'Scenario_1_baseline'
## create new gpd in which all this SOC difference output will be written
gpd_SOC_baseline_scenario_stats = gpd_baseline.copy(deep=True)
NUTS_regions = gpd_baseline.NUTS_ID.to_list()
for scenario in Scenario_files:
    scenario_name = '_'.join(Path(scenario).stem.split('_')[-3:])
    gpd_scenario = gpd.read_file(scenario)

    ## ensure the baseline and future scenario are sorted in the same way
    gpd_scenario = gpd_scenario.set_index('NUTS_ID')
    gpd_scenario = gpd_scenario.reindex(NUTS_regions)
    gpd_scenario = gpd_scenario.reset_index(drop=True)

    ### add now the desired columns to the baseline
    gpd_SOC_baseline_scenario_stats[f'SOC_mean_{Baseline_scenario_name}_[tonC/ha]'] = gpd_baseline['SOC_mean']
    gpd_SOC_baseline_scenario_stats[f'SOC_mean_{scenario_name}_[tonC/ha]'] = gpd_scenario['SOC_mean']
    gpd_SOC_baseline_scenario_stats[f'SOC_yrly_diff_{scenario_name}_[tonC/ha]'] = (gpd_SOC_baseline_scenario_stats[f'SOC_mean_{scenario_name}_[tonC/ha'] -
                                                                         gpd_SOC_baseline_scenario_stats[f'SOC_mean_{Baseline_scenario_name}_[tonC/ha]'])\
                                                                        /years_SOC_stable

    ### yearly percentual difference from reference SOC stock (absolute difference divided by SOC baseline and multiplied by 100
    gpd_SOC_baseline_scenario_stats[f'SOC_yrly_%_diff_{scenario_name}'] = (((gpd_SOC_baseline_scenario_stats[f'SOC_mean_{scenario_name}_[tonC/ha'] -
                                                                         gpd_SOC_baseline_scenario_stats[f'SOC_mean_{Baseline_scenario_name}_[tonC/ha]']) \
                                                                        /years_SOC_stable)/ \
                                                                          gpd_SOC_baseline_scenario_stats[f'SOC_mean_{Baseline_scenario_name}_[tonC/ha]']) *100


outname = 'SOC_baseline_VS_scenarios_stats_NUTS.geojson'
gpd_SOC_baseline_scenario_stats = gpd_SOC_baseline_scenario_stats.drop(columns=['SOC_mean'])
gpd_SOC_baseline_scenario_stats.to_file(os.path.join(Basefolder_AOI, outname), driver='GeoJSON')









