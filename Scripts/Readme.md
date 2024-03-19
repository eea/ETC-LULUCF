ETC/CA scripts


* Scripts --> Python scripts -which should be used to start the process (the stored py-scirpts are used for calc. of the c-squestratjion potential by NUTS)




---------------------------------------------
**BIOMASS**
* **Afforestation on living biomass**: Scripts\Biomass\LB_increase_afforestation_grassland_cropland_IPCCTier2_V20231002.py 
    - (using Natura2000, ex-wetland,peatland,... data)
    - Updated should be done in this script and the helper function:
    - https://github.com/eea/ETC-LULUCF/blob/master/src/Biomass/utils/biom_helper_functions.py)
    * All ouputs and inputs are stored on the CWS: \\cwsfileserver\projects\lulucf\f02_data\carbon_model_data\
        * input\
        * output\
    - workaround:   
        - check inputfiles
            - CLC_ACC
            - DEM
            - EnvZones
            - EU-Trees4F (ns_sdms\prob\resampled\rcp45) 2035-2065-2095
            - HighNatureValueFarmland
            - N2K
            - PEATLAND
            - WETLAND
            - NUTS
        - check used LUT:
            - JRC-yield_table (LUT_C_SEQ_AFFOR_JRC_V4.csv)
            - JRC-yield_table (LUT_C_SEQ_AFFOR_JRC_V4_GROUPED_FGS.csv)
        - check LB_increase_afforestation_grassland_cropland_IPCCTier2_V20231002.py  and /src/Biomass/utils/biom_helper_functions.py
        - START 




--------------------------------------------






