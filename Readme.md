ETC/CA LULUCF repository:

In this EEA repository the work done for the creation of a carbon sequestration potential map compatible with the LULUCF regulation will be shared.




Repo structure:
* Scripts --> Python scripts -which should be used to start the process (the stored py-scirpts are used for calc. of the c-squestratjion potential by NUTS)

* notebooks --> Showcase notebook on how to use satellite data for carbon sequestration mapping  (notebooks are not up to date!!)
* flowcharts --> Flowcharts that visualize the process used to within the notebooks for carbon sequestration mapping
* src --> in this folder and the subfolders the source code used for running the process are stored



---------------------------------------------
Main scripts:

* **Afforestation on living biomass**: Scripts\Biomass\LB_increase_afforestation_grassland_cropland_IPCCTier2_V20231002.py 
   (using Natura2000, ex-wetland,peatland,... data)
   Updated should be done in this script and the helper function: https://github.com/eea/ETC-LULUCF/blob/master/src/Biomass/utils/biom_helper_functions.py)
    * All ouputs and inputs are stored on the CWS: \\cwsfileserver\projects\lulucf\f02_data\carbon_model_data\
        * input\
        * output\
--------------------------------------------



