#!/usr/bin/env bash

export SPARK_HOME=/opt/spark3_2_0/
export PATH="$SPARK_HOME/bin:$PATH"


/data/users/Public/bontek/envs/python/python38/bin/python3.8 setup.py bdist_wheel

PYSPARK_PYTHON=/data/users/Public/bontek/envs/python/python38/bin/python3.8
${SPARK_HOME}/bin/spark-submit --driver-memory 4G --executor-memory 1g --queue default \
 --conf spark.executor.memoryOverhead=8g \
 --conf spark.driver.memoryOverhead=8g \
 --conf spark.driver.cores=4 \
 --conf spark.shuffle.service.enabled=true --conf spark.dynamicAllocation.enabled=true \
 --conf spark.yarn.appMasterEnv.PYSPARK_PYTHON=$PYSPARK_PYTHON \
 --conf spark.yarn.appMasterEnv.PYSPARK_DRIVER_PYTHON=$PYSPARK_PYTHON \
 --conf spark.executorEnv.LD_LIBRARY_PATH="/data/users/Public/bontek/envs/python/python38/lib/" \
 --conf spark.yarn.appMasterEnv.LD_LIBRARY_PATH="/data/users/Public/bontek/envs/python/python38/lib/"  \
 --conf spark.executorEnv.PYSPARK_PYTHON=$PYSPARK_PYTHON \
 --conf spark.speculation=true \
 --conf spark.yarn.submit.waitAppCompletion=true \
 --conf spark.driver.maxResultSize=16g \
 --conf spark.yarn.appMasterEnv.PROJ_LIB="/data/users/Public/bontek/envs/python/python38/share/proj/" \
 --conf spark.executorEnv.PROJ_LIB="/data/users/Public/bontek/envs/python/python38/share/proj/" \
 --queue default \
 --master yarn --deploy-mode cluster --name SOC_LUT_STRATIFICATION \
 Scripts/SOC/SOC_LUC_LUT.py

