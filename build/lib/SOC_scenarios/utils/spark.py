import os


def get_spark_sql(local=False):
    # initialise sparkContext
    sc = get_spark_context(localspark=local)
    from pyspark.sql import SQLContext
    sqlContext = SQLContext(sc)
    return sqlContext


def get_spark_context(name="CROPCAR", localspark=False,
                      pythonpath='/data/users/Public/kristofvt/'
                      'python/worldcereal/bin/python'):

    import pyspark.serializers
    import cloudpickle
    import pyspark
    from pyspark.sql import SparkSession

    pyspark.serializers.cloudpickle = cloudpickle

    if not localspark:

        # Hot fix to make keras models pickable
        # fix_pickle()

        spark = SparkSession.builder \
            .appName(name) \
            .getOrCreate()
        sc = spark.sparkContext
    else:
        os.environ['PYSPARK_PYTHON'] = pythonpath
        spark = SparkSession.builder \
            .appName(name) \
            .master('local[1]') \
            .config('spark.driver.host', '127.0.0.1') \
            .config('spark.executor.memory', '2G') \
            .config('spark.driver.memory', '2G') \
            .getOrCreate()
        sc = spark.sparkContext

    # Set log level
    sc.setLogLevel("WARN")

    return sc
