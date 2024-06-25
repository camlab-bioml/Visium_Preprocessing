# Import pythonic packages
import os
import pandas as pd

# User Specified Configurations
configfile: "config.yaml"
container: config['container']

output = 'output/' + config['v'] + '/'

sample_ids = os.listdir("data/visium/")
sce = ["Peng", "Zhou"]

# Rules
include: "pipeline/process-data.smk"

# Outputs
rule all:
    input: 
        process_data.values(),
        spots.values(),
        cell2loc_train.values(),
        cell2loc_predict.values()