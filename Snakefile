
import pandas as pd
configfile: "config.yaml"

output = 'output/' + config['v'] + '/'

samples_df = pd.read_csv(config['samples_file'])
sample_ids = samples_df.sample_name.to_list()

include: "pipeline/cellranger.smk"
include: "pipeline/process-data.smk"


rule all:
    input: 
        cellranger.values(),
        process_data.values()