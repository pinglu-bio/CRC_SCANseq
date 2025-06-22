import json
import pandas as pd
import gzip

l = []
for input_file in snakemake.input:
    df = pd.read_table(input_file)
    l.append(df)

df = pd.concat(l)
df.to_csv(snakemake.output[0], index=False)
