import json
import pandas as pd

l = []
for input_file in snakemake.input:
    with open(input_file, "r") as f:
        l.append(json.load(f))

df = pd.DataFrame.from_records(l)
df.to_csv(snakemake.output[0], index=False)
