import pandas as pd
import numpy as np

#load the local laptop atlas profile results
atlas_profile = pd.read_csv('../results/atlas_profile.csv')

#group by function and get the mean and standard error
atlas_profile = atlas_profile.groupby('function').agg(['mean', 'sem','std']).round(2)

#format the output
for index, row in atlas_profile.iterrows():
	print(f"{index}: Memory: {row['max_memory']['mean']} ±{row['max_memory']['std']}")
	print(f"{index}: Runtime: {row['runtime']['mean']} ±{row['runtime']['std']}")

print("=========================================================")

#load the atlas profile results form the HPC
atlas_profile_hpc = pd.read_csv('../results/atlas_profile_hpc.csv')

#group by function and get the mean and standard error
atlas_profile_hpc = atlas_profile_hpc.groupby('function').agg(['mean', 'sem','std']).round(2)

#format the output
for index, row in atlas_profile_hpc.iterrows():
	print(f"{index}: Memory: {row['max_memory']['mean']} ±{row['max_memory']['std']}")
	print(f"{index}: Runtime: {row['runtime']['mean']} ±{row['runtime']['std']}")
