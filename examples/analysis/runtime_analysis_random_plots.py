import scanpy as sc
from AnnSQL import AnnSQL
import time
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import gc
from scipy.stats import ttest_ind_from_stats
import scipy.stats as stats
from scipy.stats import ttest_rel

#declare the the comparisons placeholder
comparisons = None

#load all comparisons datasets like above and merge into one
for i in range(1,7):
	if comparisons is None:
		comparisons = pd.read_csv(f'../results/comparisons_random_{i}.csv')
	else:
		comparisons = pd.concat([comparisons, pd.read_csv(f'../results/comparisons_random_{i}.csv')])

#set the colors of the plots (ansql1, anndata in-mem, ansql2, anndata backed)
colors = ["#07b88e", "#a4a6a4", "#07b88e", "#a4a6a4","#FFA067"]

#plot aggregation of all filters runtime
sns.set(style="whitegrid")
plt.figure(dpi=1200)
sns.lineplot(data=comparisons, x='size', y='runtime_log', hue='method', palette=colors, errorbar="ci")
plt.ylabel("")
plt.xlabel("Library Size", fontsize=18, fontweight='bold')
plt.xticks(rotation=50)
plt.xticks(np.arange(0, 275000, 25000))
plt.xlim(0, 250000)
plt.yticks(np.log([0.01,0.1, 1, 10, 100, 1000]), [0.01,0.1, 1, 10, 100, 1000])
plt.xticks(fontsize=18, fontweight='bold')
plt.yticks(fontsize=18, fontweight='bold')	
plt.grid(False)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
# plt.legend().remove()
sns.despine(top=True)
plt.gca().spines['left'].set_linewidth(3)
plt.gca().spines['left'].set_color('black')
plt.gca().spines['bottom'].set_linewidth(3)
plt.gca().spines['bottom'].set_color('black')
plt.title("Runtime Comparisons", fontsize=22, fontweight='bold')
plt.text(72000, 6, "(6) Filters vs. Queries", fontsize=16)
plt.show()

#statistical analysis between the two types of AnnData
comparisons.groupby('method')['runtime'].describe()

# Number of comparisons (one per dataset size)
num_comparisons = len(comparisons['size'].unique())

# Adjusted alpha using Bonferroni correction
alpha_corrected = 0.05 / num_comparisons  

print("AnnSql On-Disk vs AnnData Backed (Paired t-test with Bonferroni Correction)")
for size in comparisons['size'].unique():
    annsql_backed = comparisons[(comparisons['method'] == 'AnnData Backed') & (comparisons['size'] == size)]['runtime']
    anndata_ondisk = comparisons[(comparisons['method'] == 'AnnSql On-Disk') & (comparisons['size'] == size)]['runtime']
    if len(annsql_backed) == len(anndata_ondisk):  
        t_stat, p_value = ttest_rel(annsql_backed, anndata_ondisk)
        significant = p_value < alpha_corrected
        print(f"Size: {size}, T-statistic: {t_stat}, P-value: {p_value}, "
              f"Corrected Alpha: {alpha_corrected}, Significant: {significant}")
    else:
        print(f"Size: {size}, Skipping due to unequal sample sizes.")
	
print("\nAnnSql In-Memory vs AnnData Not-Backed (Paired t-test with Bonferroni Correction)")
for size in comparisons['size'].unique():
    annsql_notbacked = comparisons[(comparisons['method'] == 'AnnData Not-Backed') & (comparisons['size'] == size)]['runtime']
    anndata_inmemory = comparisons[(comparisons['method'] == 'AnnSql In-Memory') & (comparisons['size'] == size)]['runtime']
    if len(annsql_notbacked) == len(anndata_inmemory):  
        t_stat, p_value = ttest_rel(annsql_notbacked, anndata_inmemory)
        significant = p_value < alpha_corrected
        print(f"Size: {size}, T-statistic: {t_stat}, P-value: {p_value}, "
              f"Corrected Alpha: {alpha_corrected}, Significant: {significant}")
    else:
        print(f"Size: {size}, Skipping due to unequal sample sizes.")


#Seurat vs AnnData In-Memory (Paired t-test with Bonferroni Correction)
print("\nAnnSql In-Memory vs Seurat (Paired t-test with Bonferroni Correction)")
for size in comparisons['size'].unique():
	seurat = comparisons[(comparisons['method'] == 'Seurat') & (comparisons['size'] == size)]['runtime']
	anndata_inmemory = comparisons[(comparisons['method'] == 'AnnSql In-Memory') & (comparisons['size'] == size)]['runtime']
	if len(seurat) == len(anndata_inmemory):  
		t_stat, p_value = ttest_rel(seurat, anndata_inmemory)
		significant = p_value < alpha_corrected
		print(f"Size: {size}, T-statistic: {t_stat}, P-value: {p_value}, "
			  f"Corrected Alpha: {alpha_corrected}, Significant: {significant}")
	else:
		print(f"Size: {size}, Skipping due to unequal sample sizes.")



#mean runtime for each type with a 250,000 cell library 
comparisons[comparisons['size'] == 250000].groupby('method')['runtime'].mean()

#comparision 1
plt.figure(dpi=950)
sns.lineplot(data=comparisons[comparisons['filter'] == 'Filter1'], x='size', y='runtime_log', hue='method',linewidth=4, palette=colors)
plt.ylabel("")
plt.xlabel("")
plt.yticks(np.log([0.01,0.1, 1, 10,100]), [0.01,0.1, 1, 10,100])
plt.xticks(rotation=45)
plt.legend().remove()
plt.grid(False)
plt.xticks([])

#comparision 2
plt.figure(dpi=950)
sns.lineplot(data=comparisons[comparisons['filter'] == 'Filter2'], x='size', y='runtime_log', hue='method',linewidth=4, palette=colors)
plt.ylabel("")
plt.xlabel("")
plt.yticks(np.log([0.01,0.1, 1, 10,100]), [0.01,0.1, 1, 10,100])
plt.xticks(rotation=45)
plt.legend().remove()
plt.grid(False)
plt.xticks([])

#comparision 3
plt.figure(dpi=950)
sns.lineplot(data=comparisons[comparisons['filter'] == 'Filter3'], x='size', y='runtime_log', hue='method',linewidth=4, palette=colors)
plt.ylabel("")
plt.xlabel("")
plt.yticks(np.log([0.01,0.1, 1, 10,100]), [0.01,0.1, 1, 10,100])
plt.xticks(rotation=45)
plt.legend().remove()
plt.grid(False)
plt.xticks([])

#comparision 4
plt.figure(dpi=950)
sns.lineplot(data=comparisons[comparisons['filter'] == 'Filter4'], x='size', y='runtime_log', hue='method',linewidth=4, palette=colors)
plt.ylabel("")
plt.xlabel("")
plt.yticks(np.log([0.01,0.1, 1, 10,100]), [0.01,0.1, 1, 10,100])
plt.xticks(rotation=45)
plt.legend().remove()
plt.grid(False)
plt.xticks([])

#comparision 5
plt.figure(dpi=950)
sns.lineplot(data=comparisons[comparisons['filter'] == 'Filter5'], x='size', y='runtime_log', hue='method',linewidth=4, palette=colors)
plt.ylabel("")
plt.xlabel("")
plt.yticks(np.log([0.01,0.1, 1, 10,100]), [0.01,0.1, 1, 10,100])
plt.xticks(rotation=45)
plt.legend().remove()
plt.grid(False)
plt.xticks([])

#comparision 6
plt.figure(dpi=950)
sns.lineplot(data=comparisons[comparisons['filter'] == 'Filter6'], x='size', y='runtime_log', hue='method',linewidth=4, palette=colors)
plt.ylabel("")
plt.xlabel("")
plt.yticks(np.log([0.01,0.1, 1, 10,100]), [0.01,0.1, 1, 10,100])
plt.xticks(rotation=45)
plt.legend().remove()
plt.grid(False)