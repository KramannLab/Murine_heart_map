import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sb
import matplotlib.patches as mpatches

import warnings
warnings.filterwarnings("ignore")

import matplotlib.colors as colors
gray_red = colors.LinearSegmentedColormap.from_list("grouping", ["lightgray", "red", "darkred"], N = 128)


sc.settings.verbosity = 3
sc.logging.print_version_and_date()

def correlate_means_to_gene(means, corr_gene = "EOMES"):
    import scipy.stats
    genes = means.columns.values
    cors = pd.DataFrame(index = genes, columns = ["spearman_corr", "pvalue"])
    tab = means.loc[:, [corr_gene]].values

    for gene in genes:
        tmp = scipy.stats.spearmanr(tab, means.loc[:, [gene]])   # Spearman's rho    
        cors.loc[gene, :] = tmp[0:2]
    cors.dropna(axis = 0, inplace = True)
    cors.sort_values("spearman_corr", ascending = False, inplace = True)
    
    return cors
