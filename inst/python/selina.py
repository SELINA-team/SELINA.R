from imblearn.over_sampling import SMOTE
import numpy as np
import glob
import pandas as pd
import ntpath
import datatable as dt
from functools import reduce
from collections import Counter


def read_expr(path):
    expr = dt.fread(path, header=True, sep='\t', nthreads=6)
    expr = expr.to_pandas()
    expr.index = expr.loc[:, 'Gene']
    del expr['Gene']
    expr = expr.astype(float)
    return expr
    

def delete_multiple_element(list_object, indices):
    indices = sorted(indices, reverse=True)
    for idx in indices:
        if idx < len(list_object):
            list_object.pop(idx)
    return list_object


def preprocessing(path_in, disease):
    samples = sorted(glob.glob(path_in + '/*_expr.txt'))
    metas = sorted(glob.glob(path_in + '/*_meta.txt'))
    train_sets = []
    celltypes = []
    if disease:
        diseases = []
    platforms = []
    genes = []
    print('Loading data')
    
    for i in range(len(samples)):
        train_sets.append(read_expr(samples[i]))
        meta = pd.read_csv(metas[i], sep='\t', header=0)
        celltypes.append(meta['Celltype'].to_numpy())
        platforms.append(meta['Platform'].to_numpy())
        genes.append(train_sets[i].index.to_list())
        if disease:
            diseases.append(meta['Disease'].to_numpy())

    genes = reduce(np.intersect1d, genes)

    length = len(train_sets)
    drop_index = []
    for i in range(length):
        sample_plat = np.unique(platforms[i])
        if len(sample_plat)>=2 :
            drop_index.append(i)
            for plat in sample_plat:
                index = np.where(platforms[i]==plat)[0]
                train_sets.append(train_sets[i].iloc[:,index])
                celltypes.append(celltypes[i][index])
                if disease:
                    diseases.append(diseases[i][index])
                platforms.append(platforms[i][index])
    train_sets = delete_multiple_element(train_sets,drop_index)
    celltypes = delete_multiple_element(celltypes,drop_index)
    if disease:
        diseases = delete_multiple_element(diseases,drop_index)
    platforms =delete_multiple_element(platforms,drop_index)
    
    if disease:
        length = len(train_sets)
        drop_index = []
        for i in range(length):
            disease_type = np.unique(diseases[i])
            if len(disease_type) >=2 :
                drop_index.append(i)
                for dis in disease_type:
                    index = np.where(np.array(diseases[i])==dis)[0]
                    train_sets.append(train_sets[i].iloc[:,index])
                    celltypes.append(celltypes[i][index])
                    diseases.append(diseases[i][index])
                    platforms.append(platforms[i][index])
        train_sets = delete_multiple_element(train_sets,drop_index)
        celltypes = delete_multiple_element(celltypes,drop_index)
        diseases = delete_multiple_element(diseases,drop_index)
        platforms =delete_multiple_element(platforms,drop_index)

    ct_freqs = Counter([i for item in celltypes for i in item])
    max_n = max(ct_freqs.values())
    rct_freqs = {}
    if  (max_n < 500) & (max_n>100):
        sample_n = 100
    elif max_n < 1000:
        sample_n = 500
    else:
        sample_n = 1000
    for ct,freq in ct_freqs.items():
        if freq <= sample_n:
            rct_freqs[ct] = freq

    for i in range(len(train_sets)):                                            
        sample_ct_freq = {}
        ct_freq = Counter(celltypes[i])
        if len(ct_freq)>1:
            for ct,freq in rct_freqs.items():
                if (ct in ct_freq.keys()) & (ct_freq[ct] >= 4):
                    sample_ct_freq[ct] = round(sample_n * ct_freq[ct]/freq)
            smo = SMOTE(sampling_strategy = sample_ct_freq,random_state=1, k_neighbors=3)
            train_sets[i],celltypes[i] = smo.fit_resample(train_sets[i].T,celltypes[i])
            train_sets[i] = train_sets[i].T
            platforms[i] = np.unique(platforms[i]).tolist() * train_sets[i].shape[1]
            if disease:
                diseases[i] = np.unique(diseases[i]).tolist() * train_sets[i].shape[1]
    
    celltypes = [i for item in celltypes for i in item]
    if disease:
        diseases = [i for item in diseases for i in item]
    platforms = [i for item in platforms for i in item]
    
    for i in range(len(train_sets)):
        train_sets[i] = train_sets[i].loc[genes, ]
        train_sets[i] = np.divide(train_sets[i], np.sum(train_sets[i],
                                                        axis=0)) * 10000
        train_sets[i] = np.log2(train_sets[i] + 1)
    train_data = pd.concat(train_sets, axis=1)
    if disease:
        return train_data, celltypes, diseases, platforms, genes
    else:
        return train_data, celltypes, platforms, genes