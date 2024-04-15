import numpy as np
from pycm import *

def PCA_decomposition(module_expression):
    from sklearn.decomposition import PCA
    clf = PCA(n_components=1)
    clf.fit(module_expression)
    #clf.components_[1,:]
    print('PCA dim:',len(clf.components_[0]))
    return clf.components_[0]

def Eigengene_significance(eigengene, sample_label):
    #print("to be completetd")#clf.components_
    from scipy import stats
    import math
    module_significance = stats.pearsonr(eigengene, sample_label)[0]
    return abs(module_significance)

