import numpy as np
import pandas as pd
import os

from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt


if __name__ == '__main__':
    pa = r"/Volumes/Birdbox III/ATPsynthAS_20170307/Olig/2/641/0/"
    pathList = []
    for roots, dirs, files in os.walk(pa):
        for d in dirs:
            d0 = roots+ '//' +d
            L = os.listdir(d0)
            if 'FitResults.txt' in L:
                pathList.append(d0)

    #path1
    #pathList.append("") #path2
    #pathList.append("") #path3
    for path in pathList:
        os.chdir(path)
        fit = pd.read_table('FitResults.txt')
        F = np.array(zip(fit['X'],fit['Y']))
        db = DBSCAN(eps=3, min_samples=10).fit(F)
        labels = db.labels_
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
        print('Estimated number of clusters: %d' % n_clusters_)
        
        fit['Cluster'] = labels
        fit.to_csv(path + '/' + 'Results3.csv', sep = '\t')
        
    
