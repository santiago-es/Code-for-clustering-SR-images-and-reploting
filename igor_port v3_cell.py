# -*- coding: utf-8 -*-
"""
Created on Fri Jul 28 16:18:10 2017

@author: DRW
"""
from __future__ import print_function
import numpy as np
import tifffile
import astropy
print(astropy.version.version)

from sklearn.neighbors import NearestNeighbors
from astropy.modeling.models import Gaussian2D
# from scipy.spatial.distance import pdist
# from scipy.spatial.distance import squareform
import sklearn
# from sklearn import cluster
# import itertools
import os
import pandas as pd
from pandas import Series
import copy
import skimage as sk
import scipy
import scipy.ndimage
import scipy.stats
import scipy.misc
import skimage.measure
import skimage.morphology
import timeit

Histo_precision = {}
Histo_NumLoc = {}
Histo_NNCluster = {}
Histo_NNLoc = {}
Histo_Length = {}
Histo_Eccen = {}
Histo_Density = {}


def user_input():
    """
    Input variables here. Shouldn't need to touch other functions.

    Known bugs:
        1. ///Fixed/// 20170805 - length analysis drops a cluster/FOV at some point.
        2. ///Discard///20170805 - the width=precision image may not be right. Not sure if
            the precision should be SD, FWHM or what. Currently treating as
            if it is SD.   # Eric - the all_fits_with_header.txt is supposed to input into imageJ and generate 
                                    this image using the GDSC SMLM results manager. The code is available in 
                                    'Dropbox/Shared Anlaysis/code/Generate_SR.py'. It needs to be run in ImageJ Jython platform
                                    
        3. ///Fixed/// 20170807 - Can't properly seperate the analysis (i.e. True/False) 
        4. ///Fixed/// 20170808 - fibril length analysis image is rotated
        5. ///Discard///20170808 - Summary of summary is wrong
        6. ///Fixed/// 20170809 - Memory overflown on large clusters (over 50000 localisations)
        
    Possible improvements:
        ///DONE/// Add averages of the length analysis to summary.txt
        ///DONE///Change count() to speed up the analysis.
         ///DONE/// Improve the summary part
         ///DONE////Add histogram generation
    """

    ############################################

    ############################################
    
    directory = r"H:\Dropbox (Cambridge University)\Shared Analysis\iPSC_fib_length"
    pixel_size = 131.5

    dbscan_analysis = [False, 3, 10] # epsilon, minimum_samples
    save_superres = False # Consider discard // I like it - much quicker than IJ for rough checking
    nearest_neighbor = True
    length_analysis = [True, 2, 2, 5] # threshold, gauss_sigma, cluster_size
    save_summary = True
    save_with_header = True
    save_histogram = True
    commands = [directory, pixel_size, dbscan_analysis, save_superres, 
                nearest_neighbor, length_analysis, save_summary, 
                save_with_header, save_histogram]
    master_func(commands)


def master_func(commands):  
    if commands[2][0] == True:
        print('\nCluster analysis...')
        cluster_analysis(commands[0], commands[2])
    
    print('\nAnalysing clusters...')
    file_name = 'Resultsthresh_mask_20_0.5.csv'
    paths = search_for_file(commands[0], file_name)

    summary = []
    count = 1
    for path in paths:
        print('\t{0}/{1}: {2}'.format(count, len(paths), path))
        count += 1

        file_path = os.path.join(path, file_name)
        real_locs = load_results(file_path, commands[2][1])
        extracted_data = easy_numbers(real_locs, path)

        if commands[3] == True:
            output_images(real_locs, commands[1], path)
        if commands[4] == True:             # NN analysis
            cwise_avg, cwise_std, lwise_avg, lwise_std = nearest_neighbour(real_locs, path)
            extracted_data = extracted_data + [cwise_avg, cwise_std, lwise_avg, lwise_std]
        else:
            extracted_data = extracted_data + [np.nan]*4  # Eric - somehow fill the array
        if commands[7] == True:
            headerise(path, commands[2][1])
        if commands[5][0] == True:                  # length analysis
            lavg, lsd, eavg, esd = lengths(path, commands[1], commands[5][1], commands[5][2], commands[5][3])
            extracted_data = extracted_data + [lavg, lsd, eavg, esd, extracted_data[4]/lavg]
        else:
            extracted_data = extracted_data + [np.nan]*5  # Eric - somehow fill the array
        summary.append(extracted_data)

    if commands[6] == True:
        summarise(commands[0], summary)
    if commands[8] == True:
        save_histograms(commands)

def save_histograms(command):
    directory = command[0]
    hist = os.path.join(directory, 'Histogram')
    if not os.path.isdir(hist):
        os.mkdir(hist)
    os.chdir(hist)

    pd.DataFrame(Histo_precision).to_csv('precision.csv', index = False) 
    pd.DataFrame(Histo_NumLoc).to_csv('NumLoc.csv', index = False)
    pd.DataFrame(Histo_NNCluster).to_csv('NNCluster.csv', index = False)
    pd.DataFrame(Histo_NNLoc).to_csv('NNLoc.csv', index = False)
    if command[5][0] == True: 
        for rows in Histo_NumLoc.keys():
            Histo_Density[rows] = Histo_NumLoc[rows].divide(Histo_Length[rows])
        pd.DataFrame(Histo_Length).to_csv('Length.csv', index = False)
        pd.DataFrame(Histo_Eccen).to_csv('Eccentricity.csv', index = False)
        pd.DataFrame(Histo_Density).to_csv('Density.csv', index = False)

def cluster_analysis(directory, dbscan_params):
    # Eric's dbscan
    epsilon = dbscan_params[1]
    minimum_samples = dbscan_params[2]
    path_list = search_for_file(directory, 'FitResults.txt')

    count = 1
    for path in path_list:
        print('\t{0}/{1}: {2}'.format(count, len(path_list), path))
        count += 1
        file_path = os.path.join(path, 'FitResults.txt')
        fitarray = pd.read_table(file_path) # Eric's madness.
        F = np.array(zip(fitarray['X'],fitarray['Y']))

        try:
            db = sklearn.cluster.DBSCAN(eps=epsilon, min_samples=minimum_samples).fit(F)
            labels = db.labels_

            fitarray['Cluster'] = labels
            fitarray.to_csv(os.path.join(path, 'Results'+str(epsilon)+'.csv'),
                       sep = '\t')
        except ValueError:
            pass # No localisations in FitResults, not important.


def search_for_file(path, desired_file):
    # Search a directory tree for files of a known name.
    directories = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if desired_file in name:
                directories.append(root)
    return directories


def load_results(file_path, epsilon):
    # Load results using numpy, sorry Eric
    f = np.loadtxt(open(file_path, "rb"), delimiter="\t", skiprows=1, 
                   usecols=[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17])
    real_locs = f[np.where(f[:,15] >= 0)]
    # 0 frame, 1 origx, 2 origy, 3 origvalue, 4 error, 5 noise, 6 SNR
    # 7 background, 8 signal, 9 angle, 10 x, 11 y, 12 xsd, 13 ysd, 14 precision
    # 15 cluster    
    return real_locs


def output_images(results, pixel_size, path):  # Consider discard
    # saves localisation only and width=precision images
    real_locs = results[np.where(results[:,15] >= 0)]
    path = os.path.join(path, 'SR images (Python)')
    if not os.path.exists(path):
        os.makedirs(path)
    
    # localisations only
    # just add 1 to an empty array for every localisation
    locs_only = np.zeros((4096, 4096))
    for row in real_locs:
        locs_only[int(row[11]*8)][int(row[10]*8)] += 1
    locs_only = locs_only.astype('uint16')
    tifffile.imsave(os.path.join(path, 'locs_only.tif'), locs_only)
    
    # width=precision
    #little but more complicated
    precision = np.zeros((4096, 4096))
    subpixel_size = pixel_size/8
    for row in real_locs:
        x = int(row[10]*8)
        y = int(row[11]*8)
        pp = row[14]/subpixel_size
        # astropy is used to generate the Gaussian kernel using this method
        # I am not sure if the precision in pixels (pp) is actually the SD??
        kernel = Gaussian2D(x_mean=5, y_mean=5, x_stddev=pp, y_stddev=pp)
        # the kernel is added to an 11x11 array. I think this should always
        # be big enough, but not sure
        grid = kernel(*np.mgrid[0:11, 0:11])
        # adding the grid to an array to build up the image one loc at a time
        precision[y-6:y+5, x-6:x+5] += grid
    precision = precision.astype('float32')
    #precision = precision[::-1,::-1]
    tifffile.imsave(os.path.join(path, 'prec_width.tif'), precision)


def easy_numbers(results, d):
    global Histo_precision
    # This just pulls the numbers out of Results3.txt that need no analysis
    print('\t\tChecking the abacus')
    loc_count = []
    cluster_num = int(np.amax(results[:,15]))
    easy = [d, cluster_num, np.round(np.mean(results[:,14]), decimals=3),
            np.round(np.std(results[:,14]), decimals=3)]
    Histo_precision[d] = Series(results[:,14])
    for x in xrange(0, cluster_num+1):
        loc_count.append(len(np.where(results[:,15] == x)[0]))
    Histo_NumLoc[d] = Series(loc_count)
    easy.append(np.mean(loc_count))
    easy.append(np.std(loc_count))
    return easy
    

def nearest_neighbour(real_locs, directory):
    # Nearest neighbour analysis, based on Mathew's code but much faster
    print('\t\tMeeting the neighbours')
    avg_nn_per_cluster = []
    nn_all_localisations = []
    all_loc_info = [['Clusternumber', 'xcoords', 'ycoords', 'Precision',
                     'NumNN']]
    
    for cluster_num in xrange(int(np.amax(real_locs[:,15]))+1):
        # Get x and y positions of each cluster
        cluster_details = real_locs[np.where(real_locs[:,15] == cluster_num)]
        cluster_xy = cluster_details[:, [10, 11]]
        precisions = cluster_details[:, 14]                         
        
#         if num_locs < 10000: # should test to optimise this number           
#             # use pdist
#             # this solution is fast for small clusters but too memory intensive for large ones
#             # first get all pairs of distances and convert to symmetrical array                              
#             all_distances = pdist(cluster_xy, metric='euclidean')
#             sq = squareform(all_distances)
#             
#             # find the shortest distance for each localisation
#             sq_sorted = np.sort(sq, axis=1)
#             try:
#                 thresh_distance = np.mean(sq_sorted[:,1])*5
#             except IndexError: # Not sure why this happens just yet
#                 print('\tIndex Error!')
#                 thresh_distance = 0
#         
#             # boolean array of near neighbors
#             sq2 = [(sq < thresh_distance)]
#         
#             # find the mean NN number for each localisation.
#             # -1 to remove self comparison
#             # cluster wise:
#             neighbors_per_loc = np.sum(sq2, axis=1)-1
#             avg_nn_per_cluster.append(np.mean(neighbors_per_loc))
#             # this section to save into all_localisation_information.txt
#             # neighbours_per_loc is still in order of cluster_details
#             for x in xrange(cluster_xy.shape[0]):
#                 all_loc_info.append([cluster_num, cluster_xy[x,:][0], 
#                                      cluster_xy[x,:][1], precisions[x],
#                                      neighbors_per_loc[0][x]])
#         
#             # localisation wise:
#             nn_all_localisations.append(neighbors_per_loc.tolist())
# 
#         
#         else:
        if len(cluster_xy):
            nbrs = NearestNeighbors(n_neighbors=2, algorithm='auto').fit(cluster_xy)   # Eric - change the NN analysis code to make it more time and space efficient
            distances, indices = nbrs.kneighbors(cluster_xy)
            distances = distances[:,1]
            thresh_distance = np.mean(distances)*5 # 5*NN dist average for the cluster
            
            neighbors_per_loc = []
            nbrs2 = NearestNeighbors(radius = thresh_distance, algorithm='auto').fit(cluster_xy)
            rng = nbrs2.radius_neighbors(cluster_xy)
            for i in rng[1]:
                neighbors_per_loc.append(len(i) - 1)
            neighbors_per_loc = np.array(neighbors_per_loc)

            avg_nn_per_cluster.append(np.mean(neighbors_per_loc))
            for x in xrange(cluster_xy.shape[0]):
                all_loc_info.append([cluster_num, cluster_xy[x,:][0], 
                                        cluster_xy[x,:][1], precisions[x],
                                        neighbors_per_loc[x]])
            for _ in neighbors_per_loc: nn_all_localisations.append(_)
    
    # save file for plotting later
    save_path = os.path.join(directory, 'All_localisation_information.txt')
    with open(save_path, 'w') as f:
        for row in all_loc_info:
            if len(row) != 5:
                print(row, len(row))
            f.write('{}\t{}\t{}\t{}\t{}\n'.format(*row))
    
    # un-nest all of the lists
    Histo_NNCluster[directory] =  Series(avg_nn_per_cluster)
    Histo_NNLoc[directory] = Series(nn_all_localisations)
    return np.mean(avg_nn_per_cluster), np.std(avg_nn_per_cluster), np.mean(nn_all_localisations), np.std(nn_all_localisations)


def lengths(path, pixel_size, thresh, gauss, csize):
    # Master func for the length analysis
    # This is the slowest part and could use attention
    print('\t\tGrabbing the ruler')
    spath = os.path.join(path, 'Length analysis')
    fp = 'All_Localisation_Information.txt'
    filepath = os.path.join(path, fp)
    skele, clusters = skeletonize(directory=spath, threshold=2, gauss_sigma=2,
                                  cluster_size=5, txtf=filepath)
    if skele.shape == (4096, 4096):
        final_counts, ecc = count(pixel_size, skele, clusters)
        save_data(final_counts, ecc, spath)
    Histo_Length[path] = Series(final_counts)
    Histo_Eccen[path] = Series(ecc)
    return np.mean(final_counts), np.std(final_counts), np.mean(ecc), np.std(ecc)


def skeletonize(directory, threshold, gauss_sigma, cluster_size, txtf):
    labels = {}
    arr = np.loadtxt(txtf, skiprows=1, usecols=[0,1,2])
    labelled = np.zeros((512*8, 512*8))#, dtype='uint16')
    try:
        max_val = int(np.amax(arr[:,0]+1))
    except IndexError:
        print('File contains no clusters!\n')
        return np.zeros((10, 10)), np.zeros((10, 10))
    
    for label in xrange(1, max_val + 1):    # Eric - set the range from 1 to max_value+1
        a = arr[np.where(arr[:,0] + 1 == label)]     # and then add 1 to the original label (which strats from 0)
        a = np.delete(a, np.s_[0], axis=1)
        a = a*8
        a = a.astype('int')             
        
        for xy in a:
            x = xy[0]
            y = xy[1]
            labelled[x][y] = label
            labels[(x, y)] = label
    
    clusters = labelled.astype('int')   
    clusters2 = copy.deepcopy(clusters)
    clusters3 = np.where(clusters2>0, 1, 0)
    
     
    skele1 = sk.morphology.skeletonize_3d(clusters3)
    skele1 = skele1*100
    skele1_blurred = scipy.ndimage.filters.gaussian_filter(skele1, gauss_sigma)
    skele1_blurred = skele1_blurred.astype(bool)
    skele2 = sk.morphology.skeletonize_3d(skele1_blurred)
    skele2 = skele2.astype('int16')      # Historical moment

    if not os.path.exists(directory):
        os.makedirs(directory)
    sp = os.path.join(directory, ('skele.png'))
    scipy.misc.imsave(sp, skele2.transpose())   # Eric- add transpose to rotate it back to the correct position


    ii = zip(*np.where(skele2 > 0))
    for i in ii:
        x, y = i
        flag = False
        for k in nearby(i):
            try:
                skele2[i] = labels[k]
                flag = True
                break
            except (KeyError, IndexError):
                pass
        if not flag: skele2[i] = 0     # Eric - Add nearby pixels to expand the search radius for the original
                                                            # label. Many of them are moved during the Gaussian blur
    return skele2, clusters

def nearby((x, y)):
    L = []
    for i in range(x-2, x+3):
        for j in range(y-2, y+3):
            L.append((i, j))
    return L

def count(pixel_size, skele, clusters):
    nm_lengths = []
    ecc = []
    max_val = np.amax(skele)
    props = skimage.measure.regionprops(clusters)
    
    for val in xrange(1, max_val+1):
        xy = np.where(skele==val)        
        xy = np.transpose(xy)
  
  #### Nearest neighbour method  --- same result but costs less space ####  - Eric
        length = 0
        if len(xy):
            nbrs = NearestNeighbors(radius = 1.5, algorithm='auto').fit(xy)
            rng = nbrs.radius_neighbors(xy)
            for i in rng[0]:
                length += sum(i)
            length = length/2 + 1
        else:
            length = 1
 #### 
 
 #### pdist method ####
        # dists = pdist(xy, metric='euclidean')          
        # dists = dists[dists < 2]
        # length = np.sum(dists) + 1
####

        try:
            ecc.append(props[val-1]['eccentricity'])
        except IndexError:
            pass
        nm_lengths.append(length/8.0*pixel_size)

    return nm_lengths, ecc


def save_data(final_counts, ecc, d): 
    if not os.path.exists(d):
        os.makedirs(d)  
    fi = os.path.join(d, 'lengths_and_eccentricity.txt')   # Eric - Removed '16bit'
    with open(fi, 'w') as f:
        for x1, x2 in zip(final_counts, ecc):
            f.write(str(x1) + '\t' + str(x2) + '\n')
            
            
def headerise(directory, num):
    # Slightly modified from Eric's code
    path_list = search_for_file(directory, 'Results{0}.csv'.format(num))
    for path in path_list:
        file_path = os.path.join(path, 'Results{0}.csv'.format(num))
        fit = pd.read_table(os.path.join(directory, file_path)) 
        fit = fit.drop(fit.columns[0], axis = 1)
        to_save = fit[fit.Cluster != -1]
        del to_save['Cluster']
        del to_save['Source']
        del to_save['SNR']
        to_save.to_csv(os.path.join(directory, 'All_fits_with_header_py.txt'),
                       sep = '\t',  index = False)
        header = "#Localisation Results File\n"+\
        "#FileVersion Text.D0.E0.V2\n\n"+\
        "#Name Image (LSE)\n"+\
        "#Source <gdsc.smlm.ij.IJImageSource><singleFrame>0</singleFrame><extraFrames>0</extraFrames><path>/Volumes/BIRDBOX/20161012_DNAPAINT_Tau_Fids_2/1/01.tiff</path></gdsc.smlm.ij.IJImageSource>\n"+\
        "#Bounds x0 y0 w512 h512\n"+\
        "#Calibration <gdsc.smlm.results.Calibration><nmPerPixel>131.5</nmPerPixel><gain>55.5</gain><exposureTime>50.0</exposureTime><readNoise>0.0</readNoise><bias>500.0</bias><emCCD>false</emCCD></gdsc.smlm.results.Calibration>\n"+\
        "#Configuration <gdsc.smlm.engine.FitEngineConfiguration><fitConfiguration><fitCriteria>LEAST_SQUARED_ERROR</fitCriteria><delta>1.0E-4</delta><initialAngle>0.0</initialAngle><initialSD0>2.0</initialSD0><initialSD1>2.0</initialSD1><computeDeviations>false</computeDeviations><fitSolver>LVM</fitSolver><minIterations>0</minIterations><maxIterations>20</maxIterations><significantDigits>5</significantDigits><fitFunction>CIRCULAR</fitFunction><flags>20</flags><backgroundFitting>true</backgroundFitting><notSignalFitting>false</notSignalFitting><coordinateShift>4.0</coordinateShift><signalThreshold>1665.0</signalThreshold><signalStrength>30.0</signalStrength><minPhotons>30.0</minPhotons><precisionThreshold>625.0</precisionThreshold><precisionUsingBackground>false</precisionUsingBackground><nmPerPixel>131.5</nmPerPixel><gain>55.5</gain><emCCD>false</emCCD><modelCamera>false</modelCamera><noise>0.0</noise><widthFactor>2.0</widthFactor><fitValidation>true</fitValidation><lambda>10.0</lambda><computeResiduals>false</computeResiduals><duplicateDistance>0.5</duplicateDistance><bias>500.0</bias><readNoise>0.0</readNoise><maxFunctionEvaluations>1000</maxFunctionEvaluations><searchMethod>POWELL</searchMethod><gradientLineMinimisation>false</gradientLineMinimisation><relativeThreshold>1.0E-6</relativeThreshold><absoluteThreshold>1.0E-16</absoluteThreshold></fitConfiguration><search>3.0</search><border>1.0</border><fitting>3.0</fitting><failuresLimit>10</failuresLimit><includeNeighbours>true</includeNeighbours><neighbourHeightThreshold>0.3</neighbourHeightThreshold><residualsThreshold>1.0</residualsThreshold><noiseMethod>QUICK_RESIDUALS_LEAST_MEAN_OF_SQUARES</noiseMethod><dataFilterType>SINGLE</dataFilterType><smooth><double>0.5</double></smooth><dataFilter><gdsc.smlm.engine.DataFilter>MEAN</gdsc.smlm.engine.DataFilter></dataFilter></gdsc.smlm.engine.FitEngineConfiguration>\n"+\
        "#"
        with open(os.path.join(directory,
                               'All_fits_with_header_py.txt'), 'r+') as log:
            content = log.read()
            log.seek(0)
            log.write(header+content)
    return


def summarise(directory, summary):              
    save_path = os.path.join(directory, 'summary.csv')
    summaryT = np.array(summary).transpose()
    
    to_save = pd.DataFrame(data = {'Filepath': Series(summaryT[0]),         # Use dataframe to make it more readable
    'NumCluster': Series(summaryT[1]), 
    'Precision_Ave': Series(summaryT[2]),
    'Precision_SD':Series(summaryT[3]),
    'NumLocalisation_Ave':Series(summaryT[4]),
    'NumLocalisation_SD':Series(summaryT[5]),
    'NN_cluster_Avg':Series(summaryT[6]),
    'NN_cluster_SD':Series(summaryT[7]),
    'NN_localisation_Avg':Series(summaryT[8]),
    'NN_localisation_SD':Series(summaryT[9]),
    'Length_Avg':Series(summaryT[10]),
    'Length_SD':Series(summaryT[11]),
    'Eccentricity_Avg':Series(summaryT[12]),
    'Eccentricity_SD':Series(summaryT[13]),
    'Density (Loc\\Length)':Series(summaryT[14])}, 
    columns = np.array(['Filepath', 'NumCluster', 'Precision_Ave', 'Precision_SD', 'NumLocalisation_Ave', 'NumLocalisation_SD',  'NN_cluster_Avg', 
    'NN_cluster_SD', 'NN_localisation_Avg', 'NN_localisation_SD,', 'Length_Avg', 'Length_SD',  'Eccentricity_Avg',  'Eccentricity_SD', 'Density (Loc\\Length)']))

    to_save.to_csv(save_path, index = False)

#     head = 'File,Clusters,Precision avg,Precision SD,Localisations avg,\          # Discard the summary of summary part because the 
                                                                                                                                                        # keywords depend on the filepath
#     Localisations SD,NN clusterwise_avg,NN clusterwise SD,\
#     NN localisationwise avg,NN localisationwise SD, Length avg, Length SD, \
#     Eccentricity avg, Eccentricity SD, Density (loc/length)\n'
# 
#     saver = ''
#     for i in xrange(len(summary[0])):
#         saver += '{},'
#     saver = saver[:-1]
#     saver += '\n'
#     with open(save_path, 'w') as f:
#         f.write(head)
#         for row in summary:
#             f.write(saver.format(*row))

#     dic = {}
#     for row in summary:
#         if row[0][:-1] not in dic.keys():
#             dic[row[0][:-2]] = [[row[1], row[2], row[4], row[6], row[8], row[10], row[12], row[14]]]
#         else:
#             dic[row[0][:-2]].append([row[1], row[2], row[4], row[6], row[8], row[10], row[12], row[14]])
# 
#     save_path = os.path.join(directory, 'summary of summary.csv')
#     head = 'File,Clusters avg,Precision avg,Localisations avg,NN clusterwise avg,\
#     NN localisationwise avg,Length avg,Eccentricity avg,Density (loc/length) avg,\
#     Clusters SD,Precision SD,Localisations SD,NN clusterwise SD,\
#     NN localisationwise SD,Length SD,Eccentricity SD,Density (loc/length) SD\n'
# 
#     with open(save_path, 'w') as f:
#         f.write(head)
# 
#     tosave = []
#     for k in dic.keys():
#         arr = np.array(dic[k])
#         means = np.mean(arr, axis=0)
#         std = np.std(arr, axis = 0)
#         summ = list(np.hstack((means, std)))
#         summ1 = [k]
#         summ1.extend(summ)
#         tosave.append(summ1)
#     
#     saver = ''
#     for i in xrange(len(tosave[0])):
#         saver += '{},'
#     saver += '\n'
# 
#     with open(save_path, 'a+') as f:
#         for row in tosave:
#             f.write(saver.format(*row))
    
    

if __name__ == "__main__":
    start_time = timeit.default_timer()
    user_input()
    print('\nRuntime: {} min'.format((timeit.default_timer() - start_time)/60))