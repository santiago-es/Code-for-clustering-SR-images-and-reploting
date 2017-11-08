# -*- coding: utf-8 -*-
"""
Created on Fri Jul 28 16:18:10 2017

@author: DRW
"""

from __future__ import print_function
import numpy as np
import tifffile
import astropy
from sklearn.neighbors import NearestNeighbors
from astropy.modeling.models import Gaussian2D
import sklearn
import os
import pandas as pd
from pandas import Series
import copy
import skimage as sk
import scipy
from sklearn import cluster
import scipy.ndimage
import scipy.stats
import scipy.misc
import skimage.measure
import skimage.morphology
import timeit
import warnings
import Tkinter as Tk
import scipy
import scipy.spatial
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Why do this to me, Eric?
Histo_precision = {}
Histo_NumLoc = {}
Histo_NNCluster = {} # number of neighbours cluster or localisation-wise
Histo_NNLoc = {}
Histo_Length = {}
Histo_Eccen = {}
Histo_Density = {}
Histo_NNdist = {} # nearest neighbour distance


def user_input():
    """
    Input variables here. Shouldn't need to touch other functions.

    Known bugs:
        1. ///Fixed///// 20170805 - length analysis drops a cluster/FOV at some point.
        2. ///Discard/// 20170805 - the width=precision image may not be right. Not sure if
                                    the precision should be SD, FWHM or what. Currently treating as
                                    if it is SD.   
                                    # Eric - the all_fits_with_header.txt is supposed to input into imageJ and generate 
                                    this image using the GDSC SMLM results manager. The code is available in 
                                    'Dropbox/Shared Anlaysis/code/Generate_SR.py'. It needs to be run in ImageJ Jython platform
                                    
        3. ///Fixed///// 20170807 - Can't properly seperate the analysis (i.e. True/False) 
        4. ///Fixed///// 20170808 - fibril length analysis image is rotated
        5. ///Discard/// 20170808 - Summary of summary is wrong
        6. ///Fixed///// 20170809 - Memory overflow on large clusters (over 50000 localisations)
        7. 20171020 - Need to double check the cropping is consistent

        
    Possible improvements:
        Colocalisation analysis - % colocalised & intensity of colocalised
        Remove edge clusters
        
        ///DONE/// Add ThT counter 
        ///DONE/// Add cropping function
        ///DONE/// Add averages of the length analysis to summary.txt
        ///DONE/// Change count() to speed up the analysis.
        ///DONE/// Improve the summary part
        ///DONE/// Add histogram generation
    
    something for anyone debugging - ph 116 123
    """
    
    ##### input #####
    directory = r"C:\Users\yz520\Desktop\OneDrive - University Of Cambridge\igorplay"
    pixel_size = 131.5
    
    crop = [False, (250,250), (500,500)] # Requires dbscan. (x,y) for top left and bottom right corners.
    ThT_counts = False
    dbscan_analysis = [True, 3, 10] # epsilon, minimum_samples
    save_superres = False
    nearest_neighbor = True
    length_analysis = [True, 2, 2, 5] # Requires nearest_neighbor. Threshold, gauss_sigma, cluster_size
    save_summary = True
    save_with_header = True
    save_histogram = True
    ##### /input #####

    commands = [directory, pixel_size, dbscan_analysis, save_superres, 
                nearest_neighbor, length_analysis, save_summary, 
                save_with_header, save_histogram, crop]

    if ThT_counts == True:
        file_name = 'ThT.tif'
        paths = search_for_unknown_file(commands[0], file_name)
        get_threshold(commands, paths[0], commands[9])
    else:
        root.destroy()
        commands = commands + [None]
        master_func(commands)
        

def master_func(commands):
    # ThT analysis
    if commands[10] != None:
        print('\nThT counts')
        tht_results = []
        file_name = 'ThT.tif'
        paths = search_for_unknown_file(commands[0], file_name)
        count = 1
        for path in paths:
            print('\t{0}/{1}: {2}'.format(count, len(paths), path))
            count += 1
            im = load_image(path, commands[9])
            maxima = find_maxima(im, commands[10])
            tht_results.append([path, maxima.shape[0]])
        save_tht(commands[0], tht_results)

    # SR analysis
    if commands[2][0] == True:
        # Suppress warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
        print('\nCluster analysis...')
        cluster_analysis(commands[0], commands[2], commands[9][0], commands[9][1], commands[9][2])
        
    print('\nAnalysing clusters...')
    file_name = 'DBSCAN_Results.txt'
    paths = search_for_file(commands[0], file_name)
    
    summary = []
    count = 1
    for path in paths:
        print('\t{0}/{1}: {2}'.format(count, len(paths), path))
        count += 1

        file_path = os.path.join(path, file_name)
        real_locs = load_results(file_path)

        if len(real_locs) > 0:
            extracted_data = easy_numbers(real_locs, path)

            if commands[3] == True:
                output_images(real_locs, commands[1], path)
            if commands[4] == True:             # NN analysis
                cwise_avg, cwise_std, lwise_avg, lwise_std, nndist_avg, nndist_std = nearest_neighbour(real_locs, path, commands[1])
                extracted_data = extracted_data + [cwise_avg, cwise_std, lwise_avg, lwise_std, nndist_avg, nndist_std]
            else:
                extracted_data = extracted_data + [np.nan]*6  # Eric - somehow fill the array
            if commands[7] == True:
                headerise(path)
            if commands[5][0] == True:          # length analysis
                lavg, lsd, eavg, esd = lengths(path, commands[1], commands[5][1], commands[5][2], commands[5][3])
                extracted_data = extracted_data + [lavg, lsd, eavg, esd, extracted_data[4]/lavg]
            else:
                extracted_data = extracted_data + [np.nan]*5  # Eric - somehow fill the array
        else:
            print('\t\tNo localisations')
            extracted_data = [path, 0]+[np.nan]*16
    
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
    pd.DataFrame(Histo_NNdist).to_csv('NNdist.csv', index = False)
    
    if command[5][0] == True: 
        for rows in Histo_NumLoc.keys():
            Histo_Density[rows] = Histo_NumLoc[rows].divide(Histo_Length[rows])
        pd.DataFrame(Histo_Length).to_csv('Length.csv', index = False)
        pd.DataFrame(Histo_Eccen).to_csv('Eccentricity.csv', index = False)
        pd.DataFrame(Histo_Density).to_csv('Density.csv', index = False)


def cluster_analysis(directory, dbscan_params, cropper, top_left, bottom_right):
    # suck it PANDAS
    epsilon = dbscan_params[1]
    minimum_samples = dbscan_params[2]
    path_list = search_for_file(directory, 'FitResults.txt')

    count = 1
    for path in path_list:
        print('\t{0}/{1}: {2}'.format(count, len(path_list), path))
        count += 1
        file_path = os.path.join(path, 'FitResults.txt')
        
        f = pd.read_table(file_path)
        f = f.drop(f.columns[0], axis = 1)
        # f = np.loadtxt(open(file_path, "rb"), delimiter="\t", skiprows=1, usecols=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])

        if cropper == True:
            dif = np.subtract(top_left, bottom_right)
            if (dif[0] > 0) or (dif[1] > 0) or (top_left[0] > 512) or (top_left[1] > 512):
                print('\t\tNot cropped. Check yo edge dawg')
            else:
                f = f[(f.X > top_left[0]) & (f.X< bottom_right[0]) & (f.Y > top_left[1]) & (f.Y  < bottom_right[1])]
        try:
            F = np.array(zip(f['X'],f['Y']))
            db = sklearn.cluster.DBSCAN(eps=epsilon, min_samples=minimum_samples).fit(F)
            labels = db.labels_
            f['Cluster'] = labels
            f.to_csv(os.path.join(path, 'DBSCAN_Results.txt'), index = False, sep = '\t')           
            # 20171106 - Eric: PANDAS dataframe keeps the datatype, ImageJ seems to only accept
            #  integer as Frame, X, Y coordinates
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


def search_for_unknown_file(path, desired_file):
    # Search a directory tree for files of a known name.
    paths = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if desired_file in name:
                paths.append(os.path.join(root, name))
    return paths


def load_results(file_path):
    # Load results using numpy, sorry Eric
    f = np.loadtxt(open(file_path, "rb"), delimiter="\t", skiprows=1)
    try:
        real_locs = f[np.where(f[:,15] >= 0)]
    except IndexError: # No localisations
        real_locs = np.zeros((1,16))
    
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
    #little bit more complicated
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
    # This just pulls the numbers out of DBSCAN_Results.txt that need no analysis
    print('\t\tChecking the abacus')
    loc_count = []

    cluster_num = int(np.amax(results[:,15]))
    easy = [d, cluster_num, len(results), np.round(np.mean(results[:,14]), decimals=3),
            np.round(np.std(results[:,14]), decimals=3)]
    Histo_precision[d] = Series(results[:,14])
    for x in xrange(0, cluster_num+1):
        loc_count.append(len(np.where(results[:,15] == x)[0]))
    Histo_NumLoc[d] = Series(loc_count)
    easy.append(np.mean(loc_count))
    easy.append(np.std(loc_count))
    return easy
    

def nearest_neighbour(real_locs, directory, pixel_size):
    # Nearest neighbour analysis, based on Mathew's code but much faster
    print('\t\tMeeting the neighbours')
    avg_nn_per_cluster = []
    nn_all_localisations = []
    nn_dist = []
    all_loc_info = [['Clusternumber', 'xcoords', 'ycoords', 'Precision',
                     'NumNN']]
    
    for cluster_num in xrange(int(np.amax(real_locs[:,15]))+1):
        # Get x and y positions of each cluster
        cluster_details = real_locs[np.where(real_locs[:,15] == cluster_num)]
        cluster_xy = cluster_details[:, [10, 11]]
        precisions = cluster_details[:, 14]                         
        
        try:
            nbrs = NearestNeighbors(n_neighbors=2, algorithm='auto').fit(cluster_xy)   # Eric - change the NN analysis code to make it more time and space efficient
            distances, indices = nbrs.kneighbors(cluster_xy)
            distances = distances[:,1]
            nn_dist.append(np.mean(distances)*pixel_size)
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
            
        # No localisations
        except ValueError: 
            all_loc_info.append([cluster_num, 0, 0, 0, 0])

    # save file for plotting later
    save_path = os.path.join(directory, 'All_localisation_information.txt')
    with open(save_path, 'w') as f:
        for row in all_loc_info:
            if len(row) != 5:
                print(row, len(row))
            f.write('{}\t{}\t{}\t{}\t{}\n'.format(*row))
    
    Histo_NNCluster[directory] =  Series(avg_nn_per_cluster)
    Histo_NNLoc[directory] = Series(nn_all_localisations)
    Histo_NNdist[directory] = Series(nn_dist)
    return np.mean(avg_nn_per_cluster), np.std(avg_nn_per_cluster), np.mean(nn_all_localisations), np.std(nn_all_localisations), np.mean(nn_dist), np.std(nn_dist)


def lengths(path, pixel_size, thresh, gauss, csize):
    # Master func for the length analysis
    # This is the slowest part and could use attention
    print('\t\tGrabbing the ruler')
    spath = os.path.join(path, 'Length analysis')
    fp = 'All_Localisation_Information.txt'
    filepath = os.path.join(path, fp)
    skele, clusters = skeletonize(directory=spath, threshold=2, gauss_sigma=2,
                                  cluster_size=5, txtf=filepath)
    if skele.shape:
        final_counts, ecc = count(pixel_size, skele, clusters)
        save_data(final_counts, ecc, spath)
    else:
        final_counts = [0]
        ecc = [0]
    Histo_Length[path] = Series(final_counts)
    Histo_Eccen[path] = Series(ecc)
    return np.mean(final_counts), np.std(final_counts), np.mean(ecc), np.std(ecc)


def skeletonize(directory, threshold, gauss_sigma, cluster_size, txtf):
    labels = {}
    arr = np.loadtxt(txtf, skiprows=1, usecols=[0,1,2])
    labelled = np.zeros((512*8, 512*8))
    try:
        max_val = int(np.amax(arr[:,0]+1))
    except IndexError:
        print('\t\tFile contains no clusters!')
        return np.zeros((10, 10)), np.zeros((10, 10))
    
    for label in xrange(1, max_val + 1):             # Eric - set the range from 1 to max_value+1
        a = arr[np.where(arr[:,0] + 1 == label)]     # and then add 1 to the original label (which starts from 0)
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
        xy = np.transpose(xy)       # 20171106 - Eric: Someone threw away the transpose?

        # Nearest neighbour method  --- same result but costs less space ####  - Eric
        length = 0
        nbrs = NearestNeighbors(radius = 1.5, algorithm='auto').fit(xy)
        rng = nbrs.radius_neighbors(xy)
        for i in rng[0]:
            length += sum(i)
        length = length/2 + 1

        try:
            ecc.append(props[val-1]['eccentricity'])
        except IndexError:
            pass
        nm_lengths.append(length/8.0*pixel_size)

    return nm_lengths, ecc


def save_data(final_counts, ecc, d): 
    if not os.path.exists(d):
        os.makedirs(d)  
    fi = os.path.join(d, 'lengths_and_eccentricity.txt')
    with open(fi, 'w') as f:
        for x1, x2 in zip(final_counts, ecc):
            f.write(str(x1) + '\t' + str(x2) + '\n')


def headerise(directory):
    # Slightly modified from Eric's code
    path_list = search_for_file(directory, 'DBSCAN_Results.txt')
    for path in path_list:
        file_path = os.path.join(path, 'DBSCAN_Results.txt')
        fit = pd.read_table(os.path.join(directory, file_path)) 
        to_save = fit[fit.Cluster != -1]
        del to_save['Cluster']
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
        
    # Use dataframe to make it more readable
    to_save = pd.DataFrame(data = {'Filepath': Series(summaryT[0]),         
    'NumCluster': Series(summaryT[1]),
    'Total_Locs': Series(summaryT[2]),
    'Precision_Ave': Series(summaryT[3]),
    'Precision_SD':Series(summaryT[4]),
    'NumLocalisation_Ave':Series(summaryT[5]),
    'NumLocalisation_SD':Series(summaryT[6]),
    'NN_cluster_Avg':Series(summaryT[7]),
    'NN_cluster_SD':Series(summaryT[8]),
    'NN_localisation_Avg':Series(summaryT[9]),
    'NN_localisation_SD':Series(summaryT[10]),
    'NN_distance_avg':Series(summaryT[11]),
    'NN_distance_SD':Series(summaryT[12]),
    'Length_Avg':Series(summaryT[13]),
    'Length_SD':Series(summaryT[14]),
    'Eccentricity_Avg':Series(summaryT[15]),
    'Eccentricity_SD':Series(summaryT[16]),
    'Density (Loc\\Length)':Series(summaryT[17])},
    columns = np.array(['Filepath', 'NumCluster', 'Total_Locs', 'Precision_Ave', 
                        'Precision_SD', 'NumLocalisation_Ave', 
                        'NumLocalisation_SD',  'NN_cluster_Avg', 
                        'NN_cluster_SD', 'NN_localisation_Avg', 
                        'NN_localisation_SD', 'NN_distance_avg',
                        'NN_distance_SD', 'Length_Avg', 'Length_SD',
                        'Eccentricity_Avg',  'Eccentricity_SD',
                        'Density (Loc\\Length)'])) 

    to_save.to_csv(save_path, index = False)


def save_tht(directory, results):
    with open(os.path.join(directory, 'tht_counts.txt'), 'w') as f:
        for row in results:
            f.write('{0}\t{1}\n'.format(row[0][:-4], row[1]))


def get_threshold(commands, im_path, cropper):  
    image = load_image(im_path, cropper)
    generate_window(commands, image)


def load_image(path, cropper):
    '''Returns a flattened image.'''
    im = tifffile.imread(path)
    im = np.mean(im, axis=0) # change here for sum, max etc.
    im = im.astype('uint16') # never repeat the horrors of the past
    if cropper[0] == True:
        top_left = cropper[1]
        bottom_right = cropper[2]
        dif = np.subtract(top_left, bottom_right)
        if (dif[0] > 0) or (dif[1] > 0) or (top_left[0] > 512) or (top_left[1] > 512):
            print('\t\tNot cropped. Check yo edge dawg')
        else:
            # pretty sure this is the correct cropping (?)
            im = im[top_left[1]:bottom_right[1]+1, top_left[0]:bottom_right[0]+1]
    return im


def generate_window(commands, image):
    '''Set up the layout for the GUI'''
    # need user input for the threshold, start with 1000 as default
    entry = {}
    txt = Tk.StringVar()
    txt.set('1000')
    thresh = Tk.Entry(master=root, textvariable=txt)
    thresh.grid(row=1, column=0, sticky='NSE')
    # save this to get the user defined threshold later
    entry['thresh'] = [thresh, txt]

    # button to update the maxima overlay
    prev_button = Tk.Button(master=root, text='Preview', command=lambda: fill_canvases(image, entry))
    prev_button.grid(row=1, column=1, sticky='NSW')

    # button to move forward with analysis
    ok_button = Tk.Button(master=root, text='Go', command=lambda: call_mf(commands, entry))
    ok_button.grid(row=1, column=3, columnspan=2, sticky='NSW')

    fill_canvases(image, entry)
    Tk.mainloop()


def fill_canvases(image, entry):
    '''Make two canvases, show the image in both but overlay maxima on one only'''
    for n in xrange(2):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.imshow(image)
        ax.grid(b=True, which='major', color='white',linestyle=':')
        xy = find_maxima(image, entry)
        
        # overlay maxima
        if n == 1:
            ax.plot(xy[:, 1], xy[:, 0], 'w+', markersize=5)
    
        # add plots to the canvas, dont open plots separately
        canvas = FigureCanvasTkAgg(fig, master=root)
        plt.close()            
        canvas.get_tk_widget().grid(row=0, column=n*2, columnspan=2, sticky='NSEW')


def find_maxima(data, entry):
    '''Returns a list of xy coordinates for maxima in an image.'''
    neighborhood_size = 5

    if type(entry) == str:
        threshold = int(entry)
    else:
        try:
            threshold = int(entry['thresh'][1].get())
        except ValueError:
            # maybe some smartarse put in a letter or something
            entry['thresh'][1].set('1000')
            threshold = 1000

    # find the spots
    data_max = filters.maximum_filter(data, neighborhood_size)
    maxima = (data == data_max)
    data_min = filters.minimum_filter(data, neighborhood_size)
    diff = ((data_max - data_min) > threshold)
    maxima[diff == 0] = 0
    
    # count the spots
    labeled, num_objects = ndimage.label(maxima)
    slices = ndimage.find_objects(labeled)
    x, y = [], []
    for dy,dx in slices:
        x_center = (dx.start + dx.stop - 1)/2
        x.append(int(x_center))
        y_center = (dy.start + dy.stop - 1)/2    
        y.append(int(y_center))
    
    # return yx positions of the spots
    return np.stack((y, x), axis=1)
    

def call_mf(commands, entry):
    # we now have the appropriate threshold
    thresh = entry['thresh'][1].get()
    root.destroy()
    master_func(commands + [thresh])



if __name__ == "__main__":
    start_time = timeit.default_timer()
    root = Tk.Tk()
    root.wm_title('Find maxima - ThT')
    user_input()
    print('\nRuntime: {} min'.format((timeit.default_timer() - start_time)/60))