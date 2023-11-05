#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  5 11:25:39 2023

@author: benozzo
"""
import numpy as np

def get_idx_mouse():
    idx_noHem = np.array([5,42,
                          0,2,3,37,39,40,
                          4,12,13,14,22,41,49,50,51,59,
                          6,19,20,21,43,56,57,58,
                          15,16,17,18,52,53,54,55,
                          1,7,8,9,10,11,38,44,45,46,47,48
                          ], dtype=int)
    separation = np.array([1, 4, 9, 13, 17, 23], dtype=float)*2 -1
    xtick_label = ['Vis', 'SM', 'antL', 'TEMP', 'MED', 'PF'] #['PF', 'MED', 'TEMP', 'antL', 'SM', 'Vis']
    return idx_noHem, separation, xtick_label

def get_idx_network_human():
    idx_noHem = np.array([0,1,32,33,34,
                         2,3,4,35,36,37,
                         5,6,7,8,38,39,40,41,42,
                         9,10,11,12,13,14,15,43,44,45,46,
                         16,17,18,47,48,
                         19,20,21,22,49,50,51,52,53,54,
                         23,24,25,26,27,28,29,30,31,55,56,57,58,59,60,61])
    label =  ['Vis', 'SomMot', 'DorsAttn', 'SalVenAttn', 'Limbic', 'Cont', 'Default']
    separation = np.array([5, 11, 20, 31, 36, 46, 62], dtype=float) -1
    return idx_noHem, separation, label

def func_net_ordering(mtx, species='mouse'):
    idx_network, separation, xtick_label = get_idx_mouse()
    if species == 'human':
        idx_network, separation, xtick_label = get_idx_network_human()
    mtx = mtx[idx_network, :][:, idx_network]
    return mtx, separation, xtick_label

def reducing_toNetwork(mtx, species='mouse'):
    #reduce full matrix to a network based smaller matrix
    
    mtx, separation, xtick_label = func_net_ordering(mtx, species=species)
        
    n = len(separation)
    reduced_mtx = np.zeros([n, n])
    separation = np.hstack([0, np.array(separation+1, dtype=int)])
    for i, idx_i in enumerate(separation[:-1]):
        for j, idx_j in enumerate(separation[:-1]):
            #print '(%d, %d) - (%d, %d)' % (idx_i, separation[i+1], idx_j, separation[j+1])
            tmp_mtx = mtx[idx_i : separation[i+1], idx_j : separation[j+1]]
            if i == j:
                #tmp_mtx is square => remove diag
                tmp_mtx = tmp_mtx[np.logical_not(np.identity(len(tmp_mtx)))]
            reduced_mtx[i, j] = np.mean(tmp_mtx) #np.mean(np.sum(tmp_mtx, 1)) #TO SET
    return reduced_mtx, xtick_label

def set_figure_netSep_human():
    separation = np.array([2, 5, 9, 16, 19, 23, 32, 35, 38, 43, 47, 49, 55, 62, 68, 74], dtype=float) -1  
    xtick_label = ['Vis', 'SomMot', 'DorsAttn', 'SalVenAttn', 'Limbic', 'Cont', 'Default', 'Vis', 'SomMot', 'DorsAttn', 'SalVenAttn', 'Limbic', 'Cont', 'Default', 'SubCorLF', 'SubCorRH']
    return separation, xtick_label

def set_figure_netSep_mouse():
    idx_hemSep = np.array([7, 8, 9, 10, 11, 15, 16, 17, 4, 18, 19, 22, 29,
                          0, 1, 2, 3, 30,
                          12, 13, 14, 
                          5,
                          6,
                          20, 21, 24, 25, 26, 27,
                          23,
                          28,
                          31, 32, 33,
                          34, 35,
                          36,
                          44, 45, 46, 47, 48, 52, 53, 54,
                          41, 55, 56, 59,
                          66,
                          37, 38, 39, 40, 67,
                          49, 50, 51,
                          42,
                          43,
                          57, 58, 61, 62, 63, 64,
                          60,
                          65,
                          68, 69, 70,
                          71, 72,
                          73], dtype=int)    
    separation = np.array([8, 12, 13, 18, 21, 22, 23, 29, 30, 31, 34, 36, 37], dtype=float) -1
    xtick_label = ['DMNmid', 'DMNpl', 'DMNsub', 'LCN', 'SAL', 'VIS', 'AUD', 'HPC', 'PIR', 'CTXsp', 'BF', 'THAL', 'HY']
    return idx_hemSep, separation, xtick_label    

def set_figure_areaSep_mouse():
    idx_noHem = np.array([5,
                          0,2,3,
                          4,12,13,14,22,
                          6,19,20,21,
                          15,16,17,18,
                          1,7,8,9,10,11,
                          23,
                          24,25,26,27,
                          28,
                          29,30,
                          31,32,33,
                          34,35,
                          36,
                          42,
                          37,39,40,
                          41,49,50,51,59,
                          43,56,57,58,
                          52,53,54,55,
                          38,44,45,46,47,48,
                          60,
                          61,62,63,64,
                          65,
                          66,67,
                          68,69,70,
                          71,72,
                          73], dtype=int)    
    separation = np.array([1, 4, 9, 13, 17, 23, 24, 28, 29, 31, 34, 36, 37], dtype=float) -1
    xtick_label = ['Vis', 'SM', 'antL', 'TEMP', 'MED', 'PF', 'PIR', 'HPF', 'CTXsp', 'STR', 'BF', 'THAL', 'HY']
    return idx_noHem, separation, xtick_label 