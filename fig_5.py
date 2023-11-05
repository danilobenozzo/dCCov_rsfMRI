#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 14:43:34 2022

@author: benozzo
"""
import matplotlib.pyplot as plt
from plotting_functions import figure_layout, set_size
import numpy as np 
import scipy.stats as stats
from scipy.io import loadmat, savemat

from network_separation_library import set_figure_areaSep_mouse, set_figure_netSep_human, reducing_toNetwork
from function_library import computeS
  
def plot_fig5(A, sigma2, species='mouse', labelfigure='mouse'):
       
    n_sbj = A.shape[0]
    nS = A.shape[1]
    traceS = np.zeros([n_sbj, nS])
    degreeS = np.zeros([n_sbj, nS]) #column degree
    Smean = np.zeros([nS, nS])
    idx_triu = np.triu(np.ones(nS) == 1, 1)
    for i_sbj in range(n_sbj):
        Si, Sigma, Q = computeS(A[i_sbj], sigma2[i_sbj]) #compute S, subject level
        traceS[i_sbj] = np.diag(np.dot(Si, Si.T))
        di = np.sum(Si, 0)
        degreeS[i_sbj] = di / np.sqrt(np.sum(di**2) / float(nS-1)) #zscore degree S
        Smean = Smean + (Si / np.sqrt(np.sum(Si[idx_triu]**2) / float(Si[idx_triu].shape[0] - 1))) #zscore S mean
        
    Smean = Smean / n_sbj

    #plot in/out node profile
    w = 6.2
    h = 2.1
    y_in = 2.5 #important for figure_layout
    x_in = w / h * y_in 
    fontsize_main, fontsize, linewidth_main, linewidth, markersize_main, markersize = figure_layout(x_in, y_in)       
    fig, ax = plt.subplots(figsize=(x_in, y_in))
    if species == "mouse":
        idx_noHem, separation, xtick_label = set_figure_areaSep_mouse()
        y_toplot = degreeS.copy()[:, idx_noHem]
    elif species == 'human':
        separation, xtick_label = set_figure_netSep_human()
        y_toplot = degreeS.copy()
    else:
        print('"error: species not valid')

    #test single roi
    p_rois = np.zeros(nS)
    for i_roi in range(nS):
        _, p_rois[i_roi] = stats.ttest_1samp(y_toplot[:, i_roi], popmean=0.0)
    
    from statsmodels.stats import multitest
    reject, p_rois_corrected, _, _ = multitest.multipletests(p_rois, alpha=.05, method='fdr_bh') #Benjamini/Hochberg 
    
    plt.plot([-.5, nS-.5], [0, 0], '-k', linewidth=linewidth)
    plt.plot(np.mean(y_toplot, 0), '.r', label=labelfigure, linewidth=linewidth_main)
    plt.plot(np.vstack([range(nS), range(nS)]), np.vstack([np.mean(y_toplot, 0)+np.std(y_toplot, 0), np.mean(y_toplot, 0)-np.std(y_toplot, 0)]), '-', color='r', linewidth=linewidth)
    plt.plot(np.arange(nS)[reject], np.ones(nS)[reject]*(-1.8), '*k', markersize=markersize)
    
    shift = .5
    delta_sep = np.diff(separation, prepend=-1) / 2

    plt.plot([separation+shift, separation+shift], [np.ones(separation.shape)*(-nS)-shift, np.ones(len(separation))*nS-shift], '-k', linewidth=linewidth) #vertical line
    if species == "mouse":
        plt.plot([separation+shift+nS/2, separation+shift+nS/2], [np.ones(separation.shape)*(-nS)-shift, np.ones(len(separation))*nS-shift], '-k', linewidth=linewidth) #vertical line
        plt.xticks(np.hstack([separation-delta_sep+shift, separation-delta_sep+shift + nS/2]), xtick_label + xtick_label, fontsize=6, rotation='vertical')
    elif species == "human":
        plt.plot([separation+shift, separation+shift], [np.ones(separation.shape)*(-nS)-shift, np.ones(len(separation))*nS-shift], '-k', linewidth=linewidth) #vertical line
        plt.xticks(separation-delta_sep+shift, xtick_label, fontsize=6, rotation='vertical')
    else:
        print('"error: species not valid')

    plt.xlim([-.5,nS-.5])
    plt.ylim([-2, 5])
    if species == 'human':
        plt.ylim([-2.5, 2.5])
    ax.tick_params(labelsize=fontsize)
    ax = set_size(w, h, ax)
    fig.tight_layout()
    
    label_fig = 'degreeS_norm_areaSep_%s' % (labelfigure)
    plt.savefig('figures/fig5/%s.svg' % (label_fig), format='svg')
    plt.show()
    ########################
    #reducing network
    mtx_toplot, xtick_label = reducing_toNetwork(Smean, species=species)
    mtx_toplot[mtx_toplot < 0] = 0 #from column to row: sender
    w = 2.1
    h = 2.1 #3.1
    y_in = 2.5#important for figure_layout
    x_in = w / h * y_in 
    fontsize_main, fontsize, linewidth_main, linewidth, markersize_main, markersize = figure_layout(x_in, y_in)       
    fig, ax = plt.subplots(figsize=(x_in, y_in))
    cm_bin = plt.cm.get_cmap('binary')
    ax.imshow(mtx_toplot, cmap=cm_bin)

    shift = .5
    separation = np.arange(mtx_toplot.shape[0])  
    
    ax.plot([separation+shift, separation+shift], [np.zeros(separation.shape)-shift, np.ones(len(separation))*len(separation)-shift], '-k', linewidth=linewidth) #vertical line
    delta_sep = np.diff(separation, prepend=-1) / 2.
    ax.set_xticks(separation-delta_sep+shift)
    ax.set_xticklabels(xtick_label, fontsize=fontsize, rotation='vertical')
    
    ax.plot([np.zeros(separation.shape)-shift, np.ones(len(separation))*len(separation)-shift], [separation+shift, separation+shift], '-k', linewidth=linewidth)
    ax.set_yticks(separation-delta_sep+shift)
    ax.set_yticklabels(xtick_label, fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)
    ax = set_size(w, h, ax)
    fig.tight_layout()
    label_fig = 'mtx_ntw_norm_areaSep_%s' % (labelfigure)
    plt.savefig('figures/fig5/%s.svg' % (label_fig), format='svg')
    plt.show()
    
    
    if species == "human":
        #comparison with seguin2019inferring GLOBAL
        mtx_toplot, _ = reducing_toNetwork(Smean, species=species)
        idx_triu = np.triu(np.ones(mtx_toplot.shape[0]) == 1, 1)
        dX = [10, 15, 20]
        for i_d, dX_i in enumerate(dX):
            print('d', dX_i)
            data_seguin_fig4a = loadmat('data/human_d%s_yeo7.mat' % dX_i, squeeze_me=True)
            mtx_seguin = np.array(data_seguin_fig4a['dif_yeo7'], dtype=float).T
            rho, p_val = stats.pearsonr(mtx_toplot[idx_triu], mtx_seguin[idx_triu])
            print('human coupling with structural diffusion efficiency: rho: %f, p: %f' % (rho, p_val))
    
    return