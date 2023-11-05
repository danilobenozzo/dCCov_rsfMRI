#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 14:43:34 2022

@author: benozzo
"""
import matplotlib.pyplot as plt
from scipy.io import loadmat, savemat
from scipy import stats, linalg
from plotting_functions import figure_layout, set_size
import numpy as np 
    
def compute_relation_invS_A_S(A, Sigma, S, sigma2):
    # relation between incomingA / degreeinvS and outgoingA / degreeS
    nS = A.shape[0]
    Q = (-0.5) * sigma2 * np.diag(np.ones(nS))
    S_new = Q + S
    inA = np.sum(A, 1)# - np.diagonal(A) #row sum
    outA = np.sum(A, 0)# - np.diagonal(A) #column sum
    inS = np.sum(S, 1) #row sum
    outS = np.sum(S, 0) #column sum
    degree_S_new = np.sum(S_new, 0)
    degree_invSigma = np.sum(linalg.inv(Sigma), 0)
    
    rho_in, p_val_in = stats.spearmanr(inA, degree_invSigma)
    rho_in_expected, p_val_in_expected = stats.spearmanr(inA, inS)
    inA_corrected = inA + (sigma2/2.)*degree_invSigma
    rho_in_corr, p_val_in_corr = stats.spearmanr(inA_corrected, inS)
    
    rho_out, p_val_out = stats.spearmanr(outA, outS)
    rho_out_Sigma, p_val_out_Sigma = stats.spearmanr(outA, degree_invSigma)
    outA_corrected = outA + (sigma2/2.)*degree_invSigma
    rho_out_corr, p_val_out_corr = stats.spearmanr(outA_corrected, outS)
    #lin_fit[i_sbj, 0], lin_fit[i_sbj, 1], r, p, se = stats.linregress(S_mouse[i_sbj, idx[0], idx[1], 0], S_mouse[i_sbj, idx[0], idx[1], 1])
    return rho_in, p_val_in, rho_in_expected, p_val_in_expected, rho_in_corr, p_val_in_corr, rho_out, p_val_out, rho_out_Sigma, p_val_out_Sigma, rho_out_corr, p_val_out_corr
 
def make_plot_single_sbj(A, Sigma, S, sigma2, end_label):
        
    inA = np.sum(A, 1)# - np.diagonal(A)
    degree_invSigma = np.sum(linalg.inv(Sigma), 0)
    m, q, r, p, se = stats.linregress(degree_invSigma, inA)
    print 'li fit degree invSigma and in A: m=%s, q=%s' % (m, q)
    
    # example one subj inA v invSigma
    w = 1.1
    h = 1.1
    y_in = 2.5#important for figure_layout
    x_in = w / h * y_in 
    fontsize_main, fontsize, linewidth_main, linewidth, markersize_main, markersize = figure_layout(x_in, y_in)    
    fig, ax = plt.subplots(1, 1, figsize=(x_in, y_in))
    plt.plot(degree_invSigma, inA, '.k', markersize=2)
    plt.plot(np.array([degree_invSigma.min(), degree_invSigma.max()]), np.array([degree_invSigma.min(), degree_invSigma.max()]) * m + q, '-r', linewidth=linewidth_main)
    plt.plot(np.array([degree_invSigma.min(), degree_invSigma.max()]), np.array([degree_invSigma.min(), degree_invSigma.max()]) * (-sigma2/2), linewidth=linewidth_main)
    plt.grid()
    ax.set_xlim([degree_invSigma.min() * 1.1, degree_invSigma.max() * 1.1])
    ax.set_xlabel('invSigma row sum', fontsize=fontsize)
    #ax.set_xticks(np.arange(-.2, .21, .2))
    #ax.set_xticklabels([], fontsize=fontsize)
    ax.set_ylim([inA.min() * 1.1, inA.max() * 1.1])
    ax.set_ylabel('A row sum', fontsize=fontsize)
    #ax.set_yticks(np.arange(-.2, .21, .2))
    #ax.set_yticklabels([], fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)
    ax = set_size(w, h, ax)
    fig.tight_layout()
    pwd_tosave = 'figures/fig4/'
    filename = 'inAVSdegreeinvSigma_%s' % end_label
    fig.savefig('%s%s.svg' % (pwd_tosave, filename), format='svg')

    # example one subj inA v row S
    inS = np.sum(S, 1)
    m, q, r, p, se = stats.linregress(inS, inA)
    print 'li fit degree in S and in A: m=%s, q=%s' % (m, q)
    
    w = 1.1
    h = 1.1
    y_in = 2.5#important for figure_layout
    x_in = w / h * y_in 
    fontsize_main, fontsize, linewidth_main, linewidth, markersize_main, markersize = figure_layout(x_in, y_in)    
    fig, ax = plt.subplots(1, 1, figsize=(x_in, y_in))
    plt.plot(inS, inA, '.k', markersize=2) #what we would expect when referring to incoming links
    plt.plot(np.array([inS.min(), inS.max()]), np.array([inS.min(), inS.max()]) * m + q, '-r', linewidth=linewidth_main)
    plt.plot(np.array([0, 0]), np.array([inA.min(), inA.max()]), linewidth=linewidth_main)
    plt.grid()
    ax.set_xlim([inS.min() * 1.1, inS.max() * 1.1])
    ax.set_xlabel('S row sum', fontsize=fontsize)
    #ax.set_xticks(np.arange(-.2, .21, .2))
    #ax.set_xticklabels([], fontsize=fontsize)
    ax.set_ylim([inA.min() * 1.1, inA.max() * 1.1])
    ax.set_ylabel('A row sum', fontsize=fontsize)
    #ax.set_yticks(np.arange(-.2, .21, .2))
    #ax.set_yticklabels([], fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)
    ax = set_size(w, h, ax)
    fig.tight_layout()
    pwd_tosave = 'figures/fig4/'
    filename = 'inAVSinS_%s' % end_label
    fig.savefig('%s%s.svg' % (pwd_tosave, filename), format='svg')

    # example one subj inA corrected v row S
    inA_corr = inA + (sigma2/2.) * degree_invSigma
    m, q, r, p, se = stats.linregress(inS, inA_corr)
    print 'li fit degree inS and in A corrected: m=%s, q=%s' % (m, q)
    
    w = 1.1
    h = 1.1
    y_in = 2.5#important for figure_layout
    x_in = w / h * y_in 
    fontsize_main, fontsize, linewidth_main, linewidth, markersize_main, markersize = figure_layout(x_in, y_in)    
    fig, ax = plt.subplots(1, 1, figsize=(x_in, y_in))
    plt.plot(inS, inA_corr, '.k', markersize=2) #what we would expect when referring to incoming links
    plt.plot(np.array([inS.min(), inS.max()]), np.array([inS.min(), inS.max()]) * m + q, '-r', linewidth=linewidth_main)
    plt.plot(0, 0, '.', markersize=markersize)
    plt.grid()
    ax.set_xlim([inS.min() * 1.1, inS.max() * 1.1])
    ax.set_xlabel('S row sum', fontsize=fontsize)
    #ax.set_xticks(np.arange(-.2, .21, .2))
    #ax.set_xticklabels([], fontsize=fontsize)
    ax.set_ylim([inA_corr.min() * 1.1, inA_corr.max() * 1.1])
    ax.set_ylabel('Delta row sum', fontsize=fontsize)
    #ax.set_yticks(np.arange(-.2, .21, .2))
    #ax.set_yticklabels([], fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)
    ax = set_size(w, h, ax)
    fig.tight_layout()
    pwd_tosave = 'figures/fig4/'
    filename = 'inAcorrVSinS_%s' % end_label
    fig.savefig('%s%s.svg' % (pwd_tosave, filename), format='svg')
    
    ############################
    
    nS = A.shape[0]
    Q = (-0.5) * sigma2 * np.diag(np.ones(nS))
    S_new = Q + S
    outA = np.sum(A, 0)# - np.diagonal(A) #somma colonne
    degree_S = np.sum(S, 0) #np.sum(S_new, 0)
    m, q, r, p, se = stats.linregress(degree_S, outA)
    print 'li fit degree S and out A: m=%s, q=%s' % (m, q)

    # example one subj outA v degreeS
    w = 1.1
    h = 1.1
    y_in = 2.5#important for figure_layout
    x_in = w / h * y_in 
    fontsize_main, fontsize, linewidth_main, linewidth, markersize_main, markersize = figure_layout(x_in, y_in)    
    fig, ax = plt.subplots(1, 1, figsize=(x_in, y_in))
    plt.plot(degree_S, outA, '.k', markersize=2)
    plt.plot(np.array([degree_S.min(), degree_S.max()]), np.array([degree_S.min(), degree_S.max()]) * m + q, '-r', linewidth=linewidth_main)
    plt.plot([0, 0], np.array([outA.min(), outA.max()]), linewidth=linewidth_main) #[-sigma2/2, -sigma2/2]
    plt.grid()
    ax.set_xlim([degree_S.min() * 1.1, degree_S.max() * 1.1])
    ax.set_xlabel('S columns sum', fontsize=fontsize)
    #ax.set_xticks(np.arange(-.2, .21, .2))
    #ax.set_xticklabels([], fontsize=fontsize)
    ax.set_ylim([outA.min() * 1.1, outA.max() * 1.1])
    ax.set_ylabel('A column sum', fontsize=fontsize)
    #ax.set_yticks(np.arange(-.2, .21, .2))
    #ax.set_yticklabels([], fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)
    ax = set_size(w, h, ax)
    fig.tight_layout()
    pwd_tosave = 'figures/fig4/'
    filename = 'outAVSdegreeS_%s' % end_label
    fig.savefig('%s%s.svg' % (pwd_tosave, filename), format='svg')
    
    # example one subj outA v invSigma
    m, q, r, p, se = stats.linregress(degree_invSigma, outA)
    print 'li fit degree invSigma and out A: m=%s, q=%s' % (m, q)

    w = 1.1
    h = 1.1
    y_in = 2.5#important for figure_layout
    x_in = w / h * y_in 
    fontsize_main, fontsize, linewidth_main, linewidth, markersize_main, markersize = figure_layout(x_in, y_in)    
    fig, ax = plt.subplots(1, 1, figsize=(x_in, y_in))
    plt.plot(degree_invSigma, outA, '.k', markersize=2)
    plt.plot(np.array([degree_invSigma.min(), degree_invSigma.max()]), np.array([degree_invSigma.min(), degree_invSigma.max()]) * m + q, '-r', linewidth=linewidth_main)
    plt.plot(np.array([degree_invSigma.min(), degree_invSigma.max()]), np.array([degree_invSigma.min(), degree_invSigma.max()]) * (-sigma2/2), linewidth=linewidth_main)
    plt.grid()
    ax.set_xlim([degree_invSigma.min() * 1.1, degree_invSigma.max() * 1.1])
    ax.set_xlabel('invSigma row sum', fontsize=fontsize)
    #ax.set_xticks(np.arange(-.2, .21, .2))
    #ax.set_xticklabels([], fontsize=fontsize)
    ax.set_ylim([outA.min() * 1.1, outA.max() * 1.1])
    ax.set_ylabel('A column sum', fontsize=fontsize)    
    #ax.set_yticks(np.arange(-.2, .21, .2))
    #ax.set_yticklabels([], fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)
    ax = set_size(w, h, ax)
    fig.tight_layout()
    pwd_tosave = 'figures/fig4/'
    filename = 'outAVSdegreeinvSigma_%s' % end_label
    fig.savefig('%s%s.svg' % (pwd_tosave, filename), format='svg')
    
    # example one subj outA corrected v degreeS
    outA_corr = outA + (sigma2/2.) * degree_invSigma
    m, q, r, p, se = stats.linregress(degree_S, outA_corr)
    print 'li fit degree degree S and out A corrected: m=%s, q=%s' % (m, q)
    
    w = 1.1
    h = 1.1
    y_in = 2.5#important for figure_layout
    x_in = w / h * y_in 
    fontsize_main, fontsize, linewidth_main, linewidth, markersize_main, markersize = figure_layout(x_in, y_in)    
    fig, ax = plt.subplots(1, 1, figsize=(x_in, y_in))
    plt.plot(degree_S, outA_corr, '.k', markersize=2) #what we would expect when referring to incoming links
    plt.plot(np.array([degree_S.min(), degree_S.max()]), np.array([degree_S.min(), degree_S.max()]) * m + q, '-r', linewidth=linewidth_main)
    plt.plot(0, 0, '.', markersize=markersize)
    plt.grid()
    ax.set_xlim([degree_S.min() * 1.1, degree_S.max() * 1.1])
    ax.set_xlabel('S column sum', fontsize=fontsize)
    #ax.set_xticks(np.arange(-.2, .21, .2))
    #ax.set_xticklabels([], fontsize=fontsize)
    ax.set_ylim([outA_corr.min() * 1.1, outA_corr.max() * 1.1])
    ax.set_ylabel('Delta column sum', fontsize=fontsize)
    #ax.set_yticks(np.arange(-.2, .21, .2))
    #ax.set_yticklabels([], fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)
    ax = set_size(w, h, ax)
    fig.tight_layout()
    pwd_tosave = 'figures/fig4/'
    filename = 'outAcorrVSdegreeS_%s' % end_label
    fig.savefig('%s%s.svg' % (pwd_tosave, filename), format='svg')
    
def plot_fig4():
    
    dataset = loadmat('data/data.mat', squeeze_me=True)
    
    j_sim = 0
    S_mouse = dataset['S_mouse'][:, :, :, j_sim]
    sigma2_mouse = dataset['sigma2_mouse'][:, j_sim]
    A_mouse = dataset['A_mouse'][:, :, :, j_sim]
    Sigma_mouse = dataset['Sigma_mouse'][:, :, :, j_sim]
    n_mouse = S_mouse.shape[0]
    
    rho = np.zeros([n_mouse, 6]) #in, in_expected, in_corr, out, out_sigma, out_corr
    p_val = np.zeros([n_mouse, 6]) #in, in _expected, in_corr, out, out_sigma, out_corr
    for i_sbj in range(n_mouse):
        rho[i_sbj, 0], p_val[i_sbj, 0], rho[i_sbj, 1], p_val[i_sbj, 1], rho[i_sbj, 2], p_val[i_sbj, 2], rho[i_sbj, 3], p_val[i_sbj, 3], rho[i_sbj, 4], p_val[i_sbj, 4], rho[i_sbj, 5], p_val[i_sbj, 5] = compute_relation_invS_A_S(A_mouse[i_sbj], Sigma_mouse[i_sbj], S_mouse[i_sbj], sigma2_mouse[i_sbj])
    
    i_sbj = 0
    make_plot_single_sbj(A_mouse[i_sbj], Sigma_mouse[i_sbj], S_mouse[i_sbj], sigma2_mouse[i_sbj], 'mouse_%s' % i_sbj)
        
    # to set for humans       
    S_human_lemon = dataset['S_human_lemon']
    sigma2_lemon = dataset['sigma2_human_lemon']
    A_human_lemon = dataset['A_human_lemon']
    Sigma_human_lemon = dataset['Sigma_human_lemon']
    n_human = S_human_lemon.shape[0]
    
    rho_lemon = np.zeros([n_human, 6])
    p_val_lemon = np.zeros([n_human, 6])
    for i_sbj in range(n_human):
        rho_lemon[i_sbj, 0], p_val_lemon[i_sbj, 0], rho_lemon[i_sbj, 1], p_val_lemon[i_sbj, 1], rho_lemon[i_sbj, 2], p_val_lemon[i_sbj, 2], rho_lemon[i_sbj, 3], p_val_lemon[i_sbj, 3], rho_lemon[i_sbj, 4], p_val_lemon[i_sbj, 4], rho_lemon[i_sbj, 5], p_val_lemon[i_sbj, 5] = compute_relation_invS_A_S(A_human_lemon[i_sbj], Sigma_human_lemon[i_sbj], S_human_lemon[i_sbj], sigma2_lemon[i_sbj])

    S_human_hcp = dataset['S_human_hcp']
    sigma2_hcp = dataset['sigma2_human_hcp']
    A_human_hcp = dataset['A_human_hcp']
    Sigma_human_hcp = dataset['Sigma_human_hcp']
    n_human = S_human_hcp.shape[0]
    
    rho_hcp = np.zeros([n_human, 6])
    p_val_hcp = np.zeros([n_human, 6])
    for i_sbj in range(n_human):
        rho_hcp[i_sbj, 0], p_val_hcp[i_sbj, 0], rho_hcp[i_sbj, 1], p_val_hcp[i_sbj, 1], rho_hcp[i_sbj, 2], p_val_hcp[i_sbj, 2], rho_hcp[i_sbj, 3], p_val_hcp[i_sbj, 3], rho_hcp[i_sbj, 4], p_val_hcp[i_sbj, 4], rho_hcp[i_sbj, 5], p_val_hcp[i_sbj, 5] = compute_relation_invS_A_S(A_human_hcp[i_sbj], Sigma_human_hcp[i_sbj], S_human_hcp[i_sbj], sigma2_hcp[i_sbj])
    
    ################## boxplot in 
    w = 1.5 #2.1
    h = 1.5 #2.1
    y_in = 2.5#important for figure_layout
    x_in = w / h * y_in 
    fontsize_main, fontsize, linewidth_main, linewidth, markersize_main, markersize = figure_layout(x_in, y_in)   
    flierprops = dict(marker='+', markerfacecolor='k', markersize=markersize, linestyle='none', markeredgecolor='k')
    whiskerprops = dict(color='k')
    medianprops = dict(color='k')
    capprops = dict(color='k')
    boxprops = dict(color='k')
    
    fig, ax = plt.subplots(1, 1, figsize=(x_in, y_in))
    delta_x = [-.2, 0, .2]
    delta_xi = [-.05, 0, .05]
    
    for d_i, delta_xii in enumerate(delta_xi):
        ax.boxplot(np.arctanh(rho[:, d_i]), positions=[delta_x[0] + delta_xii], patch_artist=True, widths=.05, flierprops=flierprops, whiskerprops=whiskerprops, medianprops=medianprops, capprops=capprops, boxprops=boxprops)
        ax.boxplot(np.arctanh(rho_lemon[:, d_i]), positions=[delta_x[1] + delta_xii], patch_artist=True, widths=.05, flierprops=flierprops, whiskerprops=whiskerprops, medianprops=medianprops, capprops=capprops, boxprops=boxprops)
        ax.boxplot(np.arctanh(rho_hcp[:, d_i]), positions=[delta_x[2] + delta_xii], patch_artist=True, widths=.05, flierprops=flierprops, whiskerprops=whiskerprops, medianprops=medianprops, capprops=capprops, boxprops=boxprops)
    ax.grid(axis='y')
    
    ax.set_xticks(delta_x)
    ax.set_xticklabels(['mouse', 'human', 'human'], fontsize=fontsize)
    ax.set_xlim([-.3, .3])
    ax.set_ylim([-3, 3])
    ax.tick_params(labelsize=fontsize)
    ax = set_size(w, h, ax)
    fig.tight_layout()
    pwd_tosave = 'figures/fig4/'
    filename = 'corr_inAVSdegreeinvSigma'
    fig.savefig('%s%s.svg' % (pwd_tosave, filename), format='svg')
    
    ################### boxpot out
    w = 1.5 #2.1
    h = 1.5 #2.1
    y_in = 2.5#important for figure_layout
    x_in = w / h * y_in 
    fontsize_main, fontsize, linewidth_main, linewidth, markersize_main, markersize = figure_layout(x_in, y_in)       
    fig, ax = plt.subplots(1, 1, figsize=(x_in, y_in))
    delta_x =  [-.2, 0, .2] #[-.1, 0, .1]
    delta_xi = [-.05, 0, .05] #[-.05, .05]
    
    for d_i, delta_xii in enumerate(delta_xi):        
        ax.boxplot(np.arctanh(rho[:, 3+d_i]), positions=[delta_x[0] + delta_xii], patch_artist=True, widths=.05, flierprops=flierprops, whiskerprops=whiskerprops, medianprops=medianprops, capprops=capprops, boxprops=boxprops)
        ax.boxplot(np.arctanh(rho_lemon[:, 3+d_i]), positions=[delta_x[1] + delta_xii], patch_artist=True, widths=.05, flierprops=flierprops, whiskerprops=whiskerprops, medianprops=medianprops, capprops=capprops, boxprops=boxprops)
        ax.boxplot(np.arctanh(rho_hcp[:, 3+d_i]), positions=[delta_x[2] + delta_xii], patch_artist=True, widths=.05, flierprops=flierprops, whiskerprops=whiskerprops, medianprops=medianprops, capprops=capprops, boxprops=boxprops)
    ax.grid(axis='y')
    
    ax.set_xticks(delta_x)
    ax.set_xticklabels(['mouse', 'human', 'human'], fontsize=fontsize)
    ax.set_xlim([-.3, .3]) #([-.15, .15])
    ax.set_yticks(np.arange(-1, 2.5, .5))
    ax.set_yticklabels(np.array(np.arange(-1, 2.5, .5), dtype=int), fontsize=fontsize)
    ax.set_ylim([-1, 2]) #([-.1, 2])
    ax.tick_params(labelsize=fontsize)
    ax = set_size(w, h, ax)
    fig.tight_layout()
    pwd_tosave = 'figures/fig4/'
    filename = 'corr_outAVSdegreeS'
    fig.savefig('%s%s.svg' % (pwd_tosave, filename), format='svg')    
    
    print "test rho inA-inS VS. inAcorr-inS (mouse, lemon, hcp)"
    t, pm = stats.ttest_ind(np.arctanh(rho[:,1]), np.arctanh(rho[:,2]))
    t, pl = stats.ttest_ind(np.arctanh(rho_lemon[:,1]), np.arctanh(rho_lemon[:,2]))
    t, ph = stats.ttest_ind(np.arctanh(rho_hcp[:,1]), np.arctanh(rho_hcp[:,2]))
    print "mouse:%4f, lemon:%4f, hcp:%4f" % (pm, pl, ph)
    
    print "test rho outA-outS VS. outAcorr-outS (mouse, lemon, hcp)"
    t, pm = stats.ttest_ind(np.arctanh(rho[:,3]), np.arctanh(rho[:,5]))
    t, pl = stats.ttest_ind(np.arctanh(rho_lemon[:,3]), np.arctanh(rho_lemon[:,5]))
    t, ph = stats.ttest_ind(np.arctanh(rho_hcp[:,3]), np.arctanh(rho_hcp[:,5]))
    print "mouse:%4f, lemon:%4f, hcp:%4f" % (pm, pl, ph)
    
    return rho, p_val, rho_lemon, p_val_lemon, rho_hcp, p_val_hcp