# -*- coding: utf-8 -*-
"""
"""
import scipy.stats as stats
import matplotlib.pyplot as plt
from scipy.io import loadmat, savemat
import numpy as np
import os
       
if __name__=='__main__':

    pwd_current = '%s/' % os.getcwd()
    pwd_data = '%sdata/' % pwd_current    
    if not(os.path.isdir('%sfigures' % pwd_current)):
        os.makedirs('%sfigures' % pwd_current)
    
    dataset = loadmat('%sdata.mat' % pwd_data, squeeze_me=True)
    
    #fig 4
    from fig_4 import plot_fig4
    if not(os.path.isdir('%sfigures/fig4' % pwd_current)):
        os.makedirs('%sfigures/fig4' % pwd_current)   
    print('--------- FIG 4 -------------')
    rho, p_val, rho_lemon, p_val_lemon, rho_hcp, p_val_hcp = plot_fig4()
    
    #fig 5
    from fig_5 import plot_fig5
    if not(os.path.isdir('%sfigures/fig5' % pwd_current)):
        os.makedirs('%sfigures/fig5' % pwd_current)    
    print('--------- FIG 5 -------------')
    #mouse dataset
    A_mouse = dataset['A_mouse'] #A, Areverse, Aasym    
    sigma2_mouse = dataset['sigma2_mouse']
    i_sim = 0 #original data, time-irreversible dynamics (A asym)
    plot_fig5(A_mouse[:,:,:,i_sim], sigma2_mouse[:,i_sim], species='mouse', labelfigure='mouse')

    #human MPI-LMBB dataset
    A_lemon = dataset['A_human_lemon']    
    sigma2_lemon = dataset['sigma2_human_lemon']
    plot_fig5(A_lemon, sigma2_lemon, species='human', labelfigure='mpi_lmbb')

    #human HCP dataset
    A_hcp = dataset['A_human_hcp']    
    sigma2_hcp = dataset['sigma2_human_hcp']
    plot_fig5(A_hcp, sigma2_hcp, species='human', labelfigure='hcp')    