#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  5 11:32:39 2023

@author: benozzo
"""
import numpy as np
from scipy import stats, linalg

def computeS(A, sigma2):
    #from Casti et al., 2023 (CDC2023)
    nS = A.shape[0]
    Q = sigma2 * np.diag(np.ones(nS))
    Sigma = - linalg.solve_continuous_lyapunov(A, Q)
    S = 0.5 * (np.dot(A, Sigma) - np.dot(Sigma, A.T))
    return S, Sigma, Q

def computeEntropy(S, Sigma, Q, A):
    #from Gilson et al. 2023 Physical Review
    phi = -2 * np.trace(np.linalg.multi_dot([linalg.inv(Q), A, S]))
    return phi