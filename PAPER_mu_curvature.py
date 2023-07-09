# -*- coding: utf-8 -*-
"""
Created on Sun Feb 26 15:45:32 2023

@author: jonrm
"""
from copy import deepcopy
import numpy as np

import matplotlib.pyplot as plt

from RED_analyt_sub import nDist

###############################################################################
###############################################################################
    
if __name__ == "__main__": 
    print('Start')
    
    plt.style.use('Nature.mplstyle')

    
    # Color cycle based on https://personal.sron.nl/~pault/
    colors= ['#004488', '#DDAA33', '#BB5566']
    
    mu_1_array = np.array([0.1, 0.2, 0.3])
    
    m_array = np.logspace(0, 4, 21)
    
    n_1 = 1.
    
    fSize = deepcopy(plt.rcParams.get('figure.figsize'))
    fSize = [i * 0.75 for i in fSize]
    
        
    fig, ax1 = plt.subplots(1, figsize=fSize)#(3.6, 2.25))
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylabel(r'Tree number density - trees (kg C)$^{-1}$ ha$^{-1}$')
    ax1.set_xlabel(r'Tree Mass $m$ - kg C')
    for i, mu_1 in enumerate(mu_1_array):
        n_array = nDist(m_array, n_1, mu_1)
        ax1.plot(m_array, n_array, '-', color=colors[i], label=str(mu_1))
    ax1.legend(title=r'$\mu_1$')
    plt.tight_layout()
    plt.savefig('mu_1_effect.pdf')
    plt.close()    
    