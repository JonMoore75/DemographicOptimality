# -*- coding: utf-8 -*-
"""
Created on Mon Aug 21 08:58:42 2023

@author: jonrm
"""

import numpy as np

import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.colors import LinearSegmentedColormap

from RED_analyt_sub import nu_s_CF_full_vec, m_bar, a_bar
from optima_sub1 import OptN_mu_p1_alpha, OptM_alpha_mu_p1

def mu_s_vec(alpha_array, m_s_array, mu_p1, phi_g=0.75):
    c = 1 - phi_g
    mu_p_array = mu_p1*m_s_array**c
    return np.outer(1./(1 - alpha_array), mu_p_array)

###############################################################################
###############################################################################
    
if __name__ == "__main__": 
    print('Start')
    
    ###############################
    # Niklas / MST Allometry
    phi_g = 0.75
    phi_a = 0.5
    
    c = 1 - phi_g
    x = 1./c
    
    a1 = 0.5
    
    mu_1 = 0.235 # mu1 for carbon mass from Rainfor
    alpha_typical = 0.1
    mu_p1 = mu_1 * (1 - alpha_typical)
    
    gamma = 0.05
    
    
    m_s_arr = np.logspace(-4, 0, 41)
    
    ###############################

    dM_sp_d_al = np.zeros(len(m_s_arr))
    
    for i, ms in enumerate(m_s_arr):    
        a_s = a1*ms**phi_a
        dM_sp_d_al[i] = OptM_alpha_mu_p1(ms, a_s, mu_p1, phi_g, phi_a)
  
    ###############################
    

    al_b3 = 4*(m_s_arr**0.5)*mu_p1**2/np.sqrt(  12*(m_s_arr**0.25)*mu_p1 + 3)
    lb3 = r'$\alpha \approx \dfrac{4 m_{s}^{1/2} \, \mu_{p1}^{2}}{\sqrt{  12 m_{s}^{1/4} \, \mu_{p1} + 3 }}$'
    
    al_b4 = 4*m_s_arr**0.5*mu_p1**2/np.sqrt(3)
    lb4 = r'$\alpha \approx \dfrac{4 m_{s}^{1/2} \, \mu_{p1}^{2}}{\sqrt{ 3 }}$'


    
    mu_str = r'$\mu_{p1}=$'+str(mu_p1)

    fig, ax = plt.subplots()
    ax.set_xscale('log')
    plt.plot(m_s_arr, dM_sp_d_al, 'k-', label='Actual Optima')
    plt.plot(m_s_arr, al_b3, '--', label=lb3)
    plt.plot(m_s_arr, al_b4, '--', label=lb4)
    plt.xlabel('Seed Size (kg C)')
    plt.ylabel(r'$\alpha$')
    plt.title(r'Accuracy of Approx for $\dfrac{\partial M_s}{\partial \alpha}=0$, '+mu_str)
    plt.legend()
    plt.show()
    
    fig, ax = plt.subplots()
    ax.set_xscale('log')
    plt.plot(m_s_arr, al_b3/dM_sp_d_al, '-', label=lb3)
    plt.plot(m_s_arr, al_b4/dM_sp_d_al, '-', label=lb4)
    plt.xlabel('Seed Size (kg C)')
    plt.ylabel(r'$\alpha$ Approx / Actual')
    plt.title(r'Accuracy of Approx for $\dfrac{\partial M_s}{\partial \alpha}=0$, '+mu_str)
    plt.legend()
    plt.show()