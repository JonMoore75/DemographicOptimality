# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 14:27:28 2023

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
    
    plt.style.use('Nature.mplstyle')
    
    # Paul Tol colours https://personal.sron.nl/~pault/
    colors= ['#004488', '#DDAA33', '#BB5566']
    clr_map = ['#FFFFE5', '#FFF7BC', '#FEE391', '#FEC44F', '#FB9A29',
                '#EC7014', '#CC4C02', '#993404', '#662506']
    clr_map.reverse()
    mycmap = LinearSegmentedColormap.from_list('YlOrBr', clr_map)
    mycmap.set_bad('#888888')
    mycmap.set_under('white')
    mycmap.set_over('white')
    
    
    ###############################
    # Niklas / MST Allometry
    phi_g = 0.75
    phi_a = 0.5
    
    c = 1 - phi_g
    x = 1./c
    
    a1 = 0.5
    
    Mlimit = 100.
    
    # Avoid m_s=0, alpha = 0 otherwise get divide by zero
    alpha_array = np.linspace(0.0, 0.5, 2001)[1:] 
    
    
    # m_s_fix = 0.2 #1.
    # f_rs = 1.0
    # mu_p1_array = np.linspace(0.2, 0.6, 2001)
    
    m_s_fix = 10**-4
    f_rs = 0.352*m_s_fix**0.926
    mu_p1_array = np.linspace(0.18, 0.25, 2001)
    
    a_s = a1*m_s_fix**phi_a

    ###############################
    
    
    X, Y = np.meshgrid(mu_p1_array, alpha_array)
    
    mu_s_array = mu_s_vec(alpha_array, m_s_fix, mu_p1_array, phi_g)
    coverage_array = nu_s_CF_full_vec(alpha_array, mu_s_array, phi_g, f_rs=f_rs)
    
    a_bar_array = a_bar(a1*m_s_fix**phi_a, mu_s_array, phi_g, phi_a)
    m_bar_array = m_bar(m_s_fix, mu_s_array, phi_g)
    
    N_array = coverage_array/ a_bar_array
    M_array = N_array*m_bar_array 
   
    ###############################

    maxM = np.ceil(np.nanmax(M_array))
    
    
    nullcline_alpha = np.zeros(len(mu_p1_array))
    for i, mu_p1 in enumerate(mu_p1_array):
        nullcline_alpha[i] = OptM_alpha_mu_p1(m_s_fix, a_s, mu_p1, phi_g, phi_a, f_rs=f_rs ) 

    ###############################
    
    maxN = np.ceil(np.nanmax(N_array))
    

    nullcline_mu_p1 = np.zeros(len(alpha_array))
    for i, alpha in enumerate(alpha_array):
        nullcline_mu_p1[i] = OptN_mu_p1_alpha(m_s_fix, a_s, alpha, phi_g, phi_a, f_rs=f_rs)
        
    nullcline_mu_p1[nullcline_mu_p1 > mu_p1_array[-1]] = np.nan
    nullcline_mu_p1[nullcline_mu_p1 < mu_p1_array[0]] = np.nan
      
    ###############################
    
    title = 'Fixed Seed Mass='+str(m_s_fix)
    
    cb_label = r'Total Biomass $M$ - kg C m$^{-2}$'
    xlabel = r'1 kg tree mortality : assimilate ratio $\mu_{p1}$'
    ylabel = r'Proportion of assimilate to reproduction $\alpha$'
    
    maxM = 30.#5*np.ceil(maxM/5)
    
    levels1 = np.linspace(0., min(maxM, Mlimit), 16)#NumLevels)
#    levels1d = [2., 6., 10., 20.]
    
    lab_pos=(0.075, 0.975)
    
    spacings = {'hspace':0.25, 'wspace': 0.4, 'left':0.075, 'right': 0.98,
              'top': 0.98}
    
    fig = plt.figure()
    gs = gridspec.GridSpec(2, 4, **spacings)
    

    ax1 = fig.add_subplot(gs[0,0:2])
    ax1.text(lab_pos[0], lab_pos[1], 'a', transform=ax1.transAxes, 
        fontweight='bold', va='top', ha='right', color='k',
        fontsize=plt.rcParams.get('xtick.labelsize')*1.5)
    cs = ax1.contourf(X, Y, M_array, cmap=mycmap, levels = levels1, extend='both')
    # cs2 = ax1.contour(X, Y, M_array, levels = levels1d, 
    #                   alpha=0.5, linestyles='dashed', colors='k')
    # ax1.clabel(cs2, inline=True, fontsize=5, fmt='%1.1f')
    ax1.plot(mu_p1_array, nullcline_alpha, '-', color=colors[1])
    ax1.set_ylim(0,0.5)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    cb=plt.colorbar(cs, ax=ax1)
    cb.set_label(cb_label)
 
    
    cb_label = r'Total Tree Density $N$ - m$^{-2}$'
    
    levels2 = np.linspace(0., maxN, 21)#NumLevels)
    
    ax2 = fig.add_subplot(gs[0,2:])
    ax2.text(lab_pos[0], lab_pos[1], 'b', transform=ax2.transAxes, 
    fontweight='bold', va='top', ha='right', 
    fontsize=plt.rcParams.get('xtick.labelsize')*1.5)
    cs = ax2.contourf(X, Y, N_array, cmap=mycmap, levels = levels2, extend='both')
    ax2.plot(nullcline_mu_p1, alpha_array, '-', color=colors[0])
    ax2.set_ylim(0,0.5)
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel(ylabel)
    cb=plt.colorbar(cs, ax=ax2)
    cb.set_label(cb_label) 
    
    ax3 = fig.add_subplot(gs[1,1:-1])
    ax3.text(lab_pos[0], lab_pos[1], 'c', transform=ax3.transAxes, 
    fontweight='bold', va='top', ha='right', 
    fontsize=plt.rcParams.get('xtick.labelsize')*1.5)
    ax3.plot(mu_p1_array, nullcline_alpha, '-', color=colors[1], label=r'Biomass optima')
    ax3.plot(nullcline_mu_p1, alpha_array, '-', color=colors[0], label=r'Tree density optima')
    ax3.set_ylim(0,0.5)
    ax3.set_xlabel(xlabel)
    ax3.set_ylabel(ylabel)
    ax3.legend()
    plt.savefig('Contourf.pdf')
    plt.close() 

 
    
