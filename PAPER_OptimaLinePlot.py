# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 16:15:53 2023

@author: jonrm
"""


import numpy as np
import matplotlib.pyplot as plt


from RED_analyt_sub import nu_s_CF_full, m_bar, a_bar, mu_L
from optima_sub1 import OptM_alpha_mu_p1, OptN_mu_p1_alpha

from Plotting import PlotLine, PlotData, MultiPanelPlot

def calc_M(alpha, mu_p1, m_s, a1, f_rs, phi_g, phi_a):
    # Calculate mortality:growth (mu) 
    mu_s = mu_L(mu_p1/(1-alpha), m_s, phi_g)
    
    mean_tree_mass = m_bar(m_s, mu_s, phi_g)
    tree_density = calc_N(alpha, mu_p1, m_s, a1, f_rs, phi_g, phi_a)
    
    biomass = tree_density*mean_tree_mass
    
    return biomass

def calc_N(alpha, mu_p1, m_s, a1, f_rs, phi_g, phi_a):
    # Calculate mortality:growth (mu) and crown area at seedling size
    mu_s = mu_L(mu_p1/(1-alpha), m_s, phi_g)
    a_s = a1*m_s**phi_a
    
    mean_crown_area = a_bar(a_s, mu_s, phi_g, phi_a)
    coverage = nu_s_CF_full(alpha, mu_s, phi_g, f_rs)
    
    tree_density = coverage / mean_crown_area 
    
    return tree_density

if __name__ == "__main__": 
    print('Start')
    
    plt.style.use('Nature.mplstyle')
    
    # Paul Tol colours https://personal.sron.nl/~pault/
    colors= ['#004488', '#DDAA33', '#BB5566']
    
    ###############################
    # Niklas / MST Allometry
    phi_g = 0.75
    phi_a = 0.5
    
    c = 1 - phi_g
    x = 1./c
    
    ###############################
    # Fixed parameters ( not varying with alpha)
    
    a1 = 0.5
    mu_p1 = 0.25
    gamma = 0.01
    p1 = gamma/mu_p1

    # m_s = 1.0
    # f_rs = 1.0
    # mu_ps_array = [0.2, 0.25, 0.3]
    
    m_s = 10**-4
    f_rs = 7*10**-5
    mu_ps_array = [0.02, 0.021, 0.022]
    
    ###############################
    # Panel a)
    
    lines = []
    
    a_s = a1*m_s**phi_a

    alpha_array = np.linspace(0.0, 0.5, 5001)[1:]
    
    xlabela = r'Proportion of assimilate to reproduction $\alpha$'
    ylabela = r'Total biomass density $M$ - kg C m$^{-2}$'
    lab_pos = (0.05, 1.05)
    
    for i, mu_ps in enumerate(mu_ps_array):
        mu_p1 = mu_ps/(m_s**(1-phi_g))
        M = calc_M(alpha_array, mu_p1, m_s, a1, f_rs, phi_g, phi_a)
        lab = r'$\mu_{p1}=$'+str(round(mu_p1, 3))
        lines.append(PlotLine(alpha_array, M, label=lab, kwargs={'color':colors[i]}))
        
        alphaOpt = OptM_alpha_mu_p1(m_s, a_s, mu_p1, phi_g, phi_a, f_rs=f_rs)
        MOpt = calc_M(alphaOpt, mu_p1, m_s, a1, f_rs, phi_g, phi_a)
        lines.append(PlotLine(alphaOpt, MOpt, label=None,
                    kwargs={'marker':'o', 'color':colors[i], 'markeredgewidth':0} ))

    pD1 = PlotData(lines, label='a', lab_pos = lab_pos)
    pD1.axSet.SetLabels(xlabela, ylabela)
    
    ###############################
    # Panel b)
    
    lines = []
    lines2 = []
    
    mu_ps_array = np.linspace(0.0, 0.04, 401)[1:]
    mu_p1_array = mu_ps_array/(m_s**c)
    alpha_array = np.array([0.05, 0.1, 0.15])

    xlabelb = r'Mortality : assimilate ratio for seedlings $\mu_{ps}$'
    xlabelb2 = r'Mortality : assimilate ratio for 1 kg tree $\mu_{p1}$'
    ylabelb = r'Total tree density $N$ - m$^{-2}$'
    
    for i, alpha in enumerate(alpha_array):
        mu_p1_array = mu_ps_array/(m_s**(1-phi_g))
        N = calc_N(alpha, mu_p1_array, m_s, a1, f_rs, phi_g, phi_a)
        lab = r'$\alpha$ = '+str(round(alpha, 2))
        lines.append(PlotLine(mu_ps_array, N, label=lab, kwargs={'color':colors[i]}))
        lines2.append(PlotLine(mu_p1_array, N, label=lab, kwargs={'color':colors[i]}))
        
        mu_p1Opt = OptN_mu_p1_alpha(m_s, a_s, alpha, phi_g, phi_a, f_rs=f_rs)
        NOpt = calc_N(alpha, mu_p1Opt, m_s, a1, f_rs, phi_g, phi_a)
        lines.append(PlotLine(mu_p1Opt*(m_s**c), NOpt, label=None,
                    kwargs={'marker':'o', 'color':colors[i], 'markeredgewidth':0} ))
        lines2.append(PlotLine(mu_p1Opt, NOpt, label=None,
                    kwargs={'marker':'o', 'color':colors[i], 'markeredgewidth':0} ))


    pD2 = PlotData(lines2, label='b', lab_pos = lab_pos)
    pD2.axSet.SetLabels(xlabelb2, ylabelb)

    filename = 'Optima'
    MultiPanelPlot([pD1, pD2], 2, fileN=filename, allowExpand=False)

