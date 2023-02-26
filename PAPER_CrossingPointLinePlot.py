# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 16:41:19 2023

@author: jonrm
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

from Plotting import PlotLine, PlotData, MultiPanelPlot#, PlotPanel

from RED_analyt_sub import N_s_CF, M_L, nu_s_CF_minRepro, nu_s_CF_full
from RED_analyt_sub import a_bar, m_bar#, g_bar, s_bar, s_bar_minRepro
from optima_sub1 import MIN_VAL, OptN_mu_p1_alpha, OptM_alpha_mu_p1
from lin_fit_sub import Lin_fit, GetSlopeIntercept

def N_M_mu_p1_diff(m_s, a_s, mu_p1, f_rs=1.0, m_r=None, phi_g=0.75, phi_a=0.5):
    """ Find alpha that optimises M for given mu_p1, then use that alpha to
    find the mu_p1 that optimises N.  The difference between the latter mu_p1 
    and the original one shows us how close we are to the crossing point of the 
    N and M optima"""
    
    # Find the alpha value that optimises biomass for given mu_p1
    alphaMOpt = OptM_alpha_mu_p1(m_s, a_s, mu_p1, phi_g, phi_a, m_r=m_r, f_rs=f_rs)
    
    # Find the mu_p1 value that optimises N assuming the alpha value we just found
    mu_p1NOpt = OptN_mu_p1_alpha(m_s, a_s, alphaMOpt, phi_g, phi_a, m_r=m_r, f_rs=f_rs)
    
    # This then tells us the difference in mu_p1 between the one given and the
    # one that optimses N.  At the crossing point of the optima nullclines these
    # are the same and the difference is zero
    return mu_p1 - mu_p1NOpt

def FindCrossPoint(m_s, a_s, f_rs=1.0, m_r=None, phi_g=0.75, phi_a=0.5):
    """ Find the value of mu_p1 and alpha that has both N and M optima"""
    
    def rootFunc(mu_p1):
        return N_M_mu_p1_diff(m_s, a_s, mu_p1, f_rs, m_r, phi_g, phi_a)

    # Find the root mu_p1
    res = opt.root_scalar(rootFunc, bracket=[MIN_VAL, 5.0])
    
    mu_p1 = res.root
    alpha = OptM_alpha_mu_p1(m_s, a_s, mu_p1, phi_g, phi_a, m_r=m_r, f_rs=f_rs)
    
    return mu_p1, alpha


def FindCrossPoint_f_rs(m_s, a_s, mu_1_target, m_r=None, phi_g=0.75, phi_a=0.5):
    """ Find the value of f_rs that gives both peak M and N for given mu_1 """
    def rootFunc(f_rs):
        mu_p1, alpha = FindCrossPoint(m_s, a_s, f_rs, m_r, phi_g, phi_a)
        
        return mu_1_target - mu_p1 / (1-alpha)
    
    # Abort if solution not possible, ie if need f_rs > 1
    if np.sign(rootFunc(MIN_VAL)) == np.sign(rootFunc(1.0)):
        return np.nan
    
    # Find the root f_rs
    res = opt.root_scalar(rootFunc, bracket=[MIN_VAL, 1.0])
    f_rs = res.root
    
    return f_rs

def calc_mu_s(m_s, mu_p1, alpha, phi_g=0.75):
    """ Calculate seedling mu ie mu at m_s"""
    c = 1 - phi_g
    
    mu_1 = mu_p1 / (1 - alpha)
    mu_s = mu_1*m_s**c
    return mu_s    

def calc_coverage(m_s, mu_s, alpha, f_rs, phi_g=0.75, m_r=None):
    """ Calculates the fractional coverage of tree crowns of the forest """
    if m_r is None:
        coverage = nu_s_CF_full(alpha, mu_s, phi_g, f_rs)
    else: # If have minimum tree size m_r below which there is no seed production
        coverage = nu_s_CF_minRepro(alpha, mu_s, m_s, m_r, phi_g, f_rs)
        
    return coverage

def calc_N(coverage, mu_s, a_s, phi_g=0.75, phi_a=0.5):
    """ Calculate the total forest tree density """
    mean_crown_area = a_bar(a_s, mu_s, phi_g, phi_a)
    return N_s_CF(coverage, mean_crown_area)
    
def calc_M(N, m_s, mu_s, phi_g=0.75):
    """ Calculate the total forest biomass density """
    mean_mass = m_bar(m_s, mu_s, phi_g) # Mean mass of all trees from m_s to infinity
    return M_L(N, mean_mass) 

def calc_stuff(m_s, a_s, mu_1_target, phi_g=0.75, phi_a=0.5, m_r=None):
    f_rs = FindCrossPoint_f_rs(m_s, a_s, mu_1_target, m_r=m_r, phi_g=phi_g, phi_a=phi_a)
    
    if np.isnan(f_rs):
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    
    mu_p1, alpha = FindCrossPoint(m_s, a_s, f_rs=f_rs, m_r=m_r, phi_g=phi_g, phi_a=phi_a)
    
    # Calculate diagnostics
    mu_s = calc_mu_s(m_s, mu_p1, alpha, phi_g)
    coverage = calc_coverage(m_s, mu_s, alpha, f_rs, phi_g, m_r)
    N = calc_N(coverage, mu_s, a_s, phi_g, phi_a)
    M = calc_M(N, m_s, mu_s, phi_g)
      

    return f_rs, mu_p1, alpha, coverage, N, M

def DrawPanel(xlabel, ylabel, x, y, lab, lab_pos, paneltype=None):
    lines = []
    
    lines.append(PlotLine(x, y, kwargs={'linestyle':'-'}))
    pD = PlotData(lines, label=lab, lab_pos=lab_pos)
    pD.axSet.SetLabels(xlabel, ylabel)
    if paneltype == 'loglog':
        pD.axSet.SetLogLog()
    elif paneltype == 'semilogx':
        pD.axSet.SetSemiLogX() 
        
    return pD

###############################################################################
###############################################################################
    
if __name__ == "__main__": 
    print('Start')
    
    plt.style.use('Nature.mplstyle')

    
    # Color cycle based on https://personal.sron.nl/~pault/
    colors= ['#004488', '#DDAA33', '#BB5566']
    
    
    phi_g = 0.75
    phi_a = 0.5

    c = 1 - phi_g
    x = 1./c
    
    a1 = 0.5 # Crown area of seedling of mass 1 kg C (guesstimate)
    
    mu_1_target = 0.235 # mu1 from Rainfor
    
    m_s_array = np.logspace(-4, 0, 9)#np.array([0.001, 0.01, 0.1, 1.])

    
    # Min reproductive size, plants smaller than this produce no seeds
    m_r_array=np.logspace(1, 3, 5)
    
    # Fill arrays with nans initially, with dimensions same as seed size array
    f_rs_array = np.full_like(m_s_array, np.nan)
    M_array = np.full_like(m_s_array, np.nan)
    N_array = np.full_like(m_s_array, np.nan)
    nu_array = np.full_like(m_s_array, np.nan)
    mu_p1_array = np.full_like(m_s_array, np.nan)
    alpha_array =  np.full_like(m_s_array, np.nan)

    
    # Versions with minimum reproductive size
    f_rs_array_minRepro = np.full((len(m_s_array), len(m_r_array)), np.nan)
    M_array_minRepro = np.full_like(f_rs_array_minRepro, np.nan)
    N_array_minRepro = np.full_like(f_rs_array_minRepro, np.nan)
    nu_array_minRepro = np.full_like(f_rs_array_minRepro, np.nan)
    mu_p1_array_minRepro = np.full_like(f_rs_array_minRepro, np.nan)
    alpha_array_minRepro =  np.full_like(f_rs_array_minRepro, np.nan)
    
    for i, m_s in enumerate(m_s_array):
        print('m_s='+str(m_s))
        a_s = a1*m_s**phi_a
        
        out = calc_stuff(m_s, a_s, mu_1_target, phi_g, phi_a)
              
        # Assign results to arrays
        f_rs_array[i], mu_p1_array[i], alpha_array[i], \
        nu_array[i], N_array[i], M_array[i] = out


    ###############################
    # Plots
    
    xlabel = r'Seed Mass $m_s$'
    lab_pos = (0.1, 0.9)
    
    panels = []
    
    pD = DrawPanel(xlabel, 'Proportion of reproductive assimilate\n surviving to seed germination '+r'$f_{rs}$',
                    m_s_array, f_rs_array, 'a', lab_pos, 'loglog')
    panels.append(pD)
    
    pD = DrawPanel(xlabel, 'Proportion assimilate allocated to\n'+r'reproduction $\alpha$', m_s_array, alpha_array, 'b', lab_pos, 'semilogx')
    panels.append(pD)   

    pD = DrawPanel(xlabel, r'1 kg tree mortality : assimilate ratio $\mu_{p1}$', m_s_array, mu_p1_array, 'c', lab_pos, 'semilogx')
    panels.append(pD)    

    pD = DrawPanel(xlabel, r'Total Biomass $M$ - kg C m$^{-2}$', m_s_array, M_array, 'd', lab_pos, 'semilogx')
    panels.append(pD)
    
    pD = DrawPanel(xlabel, r'Total Tree Density $N$ - m$^{-2}$', m_s_array, N_array, 'e', lab_pos, 'semilogx')
    panels.append(pD)
    
    pD = DrawPanel(xlabel, r'Total Crown Coverage', m_s_array, nu_array, 'f', lab_pos, 'semilogx')
    panels.append(pD)  
    
    ###############################
    
    spacings = {'hspace':0.25, 'wspace': 0.35, 'left':0.08, 'right': 0.995,
              'top': 0.99}
    
    filename = 'DO_results'
    MultiPanelPlot(panels, 3, spacings=spacings, fileN=filename, allowExpand=False)

    
    ###############################
    
    # Find slope of f_rs vs seed size line using lin regression
    y = np.log10(f_rs_array)
    x = np.log10(m_s_array)
    res = Lin_fit(x, y)
    slope, intercept, _ , _ = GetSlopeIntercept(res)
    k = 10**intercept
    
    fig, ax1 = plt.subplots(1)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.plot(m_s_array, f_rs_array, 'ko', label='Theory')  
    ax1.plot(m_s_array, k*m_s_array**slope, 'r:', label='Regression fit')
    ax1.set_ylabel(r'$f_{rs}$ proportion of reproductive assimilate surviving to seed germination')
    ax1.set_xlabel(r'Seed Mass $m_s$ - kg C')
    ax1.legend()
    plt.savefig('f_rs_regression.pdf')
    plt.close()    
    