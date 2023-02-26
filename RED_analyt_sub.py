# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 15:37:17 2020

@author: jonrm
"""

import numpy as np
import scipy.special as sps

###############################################################################
    
def gam_uinc_sp(x,z):
    """ Calculates the upper incomplete gamma function """ 
    return sps.gamma(x)*sps.gammaincc(x, z)

def F(x, z):
    """ Evaluates the common forms seen in these equations of
        exp(z)*Gamma(x, z) / z^(x-1)
        
        This evaluates to a short finite series of terms when x an integer
        
        Usually z = mu_s / (1 - phi_g)
        
        We prefer using the series form when possible as better behaved but will 
        automatically use the gamma function version when not using MST
    """
    if np.isclose(x, 3):
        return 1 + 2/z + 2/z**2
    if np.isclose(x, 4):
        return 1 + 3/z + 6/z**2 + 6/z**3
    if np.isclose(x, 5):
        return 1 + 4/z + 12/z**2 + 24/z**3 + 24/z**4
    print('Gamma evaluation with (x,z)='+str((x,z)))
    return np.exp(z)*gam_uinc_sp(x, z) / ( z**(x - 1))
        

###############################################################################
# RED Analytical Equilubrium Equations
# Common parameters:
#
# gamma - tree mortality rate (per yr)
# m_s - mass of tree of seddling size (kg C)
# g_s - tree growth rate of tree of seeding size (kg C per yr)
# a_s - tree crown area tree of seeding size (m^2)
# mu_s - ratio of mortality to growth for tree of seeding size = gamma * m_s / g_s
# g_1 - tree growth rate of 1 kg C tree (kg C per yr)
# a_1 - tree crown area tree of 1 kg C tree (m^2)
# mu_1 - ratio of mortality to growth for 1 kg C tree = gamma / g_1
# phi_g - tree growth scaling power (power law) 
# phi_a - tree crown area scaling power (power law) 
# x = 1 / ( 1 - phi_g )
# y = x * mu_s


def nDist(m, n_L, mu_1, phi_g = 0.75):
    """ Equilibrium Size Distribution - Weibull Distribution"""
    c = 1 - phi_g
    x = 1/c
    return n_L*np.exp(x*mu_1*(m[0]**c - m**c))*(m[0]/m)**phi_g

def a_L(a_1, m_L, phi_a=0.5):
    """ Crown area of tree of mass m_L """
    return a_1*m_L**phi_a

def g_L(g_1, m_L, phi_g=0.75):
    """ Growth rate of tree of mass m_L """
    return g_1*m_L**phi_g

# Gridbox variables, ie properties of whole forest of all tree sizes.  
# Subscript L denotes this is evaluated from any given size m_L to infinity
# eg To calculate at the seed size m_s then use n_s, g_s, N_s etc

def N_L(n_L, _g_L, gamma):
    """ Equilibrium Total Tree Density """
    return n_L*_g_L/gamma

def nu_L(_N_L, _a_bar):
    """ Equilibrium Fractional Coverage - Proportion of ground covered by 
    tree crowns """
    return _N_L*_a_bar

def M_L(_N_L, _m_bar):
    """ Equilibrium Total Biomass Density """
    return _N_L*_m_bar

def G_L(_N_L, _g_bar):
    """ Equilibrium Total Growth Density """
    return _N_L*_g_bar


def mu_L(mu_1, m_L, phi_g = 0.75 ):
    """ Calculate ratio of mortality to growth rate at tree size m_L """
    c = 1 - phi_g
    return mu_1*m_L**c

# Mean tree variables, ie gridbox total divided by number trees (all sizes)
# for tree distribution from size m_L to infinity

def a_bar(_a_L, mu_L, phi_g = 0.75, phi_a=0.5):
    """ Mean treee crown area """
    x = 1/( 1 - phi_g )
    y = x*mu_L

    return _a_L*F(phi_a*x+1, y) 

def g_bar(_g_L, mu_L, phi_g = 0.75):
    """ Mean tree growth rate """
    x = 1/( 1 - phi_g )
    y = x*mu_L

    return _g_L*F(x, y)

def m_bar(m_L, mu_L, phi_g = 0.75):
    """ Mean tree mass """
    x = 1/( 1 - phi_g )
    y = x*mu_L

    return m_L*F(x+1, y) 

def h_bar(_h_L, mu_L, phi_g=0.75, phi_h=0.25):
    """ Mean tree height """
    x = 1/( 1 - phi_g )
    y = x*mu_L    

    return _h_L*F(phi_h*x+1, y) 

def ML_s(_ml_L, _N_L, mu_L, phi_g=0.75, phi_ml=0.75):
    """ Equilibrium Total Leaf Carbon Mass Density.  This will give LAI in terms 
    of gridbox area if multiply by sigma_l the specific leaf carbon density"""
    x = 1/( 1 - phi_g )
    y = x*mu_L

    return _ml_L*_N_L*F(phi_ml*x+1, y)

def Lnu_s(_l_L, mu_L, phi_g=0.75, phi_l=0.25):
    """ Equilibrium LAI for covered area. ie Leaf area per unit area covered 
    by tree crowns """
    x = 1/( 1 - phi_g )
    y = x*mu_L 

    return _l_L*F(phi_l*x+1, y)


def s_bar(_g_s, m_s, mu_s, alpha, phi_g = 0.75, f_rs = 1.0):
    """ Mean tree seed production rate 
        alpha - is proportion assimilate to reproduction, f_rs is fraction
        of alpha going to seeds with remainder lost due to reproductive costs.
    """
    return f_rs*alpha*g_bar(_g_s, mu_s, phi_g)/(m_s*(1 - alpha))

def s_bar_minRepro(_g_s, m_s, m_r, mu_s, alpha, phi_g = 0.75, f_rs = 1.0):
    """ Mean tree seed production rate 
        alpha - is proportion assimilate to reproduction, f_rs is fraction
        of alpha going to seeds with remainder lost due to reproductive costs.
        min reproductive size m_r > seedling size m_s"""
    x = 1/( 1 - phi_g )
    g_r = _g_s*(m_r/m_s)**phi_g
    mu_r = mu_s*(m_r/m_s)**(1- phi_g)
    return s_bar(g_r, m_s, mu_r, alpha, phi_g, f_rs)*np.exp(x*(mu_s - mu_r))



###############################################################################
# RED Closed Form Equilibrium Equations - when a fixed proportion (alpha) of NPP
# is allocated to reproduction

def nu_s_CF(_gamma, _s_bar, mu_seedpool=0.0, clamp=True):  
    """ Equilibrium Fractional Coverage - Proportion of ground covered by 
    tree crowns """
    f_seedpool = 1/(1. + mu_seedpool)
    res = 1 - _gamma / (_s_bar * f_seedpool) 
    
    if clamp:
        res = np.clip(res, 0, 1)
    
    return res

def nu_s_CF_full(alpha, mu_s, phi_g = 0.75, f_rs = 1.0, mu_seedpool=0.0, clamp=True):
    """ Version that doesn't require gamma and g to be be known directly only
    their ratio """
    f_seedpool = 1/(1. + mu_seedpool)
    
    res = mu_s / s_bar(1., 1., mu_s, alpha, phi_g, f_rs)
    
    res =  1 - res / f_seedpool
     
    if clamp:
        res = np.clip(res, 0, 1)
    
    return res 

def nu_s_CF_minRepro(alpha, mu_s, m_s, m_r, phi_g = 0.75, f_rs = 1.0, mu_seedpool=0.0, clamp=True):
    """ Equilibrium Fractional Coverage - Proportion of ground covered by 
    tree crowns 
    for when there is a min reproductive size m_r > seedling size m_s""" 
    
    f_seedpool = 1/(1. + mu_seedpool)
    
    res = mu_s / (m_s*s_bar_minRepro(1., m_s, m_r, mu_s, alpha, phi_g, f_rs))
    
    
    res = 1 - res / f_seedpool

    if clamp:
        res = np.clip(res, 0, 1)
    
    return res

def nu_s_CF_full_vec(alpha, mu_s, phi_g = 0.75, f_rs = 1.0, mu_seedpool=0.0, clamp=True):
    """ Version of nu_s_CF_full for when both alpha and mu_s are vectors/arrays """
    f_seedpool = 1/(1. + mu_seedpool) # Loss ratio in seed pool
    
    res_al = (1 - alpha) / alpha
    res_mu = mu_s /  (f_rs*g_bar(1.0, mu_s, phi_g) * f_seedpool)
    
    if res_mu.ndim > 1: # if res_mu already 2d, make sure res_al in right shape 
        # [:,np.newaxis] adds new axis to array converts 1d array to 2d with 1 column etc
        # eg (10,) to (10,1)
        res_al = res_al[:,np.newaxis]
        res = 1 - res_al*res_mu 
    else:
        res =  1 - np.outer(res_al, res_mu) 
     
    if clamp:
        res = np.clip(res, 0, 1)
    
    return res

# Tree density

def N_s_CF(coverage, mean_crown_area):
    """ Equilibrium Total Tree Density is the """
    return coverage / mean_crown_area
    


    
