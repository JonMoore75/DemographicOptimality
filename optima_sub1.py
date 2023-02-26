# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 15:37:17 2020

@author: jonrm
"""

import numpy as np
import scipy.optimize as opt

from RED_analyt_sub import N_s_CF, M_L, nu_s_CF_minRepro, nu_s_CF_full, a_bar, m_bar

# Set max and min mu_p1 bounds for optima finding
MIN_VAL = np.finfo(float).eps
MAX_mu_VAL = 6.
MAX_alpha_VAL = 0.9999
MAX_ms_VAL = 10.
                        
def OptN_mu_p1_alpha(m_s, a_s, alpha, phi_g=0.75, phi_a=0.5, m_r=None, mu_seed=0.0, f_rs=1.0):
    """ Finds mu_p1 that optimises Total Tree Density N for const alpha """
    c = 1 - phi_g
    
    if m_r is None:
        m_r = m_s
    
    def FitFunc(mu_p1):
        mu_1 = mu_p1/(1-alpha)
        mu_s = mu_1*m_s**c
    
        coverage = nu_s_CF_minRepro(alpha, mu_s, m_s, m_r, phi_g, 
                                    mu_seedpool=mu_seed, f_rs=f_rs, clamp=False)
        mean_crown_area = a_bar(a_s, mu_s, phi_g, phi_a)
        
        return -N_s_CF(coverage, mean_crown_area)
    
    #res = opt.minimize_scalar(FitFunc, method='brent', bracket=(MIN_VAL, MAX_mu_VAL))
    res = opt.minimize_scalar(FitFunc, method='bounded', bounds=(MIN_VAL, MAX_mu_VAL))
    
    if res.success:
        return res.x
    print('No mu_p1 soln alpha='+str(alpha)+', m_r='+str(m_r), res)
    return np.nan


def OptN_mu1_alpha(m_s, a_s, alpha, phi_g=0.75, phi_a=0.5, m_r=None, mu_seed=0.0, f_rs=1.0):
    """ Finds mu_1 that optimises Total Tree Density N for const alpha """
    c = 1 - phi_g
    
    if m_r is None:
        m_r = m_s
    
    def FitFunc(mu_1):
        mu_s = mu_1*m_s**c
        
         
        coverage = nu_s_CF_minRepro(alpha, mu_s, m_s, m_r, phi_g, 
                                    mu_seedpool=mu_seed, f_rs=f_rs, clamp=False)
        mean_crown_area = a_bar(a_s, mu_s, phi_g, phi_a)
        
        return -N_s_CF(coverage, mean_crown_area)
    
    #res = opt.minimize_scalar(FitFunc, method='brent', bracket=(MIN_VAL, MAX_mu_VAL))
    res = opt.minimize_scalar(FitFunc, method='bounded', bounds=(MIN_VAL, MAX_mu_VAL))
    
    if res.success:
        return res.x
    print('No mu_1 soln alpha='+str(alpha)+', m_r='+str(m_r), res)
    return np.nan 

def OptN_alpha_mu_p1(m_s, a_s, mu_p1, phi_g=0.75, phi_a=0.5, m_r=None, mu_seed=0.0, f_rs=1.0):
    """ Finds alpha that optimises Total Tree Density for given mu_p1"""
    c = 1 - phi_g
    
    if m_r is None:
        m_r = m_s
    
    def FitFunc(alpha):
        mu_s = mu_p1*m_s**c/(1-alpha)
         
        coverage = nu_s_CF_minRepro(alpha, mu_s, m_s, m_r, phi_g, 
                                    mu_seedpool=mu_seed, f_rs=f_rs, clamp=False)
        mean_crown_area = a_bar(a_s, mu_s, phi_g, phi_a)
        
        return -N_s_CF(coverage, mean_crown_area)
    
    res = opt.minimize_scalar(FitFunc, method='bounded', bounds=(MIN_VAL, MAX_alpha_VAL))
    
    if res.success:
        return res.x.item()
    
    print('No alpha soln mu_p1='+str(mu_p1)+', m_r='+str(m_r), res)
    return np.nan

def OptN_alpha_mu_1(mu_ps, a1, mu_1, phi_g=0.75, phi_a=0.5, mu_seed=0.0, f_rs=1.0):
    """ Finds alpha that optimises Total Tree Density for given mu1"""
    x = 1/(1 - phi_g)
    
    def FitFunc(alpha):
        mu_s = mu_ps/(1 - alpha)
        m_s = ( mu_s / mu_1 )**x
        a_s = a1*m_s**phi_a
         
        coverage = nu_s_CF_full(alpha, mu_s, phi_g, mu_seedpool=mu_seed, f_rs=f_rs,
                                clamp=False)
        mean_crown_area = a_bar(a_s, mu_s, phi_g, phi_a)
        
        return -N_s_CF(coverage, mean_crown_area)
    
    res = opt.minimize_scalar(FitFunc, method='bounded', bounds=(MIN_VAL, MAX_alpha_VAL))
    
    if res.success:
        return res.x.item()
    
    print('No alpha soln mu_1='+str(mu_1), res)
    return np.nan

def OptM_alpha_mu_p1(m_s, a_s, mu_p1, phi_g=0.75, phi_a=0.5, m_r=None, f_rs=1.0, mu_seed=0.0):
    """ Finds alpha that optimises Biomass for constant mu_p1"""
    c = 1 - phi_g
    
    if m_r is None:
        m_r = m_s
    
    def FitFunc(alpha):
        mu_1 = mu_p1/(1-alpha)
        mu_s = mu_1*m_s**c

            
        mean_crown_area = a_bar(a_s, mu_s, phi_g, phi_a)
        mean_mass = m_bar(m_s, mu_s, phi_g)
        
        coverage = nu_s_CF_minRepro(alpha, mu_s, m_s, m_r, phi_g, 
                                    mu_seedpool=mu_seed, f_rs=f_rs, clamp=False)
        N = N_s_CF(coverage, mean_crown_area)
        
        return -M_L(N, mean_mass)   
    
    res = opt.minimize_scalar(FitFunc, method='bounded', bounds=(MIN_VAL, MAX_alpha_VAL))
    if res.success:
        return res.x.item()
    print('No alpha soln mu_p1='+str(mu_p1)+', m_r='+str(m_r), res)
    return np.nan 

def OptM_alpha_mu_p1_2(mu_s, mu_p1, a1, phi_g=0.75, phi_a=0.5, m_r=None, f_rs=1.0, mu_seed=0.0):
    """ Finds alpha that optimises Biomass for constant mu_p1"""
    x = 1./(1 - phi_g)
    
    def FitFunc(alpha):
        mu_1 = mu_p1/(1-alpha)
        m_s = (mu_s / mu_1)**x
        a_s = a1*m_s**phi_a

        if m_r is None:
            m_rr = m_s
        else:
            m_rr = m_r
            
        mean_crown_area = a_bar(a_s, mu_s, phi_g, phi_a)
        mean_mass = m_bar(m_s, mu_s, phi_g)
        
        coverage = nu_s_CF_minRepro(alpha, mu_s, m_s, m_rr, phi_g, 
                                    mu_seedpool=mu_seed, f_rs=f_rs, clamp=False)
        N = N_s_CF(coverage, mean_crown_area)
        
        return -M_L(N, mean_mass)   
    
    res = opt.minimize_scalar(FitFunc, method='bounded', bounds=(MIN_VAL, MAX_alpha_VAL))
    if res.success:
        return res.x.item()
    print('No alpha soln mu_p1='+str(mu_p1)+', m_r='+str(m_r), res)
    return np.nan


def OptM_ms_alpha(alpha, a_1, mu_p1, phi_g=0.75, phi_a=0.5, m_r=None, mu_seed=0.0, f_rs=1.0):
    """ Finds seed mass that optimises Biomass for a given mu_p1 and alpha """
    c = 1 - phi_g
    mu_1 = mu_p1/(1-alpha)
    
    
    def FitFunc(m_s):
        _a_s = a_1*m_s**phi_a
        mu_s = mu_1*m_s**c
        
        if m_r is None:
            m_rr = m_s
        else:
            m_rr = m_r
        
        coverage = nu_s_CF_minRepro(alpha, mu_s, m_s, m_rr, phi_g, 
                                    mu_seedpool=mu_seed, f_rs=f_rs, clamp=False)
        mean_crown_area = a_bar(_a_s, mu_s, phi_g, phi_a)
        mean_mass = m_bar(m_s, mu_s, phi_g)
        
        return -(coverage*mean_mass / mean_crown_area)
    
    res = opt.minimize_scalar(FitFunc, method='bounded', bounds=(MIN_VAL, MAX_ms_VAL))
    
    if res.success:
        return res.x.item()
    
    print('No m_s soln mu_p1='+str(mu_p1)+', m_r='+str(m_r), res)
    return np.nan 

def OptM_mu_ps_alpha(alpha, a_1, mu_1, phi_g=0.75, phi_a=0.5, m_r=None, mu_seed=0.0, f_rs=1.0):
    """ Finds seed mass that optimises Biomass for a given mu_p1 and alpha """
    x = 1./(1 - phi_g)
    
    def FitFunc(mu_ps):
        mu_s = mu_ps/(1-alpha)
        m_s = (mu_s/mu_1)**x
        a_s = a_1*m_s**phi_a
        
        if m_r is None:
            m_rr = m_s
        else:
            m_rr = m_r
        
        coverage = nu_s_CF_minRepro(alpha, mu_s, m_s, m_rr, phi_g, 
                                    mu_seedpool=mu_seed, f_rs=f_rs, clamp=False)
        mean_crown_area = a_bar(a_s, mu_s, phi_g, phi_a)
        mean_mass = m_bar(m_s, mu_s, phi_g)
        
        return -(coverage*mean_mass / mean_crown_area)
    
    res = opt.minimize_scalar(FitFunc, method='bounded', bounds=(MIN_VAL, MAX_ms_VAL))
    
    if res.success:
        return res.x.item()
    
    print('No mu_ps soln mu_1='+str(mu_1)+', m_r='+str(m_r), res)
    return np.nan


def Opt_nu_alpha_mu_ps(mu_ps, m_s, phi_g=0.75, phi_a=0.5, m_r=None, \
                       f_rs=1.0, mu_seed=0.0):
    """ Finds alpha that optimises coverage for constant mu_ps and m_s"""
    
    if m_r is None:
        m_r = m_s
    
    def FitFunc(alpha):
        mu_s = mu_ps/(1-alpha)
        
        coverage = nu_s_CF_minRepro(alpha, mu_s, m_s, m_r, phi_g, 
                                    mu_seedpool=mu_seed, f_rs=f_rs, clamp=False)
        
        return -coverage   
    
    res = opt.minimize_scalar(FitFunc, method='bounded', bounds=(MIN_VAL, MAX_alpha_VAL))
    if res.success:
        return res.x.item()
    print('No alpha soln mu_ps='+str(mu_ps)+', m_r='+str(m_r), res)
    return np.nan 

def Opt_nu_alpha_mu_ps2(mu_ps, mu_1, phi_g=0.75, phi_a=0.5, m_r=None, \
                       f_rs=1.0, mu_seed=0.0):
    """ Finds alpha that optimises coverage for constant mu_ps and mu_1"""
    x = 1./(1-phi_g)
  
    def FitFunc(alpha):
        mu_s = mu_ps/(1-alpha)
        m_s = (mu_s/mu_1)**x
        
        if m_r is None:
            m_rr = m_s
        
        coverage = nu_s_CF_minRepro(alpha, mu_s, m_s, m_rr, phi_g, 
                                    mu_seedpool=mu_seed, f_rs=f_rs, clamp=False)
        
        return -coverage   
    
    res = opt.minimize_scalar(FitFunc, method='bounded', bounds=(MIN_VAL, MAX_alpha_VAL))
    if res.success:
        return res.x.item()
    print('No alpha soln mu_ps='+str(mu_ps)+', m_r='+str(m_r), res)
    return np.nan 