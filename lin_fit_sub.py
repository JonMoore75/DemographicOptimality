# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 12:03:34 2022

@author: jonrm
"""
import numpy as np
import statsmodels.api as sm

def Lin_fit(x, y):
    # Make sure all infinities convert to NaN so are ignored by fitting routine
    y[ y == np.Inf] = np.NaN
    y[ y == -np.Inf] = np.NaN
    
    #  Then remove the NaNs
    valid = ~np.isnan(y)
    y = y[valid]
    x = x[valid]
    
    # Do regression first forming the matrix X
    X = sm.add_constant(x)
    model = sm.OLS(y,X,missing='drop')
    results = model.fit()
    
    return results

def GetSlopeIntercept(results):
    slope_err, intercept_err = results.bse[1], results.bse[0]
    slope, intercept= results.params[1], results.params[0] 

    return slope, intercept, slope_err, intercept_err