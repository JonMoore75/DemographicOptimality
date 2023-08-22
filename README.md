# Demographic Optimality

This repository is the python code used for the paper Demographic optimality in forests

The various files correspond to the following

Red_analyt_sub.py - The RED equilibrium equations based on DET  
optima_sub1.py - Various routines for finding optima of forest level-metrics in terms of parameters  
lin_fit_sub.py - Basic linear regression  
plotting.py - routines to help simplify the code for the plots, particulary multi-panel plots, uses matplotlib  

PAPER_common_mu.py - Code for figure 1, showing the similar size distributions of Amazon and US forests
PAPER_mu_curvature.py - Code for figure 2, showing the effect mu_1 has on the curvature of RED forest size-distributions  
PAPER_OptimaLinePlot.py - Code for figure 4, showing the optima in terms of a line plot for both biomass and tree density  
PAPER_OptimaContourNullclines.py - Code for figure 5, showing the optima in terms of alpha and mu_p1 as contours and a panel showing the crossing point of the optima  
PAPER_CrossingPointLinePlot.py - Code for figure 6, showing the result when the crossing point of the optima is matched to the US/Amazon mu_1=0.235 as a function of seed size  

SUPPLEMENT_BiomassNullclineAccuracy.py - Code for figure in supplement showing the accuray of the biomass nullcline / optima approximation equation
