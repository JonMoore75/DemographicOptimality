# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 15:33:29 2023

@author: jonrm
"""

import pickle
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from dataclasses import dataclass, field

@dataclass
class PlotLine:
    x_data: np.array
    y_data: np.array
    fmt:str = None
    label: str = None
    color: str = None
    
def SetMaxXLim2Nearest100(ax, xvalues):
    """ Rounds up the maximum axis DBH value to nearest 100cm  """
    xmax = (int(max(xvalues)/100.)+1)*100 
    xmin = ax.get_xlim()[0]
    ax.set_xlim(xmin, xmax)
    
def LabelFmt(y, pos): 
    # Find the number of decimal places required, is 0 for numbers >=1
    decimalplaces = int(np.ceil(np.maximum(-np.log10(y), 0)))     
    # Insert that number into a format string
    formatstring = '{:.{prec}f}'.format(y, prec=decimalplaces)
    # Return the formatted tick label
    return formatstring#.format(y)

def UseScalarLogAxisLabels(axis, minor=True):
    """ Changes log axis labels from 10^2 style to scalar value e.g. 100 """
    axis.set_major_formatter(ticker.ScalarFormatter())
    axis.set_major_formatter(ticker.FuncFormatter( LabelFmt ))
    
    if minor:
        axis.set_minor_formatter(ticker.ScalarFormatter())
        axis.set_minor_formatter(ticker.FuncFormatter( LabelFmt ))
        
def LabelEveryOtherMinorTickLabel(axis):
    """ Adds labels to minor tick labels """
    for label in axis.get_ticklabels(minor=True)[1::2]:
        label.set_visible(False)
    
if __name__ == "__main__": 
    print('Start')

    mpl.rcParams['figure.dpi'] = 200
    
    #plt.style.use(thisfilepath+'/MatPlotLibStyles/StevePlot.mplstyle')
    #plt.style.use(thisfilepath+'/MatPlotLibStyles/custom.mplstyle')
    #plt.style.use(thisfilepath+'/MatPlotLibStyles/Nature.mplstyle')
    
    with open('US_Rainfor_Data.pkl', 'rb') as plot_file:
        lines = pickle.load(plot_file)
        
    # lines is a list of PlotLine classes
    # Each PlotLine class has x_data, y_data, fmt, label and color
    # Lines should be in order USDA Data, USDA fit, Rainfor Data, Rainfor Fit
    # So for example lines[0] is USDA Data and lines[3] is Rainfor Fit
    # fmt is the matplotlib plot fmt e.g. 'o' or '-'
    # color is the line color
    # label is the label for the legend
    
    lines[2].color = 'g'
    lines[3].color = 'm'
        
    fig, ax = plt.subplots(1)  

    # log-log axes
    ax.set_yscale('log')
    ax.set_xscale('log')
    
    for line in lines:
        ax.plot(line.x_data, line.y_data, line.fmt, color=line.color, label=line.label)
    
    ax.set_xlabel('Trunk Diameter (cm)')
    ax.set_ylabel('Number Density ( Trees cm$^{-1}$ ha$^{-1}$ )')
    
    SetMaxXLim2Nearest100(ax, xvalues = lines[0].x_data)
    UseScalarLogAxisLabels(ax.xaxis)
    LabelEveryOtherMinorTickLabel(ax.xaxis)
    ax.legend()
    plt.tight_layout()
    plt.show()
    
    