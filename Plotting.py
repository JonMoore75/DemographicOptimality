# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 12:57:00 2023

@author: jonrm
"""
import numpy as np

from copy import deepcopy
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.ticker as ticker
import matplotlib.transforms as transforms
from dataclasses import dataclass, field


@dataclass
class axesSettings:
    ytickOn: bool = True
    ylabelOn: bool = True
    xtickOn: bool = True
    xlabelOn: bool = True
    legendOn: bool = True
    ScalarXTicks: bool = False
    xLog: bool = False
    yLog: bool = False
    xlabel: str = ''
    ylabel: str = ''
        
    def SetLabels(self, xlabel, ylabel):
        self.xlabel = xlabel
        self.ylabel = ylabel
        
    def SetLogLog(self):
        self.xLog = True
        self.yLog = True 
        
    def SetSemiLogX(self):
        self.xLog = True
        self.yLog = False

def LabelFmt(y, pos): 
    # Find the number of decimal places required
    if np.isclose(y, 0.0):
        decimalplaces = 1
    else:
        decimalplaces = int(np.ceil(np.maximum(-np.log10(abs(y)), 0))) 
    # Insert that number into a format string
    formatstring = '{:.{prec}f}'.format(y, prec=decimalplaces)
    # Return the formatted tick label
    return formatstring

def UseScalarLogAxisLabels(axis, minor=True):
    axis.set_major_formatter(ticker.ScalarFormatter())
    axis.set_major_formatter(ticker.FuncFormatter( LabelFmt ))
    
    if minor:
        axis.set_minor_formatter(ticker.ScalarFormatter())
        axis.set_minor_formatter(ticker.FuncFormatter( LabelFmt ))
        
class PlotLine:
    def __init__(self, x_data, y_data, label=None, kwargs={}):
        self.x_data = x_data
        self.y_data = y_data
        self.label = label
        self.kwargs = kwargs
        
@dataclass
class MarkerLine:
    """ This is a line added to plot vertically or horizontally """
    isVert: bool
    value: float
    label: str = ''
    style: str = ':'

@dataclass      
class PlotData:
    lines: [PlotLine] = field(default_factory=list)
    axlines: [MarkerLine] = field(default_factory=list)
    label: str = ''
    lab_pos: (float, float) =(0.05, 1.025)
    ylim: (float, float) = None
    xlim: (float, float) = None
    axSet: axesSettings = None 
    
    def __post_init__(self):
        if self.axSet is None:
            self.axSet = axesSettings()
        
# class PlotData:
#     def __init__(self, lines=[], axlines=[], label='', 
#                  lab_pos=(0.05, 1.025), axSet=None):
#         if axSet is None:
#             axSet = axesSettings()
#         self.lines   = lines
#         self.label   = label
#         self.lab_pos = lab_pos
#         self.axSet   = axSet
#         self.axlines = axlines # Type MarkerLine
#         self.ylim = None
#         self.xlim = None

class AxLimits:
    """ Used to record the limits of all axes in multipanel plots and 
    then set them to common limits"""
    def __init__(self):
        self.axList = []
        self.xminList = []
        self.xmaxList = []
        self.yminList = []
        self.ymaxList = []
        
    def AddLim(self, ax, xlim, ylim):
        self.axList.append(ax)
        self.xminList.append(xlim[0])
        self.xmaxList.append(xlim[1])
        self.yminList.append(ylim[0])
        self.ymaxList.append(ylim[1]) 
        
    def GetCommonXLimits(self):
        return min(self.xminList), max(self.xmaxList)
    
    def GetCommonYLimits(self):
        return min(self.yminList), max(self.ymaxList)
    
    def SetAlltoCommonLimits(self):
        for ax in self.axList:
            ax.set_xlim(*self.GetCommonXLimits())
            ax.set_ylim(*self.GetCommonYLimits())
            
def EnactAxesSettings(ax, axSet, verbose=False):
    if axSet.xlabelOn:
        ax.set_xlabel(axSet.xlabel)
    if axSet.ylabelOn:
        ax.set_ylabel(axSet.ylabel)

    if not axSet.xlabelOn:
        ax.set_xlabel('')
        if verbose: print('No xlab')    
    if not axSet.ylabelOn:
        ax.set_ylabel('')
        if verbose: print('No ylab')
    if not axSet.ytickOn:
        ax.set_yticklabels([])
        if verbose: print('No ytick')
        ax.yaxis.set_major_formatter(ticker.NullFormatter())
    if not axSet.xtickOn:
        ax.set_xticklabels([])
        ax.xaxis.set_major_formatter(ticker.NullFormatter())
        if verbose: print('No xtick')
    elif axSet.ScalarXTicks:
        UseScalarLogAxisLabels(ax.xaxis)
        # Display every other minor tick label
        for label in ax.xaxis.get_ticklabels(minor=True)[1::2]:
            label.set_visible(False)

        
def SetMaxXLim2Nearest100(ax, xvalues):
    xmax = (int(max(xvalues)/100.)+1)*100 # Round up to next 100cm DBH
    xmin = ax.get_xlim()[0]
    ax.set_xlim(xmin, xmax)
    
def PlotPanel(ax, pD, commonLegend=False, verbose=False):
    """ Plots both axes in log10 scale 
        pD is a PlotData object containing the data """
    
    # log axes
    if pD.axSet.yLog:
        ax.set_yscale('log')
    if pD.axSet.xLog:
        ax.set_xscale('log')
    
    # Plot label
    ax.text(pD.lab_pos[0], pD.lab_pos[1], pD.label, transform=ax.transAxes, 
            fontweight='bold', va='top', ha='right', 
            fontsize=plt.rcParams.get('xtick.labelsize')*1.5)

    legCount = 0       
    
    for line in pD.lines:
        ax.plot(line.x_data, line.y_data, label=line.label, **line.kwargs)
        if line.label is not None:
            legCount += 1
            
    
    if len(pD.axlines) > 0:
        # the x coords of this transformation are data, and the
        # y coord are axes
        trans = transforms.blended_transform_factory(
                ax.transData, ax.transAxes)
        for axline in pD.axlines:
            if axline.isVert:
                ax.axvline(axline.value, linestyle=axline.style)
                x = axline.value
                ax.text(x, 0.1, axline.label, transform=trans, \
                        rotation="vertical", va='bottom', ha='left', 
                        fontsize=1.2*plt.rcParams.get('xtick.labelsize'))

            else:
                ax.axhline(axline.value, linestyle=axline.style)
    
    if not commonLegend and legCount > 0:
        ax.legend(loc=0)
    
    #SetMaxXLim2Nearest100(ax, xvalues=pD.lines[0].x_data)
            
    EnactAxesSettings(ax, pD.axSet, verbose)
    
def MultiPanelPlot(plotList, NumPanelCol, fileN=None, commonLegend = False,
                   commonXScale=False, commonAxisLabels = False, 
                   allowExpand = True, spacings=None):
    """ Plots multipanel plot
        plotList is expected to be list of PlotData items """
    NumPanels = len(plotList)
    NumPanelRow = int(np.ceil(float(NumPanels)/NumPanelCol))
    
    axLims = AxLimits()
    
    fSize = deepcopy(plt.rcParams.get('figure.figsize'))
    if allowExpand:
        fSize[0] = fSize[0]*(1 + 0.4*NumPanelCol)
        fSize[1] = fSize[1]*(1 + 0.4*NumPanelRow)
    
    if spacings is None:
        spacings = {'hspace':0.075, 'wspace': 0.15, 'left':0.075, 'right': 0.95,
                 'top': 0.95}
       
        if commonAxisLabels:
            spacings = {'hspace':0.075, 'wspace': 0.075, 'left':0.1, 'right': 0.95,
                      'top': 0.95}        
        
    kwargs = spacings
    
    fig = plt.figure(figsize=fSize)
    gs = gridspec.GridSpec(NumPanelRow, NumPanelCol, **kwargs)
    
    for i in range(0,NumPanelRow):
        for j in range(0,NumPanelCol):
            plotNum = i*NumPanelCol + j
            
            if plotNum >= len(plotList):
                break
            gp = gs[i,j]
            ax = fig.add_subplot(gp)
            
            pD = plotList[plotNum]
            
            if commonAxisLabels:
                # Switch off axes labels for individual panels and remove tick 
                # labels except for xticks on bottom row and yticks on left column
                pD.axSet.ylabelOn = False
                pD.axSet.xlabelOn = False
                if j > 0:
                    pD.axSet.ytickOn = False
                if i < NumPanelRow-1:
                    pD.axSet.xtickOn = False
                    
            if pD.ylim is not None:
                ax.set_ylim(pD.ylim)
            if pD.xlim is not None:
                ax.set_xlim(pD.xlim)
                
#            ax.tick_params(axis='both', which='minor', labelsize=10)
#            ax.tick_params(axis='both', which='major', labelsize=10)
                
            PlotPanel(ax, pD, commonLegend, verbose=True)
            
            if commonLegend:
                # Display legend in top right plot
                if j == NumPanelCol-1 and i == 0:
                    ax.legend(loc=0)
                
            # Record each panels y and x limits
            axLims.AddLim(ax, ax.get_xlim(), ax.get_ylim())
                
                
    # Set all axes limits the same  
    if commonXScale:      
        axLims.SetAlltoCommonLimits()
        
    if commonAxisLabels:
        pD = plotList[0]
        fig.text(0.025, 0.5, pD.axSet.ylabel, rotation="vertical", va="center", 
                  fontsize=1.4*plt.rcParams.get('ytick.labelsize'))
        fig.text(0.5, 0.05, pD.axSet.xlabel, ha="center", 
                  fontsize=1.4*plt.rcParams.get('xtick.labelsize'))
    PlotorShow_PDFPNG(fileN)
    
def PlotorShow_PDFPNG(fileN, plotpath=''):
    if fileN is None:
        plt.show()
    else:
        print('Saving to: '+plotpath+fileN+'.pdf')
        plt.savefig(plotpath+fileN+'.pdf')
        print('Saving to: '+plotpath+fileN+'.png')
        plt.savefig(plotpath+fileN+'.png')
        plt.close() 
###############################################################################
###############################################################################
    
if __name__ == "__main__": 
    pass