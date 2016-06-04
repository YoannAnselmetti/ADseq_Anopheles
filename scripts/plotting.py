__author__ = "Cedric Chauve"
__date__ = "May 2016"

# Plotting routines

import sys, math, numpy as np
import matplotlib, matplotlib.pyplot as plt

# Plotting a distribution of values per species in stacked bars (c = floating color modifier)
def plot_scores_distribution_per_species(T,nb_stacks,c,xlabel,title,h,show):
    ivalues=range(0,nb_stacks)
    bars={}
    bars[0]=[]
    for species in T.keys():
        bars[0].append(0.0)
    for t in ivalues:
        bars[t+1]=[]
        for species in T.keys():
            bars[t+1].append(T[species][t])
 
    colors = plt.cm.BuPu(np.linspace(0, 0.5, len(ivalues)))
    ind    = np.arange(len(T.keys()))
    lft    = sum([np.array(bars[0])])
    height = h

    for t in ivalues:
        plt.barh(ind, np.array(bars[t+1]), height=height, color=c*colors[t], left=lft)
        lft+=sum([np.array(bars[t+1])])
    plt.title(title)
    plt.yticks(ind+height/2,  T.keys())
    plt.xlabel(xlabel)    
    if show:
        plt.show()
