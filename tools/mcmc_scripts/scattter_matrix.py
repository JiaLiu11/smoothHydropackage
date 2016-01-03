#!/usr/bin/env python

# filename: scatter_matrix.py
# Purpose: plot density, 2-2 correlation and heat map of highly dimensional 
#          parameters, in analog with scatterMatrix in GGally of R.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from os import path
from scipy.stats import gaussian_kde
import csv

# run constant
init_model = "MCGlb_1" # kln, glb
infile_name = path.join(init_model,'mcmc_50000.csv')
outfile_name= 'params_dist_%s.png'%init_model


if init_model.split('_')[0] == "MCKLN":
    params_limit = np.array([[0.1, 2.0],[0.08, 0.28],[0, 3]])
else:
    params_limit = np.array([[0.1, 2.9],[0.0, 0.21],[0, 3]])
params_name = [r"$\tau_s$ (fm/c)", r"$\eta/s$", r"bulk norm"]

# plot constants
cmap = plt.cm.get_cmap("jet")
plotfontsize = 16
bin_number = 50

def plotDensity(x_array, x_l, x_h, gaussian_width, bins, plot_object, **kwargs):
    """
    plot out the distribution of one parameter 
    """
    # create gaussian density object
    density = gaussian_kde(x_array)
    # adjust the width of gaussian
    density.covariance_factor = lambda : gaussian_width
    density._compute_covariance()
    # prepare data to plot
    xs = np.linspace(x_l, x_h, 200)
    hplot = plot_object.plot(xs, density(xs), **kwargs)
    plot_object.hist(x_array,range=(x_l,x_h),bins=bins,normed=True, alpha=0.7)
    return hplot


def plotHeatmap(x_array, y_array, bin_num, plot_cmap, plot_object, **kwargs):
    """
    plot 2D histogram shows the joint distribution of two parameters
    """
    hplot = plot_object.hist2d(x_array, y_array, 
        bins=bin_num, cmap=plot_cmap, vmin=0, vmax=3e2, **kwargs)
    return hplot

def plotScatter(x_array, y_array, xy_boundary, plot_object, **kwargs):
    """
    scatter 2D plot 
    """
    m, b=np.polyfit(x_array, y_array,1)
    plot_object.plot(x_array, y_array, 'b.')
    plot_object.plot(x_array, m*x_array+b, 'r-', linewidth=2)
    x_l_bd, x_r_bd, y_l_bd, y_r_bd = xy_boundary
    plot_object.set_xlim(x_l_bd, x_r_bd)
    plot_object.set_ylim(y_l_bd, y_r_bd)

# load in data
data_df = pd.read_csv(infile_name, sep=',',header=True)
data_selected = data_df.values

# attempt panel plot
fig, axList = plt.subplots(3, 3, sharex = False, sharey = False)
for irow in range(3):
    for icol in range(3):
        x_data = data_selected[:,icol]
        y_data = data_selected[:,irow]
        x_left_bd = np.min(x_data) #params_limit[icol,0]
        x_right_bd= np.max(x_data) #params_limit[icol,1]
        y_left_bd = np.min(y_data) #params_limit[irow,0]
        y_right_bd= np.max(y_data) #params_limit[irow,1]
        xy_boundary = [x_left_bd, x_right_bd, y_left_bd, y_right_bd]
        # diagonal: density plot
        if(irow==icol):
            plot_data = data_selected[:,irow], 
            plotDensity(plot_data, x_left_bd, x_right_bd,
                0.15, 24, axList[irow, icol], 
                color="red", linewidth=2)
        # lower: heat map
        elif(irow>icol):
            plotHeatmap(x_data, y_data, 
              bin_number, cmap, axList[irow, icol])
            # plotScatter(x_data, y_data, xy_boundary, axList[irow, icol])
        # upper: scatter and correlation
        elif(irow<icol):
            plotScatter(x_data, y_data, xy_boundary, axList[irow, icol])

# hide tick labels
plt.setp([a.get_xticklabels() for a in axList[0,:]], visible=False)
plt.setp([a.get_xticklabels() for a in axList[1,:]], visible=False)
plt.setp([a.get_yticklabels() for a in axList[:,1]], visible=False)

#plt.setp([a.get_yticklabels() for a in axList[:,2]], visible=False)
plt.setp(axList[0,0].get_yticklabels(), visible=False)
plt.setp(axList[2,2].get_yticklabels(), visible=False)

# # adjust xticks
# axList[2,0].xaxis.set_ticks(np.linspace(params_limit[0,0],params_limit[0,1],5))
# axList[0,1].xaxis.set_ticks(np.linspace(params_limit[1,0],params_limit[1,1],6))
# axList[2,1].xaxis.set_ticks(np.linspace(params_limit[1,0],params_limit[1,1],6))
# axList[2,2].xaxis.set_ticks(np.linspace(params_limit[2,0],params_limit[2,1],8))

# # adjust yticks
# axList[1,0].yaxis.set_ticks(np.linspace(params_limit[1,0],params_limit[1,1],6))
# axList[2,0].yaxis.set_ticks(np.linspace(params_limit[2,0],params_limit[2,1],8))
# axList[0,2].yaxis.set_ticks(np.linspace(params_limit[0,0],params_limit[0,1],5))
# axList[1,2].yaxis.set_ticks(np.linspace(params_limit[1,0],params_limit[1,1],6))
axList[0,2].yaxis.set_ticks_position('right')
axList[1,2].yaxis.set_ticks_position('right')
# assign axis labels and allign y labels
i=0
for a in axList[2,:]:
    a.set_xlabel(params_name[i],
        {'fontsize': plotfontsize+6},labelpad=15) 
    a.tick_params('both',length=6,width=1,which='major',labelsize=plotfontsize-5.5)    
    i +=1
i=0
for a in axList[:,0]:
    a.set_ylabel(params_name[i], 
        {'fontsize': plotfontsize+6},labelpad=-5)    
    a.tick_params('both',length=6,width=1,which='major',labelsize=plotfontsize-5.5)    
    a.yaxis.set_label_coords(-0.2, 0.5)
    i +=1
for a in axList[0:2,2]:
    a.tick_params('both',length=6,width=1,which='major',labelsize=plotfontsize-5.5)

fig.subplots_adjust(left = 0.10, bottom = 0.10,right = 0.92, 
                    top = 0.98, hspace=0.15, wspace=0.15)
fig.set_size_inches(10,8)


plt.savefig(outfile_name, dpi=200)
#plt.show()