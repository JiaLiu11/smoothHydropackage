#!/usr/bin/env python

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys

# set plotting parameters
fontsize = 16
sns.set(style="white", color_codes=True)
sns.set_style("ticks", {"xtick.major.size": 4, "ytick.major.size": 4})
rc={'font.size': fontsize, 'axes.labelsize': fontsize, 'legend.fontsize': fontsize, 
    'axes.titlesize': fontsize, 'xtick.labelsize': fontsize-2, 'ytick.labelsize': fontsize-2}
sns.set_context(rc=rc)
#sns.set_style("darkgrid", {"axes.facecolor": ".9"})

# set parameter ranges
taus_bd  = (0.1, 3)
etas_bd  = (0, 0.21)
tsw_bd   = (130, 170)
bNorm_bd = (0,4)
taus_bins = np.arange(taus_bd[0], taus_bd[1]+0.1, 0.1)
etas_bins = np.arange(etas_bd[0], etas_bd[1]+0.01, 0.01)
tsw_bins  = np.arange(tsw_bd[0], tsw_bd[1]+2.5, 2.5)
bNorm_bins= np.arange(bNorm_bd[0], bNorm_bd[1]+0.05, 0.05)

# load data
params_table = np.loadtxt('params_list_glb_4params.dat')
tau_s = params_table[:,0]
eta_s = params_table[:,1]
t_sw  = params_table[:,2]
visbulkNorm = params_table[:,3]


# # start to plot tau/s v.s. eta/s
# g = sns.JointGrid(tau_s, eta_s,
#     xlim=taus_bd, ylim=etas_bd)
# g.ax_marg_x.hist(tau_s, bins=taus_bins, 
#     alpha=0.7, edgecolor="white")
# g.ax_marg_y.hist(eta_s, bins=etas_bins, 
#     alpha=0.7, orientation="horizontal", edgecolor="white")
# g.plot_joint(plt.scatter, 
#     s=25,  edgecolor="white", cmap='Blues')
# g.set_axis_labels(r'$\tau_s$',r'$\eta/s$', fontsize=fontsize+4)
# g.fig.savefig("taus_etas_glb.pdf", bbox_inches="tight")
# g.fig.clear()

# # start to plot tau/s v.s. bulk norm
# g = sns.JointGrid(tau_s, visbulkNorm,
#     xlim=taus_bd, ylim=bNorm_bd)
# g.ax_marg_x.hist(tau_s, bins=taus_bins, 
#     alpha=0.7, edgecolor="white")
# g.ax_marg_y.hist(visbulkNorm, bins=bNorm_bins, 
#     alpha=0.7, orientation="horizontal", edgecolor="white")
# g.plot_joint(plt.scatter, 
#     s=25,  edgecolor="white", cmap='Blues')
# g.set_axis_labels(r'$\tau_s$',r'bulk norm.', fontsize=fontsize+4)
# g.fig.savefig("taus_bulkNorm_glb.pdf", bbox_inches="tight")
# g.fig.clear()

# # start to plot eta_s v.s. bulk norm
# g = sns.JointGrid(eta_s, visbulkNorm,
#     xlim=etas_bd, ylim=bNorm_bd)
# g.ax_marg_x.hist(eta_s, bins=etas_bins, 
#     alpha=0.7, edgecolor="white")
# g.ax_marg_y.hist(visbulkNorm, bins=bNorm_bins, 
#     alpha=0.7, orientation="horizontal", edgecolor="white")
# g.plot_joint(plt.scatter, 
#     s=25,  edgecolor="white", cmap='Blues')
# g.set_axis_labels(r'$\eta/s$',r'bulk norm.', fontsize=fontsize+4)
# g.fig.savefig("etas_bulkNorm_glb.pdf", bbox_inches="tight")

# start to plot eta_s v.s. t_sw
g = sns.JointGrid(eta_s, t_sw,
    xlim=etas_bd, ylim=tsw_bd)
g.ax_marg_x.hist(eta_s, bins=etas_bins, 
    alpha=0.7, edgecolor="white")
g.ax_marg_y.hist(t_sw, bins=tsw_bins, 
    alpha=0.7, orientation="horizontal", edgecolor="white")
g.plot_joint(plt.scatter, 
    s=25,  edgecolor="white", cmap='Blues')
g.set_axis_labels(r'$\eta/s$',r'$T_{sw}$', fontsize=fontsize+4)
g.fig.savefig("etas_tsw_glb.pdf", bbox_inches="tight")

# start to plot eta_s v.s. t_sw
g = sns.JointGrid(visbulkNorm, t_sw,
    xlim=bNorm_bd, ylim=tsw_bd)
g.ax_marg_x.hist(visbulkNorm, bins=bNorm_bins, 
    alpha=0.7, edgecolor="white")
g.ax_marg_y.hist(t_sw, bins=tsw_bins, 
    alpha=0.7, orientation="horizontal", edgecolor="white")
g.plot_joint(plt.scatter, 
    s=25,  edgecolor="white", cmap='Blues')
g.set_axis_labels(r'bulk norm.',r'$T_{sw}$', fontsize=fontsize+4)
g.fig.savefig("bulkNorm_tsw_glb.pdf", bbox_inches="tight")
