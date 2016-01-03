#!/usr/bin/env python

# filename: chiSqAnalysis_Glb.py
# Purpose: plot density, 2-2 correlation and heat map of highly dimensional 
#          parameters, in analog with scatterMatrix in GGally of R.

import numpy as np
import matplotlib.pyplot as plt

# run constant
init_model = "hasfs_glb" # kln, glb
inputFile_name = "%s_tet_chiSq_4params.dat"%init_model
chiSq_criteria = 90 # 5 DOF: v2, v3, pion <pT>, proton <pT>, kaon <pT>
taus_bd = [0.1, 3]
etas_bd = [0,  0.21]
tdec_bd = [0.13, 0.17]
visbulknorm_bd = [0,4]
chisq_bd= [0, 850]
exp_factor = 1e5
init_model_label = "Glb w/ pre-eq." if(init_model=="hasfs_glb") else "Glb w/o pre-eq."

dof=3

# load in data
raw_data = np.genfromtxt(inputFile_name, missing_values='NA')
# unpack data
chisq_array= raw_data[:, 4]

# select according to chi squared criteria
data_selected = raw_data[raw_data[:,-1]<=chiSq_criteria,:] 
taus_selected_array = data_selected[:, 0]
etas_selected_array = data_selected[:, 1]
visbulknorm_selected_array = data_selected[:, 3]
chisq_selected_array= data_selected[:, 4]

# plot constants
plotfontsize = 18
plot_label_list=['(a)', '(b)', '(c)', '(d)']

def plotDensity(x_array, x_bd, gaussian_width, binwidth, plot_object, **kwargs):
	"""
	plot out the distribution of one parameter 
	"""
	#from scipy.stats import gaussian_kde
	# create gaussian density object
	#density = gaussian_kde(x_array)
	# adjust the width of gaussian
	#density.covariance_factor = lambda : gaussian_width
	#density._compute_covariance()
	# prepare data to plot
	x_l, x_h = x_bd
	#xs = np.linspace(x_l, x_h, 200)
	#hplot = plot_object.plot(xs, density(xs), **kwargs)
	plot_object.hist(x_array, bins=np.arange(x_l, x_h+binwidth, binwidth),
		color='red')
	#return hplot


def plotScatter(x_array, y_array, weight_array, xy_boundary, plot_object, **kwargs):
	"""
	Input: weight_array is the size of scattered symbol
	scatter 2D plot with size identifying the size of chi Square.
	Larger size means smaller chi square
	"""

	plot_object.scatter(x_array, y_array, 
		edgecolors='w', s=weight_array, facecolors='k')#, alpha=0.7)
#		c='black', s=weight_array, edgecolors='none', alpha=0.7)
	# set boundary limit
	x_l_bd, x_r_bd, y_l_bd, y_r_bd = xy_boundary
	plot_object.set_xlim(x_l_bd, x_r_bd)
	plot_object.set_ylim(y_l_bd, y_r_bd)

def genWeight(chisq_array, scale_factor):
	"""
	Custom function which generate the weight from chi square:
	smaller the chi square, larger the size.
	"""
	temp = 1.0/chisq_array
#	temp_scaled = (temp-np.min(temp))/np.std(temp)*scale_factor
	temp_scaled = scale_factor*(2.0/(1+chisq_array/dof))**4
	return temp_scaled


# attempt panel plot
fig, axList = plt.subplots(2, 2, sharex = False, sharey = False)

# plot distribution
plotDensity(chisq_array, chisq_bd, 0.2, 50, axList[0,0])
axList[0,0].set_xlabel(r'$\chi^2$', fontsize=plotfontsize+8)
axList[0,0].set_ylabel(r'count', fontsize=plotfontsize+8)
axList[0,0].set_ylim(0,300)
axList[0,0].text(0.26, 0.85, init_model_label,
    fontsize=plotfontsize+4,
    transform = axList[0,0].transAxes)
axList[0,0].set_xticks(np.array([0, 200, 400, 600, 800]))

# plot visbulknorm vs. eta/s
chisq_weight = genWeight(chisq_selected_array, exp_factor)
plotScatter(etas_selected_array, visbulknorm_selected_array, chisq_weight, 
	etas_bd+visbulknorm_bd, axList[0,1])
axList[0,1].set_xlabel(r'$\eta/s$', fontsize=plotfontsize+8)
axList[0,1].set_ylabel(r'bulk norm', fontsize=plotfontsize+8)
#axList[0,1].set_yticks(np.array([100, 120, 140, 160]))
axList[0,1].set_xticks(np.linspace(0, 0.25, 6))
axList[0,1].set_xticklabels(['%g'%i for i in np.linspace(0, 0.25, 6)])

# plot taus vs. eta/s
plotScatter(taus_selected_array, etas_selected_array, chisq_weight, 
	taus_bd+etas_bd, axList[1,0])
axList[1,0].set_xlabel(r'$\tau_s$ (fm)', fontsize=plotfontsize+8)
axList[1,0].set_ylabel(r'$\eta/s$', fontsize=plotfontsize+8)
axList[1,0].set_yticks(np.linspace(0, 0.25, 6))
axList[1,0].set_yticklabels(['%g'%i for i in np.linspace(0, 0.25, 6)])

# plot taus  vs. visbulknorm
plotScatter(taus_selected_array, visbulknorm_selected_array, chisq_weight, 
	taus_bd+visbulknorm_bd, axList[1,1])
axList[1,1].set_xlabel(r'$\tau_s$ (fm)', fontsize=plotfontsize+8)
axList[1,1].set_ylabel(r'bulk norm', fontsize=plotfontsize+8)
#axList[1,1].set_yticks(np.array([100, 120, 140, 160]))
# adjust tick label size
for i in [0,1]:
	for j in [0,1]:
		axList[i,j].tick_params('both',length=6,width=1,which='major',labelsize=plotfontsize+2)

i=0
for a in axList[0,:]:
	a.text(0.9, 0.05, plot_label_list[i],
		    fontsize=plotfontsize,
    		transform = a.transAxes)
	i+=1
for a in axList[1,:]:
	a.text(0.9, 0.05, plot_label_list[i],
		    fontsize=plotfontsize,
    		transform = a.transAxes)
	i+=1

# adjust figure
fig.subplots_adjust(left = 0.11, bottom = 0.1,right = 0.97, 
                    top = 0.98, hspace=0.3, wspace=0.3)
fig.set_size_inches(10,8)


plt.savefig('chiSq_analysis_%s_4params_1.pdf'%init_model)
plt.show()