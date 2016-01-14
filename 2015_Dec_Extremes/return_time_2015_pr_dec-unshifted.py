#!/usr/local/bin/python2.7
#############################################################################
# Program : return_time_2014_bootstrap.py
# Author  : Peter Uhe, based on script by Sarah Sparrow (assisted by nathalie Schaller)
# Date	: 03/07/2015
# Purpose : Plot return time periods for 2015 all forcing, climatology and natural data

#############################################################################

import sys
import os
import numpy as np
import math
import datetime
import fnmatch
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import glob
from scipy.io.netcdf import netcdf_file
import matplotlib.animation as animat
import random
import time
import pickle
sys.path.append('../CADrought2015')
sys.path.append('../wah_eu_2014')
from sort_umids import choose_tasknames1
from return_time_plot import *
from region_returntimes import ret_time_bootstrap_stats

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

pkl_dir='/home/cenv0437/cpdn_analysis/2015_Dec_Extremes/pkl'

def return_time_pr(fig,datad,mecs,cols,region,fname_out=False):


	#----------------------------------------------------------------------
	# Set the plot fint size
	#---------------------------------------------------------------------
	font = {'family' : 'sans-serif','size'   : 20}
	matplotlib.rc('font', **font)


	### Set up plot
	plt.clf()
	fig.set_size_inches(8,8)
	ax = fig.add_subplot(1,1,1)

	for label,data in datad.iteritems():
		
		# Return Times
		y_data, x_data = calc_return_times(data,  direction="descending",  period=1)
		l1=ax.semilogx(x_data,y_data, marker='o',markersize=5,
				   linestyle='None',mec=mecs[label],mfc=cols[label],
				   color=cols[label],fillstyle='full',
				   label=label,zorder=5)
		# Confidence Interval	   
		conf = calc_return_time_confidences(data,direction="descending",bsn=1e3)  
		conf_5=conf[0,:].squeeze()
		conf_95=conf[1,:].squeeze()
		cl1=ax.fill_between(x_data,conf_5,conf_95,facecolor=mecs[label],edgecolor=mecs[label],alpha=0.3,linewidth=1.5,zorder=4)

	
	# Plotting Stuff
	fig.subplots_adjust(bottom=0.15)
	ax.set_ylabel("Precipitation threshold (mm/day)",fontsize=16)
	ax.set_xlabel("Chance of average precipitation greater than threshold",fontsize=16)
	plt.setp(ax.get_xticklabels(),fontsize=12)
	plt.setp(ax.get_yticklabels(),fontsize=12)
	ax.set_title(region+" Dec 2015 Precipitation\n")
	
#	ylims={'California':(4,9),'Oregon':(0,5),'Washington':(0,5)}
#	ax.set_ylim(0,10) 
	ax.set_xlim(1,1e3)
	labels=['','','1/10','1/100','1/1000']
	ax.set_xticklabels(labels)
	plt.legend(loc='lower right',numpoints=1)
	
	
	# Save figure if out_dir is defined
	if fname_out:
#		if not os.path.exists(out_dir):
#			os.mkdir(out_dir)
#		fname_out=os.path.join(out_dir,"return_time_pr_"+str(nsims).zfill(5)++"_"+region+".png")
		fig.savefig(fname_out,dpi=28.75*2)

def region_mean(data,region,lat_coord,lon_coord,lsm=None):

	# Land sea mask
	if lsm==None:
		lsm = netcdf_file('/home/cenv0437/scratch/data_from_ouce/land_mask_eu50.nc').variables['lsm'][0,0,:]+1

	# Central England
	if region=='Central England':
		region_mask=~np.logical_and(np.logical_and(lat_coord>51 ,lat_coord<54) , lon_coord>357)
	# Northern England
	if region=='Northern England':
		region_mask=~np.logical_and(np.logical_and(lat_coord>54 ,lat_coord<57) , np.logical_or(lon_coord>354,lon_coord<2))

	region_mask=np.logical_or(region_mask,lsm)
	
	region_data=np.ma.masked_where(region_mask,data)
#	region_data=np.ma.masked_values(region_data,-1.073742e+09) # mask missing values
	return region_data.mean()

# Data with time array as first dimension
def region_mean2(region,data,lat_coord,lon_coord,lsm=None):
	nt,ny,nx=data.shape
	mean_data=np.zeros([nt])

	# Land sea mask
	if lsm==None:
		lsm = netcdf_file('/home/cenv0437/scratch/data_from_ouce/land_mask_eu50.nc').variables['lsm'][0,0,:]+1

	if region=='Central England':
	# Central England
		region_mask=~np.logical_and(np.logical_and(lat_coord>51 ,lat_coord<54) , lon_coord>357)
	
	# Northern England
	if region=='Northern England':
		region_mask=~np.logical_and(np.logical_and(lat_coord>54 ,lat_coord<57) , np.logical_or(lon_coord>354,lon_coord<2))
		
	region_mask=np.logical_or(region_mask,lsm)
	
	for t in range(nt):
		mean_data[t]=np.ma.masked_where(region_mask,data[t,:]).mean()
#	region_data=np.ma.masked_values(region_data,-1.073742e+09) # mask missing values
	return mean_data

#Main controling function
def main():

	region='Central England'

	# Model climatology
	lsm = netcdf_file('/home/cenv0437/scratch/data_from_ouce/land_mask_eu50.nc').variables['lsm'][0,0,4:-7,4:-4]+1

	# Load template region file:
	fname_grid='/home/cenv0437/scratch/data_from_ouce/hadam3p_eu_template.nc'
	f_grid=netcdf_file(fname_grid)
	lat_coords=f_grid.variables['global_latitude0'][:]
	lon_coords=f_grid.variables['global_longitude0'][:]
	regional_ny,regional_nx=lon_coords.shape


	months=['2015-12','2016-01']
		
	f_pickle=pkl_dir+'/pr_results_dec-unshifted_'+region+'.pkl'
	if os.path.exists(f_pickle):
		fdata=open(f_pickle,'rb')
		data_hist=pickle.load(fdata)
		data_nat=pickle.load(fdata)
		data_clim=pickle.load(fdata)
		fdata.close()
	else:
		data_hist={}
		data_nat={}
		data_clim={}
		
	
# 				for i,monthstr in enumerate(months): # Jan to Oct
# 					if i==0:
# 						year=year0
# 					else:
# 						year=year0+1
# 					fname2=fname[:-10]+str(year)+'-'+monthstr[-2:]+'.nc'
# 					print os.path.basename(fname2)
# 					tmp[i*30:(i+1)*30]=netcdf_file(fname2,'r').variables['field90'][:,0,:]*60*60*24 # conver to mm/day	
				
	# Actual 
# 	fnames_data_hist=glob.glob('/home/cenv0437/scratch/hadam3p_eu/batch266/ga.pd/field90/field90_????_2015-11.nc')
# 	for i,fname in enumerate(fnames_data_hist):
# 		umid=os.path.basename(fname).split('_')[1]
# 		try:
# 			if umid not in data_hist.keys():
# 				tmp=np.zeros([60])
# 				for i,monthstr in enumerate(months):
# 					fname2=fname[:-10]+monthstr+'.nc'
# 					print os.path.basename(fname2)
# 					tmp[i*30:(i+1)*30]=region_mean2(region,netcdf_file(fname2,'r').variables['field90'][:,0,:],lat_coords,lon_coords)*60*60*24
# 				data_hist[umid]=tmp[:30].mean(0)
# 		except Exception,e:
# 			print "Error loading files",e
# 			continue
# 			
# 				
# 	# Natural 
# 	fnames_data_nat=glob.glob('/home/cenv0437/scratch/hadam3p_eu/batch267/ga.pd/field90/field90_????_2015-11.nc')
# 	for i,fname in enumerate(fnames_data_nat):
# 		umid=os.path.basename(fname).split('_')[1]
# 		try:
# 			if umid not in data_nat.keys():
# 				tmp=np.zeros([60])
# 				for i,monthstr in enumerate(months):
# 					fname2=fname[:-10]+monthstr+'.nc'
# 					print os.path.basename(fname2)
# 					tmp[i*30:(i+1)*30]=region_mean2(region,netcdf_file(fname2,'r').variables['field90'][:,0,:],lat_coords,lon_coords)*60*60*24
# 				data_nat[umid]=tmp[:30].mean(0)
# 		except Exception,e:
# 			print "Error loading files",e
# 			continue
# 			
# 	# Clim (all months)
# 	fnames_data_clim=glob.glob('/home/cenv0437/scratch/hadam3p_eu/batch211/ga.pd/field90/field90_????_????-12.nc')
# 	for fname in fnames_data_clim:
# 		umid=os.path.basename(fname).split('_')[1]
# 		year0=int(os.path.basename(fname).split('_')[2][:4])
# 		try:
# 			if umid not in data_clim.keys():
# 				tmp=np.zeros([60])
# 				for i,monthstr in enumerate(months):
#  					if i==0:
#  						year=year0
#  					else:
#  						year=year0+1
#  					fname2=fname[:-10]+str(year)+'-'+monthstr[-2:]+'.nc'
# 					print os.path.basename(fname2)
# 					tmp[i*30:(i+1)*30]=region_mean2(region,netcdf_file(fname2,'r').variables['field90'][:,0,:],lat_coords,lon_coords)*60*60*24
# 				data_clim[umid]=tmp[:30].mean(0)		
# 		except Exception,e:
# 			print "Error loading files",e
# 			continue


	# Save data for reuse later
	fdata=open(f_pickle,'wb')
	pickle.dump(data_hist,fdata,-1)
	pickle.dump(data_nat,fdata,-1)
	pickle.dump(data_clim,fdata,-1)
	fdata.close()
	
	mecs={}
	cols={}
	
	# Define color scheme:
	cols['2015 Actual']='orange'
	mecs['2015 Actual']='darkorange'
	
	cols['2015 Natural']="MediumSeaGreen"
	mecs['2015 Natural']="DarkGreen"
	
	cols['1985-2014 Climatology']='RoyalBlue'
	mecs['1985-2014 Climatology']='MediumBlue'
	
	cols['NaturalFix']='Purple'
	mecs['NaturalFix']='DarkMagenta'


	alpha=0.07
			
	
	
	datad={}

	num_ensembles=len(data_hist.values())
	print 'ensembles used hist',num_ensembles
	data_hist=np.array(data_hist.values())
	print data_hist.min(),data_hist.max()
	datad['2015 Actual']=data_hist
	
	num_ensembles=len(data_nat.values())
	print 'ensembles used nat',num_ensembles
	data_nat=np.array(data_nat.values())
	print data_nat.min(),data_nat.max()
	datad['2015 Natural']=data_nat
	
	num_ensembles=len(data_clim.values())
	print 'ensembles used clim',num_ensembles
	datad['1985-2014 Climatology']=np.array(data_clim.values())
	
	#TODO  Last thing to do is bias correct


	n_hist=len(datad['2015 Actual'])
	n_nat=len(datad['2015 Natural'])
	n_clim=len(datad['1985-2014 Climatology'])

	# Use 1 in 100 year event as reference
	nyears=100 
	
	# Threshold relative to actual ensemble
	tmp=datad['2015 Actual']
	tmp.sort()
	thresh=tmp[-(n_hist/nyears)]
	
	# Observed 2015 threshold
#			thresh=5.0

	condition=datad['2015 Actual']>thresh
	count_hist=condition.sum(0)/(n_hist*1.0)
	return_hist=ret_time_bootstrap_stats(datad['2015 Actual'],thresh)
	
	condition=datad['2015 Natural']>thresh
	count_nat=condition.sum(0)/(n_nat*1.0)
	return_nat=ret_time_bootstrap_stats(datad['2015 Natural'],thresh)
	
	condition=datad['1985-2014 Climatology']>thresh
	count_clim=condition.sum(0)/(n_clim*1.0)
	return_clim=ret_time_bootstrap_stats(datad['1985-2014 Climatology'],thresh)
	
	ratio_anth=(count_hist)/(count_nat)
	ratio_anom=(count_hist)/count_clim
	
	print "1 in "+str(nyears)+" event"
#			print '2015 event'
	print "Threshold",thresh
	print "Return times, actual, natural, climatology"
	print 1./count_hist,1./count_nat,1./count_clim
	print return_hist,return_nat,return_clim
	print "Relative likelihood anthropogenic",ratio_anth
	print "Relative likelihood 2015 anomaly", ratio_anom



	
	fout='/home/cenv0437/cpdn_analysis/2015_Dec_Extremes/figs/pr_return_dec-unshifted_'+region+'_uncorrected'
	# Make still image
	f_ext='.png'
	return_time_pr(plt.figure(),datad,mecs,cols,region,fname_out=fout+f_ext)
	print "saved as",fout+f_ext

	print 'Finished!'

#Washerboard function that allows main() to run on running this file
if __name__=="__main__":
  main()




