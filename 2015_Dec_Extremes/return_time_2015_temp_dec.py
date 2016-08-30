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

def return_time_temp(fig,datad,mecs,cols,region,fname_out=False):


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

	if region=='Central England':
		plt.plot([1,10000],[9.7,9.7],'m',label='Observed')
	
	# Plotting Stuff
	fig.subplots_adjust(bottom=0.15)
	ax.set_ylabel("Temperature threshold (degrees Celsius)",fontsize=16)
	ax.set_xlabel("Chance of average temperature greater than threshold",fontsize=16)
	plt.setp(ax.get_xticklabels(),fontsize=12)
	plt.setp(ax.get_yticklabels(),fontsize=12)
	ax.set_title(region+" Dec 2015 Temperature\n")
	
#	ylims={'California':(4,9),'Oregon':(0,5),'Washington':(0,5)}
	ax.set_ylim(5,10) 
	ax.set_xlim(10,1e3)
	labels=['','1/10','1/100','1/1000']
	ax.set_xticklabels(labels)
	plt.legend(loc='best',numpoints=1,frameon=False)
	
	
	# Save figure if out_dir is defined
	if fname_out:
#		if not os.path.exists(out_dir):
#			os.mkdir(out_dir)
#		fname_out=os.path.join(out_dir,"return_time_temp_"+str(nsims).zfill(5)++"_"+region+".png")
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


	monthstr='2015-12'
		
	f_pickle=pkl_dir+'/temp_results_dec_'+region+'.pkl'
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
		
		
				
	# Actual 
# 	fnames_data_hist=glob.glob('/home/cenv0437/scratch/hadam3p_eu/batch266/ga.pe/field16/field16_????_2015-11.nc')
# 	for i,fname in enumerate(fnames_data_hist):
# 		umid=os.path.basename(fname).split('_')[1]
# 		try:
# 			if umid not in data_hist.keys():
# 				tmp=np.zeros([2])
# 				fname2=fname[:-10]+monthstr+'.nc'
# 				print os.path.basename(fname2)
# 				tmp0=netcdf_file(fname2,'r').variables['field16'][0,0,:]
# 				data_hist[umid]=region_mean(tmp0,region,lat_coords,lon_coords)-273.15
# 		except Exception,e:
# 			print "Error loading files",e
# 			continue
# 			
# 				
# 	# Natural 
# 	fnames_data_nat=glob.glob('/home/cenv0437/scratch/hadam3p_eu/batch267/ga.pe/field16/field16_????_2015-11.nc')
# 	for i,fname in enumerate(fnames_data_nat):
# 		umid=os.path.basename(fname).split('_')[1]
# 		try:
# 			if umid not in data_nat.keys():
# 				fname2=fname[:-10]+monthstr+'.nc'
# 				print os.path.basename(fname2)
# 				tmp0=netcdf_file(fname2,'r').variables['field16'][0,0,:]
# 				data_nat[umid]=region_mean(tmp0,region,lat_coords,lon_coords)-273.15
# 		except Exception,e:
# 			print "Error loading files",e
# 			continue
# 			
# 	# Clim (all months)
# 	fnames_data_clim=glob.glob('/home/cenv0437/scratch/hadam3p_eu/batch211/ga.pe/field16/field16_????_????-12.nc')
# 	for fname in fnames_data_clim:
# 		umid=os.path.basename(fname).split('_')[1]
# 		try:
# 			if umid not in data_clim.keys():
# 				print os.path.basename(fname)
# 				tmp0=netcdf_file(fname,'r').variables['field16'][0,0,:]
# 				data_clim[umid]=region_mean(tmp0,region,lat_coords,lon_coords)-273.15
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
        cols['2015 Actual']='red'
        mecs['2015 Actual']='darkred'

        cols['2015 Natural']="RoyalBlue"
        mecs['2015 Natural']="MediumBlue"

        cols['1985-2013 Climatology']='orange'
        mecs['1985-2013 Climatology']='darkorange'

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
	datad['1985-2013 Climatology']=np.array(data_clim.values())
	
	#TODO  Last thing to do is bias correct
	# ERAI Climatology (from climate explorer)
	obs={'Central England':4.82657,'Northern England':4.16404}
	clim=datad['1985-2013 Climatology'].mean()
	bias=clim-obs[region]
	print 'obs,clim,bias'
	print obs[region],clim,bias
	
	for key,data in datad.iteritems():
		datad[key]=data-bias
		
	print 'adjusted clim:',datad['1985-2013 Climatology'].mean()

	n_hist=len(datad['2015 Actual'])
	n_nat=len(datad['2015 Natural'])
	n_clim=len(datad['1985-2013 Climatology'])

	# Use 1 in 100 year event as reference
	nyears=100 
	
	# Threshold relative to actual ensemble
	tmp=datad['2015 Actual']
	tmp.sort()
	thresh=tmp[-(n_hist/nyears)]
	
	# Observed 2015 threshold
#	thresh=9.7

	condition=datad['2015 Actual']>thresh
	count_hist=condition.sum(0)/(n_hist*1.0)
	return_hist=ret_time_bootstrap_stats(datad['2015 Actual'],thresh)
	
	condition=datad['2015 Natural']>thresh
	count_nat=condition.sum(0)/(n_nat*1.0)
	return_nat=ret_time_bootstrap_stats(datad['2015 Natural'],thresh)
	
	condition=datad['1985-2013 Climatology']>thresh
	count_clim=condition.sum(0)/(n_clim*1.0)
	return_clim=ret_time_bootstrap_stats(datad['1985-2013 Climatology'],thresh)
	
	ratio_anth=(count_hist)/(count_nat)
	ratio_anom=(count_hist)/count_clim
	
	print "1 in "+str(nyears)+" event"
#print '2015 event'
	print "Threshold",thresh
	print "Return times, actual, natural, climatology"
	print 1./count_hist,1./count_nat,1./count_clim
	print return_hist,return_nat,return_clim
	print "Relative likelihood anthropogenic",ratio_anth
	print "Relative likelihood 2015 anomaly", ratio_anom



	
	fout='/home/cenv0437/cpdn_analysis/2015_Dec_Extremes/figs/temp_return_dec_'+region+'_corrected'
	# Make still image
	f_ext='.png'
	return_time_temp(plt.figure(),datad,mecs,cols,region,fname_out=fout+f_ext)
#	print "saved as",fout+f_ext

	print 'Finished!'

#Washerboard function that allows main() to run on running this file
if __name__=="__main__":
  main()




