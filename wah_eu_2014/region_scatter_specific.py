#!/usr/local/bin/python2.7
# NOTE python2.7 must be used for this script

import glob,os
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle
from mpl_toolkits.basemap import Basemap,cm
from region_returntimes import get_region_returntimes2, remap_p5deg_to_rotated, find_closest
from scipy.io.netcdf import netcdf_file

pkl_dir='/gpfs/projects/cpdn/scratch/cenv0437/pickle_data/'

if __name__=='__main__':

# Load model data from pickle files
	fdata=open(pkl_dir+'historical_data.pkl','rb')
	historical_data=pickle.load(fdata)
	fdata.close()
	fdata=open(pkl_dir+'natural_data.pkl','rb')
	natural_data=pickle.load(fdata)
	fdata.close()
	fdata=open(pkl_dir+'clim_hist_data.pkl','rb')
	clim_hist_data=pickle.load(fdata)
	fdata.close()
	fdata=open(pkl_dir+'clim_nat_data.pkl','rb')
	clim_nat_data=pickle.load(fdata)
	fdata.close()
	print 'loaded data from pkl files'

	obs_p5deg='/home/cenv0437/scratch/data_from_ouce/CRU_TS_Absolute_plus_BEST_anomaly_201312_201412.nc'
	obsname='BEST'
	#obs_p5deg='/home/cenv0437/scratch/data_from_ouce/CRU_TS_Absolute_plus_CRU4_anomaly_201312_201412.nc'
#	obsname='CRU4'

	#obs_p5deg='/home/cenv0437/scratch/data_from_ouce/CRU_TS_Absolute_plus_GISS_anomaly_201312_201412.nc'
#	obsname='GISS'



	# Some specific regions:
#	regions=[(51.2, 13.8)]
	regions=[(62.561671, 13.958159)]
	regions=[(51.749181, -1.281326)]

# Cacluate return times for many regions

	print historical_data.shape	
	nens,nlat,nlon=historical_data.shape
	hist=np.zeros([len(regions),50])
	nat=np.zeros([len(regions),50])
	clim_hist=np.zeros([len(regions),50])
	clim_nat=np.zeros([len(regions),50])
	lp=np.zeros([len(regions),50])
	hist_var=np.zeros([len(regions),50]) #variance of hist ensemble 
	nat_var=np.zeros([len(regions),50]) # variance of nat ensemble
	

	# Regrid obs to rotated regional grid
	rot_template='/home/cenv0437/scratch/data_from_ouce/hadam3p_eu_z4ao_2013_1_009238311_0_tasmean.nc'
	f_rot=netcdf_file(rot_template,'r')
	lat_coord=f_rot.variables['global_latitude0'][4:-7,4:-4]
	lon_coord=f_rot.variables['global_longitude0'][4:-7,4:-4]
	f_rot.close()	
	obs=remap_p5deg_to_rotated(obs_p5deg,rot_template)[:,:]
		

	for reg,rval in enumerate(regions):
		j0,i0=find_closest(rval[0],rval[1],lat_coord,lon_coord)
		print 'centrepoint:',lat_coord[j0,i0],lon_coord[j0,i0]
		for n in range(50): #width/height of domain
			print n,
			try:
					vals=get_region_returntimes2(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,i0-n,i0+n+1,j0-n,j0+n+1)
			except Exception,e:
				print e
				break
			if vals[4]>0: # If there are any land points to analyse
				hist[reg,n]=vals[0]
				nat[reg,n]=vals[1]
				clim_hist[reg,n]=vals[2]
				clim_nat[reg,n]=vals[3]
				lp[reg,n]=vals[4]
				hist_var[reg,n]=vals[5]
				nat_var[reg,n]=vals[6]

###############################################
# Plotting	

#First get rid of annoying Mac hidden files
	delfiles=glob.glob('scatter_figs/._*.png')
	for f in delfiles:
		os.remove(f)
		
# Some data manipulation...
	lp=lp.flatten()
	hist=hist.flatten()
	nat=nat.flatten()
	clim_hist=clim_hist.flatten()
	clim_nat=clim_nat.flatten()
	

	
	condition=np.logical_or(np.logical_not(np.isfinite(hist)),hist>1000)
	hist=np.ma.masked_where(condition,hist)
	

#	condition=np.logical_or(np.logical_not(np.isfinite(nat)),nat>4000)
	print (~np.isfinite(nat)).sum()
	nat=np.ma.masked_where(nat==np.NaN,nat) # mask out any possible nan's
	print (~np.isfinite(nat)).sum()
	nat[~np.isfinite(nat)]=natural_data.shape[0]	
	print (~np.isfinite(nat)).sum()
	print nat.min(),nat.max()
	condition=np.logical_or(np.logical_not(np.isfinite(clim_hist)),clim_hist>4000)
	clim_hist=np.ma.masked_where(condition,clim_hist)
	
#	condition=np.logical_or(np.logical_not(np.isfinite(clim_nat)),clim_nat>4000)

	print (~np.isfinite(clim_nat)).sum()
	clim_nat=np.ma.masked_where(clim_nat==np.NaN,clim_nat) # mask out any possible nan's
	print (~np.isfinite(clim_nat)).sum()
	clim_nat[~np.isfinite(clim_nat)]=clim_nat_data.shape[0]
	print (~np.isfinite(clim_nat)).sum()
	print clim_nat.min(),clim_nat.max()


	print '###########################################################\n'

	plt.figure(1)
	plt.plot(lp**0.5*50.,1/nat.flatten(),'.',label='nat')
#	plt.plot(lp,clim_hist.flatten(),'g.',label='clim_hist')
#	plt.plot(lp,hist.flatten(),'.b',label='hist')
	plt.semilogy(lp**0.5*50.,1/clim_nat.flatten(),'.',label='clim_nat')
#	plt.title('p0 vs size of region')
#	plt.legend()
	plt.ylim([0.001,1])
#	plt.xlim([0,5000])
	plt.xlabel('Scale of region (km)')
	plt.ylabel('p0')
	plt.savefig('figure_nat_specific.png')
	
	plt.plot(lp**0.5*50.,1/hist.flatten(),'.',label='hist')
#	plt.plot(lp,clim_hist.flatten(),'g.',label='clim_hist')
#	plt.plot(lp,hist.flatten(),'.b',label='hist')
	plt.semilogy(lp**0.5*50.,1/clim_hist.flatten(),'.',label='clim_hist')
	plt.title('Probablility for regions centred on: '+str((lat_coord[j0,i0],lon_coord[j0,i0])))
	plt.legend(loc='best')
	plt.ylim([0.001,1])
#	plt.xlim([0,5000])
	plt.xlabel('Scale of region (km)')
	plt.ylabel('Probability of event')
	plt.savefig('figure_hist_specific.png')
	
	plt.figure(2)
	plt.plot(lp**0.5*50.,(1-hist/nat),'.',label='2014')
#	plt.plot(lp,clim_hist.flatten(),'g.',label='clim_hist')
#	plt.plot(lp,hist.flatten(),'.b',label='hist')
	plt.plot(lp**0.5*50.,(1-clim_hist/clim_nat),'.',label='clim')
	plt.title('far vs size of region')
	plt.legend()
#	plt.ylim([0,3000])
#	plt.xlim([0,5000])
	plt.xlabel('Scale of region (km)')
	plt.ylabel('FAR')
	plt.savefig('figure_far_specific.png')
	
	