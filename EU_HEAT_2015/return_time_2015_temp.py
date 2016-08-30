#############################################################################
# Program : return_time_2014_bootstrap.py
# Author  : Sarah Sparrow (assisted by nathalie Schaller)
# Date	: 05/03/2014
# Purpose : Plot return time periods for 2014 all forcing and natural data
# Updates : 02/04/2014 Updated from return_time_2014.py to include bootstrap 
#		   error bars
#		   25/04/2014 Updated to include CMIP5 model bootstraps
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
import pickle
from scipy.stats import gumbel_r

EU = os.path.expanduser
sys.path.append('../CADrought2015')
from return_time_plot import *

sys.path.append('../wah_eu_2014')
from sort_umids import choose_tasknames1,choose_mask1
from region_returntimes import plot_region,find_closest,find_closest_1d

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

pkl_dir='/gpfs/projects/cpdn/scratch/cenv0437/pickle_data/'


def return_time_temp(hist_data,nat_data,days_aved,average,reg):

	print hist_data.shape
	print nat_data.shape
	if len(hist_data.shape)>1:
		period=hist_data.shape[1]
	else:
		period=1
		
#	hist_data=hist_data.flatten()
#	nat_data=nat_data.flatten()

	#----------------------------------------------------------------------
	# Set the plot fint size
	#---------------------------------------------------------------------
	font = {'family' : 'sans-serif',
			'size'   : 20}

	matplotlib.rc('font', **font)

	#--------------------------------------------------------------------------
	# Define model names for files and type
	#--------------------------------------------------------------------------
	types=["all_forcings","natural"]
	
#	season_col_allf="RoyalBlue" 
#	allf_mec="MediumBlue"
	season_col_allf='orange'
	allf_mec='darkorange'
	
	season_col_nat="MediumSeaGreen"
	nat_mec="DarkGreen"
	
	season_col_clim='RoyalBlue'
	clim_mec='MediumBlue'

	fig = plt.figure()
	fig.set_size_inches(8,8)
	ax = fig.add_subplot(1,1,1)
	fig.subplots_adjust(bottom=0.15)
	ax.set_ylabel("Temperature threshold (degrees Celcius)",fontsize=16)
	ax.set_xlabel("Chance of "+str(days_aved)+" day temperature greater than threshold",fontsize=16)
	plt.setp(ax.get_xticklabels(),fontsize=12)
	plt.setp(ax.get_yticklabels(),fontsize=12)


	# Hist data
	y_data_all, x_data_all = calc_return_times(hist_data,direction="descending",period=period)
	print 'histy',y_data_all.min(),y_data_all.max()
	print 'histx',x_data_all.min(),x_data_all.max()
	conf_all = calc_return_time_confidences(hist_data,direction="descending",bsn=1e2)  
	l1=ax.semilogx(x_data_all,y_data_all, marker='o',markersize=.1,
				   linestyle='None',mec=allf_mec,mfc=season_col_allf,
				   color=season_col_allf,fillstyle='full',
				   label="Actual",zorder=5)
	conf_all_5=conf_all[0,:].squeeze()
	conf_all_95=conf_all[1,:].squeeze()
	cl1=ax.fill_between(x_data_all,conf_all_5,conf_all_95,facecolor=allf_mec,edgecolor=allf_mec,alpha=0.3,linewidth=1.5,zorder=4)


	# Nat data
	y_data_all, x_data_all = calc_return_times(nat_data,  direction="descending",  period=period)
	conf_all = calc_return_time_confidences(nat_data,direction="descending",bsn=1e2)  
	l1=ax.semilogx(x_data_all,y_data_all, marker='o',markersize=.1,
				   linestyle='None',mec=nat_mec,mfc=season_col_nat,
				   color=season_col_nat,fillstyle='full',
				   label="Natural",zorder=5)
	conf_all_5=conf_all[0,:].squeeze()
	conf_all_95=conf_all[1,:].squeeze()
	cl1=ax.fill_between(x_data_all,conf_all_5,conf_all_95,facecolor=nat_mec,edgecolor=nat_mec,alpha=0.3,linewidth=1.5,zorder=4)
	

	nens_hist=hist_data.shape[0]
	nens_nat=nat_data.shape[0]
	nsims=nens_hist +nens_nat # Plus number from other simulations
	ax.set_title(reg+" JJA 2014 "+str(days_aved)+" day ave of "+average+" Temp")

	ax.set_ylim(20,50) #0,450
	ax.set_xlim(1,1e3) 
	labels=['','','1/10','1/100','1/1000']
	ax.set_xticklabels(labels)

	plt.legend(loc='lower right',markerscale=50,numpoints=1)
	fname_out="return_time_"+reg+"_"+str(days_aved)+"day"+average+".png"
#	fig.savefig(fname_out,dpi=28.75*2)
	if not os.path.exists('plots/'+average):
		os.mkdir('plots/'+average)
	fig.savefig('plots/'+average+'/'+fname_out)


# Load ensemble of files into a single array
def load_timeseries(filenames,region,months,bias):
	data=np.ma.zeros([len(filenames),len(months)*30])
	tmp=np.zeros([len(months)*30])
	i=0
	[j_s,j_e,i_s,i_e]=region
	for f in filenames:
#		print i,os.path.basename(f)
		try:
			f2=f.replace('field16','field16_1')
			for j,monthstr in enumerate(months):
				f_month=f[:-10]+monthstr+'.nc'
				f2_month=f2[:-10]+monthstr+'.nc'
				var1=netcdf_file(f_month,'r').variables['field16'][:,0,4:-7,4:-4]
				var1=np.ma.masked_values(var1,-1.07374e+09)
				var2=netcdf_file(f2_month,'r').variables['field16_1'][:,0,4:-7,4:-4]
				var2=np.ma.masked_values(var2,-1.07374e+09)
				tmp[j*30:(j+1)*30]=((var1+var2)/2.-bias)[:,j_s:j_e,i_s:i_e].mean(1).mean(1)
			
			if tmp.max()>350.0 or tmp.min()<170 or not np.all(np.isfinite(tmp)):
				print 'error: wierd vals',f
				continue
			else:
#				print tmp.min(),tmp.max()
				data[i,:]=tmp
				i=i+1
		except:
			print 'Error, cannot load files',f
			raise
			#continue
	return data[:i,:]
	
# Load ensemble of files into a single array
def load_timeseries2(filenames,field,regions,months,bias):
	data=np.ma.zeros([len(filenames),len(months)*30,len(regions)])
	tmp=np.zeros([len(months)*30,len(regions)])
	i=0
	for f in filenames:
#		print i,os.path.basename(f)
		try:
			for j,monthstr in enumerate(months):
				f_month=f[:-10]+monthstr+'.nc'
				# Load variable
				var1=netcdf_file(f_month,'r').variables[field][:,0,4:-7,4:-4]
				# Mask out missing values
				var1=np.ma.masked_values(var1,-1.07374e+09)
				for k,region in enumerate(regions):
					x=np.array(region).shape
					if len(x)>1: # Mask for region
						for t in range(30):
							tmp[j*30+t,k]=np.ma.masked_where(~region,var1[t,:]).mean()
					elif x[0]==4: # Box of grid points
						tmp[j*30:(j+1)*30,k]=var1[:,region[0]:region[1],region[2]:region[3]].mean(1).mean(1)
					elif x[0]==2: # Single point
						tmp[j*30:(j+1)*30,k]=var1[:,region[0],region[1]]
						
			
			if tmp.max()>350.0 or tmp.min()<170 or not np.all(np.isfinite(tmp)):
				print 'error: wierd vals',f
				continue
			else:
#				print tmp.min(),tmp.max()
				data[i,:]=tmp
				i=i+1
		except IndexError,e:
			print 'Index Error file',f
			print e
			continue
		except IOError, e:
			print 'io Error file',f
			print e
                        continue
		except:
			print 'Error, cannot load file',f
			raise
	return data[:i,:]

def replace_field(filename):
	return filename.replace('field16','field16_1')

#Main controling function
def main():



###############  Model climatological bias
	bias_file='/home/cenv0437/scratch/data_from_ouce/EU_JJA_temp-cruts_0.44_mean.nc'
	tmp=netcdf_file(bias_file).variables['temp_bias']
	bias=np.ma.masked_values(tmp[0,0,:,:],tmp.missing_value)

#	print 'bias mask',bias.mask.sum()
	
############ Get rotated grid info
	rot_template='/home/cenv0437/scratch/data_from_ouce/hadam3p_eu_z4ao_2013_1_009238311_0_tasmean.nc'
	f_rot=netcdf_file(rot_template,'r')
	lat_coord=f_rot.variables['global_latitude0'][4:-7,4:-4]
	lon_coord=f_rot.variables['global_longitude0'][4:-7,4:-4]
	f_rot.close()

	# Select region:
	reg_latlon=np.logical_and ( np.logical_or(lon_coord>=355.,lon_coord<20) , np.logical_and(42<=lat_coord,lat_coord<=55.) )
	region_grid=[42,69,37,64]
	plot_region(reg_latlon,lat_coord,lon_coord,0,-1,0,-1,'Lon_coords')

	region_names=['EU-latlon','EU-grid','de Bilt','Oxford','Paris','Madrid','Berlin','Zurich','Lyon','Mannheim']
	# Missed Lyon: [45.76,4.83]
	obsdata_hand=[32.5,27.5,34.5,39.3,31.8,30.8,37.4]
	
	points_latlon=[[52.1,5.2],[51.75,-1.3],[48.9,2.3],[40.4,-3.7],[52.5,13.4],[47.4,8.4],[45.76,4.83],[49.5, 8.46]]
	regions=[reg_latlon,region_grid]
	for coord in points_latlon:
		regions.append(find_closest(coord[0],coord[1],lat_coord,lon_coord))
	for i,region in enumerate(regions[2:]):
		plot_region(lat_coord,lat_coord,lon_coord,region[0],region[0]+1,region[1],region[1]+1,region_names[i+2])
	months=['2014-06','2014-07','2014-08']
	

#################  Model Data:

	read_data=False
	if read_data or not os.path.exists(pkl_dir+'EU_heat_max.pkl'):

		infiles=glob.glob('/home/cenv0437/scratch/hadam3p_eu/batch100/ga.pd/field16/field16_????_2014-06.nc')
		# Have to choose ensemble members by umid
		historical_files=choose_tasknames1(infiles,'z200','z2gn')+choose_tasknames1(infiles,'z4c0','z4sn')
		natural_files=[x for x in infiles if x not in historical_files]

		print 'len historical batch100',len(historical_files)
		historical_files=historical_files+glob.glob('/home/cenv0437/scratch/hadam3p_eu/batch166/ga.pd/field16/field16_????_2014-06.nc') #addition from batch 166

	
#################################################
# Load Data	

	# Load historical data into single array and bias correct
		print '2014 historical files:',len(historical_files)
		hist_data_max=load_timeseries2(historical_files,'field16',regions,months,bias)
		print 'loaded all forcings 2014 data... max'
		historical_files_min=map(replace_field,historical_files)
		hist_data_min=load_timeseries2(historical_files_min,'field16_1',regions,months,bias)
		print 'all forcings 2014 data... min'
		
	# Load natural simulation data into single array and bias correct
		print '2014 natural files:',len(natural_files)
		nat_data_max= load_timeseries2(natural_files,'field16',regions,months,bias)
		print 'loaded natural 2014 data... max'
		natural_files_min=map(replace_field,natural_files)
		nat_data_min=load_timeseries2(natural_files_min,'field16_1',regions,months,bias)
		print 'natural 2014 data... min'

		
#		# Write data to pickle
#		fdata=open(pkl_dir+'EU_heat.pkl','wb')
#		pickle.dump(hist_data,fdata,-1)
#		pickle.dump(nat_data,fdata,-1)
#		fdata.close()
		
		# Write data to pickle
		fdata=open(pkl_dir+'EU_heat_max.pkl','wb')
		pickle.dump(hist_data_max,fdata,-1)
		pickle.dump(nat_data_max,fdata,-1)
		fdata.close()

		# Write data to pickle
		fdata=open(pkl_dir+'EU_heat_min.pkl','wb')
		pickle.dump(hist_data_min,fdata,-1)
		pickle.dump(nat_data_min,fdata,-1)
		fdata.close()


	else: #load from pickle
		fdata=open(pkl_dir+'EU_heat_max.pkl','rb')
		hist_data_max=pickle.load(fdata)
		nat_data_max=pickle.load(fdata)
		fdata.close()
		
		fdata=open(pkl_dir+'EU_heat_min.pkl','rb')
		hist_data_min=pickle.load(fdata)
		nat_data_min=pickle.load(fdata)
		fdata.close()
		
		print 'loaded data from pkl files'

	#################  Load obs / forecast

	fname_obs='/home/cenv0437/scratch/data_from_ouce/erai_tmax_daily_e_12_max_25_3v_july7.nc'
	#fname_obs='/home/cenv0437/scratch/data_from_ouce/erai_tmax_daily_e_12_max_5_3v.nc'
#	fname_obs='/home/cenv0437/scratch/data_from_ouce/erai_tmax_daily_e_12_mean_25.nc'
	f_obs=netcdf_file(fname_obs,'r')
	obs=np.ma.masked_values(f_obs.variables['max_tmax'][:],3.e33)-273.15
	print obs.shape
	print obs[-12:,53,8]
	obslat=f_obs.variables['lat'][:]
	obslon=f_obs.variables['lon'][:]
	obs_indices=[]
	for i,latlon in enumerate(points_latlon):
		tmp=find_closest_1d(latlon[0],latlon[1]%360,obslat,obslon)
		obs_indices.append(tmp)
		print region_names[i+2],latlon,tmp[0],tmp[1],obslat[tmp[0]],obslon[tmp[1]],obs[436:439,tmp[0],tmp[1]]

	
	for region in range(len(region_names[2:])):
#	region=4

		average='season_max'
		print region_names[region+2]
		xcoord=regions[region+2][0]
		ycoord=regions[region+2][1]
		hist_data=hist_data_max[:,:,region+2]
		nat_data=nat_data_max[:,:,region+2]
	
	#	average='min'
	#	hist_data=hist_data_min[:,:,region] 
	#	nat_data=nat_data_min[:,:,region]
	
	
#		print "Data shape:"
#		print 'hist',hist_data.shape
#		print 'nat',nat_data.shape
	
#		print "Data range"	
#		print 'hist',hist_data.min(),hist_data.mean(),hist_data.max()
#		print 'nat',nat_data.min(),nat_data.mean(),nat_data.max()


	# 1 day ave
#		print 'One day'
#		days_aved=1
#		hist_data_xday=hist_data - 273.15
#		nat_data_xday=nat_data - 273.15

#		print 'hist_xday',hist_data_xday.min(),hist_data_xday.mean(),hist_data_xday.max()
#		print 'nat_xday',nat_data_xday.min(),nat_data_xday.mean(),nat_data_xday.max()
#		# Seasonmax
#		return_time_temp(hist_data_xday.max(1),nat_data_xday.max(1),days_aved,average,region_names[region])
		# All data
#		return_time_temp(hist_data_xday,nat_data_xday,days_aved,average,region_names[region])

	# 3 day ave
		average='seasonmax_tmax3'
		days_aved=3
		hist_data=hist_data_max[:,:,region]
		nat_data=nat_data_max[:,:,region]
		hist_data_3day=(hist_data[:,:-2]+hist_data[:,1:-1]+hist_data[:,2:])/3.-273.15
		nat_data_3day=(nat_data[:,:-2]+nat_data[:,1:-1]+nat_data[:,2:])/3.-273.15

		ind=obs_indices[region]
		obs_data=obs[437:439,ind[0],ind[1]].max()
		
#		bias: Correct the jja ave of the monthly max of the three daily mean Tmax. 
		obs_jj=(obs[12*11+5:12*31+5:12,ind[0],ind[1]]+obs[12*11+6:12*31+6:12,ind[0],ind[1]]+obs[12*11+7:12*31+7:12,ind[0],ind[1]])/3.
		hist_jjmean=(hist_data_3day[:,0:30].max(1)+hist_data_3day[:,30:60].max(1) +hist_data_3day[:,60:].max(1))/3.
		hist_data_bias=hist_jjmean.mean(0)-obs_jj.mean()
		print hist_data_bias
		
		
		hist_data_3day=hist_data_3day.max(1)-hist_data_bias
		nat_data_3day=nat_data_3day.max(1)-hist_data_bias
		
		# Seasonmax
		return_time_temp(hist_data_3day,nat_data_3day,days_aved,average,region_names[region+2])
		# All data
#		return_time_temp(hist_data_xday,nat_data_xday,days_aved,average,region_names[region])

	# 5 day ave
#		print 'Five day'
#		days_aved=5
#		hist_data_xday=(hist_data[:,:-4]+hist_data[:,1:-3]+hist_data[:,2:-2]+hist_data[:,3:-1]+hist_data[:,4:])/5. -273.15
#		nat_data_xday=(nat_data[:,:-4]+nat_data[:,1:-3]+nat_data[:,2:-2]+nat_data[:,3:-1]+nat_data[:,4:])/5. - 273.15
	
#		print 'hist_xday',hist_data_xday.min(),hist_data_xday.mean(),hist_data_xday.max()
#		print 'nat_xday',nat_data_xday.min(),nat_data_xday.mean(),nat_data_xday.max()
		# Seasonmax
#		return_time_temp(hist_data_xday.max(1),nat_data_xday.max(1),days_aved,average,region_names[region])
		# All data
#		return_time_temp(hist_data_xday,nat_data_xday,days_aved,average,region_names[region])
	
	
	print 'Finished!'

#Washerboard function that allows main() to run on running this file
if __name__=="__main__":
  main()




