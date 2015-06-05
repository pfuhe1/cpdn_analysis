import numpy as np
from netcdf_file import netcdf_file
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle
from mpl_toolkits.basemap import Basemap,cm	
from region_returntimes import get_region_returntimes, remap_p5deg_to_rotated

if __name__=='__main__':
	
#######################################
# Geert Jan data
	# Open file:
	fname='../../scratch/data_from_ouce/m23134.nc'
	f_gy=netcdf_file(fname,'r')
	
# Standard deviation
	stdev=f_gy.variables['sd'][0,0,:]
	stdev=np.ma.masked_values(stdev,3.e33)
	print stdev.min(),stdev.max()

	lat=f_gy.variables['lat'][:]
	lon=f_gy.variables['lon'][:]
	ny=lat.shape[0]
	nx=lon.shape[0]
	lon2,lat2=np.meshgrid(lon,lat)


#######################################
# Andrew King data

#	fname='../../scratch/data_from_ouce/FAR_bestestimate_2014thresh.nc'
	fname='../../scratch/data_from_ouce/CMIP5_noise_detrend.nc'
#	fname='../../scratch/data_from_ouce/FAR_bestestimate.nc'
	f_ak=netcdf_file(fname,'r')
	
# FAR
	noise_ak=f_ak.variables['noise'][:]
#	var=np.ma.masked_values(ratio,1.e20)
#	noise_ak=np.ma.masked_values(noise_ak,-999.)
#	noise_ak=np.ma.masked_values(noise_ak,0.0)
#	noise_ak=np.ma.masked_where(~np.isfinite(noise_ak),noise_ak)
#	far_ak[~np.isfinite(far_ak)]=1.0
	
	print np.nanmin(noise_ak),np.nanmax(noise_ak)
	print 'FAR_CMIP5 mean:',noise_ak.mean()
	
	lat=np.append(f_ak.variables['latitude'][:],[90])
	lon=np.append(f_ak.variables['longitude'][:],[360])
	ny=lat.shape[0]
	nx=lon.shape[0]
	lon_ak,lat_ak=np.meshgrid(lon,lat)


#######################################
# WAH data

	pkl_dir='/gpfs/projects/cpdn/scratch/cenv0437/pickle_data/'
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
	print 'hist2014',historical_data.shape[0]
	print 'nat2014',natural_data.shape[0]
	print 'histClim',clim_hist_data.shape[0]
	print 'natClim',clim_nat_data.shape[0]
	
	obs_p5deg= '/home/cenv0437/scratch/data_from_ouce/CRU_TS_Absolute_plus_BEST_anomaly_201312_201412.nc'
	obsname='BEST'
	rot_template= '/home/cenv0437/scratch/data_from_ouce/hadam3p_eu_z4ao_2013_1_009238311_0_tasmean.nc'
	f_rot=netcdf_file(rot_template,'r')
	lat_coord=f_rot.variables['global_latitude0'][4:-7,4:-4]
	lon_coord=f_rot.variables['global_longitude0'][4:-7,4:-4]
	f_rot.close()	
	obs=remap_p5deg_to_rotated(obs_p5deg,rot_template)[:,:]
	


	# FAR
	noise2014=natural_data.std(0) #/(hist_mean.mean()-nat_mean.mean()) # Areas where count_hist =0 -> NaN
	print ' Noise 2014 mean:',noise2014.mean()
	noiseclim=clim_nat_data.std(0) #/(histclim_mean.mean()-natclim_mean.mean()) # Areas where count_hist =0 -> NaN
	print 'Noise Clim mean:',noiseclim.mean()

#######################################
# Plotting stuff
	
	#Copy projection region from previous set up
	lon0= 17.3525972366
	lat0= 50.1672
	llcrnrlon= 348.619
	llcrnrlat= 23.1428909302
	urcrnrlon= 63.4356536865
	urcrnrlat= 66.9738
	m = Basemap(projection='stere',lon_0=lon0,lat_0=lat0, llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat)
	plt.set_cmap('coolwarm')
	
#	m = Basemap(projection='stere',lon_0=lon2[ny/2,nx/2]+3.0,lat_0=lat2[ny/2,nx/2],resolution='c',llcrnrlon=lon2[0,0],llcrnrlat=lat2[0,0]-5.0,urcrnrlon=lon2[-1,-1]+10.0,urcrnrlat=lat2[-1,-1]+5.0)

################################
# Create basmap mappings for different grids

 	x,y=m(lon2,lat2) # Geert Jan grid
	x2,y2=m(lon_ak,lat_ak) # CMIP5 grid
	x3,y3=m(lon_coord[:],lat_coord[:])  # WAH grid


###############################
	
	plt.set_cmap('coolwarm')		
	circles=[30,40,50,60,70]
	meridians=[-40,-30,-20,-10,0,10,20,30,40,50,60]
	axes_list=[]

	plt.figure(6)
	folder='combined_figs/'

	plt.subplot(221)
	axes_list.append(plt.gca())
	plt.title('a)')
	c=m.pcolor(x,y,stdev,vmin=0.6,vmax=1.)
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	

	plt.subplot(222)
	axes_list.append(plt.gca())
	plt.title('b)')
	c=m.pcolor(x2,y2,noise_ak,vmin=0.,vmax=2)#,np.arange(0.6,1.05,.05),extend='both')
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)

	
	plt.subplot(223)
	axes_list.append(plt.gca())
	plt.title('c)')
	c=m.pcolor(x3,y3,noise2014,vmin=0.,vmax=2)
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	
	plt.subplot(224)
	axes_list.append(plt.gca())
	plt.title('d)')
	c=m.pcolor(x3,y3,noiseclim,vmin=0.,vmax=2)
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	
	
	# Finish up plot
	plt.tight_layout()
	plt.colorbar(c,orientation='vertical',ax=axes_list,aspect=40,pad=0.05,extend='min')

	plt.savefig(folder+'noise_combined.png')
		
