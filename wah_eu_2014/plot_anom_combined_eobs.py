import numpy as np
from netcdf_file import netcdf_file
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle
from mpl_toolkits.basemap import Basemap,cm	
from region_returntimes import get_region_returntimes, remap_p5deg_to_rotated2
from datetime import datetime

if __name__=='__main__':
	
#######################################
# Geert Jan data
	# Open file:
#	fname='../../scratch/data_from_ouce/trends_berkeley_tavg_2986.nc'
	fname='../../scratch/data_from_ouce/s21249_1900_2013_anom_2014_ann_nc3.nc'
	f_gy=netcdf_file(fname,'r')
	
# Anom
	anom_gy=f_gy.variables['tg'][0,0,:]
#	var=np.ma.masked_values(ratio,1.e20)
	anom_gy=np.ma.masked_values(anom_gy,3.e33)
	print anom_gy.min(),anom_gy.max()


	lat=f_gy.variables['lat'][:]
	lon=f_gy.variables['lon'][:]
	ny=lat.shape[0]
	nx=lon.shape[0]
	lon2,lat2=np.meshgrid(lon,lat)


#######################################
# Andrew King data

	fname='../../scratch/data_from_ouce/CMIP5_Obs2014_DifferenceV2.nc'
#	fname='../../scratch/data_from_ouce/FAR_bestestimate.nc'
	f_ak=netcdf_file(fname,'r')
	
# FAR
	anom_ak=f_ak.variables['Obs-Models'][:]
#	anom_ak=np.ma.masked_values(anom_ak,1.e20)
	anom_ak=np.ma.masked_values(anom_ak,-99.)
#	trend_ak=np.ma.masked_values(trend_ak,0.0)
#	trend_ak=np.ma.masked_where(~np.isfinite(trend_ak),trend_ak)
#	far_ak[~np.isfinite(far_ak)]=1.0
	
	print np.nanmin(anom_ak),np.nanmax(anom_ak)
	print 'FAR_CMIP5 mean:',anom_ak.mean()
	
	lat=np.append(f_ak.variables['latitude'][:],[90])
	lon=np.append(f_ak.variables['longitude'][:],[360])
	ny=lat.shape[0]
	nx=lon.shape[0]
	lon_ak,lat_ak=np.meshgrid(lon,lat)


#######################################
# WAH data

	pkl_dir='/gpfs/projects/cpdn/scratch/cenv0437/pickle_data/'
	fdata=open(pkl_dir+'historical_data2.pkl','rb')
	historical_data=pickle.load(fdata)
	fdata.close()
	fdata=open(pkl_dir+'natural_data2.pkl','rb')
	natural_data=pickle.load(fdata)
	fdata.close()
	fdata=open(pkl_dir+'clim_hist_data2.pkl','rb')
	clim_hist_data=pickle.load(fdata)
	fdata.close()
	fdata=open(pkl_dir+'clim_nat_data2.pkl','rb')
	clim_nat_data=pickle.load(fdata)
	fdata.close()
	print 'loaded data from pkl files'
	print 'hist2014',historical_data.shape[0]
	print 'nat2014',natural_data.shape[0]
	print 'histClim',clim_hist_data.shape[0]
	print 'natClim',clim_nat_data.shape[0]
	

	
	# Model climatology relative to 1961-1990
	clim_file=netcdf_file('/home/cenv0437/tas_clim_1960-1990.nc','r')
	clim=clim_file.variables['tas'][0,0,:]
	lat_coord=clim_file.variables['global_latitude0'][:]
	lon_coord=clim_file.variables['global_longitude0'][:]
	clim_file.close()
	
	
	obs_file=netcdf_file('../../scratch/data_from_ouce/tg_0.50deg_reg_v11.0.nc','r')
	obs_var=obs_file.variables['tg']
	obs_fillvalue=obs_var._FillValue
	obs_scale_factor=obs_var.scale_factor
	# Calculate anomaly of observations relative to reference period 1961-1990
	ref_start=(datetime(1961,1,1)-datetime(1950,1,1)).days
	ref_end=(datetime(1990,12,31)-datetime(1950,1,1)).days
	anom_start=(datetime(2014,1,1)-datetime(1950,1,1)).days
	anom_end=(datetime(2014,12,31)-datetime(1950,1,1)).days
	print 'ref',ref_start,ref_end,'anom',anom_start,anom_end
	ref=(np.ma.masked_values(obs_var[ref_start:ref_end,:,:],obs_fillvalue)*obs_scale_factor).mean(0)
	anom=(np.ma.masked_values(obs_var[anom_start:anom_end,:,:],obs_fillvalue)*obs_scale_factor).mean(0)
	obs=anom-ref
	# Regrid obs
	rot_template= '/home/cenv0437/tas_clim_1960-1990.nc'
	obs=remap_p5deg_to_rotated2(obs,rot_template,obs_fillvalue)
	
	hist_mean=historical_data.mean(0)
	nat_mean=natural_data.mean(0)
	histclim_mean=clim_hist_data.mean(0)
	natclim_mean=clim_nat_data.mean(0)	
	anom2014_wah=obs-(hist_mean-clim)
	anomClim_wah=obs-(histclim_mean-clim)
	# Mask
	anom2014_wah=np.ma.masked_where(obs.mask,anom2014_wah)
	anomClim_wah=np.ma.masked_where(obs.mask,anomClim_wah)

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
	vmin=-2
	vmax=2

	plt.subplot(221)
	axes_list.append(plt.gca())
	plt.title('a)')
	c=m.pcolor(x,y,anom_gy,vmin=vmin,vmax=vmax)
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	

	plt.subplot(222)
	axes_list.append(plt.gca())
	plt.title('b)')
	c=m.pcolor(x2,y2,anom_ak,vmin=vmin,vmax=vmax)
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)

	
	plt.subplot(223)
	axes_list.append(plt.gca())
	plt.title('c)')
	c=m.pcolor(x3,y3,anom2014_wah,vmin=vmin,vmax=vmax)
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	
	plt.subplot(224)
	axes_list.append(plt.gca())
	plt.title('d)')
	c=m.pcolor(x3,y3,anomClim_wah,vmin=vmin,vmax=vmax)
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	
	
	# Finish up plot
	plt.tight_layout()
	plt.colorbar(c,orientation='vertical',ax=axes_list,aspect=40,pad=0.05,extend='min')

	plt.savefig(folder+'anom_combined_eobs.png')
	
	# Check obs anomaly
	plt.figure()
	plt.title('EOBS anomaly for 2014 vs 1961-1990')
	c=m.pcolor(x3,y3,obs,vmin=0,vmax=3)
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	plt.colorbar()
	plt.savefig('EOBS_check.png')
