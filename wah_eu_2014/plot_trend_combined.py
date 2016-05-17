import numpy as np
from netcdf_file import netcdf_file
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle
from mpl_toolkits.basemap import Basemap,cm	
from region_returntimes import get_region_returntimes, remap_p5deg_to_rotated

if __name__=='__main__':
	
	obsname='BEST'
	
#######################################
# Geert Jan data
	# Open file:
#	fname='../../scratch/data_from_ouce/trends_berkeley_tavg_2986.nc'
#	fname='../../scratch/data_from_ouce/g20150603_0937_20151.nc'

	if obsname=='BEST':
		fname='../../scratch/data_from_ouce/geertjan_trend_berkeley.nc'
		fname='../../scratch/data_from_ouce/g20160512_1836_2290_nc3.nc'
	else:
		fname='../../scratch/data_from_ouce/geertjan_trend_crutem.nc'
		
	f_gy=netcdf_file(fname,'r')
# Regression
	regr=f_gy.variables['regr'][0,:] # (convert to absolute temp by multiplying by global warming)
#	var=np.ma.masked_values(ratio,1.e20)
	regr=np.ma.masked_values(regr,3.e33)*0.95
	print regr.min(),regr.max()


	lat=f_gy.variables['lat'][:]
	lon=f_gy.variables['lon'][:]
	ny=lat.shape[0]
	nx=lon.shape[0]
	lon2,lat2=np.meshgrid(lon,lat)

#######################################
# Andrew King data

	fname='../../scratch/data_from_ouce/CMIP5_signal.nc'
	fname='../../scratch/data_from_ouce/CMIP5_signal_v2.nc'
#	fname='../../scratch/data_from_ouce/FAR_bestestimate.nc'
	f_ak=netcdf_file(fname,'r')
	
# FAR
	trend_ak=f_ak.variables['signal'][:]
#	var=np.ma.masked_values(ratio,1.e20)
#	trend_ak=np.ma.masked_values(trend_ak,-99.)
#	trend_ak=np.ma.masked_values(trend_ak,0.0)
#	trend_ak=np.ma.masked_where(~np.isfinite(trend_ak),trend_ak)
#	far_ak[~np.isfinite(far_ak)]=1.0
	
	print np.nanmin(trend_ak),np.nanmax(trend_ak)
	print 'FAR_CMIP5 mean:',trend_ak.mean()
	
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
	
#	rot_template= '/home/cenv0437/scratch/data_from_ouce/hadam3p_eu_z4ao_2013_1_009238311_0_tasmean.nc'
	rot_template='/home/cenv0437/scratch/data_from_ouce/hadam3p_eu_template.nc'
	f_rot=netcdf_file(rot_template,'r')
	lat_coord=f_rot.variables['global_latitude0'][4:-7,4:-4]
	lon_coord=f_rot.variables['global_longitude0'][4:-7,4:-4]
	f_rot.close()	
	
	hist_mean=historical_data.mean(0)
	nat_mean=natural_data.mean(0)
	histclim_mean=clim_hist_data.mean(0)
	natclim_mean=clim_nat_data.mean(0)	

	# FAR
	trend2014=(hist_mean-nat_mean) #/(hist_mean.mean()-nat_mean.mean()) # Areas where count_hist =0 -> NaN
	print 'FAR2014 mean:',trend2014.mean()
	trendclim=(histclim_mean-natclim_mean) #/(histclim_mean.mean()-natclim_mean.mean()) # Areas where count_hist =0 -> NaN
	print 'FARClim mean:',trendclim.mean()

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
#	plt.set_cmap('coolwarm')
	plt.set_cmap('OrRd')
	
#	m = Basemap(projection='stere',lon_0=lon2[ny/2,nx/2]+3.0,lat_0=lat2[ny/2,nx/2],resolution='c',llcrnrlon=lon2[0,0],llcrnrlat=lat2[0,0]-5.0,urcrnrlon=lon2[-1,-1]+10.0,urcrnrlat=lat2[-1,-1]+5.0)

################################
# Create basmap mappings for different grids

 	x,y=m(lon2,lat2) # Geert Jan grid
	x2,y2=m(lon_ak,lat_ak) # CMIP5 grid
	x3,y3=m(lon_coord[:],lat_coord[:])  # WAH grid

	print "trend:",(regr).mean(),trend_ak.mean(),trend2014.mean(),trendclim.mean()

###############################
	
	circles=[30,40,50,60,70]
	meridians=[-40,-30,-20,-10,0,10,20,30,40,50,60]
	axes_list=[]

	plt.figure(6)
	folder='combined_figs/'
	vmin=0
	vmax=2

	plt.subplot(221)
	axes_list.append(plt.gca())
	plt.title('a)')
	c=m.pcolor(x,y,regr,vmin=vmin,vmax=vmax)
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	

	plt.subplot(222)
	axes_list.append(plt.gca())
	plt.title('b)')
	c=m.pcolor(x2,y2,trend_ak,vmin=vmin,vmax=vmax)
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)

	
	plt.subplot(223)
	axes_list.append(plt.gca())
	plt.title('c)')
	c=m.pcolor(x3,y3,trend2014,vmin=vmin,vmax=vmax)
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	
	plt.subplot(224)
	axes_list.append(plt.gca())
	plt.title('d)')
	c=m.pcolor(x3,y3,trendclim,vmin=vmin,vmax=vmax)
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	
	
	# Finish up plot
	plt.tight_layout()
	plt.colorbar(c,orientation='vertical',ax=axes_list,aspect=40,pad=0.05,extend='both')

	plt.savefig(folder+'trend_combined_'+obsname+'_v4.png')
		
