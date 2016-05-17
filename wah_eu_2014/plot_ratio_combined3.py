import numpy as np
from netcdf_file import netcdf_file
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle
from mpl_toolkits.basemap import Basemap,cm	
from region_returntimes import get_region_returntimes, remap_p5deg_to_rotated,remap_1deg_to_rotated,ave_bias,find_closest_1d
from matplotlib.colors import LogNorm

if __name__=='__main__':
	
#######################################
# Geert Jan data
	# Open file:
#	fname='../../scratch/data_from_ouce/trends_berkeley_tavg_2986.nc'
#	fname='../../scratch/data_from_ouce/trends_berkeley_tavg_14560.nc'
	fname='../../scratch/data_from_ouce/trends_berkeley_tavg_1_mean_80.oldenborgh@knmi.nl_20685_nc3.nc'
	f_gy=netcdf_file(fname,'r')
	
# Ratio
	ratio=f_gy.variables['ratio'][0,:]
#	var=np.ma.masked_values(ratio,1.e20)
	ratio=np.ma.masked_values(ratio,3.e33)
	print ratio.min(),ratio.max()
# LORatio
	loratio=f_gy.variables['loratio'][0,:]
#	loratio=np.ma.masked_values(loratio,1.e20)
	loratio=np.ma.masked_values(loratio,3.e33)
	print loratio.min(),loratio.max()
	

	lat=f_gy.variables['lat'][:]
	lon=f_gy.variables['lon'][:]
	ny=lat.shape[0]
	nx=lon.shape[0]
	lon2,lat2=np.meshgrid(lon,lat)


#######################################
# Andrew King data

#	fname='../../scratch/data_from_ouce/FAR_bestestimate_2014thresh.nc'
	fname='../../scratch/data_from_ouce/FAR10p_2014thresh_v2.nc'
#	fname='../../scratch/data_from_ouce/FAR_bestestimate.nc'
	f_ak=netcdf_file(fname,'r')
	
# FAR
	far_ak=f_ak.variables['FAR'][:]
#	var=np.ma.masked_values(ratio,1.e20)
	far_ak=np.ma.masked_values(far_ak,-999.)
	far_ak=np.ma.masked_values(far_ak,0.0)
#	far_ak[far_ak==-999.]=0.9
	#far_ak=np.ma.masked_where(~np.isfinite(far_ak),far_ak)
	ratio_ak= 1./(1.-far_ak+1.e-16)
#	far_ak[~np.isfinite(far_ak)]=1.0

# Can't remember why I would have done this...	
#	lat=np.append(f_ak.variables['latitude'][:],[90])
#	lon=np.append(f_ak.variables['longitude'][:],[360])

	lat=f_ak.variables['latitude'][:]
	lon=f_ak.variables['longitude'][:]
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
	
	# Get old bias

	bias_files='/home/cenv0437/scratch/data_from_ouce/EU_???_temp-cruts_0.44_mean.nc'
	bias = ave_bias(bias_files)[:,:]

	# Model climatology
	clim_file=netcdf_file('/home/cenv0437/tas_clim_1960-1990.nc','r')
	clim=clim_file.variables['tas'][0,0,:]
	lat_coord=clim_file.variables['global_latitude0'][:]
	lon_coord=clim_file.variables['global_longitude0'][:]
	clim_file.close()

	# OBS anom

	obs_file=netcdf_file('TAVG_LatLong1_1961_1991_anom_nc3.nc','r')
	obs=obs_file.variables['temperature'][-12:,:,:].mean(0)
	obs_fillvalue=obs_file.variables['temperature']._FillValue
	obslat=obs_file.variables['latitude'][:]
	obslon=obs_file.variables['longitude'][:]
	# Regrid obs
	obs=remap_1deg_to_rotated(obs,'/home/cenv0437/tas_clim_1960-1990.nc',obs_fillvalue)
	
	# Use bias field to mask away ocean
	clim=np.ma.masked_where(bias.mask,clim)
	
	condition=historical_data-clim>obs
	count_hist=condition.sum(0)/(historical_data.shape[0]*1.0)
	
	condition=natural_data-clim>obs
	count_nat=condition.sum(0)/(natural_data.shape[0]*1.0)

	condition=clim_hist_data-clim>obs
	count_clim_hist=condition.sum(0)/(clim_hist_data.shape[0]*1.0)
	
	condition=clim_nat_data-clim>obs
	count_clim_nat=condition.sum(0)/(clim_nat_data.shape[0]*1.0)	

	# ratio
	ratio2014=((count_hist)/(count_nat+1.e-16)) # Areas where count_hist =0 -> NaN
	ratioclim=((count_clim_hist)/(count_clim_nat+1.e-16)) # Areas where count_hist =0 -> NaN

#####################################

# Print out specific regions:
        regions=[(51.2, 13.8),(62.561671, 13.958159),(51.749181, -1.281326),(55.685268, 37.274331),(43.959139, 20.910196),(40.402960, -3.732137)]
        rnames=['Germany','Sweden','England','Russia','Serbia','Spain']
	
	for i,reg in enumerate(regions):
		print rnames[i]
		jp,ip=find_closest_1d(reg[0],reg[1],lat,lon)
		print ratio_ak[jp,ip]


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

	x,y=m(lon2,lat2)
	x2,y2=m(lon_ak,lat_ak)
	x3,y3=m(lon_coord[:],lat_coord[:])


###############################
	
		
	circles=[30,40,50,60,70]
	meridians=[-40,-30,-20,-10,0,10,20,30,40,50,60]
	axes_list=[]

	plt.figure(6)
	folder='combined_figs/'

	plt.subplot(221)
	axes_list.append(plt.gca())
	plt.title('a)')
	c=m.pcolor(x,y,loratio,norm=LogNorm(vmin=1., vmax=1.e3))
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	

	plt.subplot(222)
	axes_list.append(plt.gca())
	plt.title('b)')
	c=m.pcolor(x2,y2,ratio_ak,norm=LogNorm(vmin=1., vmax=1.e3))#,np.arange(0.6,1.05,.05),extend='both')
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)

	
	plt.subplot(223)
	axes_list.append(plt.gca())
	plt.title('c)')
	c=m.pcolor(x3,y3,ratio2014,norm=LogNorm(vmin=1., vmax=1.e3))
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	
	plt.subplot(224)
	axes_list.append(plt.gca())
	plt.title('d)')
	c=m.pcolor(x3,y3,ratioclim,norm=LogNorm(vmin=1., vmax=1.e3))
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	
	
	# Finish up plot
	plt.tight_layout()
#	c.cmap.set_under('b')
	cbar=plt.colorbar(c,orientation='vertical',ax=axes_list,aspect=40,pad=0.05,extend='max')
	cbar.ax.set_yticklabels(['1', '10', '100','1000'])

	plt.savefig(folder+'ratio_combined3.png')
		
