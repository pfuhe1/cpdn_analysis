import numpy as np
from netcdf_file import netcdf_file
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle
from mpl_toolkits.basemap import Basemap,cm	

if __name__=='__main__':


	# Open file:
#	fname='../../scratch/data_from_ouce/trends_berkeley_tavg_2986.nc'
	fname='../../scratch/data_from_ouce/trends_berkeley_tavg_14560.nc'
	f_gy=netcdf_file(fname,'r')
	
# Ratio
	ratio=f_gy.variables['ratio'][0,:]
#	var=np.ma.masked_values(ratio,1.e20)
	ratio=np.ma.masked_values(ratio,3.e33)
	print ratio.min(),ratio.max()

# LORatio
	loratio=f_gy.variables['loratio'][0,:]
#	var=np.ma.masked_values(loratio,1.e20)
	loratio=np.ma.masked_values(loratio,3.e33)
	print loratio.min(),loratio.max()
	
	lat=f_gy.variables['lat'][:]
	lon=f_gy.variables['lon'][:]
	ny=lat.shape[0]
	nx=lon.shape[0]
	lon2,lat2=np.meshgrid(lon,lat)

	
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
	x,y=m(lon2[:],lat2[:])
	
	plt.figure(6)
	folder='tmp_figs/'

#	plt.title('FAR = '+'{:.2f}'.format(far.mean()))
#	c=m.contourf(x,y,np.log(var))#,np.arange(0.0,1.05,.05),extend='both')
	c=m.contourf(x,y,1-1/ratio,np.arange(0.4,1.05,.05),extend='both')
	m.drawcoastlines()
	m.drawcountries()
#	m.drawparallels(circles)
#	m.drawmeridians(meridians)
	# Make an axis for the colorbar on the right side
	plt.colorbar(c,orientation='horizontal',aspect=40,pad=0.05)
	plt.savefig(folder+'geert_jan_far.png')	
	
	plt.figure(2)
#	plt.title('FAR = '+'{:.2f}'.format(far.mean()))
#	c=m.contourf(x,y,np.log(var))#,np.arange(0.0,1.05,.05),extend='both')
	c=m.contourf(x,y,1-1/loratio,np.arange(0.4,1.05,.05),extend='both')
	m.drawcoastlines()
	m.drawcountries()
#	m.drawparallels(circles)
#	m.drawmeridians(meridians)
	# Make an axis for the colorbar on the right side
	plt.colorbar(c,orientation='horizontal',aspect=40,pad=0.05)
	plt.savefig(folder+'geert_jan_farlo.png')	
		
		
