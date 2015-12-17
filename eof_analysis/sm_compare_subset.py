# Eof analysis of soil moisture from dumps
#soil levels
#10,30,60,200  (thickness, cm)
#0,10,40,100,300 (boundary cm)

import sys,os
import numpy as np
import glob
from scipy.io import netcdf_file
#eofs installed to /gpfs/home/cenv0437/.local/lib/python2.7/site-packages/
sys.path.append('/gpfs/home/cenv0437/.local/lib/python2.7/site-packages/')
from eofs.standard import Eof
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap,cm

# Load ensemble of files into a single array
def load_regional(filenames):
	data=np.ma.zeros([len(filenames),119-12,122-12])
	i=0
	for f in filenames:
		try:
			print i
			tmp=netcdf_file(f,'r').variables['sm'][:]
			if not np.all(np.isfinite(tmp)):
				print 'error: wierd vals',f
				continue
			else:
				# Integral over top 1m of soil
				data[i,:]=tmp[0,0,6:-6,6:-6]*.1+tmp[0,1,6:-6,6:-6]*.3+tmp[0,2,6:-6,6:-6]*.6 
				i=i+1
		except Exception,e:
			print 'Error, cannot load files',f,e
			continue
	return data[:i,:]
	
	
# Load ensemble of files into a single array
def load_global(filenames):
	data=np.ma.zeros([len(filenames),145,192])
	i=0
	for f in filenames:
		try:
			tmp=netcdf_file(f,'r').variables['sm']
			if not np.all(np.isfinite(tmp)):
				print 'error: wierd vals',f
				continue
			else:
				# Integral over top 1m of soil
				data[i,:]=tmp[0,0,6:-6,6:-6]*.1+tmp[0,0,6:-6,6:-6]*.3+tmp[0,0,6:-6,6:-6]*.6 
				i=i+1
		except:
			print 'Error, cannot load files',f
			continue
	return data[:i,:]
	
def plot_region_pnw(data,lat_coord,lon_coord,lat_nw,lat_se,lon_nw,lon_se,title,varfrac):
	plt.clf()
	circles=[30,40,50,60,70]
	meridians=[-130,-120,-110]
	ny,nx=lat_coord.shape
	m = Basemap(projection='stere',lon_0=lon_coord[ny/2,nx/2],lat_0=lat_coord[ny/2,nx/2],resolution='l',llcrnrlon=lon_coord[-1,0],llcrnrlat=lat_coord[-1,0],urcrnrlon=lon_coord[0,-1],urcrnrlat=lat_coord[0,-1])
	x,y=m(lon_coord[:],lat_coord[:])
	
	plt.title('Plotting data for: '+title+'\nvariance:'+str(varfrac))
	try:
		max=np.absolute(data[lat_nw:lat_se,lon_nw:lon_se]).max()
		m.pcolor(x[lat_nw:lat_se,lon_nw:lon_se],y[lat_nw:lat_se,lon_nw:lon_se],data[lat_nw:lat_se,lon_nw:lon_se],vmin=-max,vmax=max)
		plt.colorbar()
	except:
		raise
		m.plot(x[lat_nw:lat_se,lon_nw:lon_se],y[lat_nw:lat_se,lon_nw:lon_se],'r.') # Just put a dot at the x y point
	m.drawcoastlines()
	m.drawcountries()
	y=m.drawstates()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	try:
		os.remove(output_dir+'/._'+title+'.png')
	except:
		pass
	plt.savefig(output_dir+'/'+title+'.png')
	
if __name__=='__main__':
	
		# Paths
		umids_chosen='/home/cenv0437/cpdn_analysis/eof_analysis/eof_plots_nat/umids_chosen.txt'	
		output_dir='eof_plots_nat'
		in_dir='/gpfs/projects/cpdn/scratch/cenv0437/dump_sm/hadam3p_eu/batch221/nat/'

		if not os.path.exists(output_dir):
			os.mkdir(output_dir)
		filenames=glob.glob(in_dir+'region_*.nc')
#		for batchid in range(30,42):
#			filenames.extend(glob.glob('/gpfs/projects/cpdn/scratch/cenv0437/dump_sm/hadam3p_pnw/batch181/sub-batch'+str(batchid)+'/region_*.nc'))
		fexample=netcdf_file('/home/cenv0437/scratch/batch_100/hadam3p_eu_z2ht_2013_1_009208650_0_tasmean.nc','r')
		lat_coord=fexample.variables['global_latitude0'][6:-6,6:-6]
		lon_coord=fexample.variables['global_longitude0'][6:-6,6:-6]

		print "\nWhole ensemble..."	

#		print "\n Subset of 200..."
#		np.random.shuffle(filenames)
#		filenames=filenames[:200]


		data=load_regional(filenames)
		data=np.ma.masked_values(data,2.e+20)
		print "data loaded",data.shape
	
		# Set up info
		plt.set_cmap('RdBu')
		neofs=5
		nens=data.shape[0]
		nwanted=57

		datamean=data[:].mean(0)
		sel_points=[[40,30],[59,60],[80,30],[30,75],[90,55]]
		for x,y in sel_points:
			datamean[y,x]=-200.
		plot_region_pnw(datamean,lat_coord,lon_coord,0,-1,0,-1,'mean',data.mean())

		print "\nSubsampled ensemble..."
		fnames_subset=[]
		for umid in open(umids_chosen):
			fnames_subset.append(in_dir+'region_'+umid.strip()+'_sm.nc')
		data_subset=load_regional(fnames_subset)
		data_subset=np.ma.masked_values(data_subset,2.e20)

		for i,pt in enumerate(sel_points):
			plt.clf()
			#plt.subplot(121)
			plt.boxplot([data[:,pt[1],pt[0]],data_subset[:,pt[1],pt[0]]],whis=[5,95])
			#plt.subplot(122)
			#plt.boxplot(data_subset[:,pt[1],pt[0]])
			plt.savefig(output_dir+'/boxplots_'+str(i)+'.png')
		

