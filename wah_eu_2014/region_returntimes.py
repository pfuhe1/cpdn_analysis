#!/usr/local/python2.7
# NOTE python2.7 must be used for this script

import glob,os
from sort_umids import sort_tasknames,choose_mask,choose_tasknames
import numpy as np
from netcdf_file import netcdf_file
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle
from mpl_toolkits.basemap import Basemap,cm

# Takes 4 seasonal bias files and calculates a yearly average
def ave_bias(bias_files):
	data=np.ma.zeros([108,114])
	i=0
	for f in glob.glob(bias_files):
		i=i+1
		tmp=netcdf_file(f).variables['temp_bias']
		tmp=np.ma.masked_values(tmp[0,0,:,:],tmp.missing_value)
		data=data+tmp
	if i==4:
		return data/4.0
	else:
		raise Exception('Incorrect number of seasonal bias files')

# Load ensemble of files into a single array
def load_ensemble(filenames):
	data=np.ma.zeros([len(filenames),108,114])
	for i,f in enumerate(filenames):
		tmp=netcdf_file(f,'r').variables['tas'][0,0,4:-7,4:-4] # Cut off border points
		if tmp.max()>350.0 or tmp.min<170 or not np.all(np.isfinite(tmp)):
			print 'error: wierd vals',f
			data[i,:]=np.ma.masked # Mask out this whole ensemble member
		else:
			data[i,:]=tmp
	return data
	
def remap_N96_to_rotated(N96_file,rot_template):
	f_rot=netcdf_file(rot_template,'r')
	g_lat=f_rot.variables['global_latitude0'][4:-7,4:-4]
	g_lon=f_rot.variables['global_longitude0'][4:-7,4:-4]
	N96_data=netcdf_file(N96_file,'r').variables['tasmean'][0,0,:]
	print N96_data.shape
	dataout=np.zeros([108,114])
	for i in range(114):
		for j in range(108):
			# Calculate nearest neighbor indices 
			index_lat=int(((g_lat[j,i]+90)/1.25)+0.5)
			index_lon=int((g_lon[j,i]/1.875)+0.5)%192
			dataout[j,i]=N96_data[index_lat,index_lon]
	dataout=np.ma.masked_values(dataout,0.0)
	return dataout		
	
def remap_p5deg_to_rotated(p5deg_file,rot_template):
	f_rot=netcdf_file(rot_template,'r')
	g_lat=f_rot.variables['global_latitude0'][4:-7,4:-4]
	g_lon=f_rot.variables['global_longitude0'][4:-7,4:-4]
	p5deg_nc=netcdf_file(p5deg_file,'r').variables['tasmean']
	p5deg_data=p5deg_nc[0,:]
	print p5deg_data.shape
	dataout=np.zeros([108,114])
	for i in range(114):
		for j in range(108):
			# Calculate nearest neighbor indices 
			index_lat=int(((g_lat[j,i]+89.75)/0.5)+0.5) 
			index_lon=int((g_lon[j,i]/0.5)+0.5)%720
			dataout[j,i]=p5deg_data[index_lat,index_lon]
	dataout=np.ma.masked_values(dataout,p5deg_nc._FillValue)
	return dataout


# Create netcdf file on rotated grid, removing the 'sponge layer' at the edge of the domain
def create_netcdf_stripped(template,data,outname,lat_start=4,lat_end=7,lon_start=4,lon_end=4):
	# create outfile object
	outfile=netcdf_file(outname,'w')
	
	# Create dimensions copied from template file
	temp=template.variables['tas']
	for dim in temp.dimensions:
		if dim=='time1': 
			leng=1
			outfile.createDimension(dim,leng)
			outfile.createVariable(dim,'f',(dim,))
			outfile.variables[dim][:]=template.variables[dim][:]
			outfile.variables[dim].__setattr__('axis','T')
		elif dim[:3]=='lat': 
			leng=int(template.dimensions[dim])-lat_start-lat_end
			outfile.createDimension(dim,leng)
			outfile.createVariable(dim,'f',(dim,))
			outfile.variables[dim][:]=template.variables[dim][lat_start:-lat_end]
			outfile.variables[dim].__setattr__('axis','Y')
		elif dim[:3]=='lon' : 
			leng=int(template.dimensions[dim])-lon_start-lon_end
			outfile.createDimension(dim,leng)
			outfile.createVariable(dim,'f',(dim,))
			outfile.variables[dim][:]=template.variables[dim][lon_start:-lon_end]
			outfile.variables[dim].__setattr__('axis','X')
		else: 
			leng=int(template.dimensions[dim])
			outfile.createDimension(dim,leng)
			outfile.createVariable(dim,'f',(dim,))
			outfile.variables[dim][:]=template.variables[dim][:]
			
		for att in template.variables[dim]._attributes:
			outfile.variables[dim].__setattr__(att,template.variables[dim].__getattribute__(att))
		
		if dim=='time1':	
			outfile.variables[dim].__setattr__('calendar','gregorian') # Hack so Karsten can read these files in grads
	# Rotated pole variable
	outfile.createVariable('rotated_pole0','c',())
	for att in template.variables['rotated_pole0']._attributes:
			outfile.variables['rotated_pole0'].__setattr__(att,template.variables['rotated_pole0'].__getattribute__(att))
	
	# Global latitude
	outfile.createVariable('global_latitude0','f',template.variables['global_latitude0'].dimensions)
	outfile.variables['global_latitude0'][:]=template.variables['global_latitude0'][lat_start:-lat_end,lon_start:-lon_end]
	for att in template.variables['global_latitude0']._attributes:
			outfile.variables['global_latitude0'].__setattr__(att,template.variables['global_latitude0'].__getattribute__(att))
			
	# Global longitude
	outfile.createVariable('global_longitude0','f',template.variables['global_longitude0'].dimensions)
	outfile.variables['global_longitude0'][:]=template.variables['global_longitude0'][lat_start:-lat_end,lon_start:-lon_end]
	for att in template.variables['global_longitude0']._attributes:
			outfile.variables['global_longitude0'].__setattr__(att,template.variables['global_longitude0'].__getattribute__(att))
	
	# Create data variable (named tas)
	outfile.createVariable('tas','f',temp.dimensions)
	outfile.variables['tas'].assignValue(data.filled(temp._FillValue))
	for att in temp._attributes:
		outfile.variables['tas'].__setattr__(att,temp.__getattribute__(att))
	
	outfile.close()

###############

def land_points(obs):
	#number of land points
	return np.logical_not(obs.mask).sum()

# Calculate the return time of the mean of the data being greater than an observational dataset
def ret_time(data,obs):

	ensembles=data.shape[0]*1.0 # First dimension of data is for ensembles. 
	try:
		return ensembles/(data.mean(1).mean(1)>obs.mean(0).mean(0)).sum(0)
	except:
		return np.inf

#################

def print_model_returntimes(natural_data,obs):
# CCSM4	
	model_mask=choose_mask(natural_files,'z2go','z2m7')+choose_mask(natural_files,'z4so','z5kf')
	print sum(model_mask)
	ret_model=ret_time(natural_data[model_mask,:,:],obs)
	print '2014, CCSM4:', ret_model
	
# Canesm
	model_mask=choose_mask(natural_files,'z32w','z38f')+choose_mask(natural_files,'z7vs','z8nj')
	print sum(model_mask)
	ret_model=ret_time(natural_data[model_mask,:,:],obs)
	print '2014, Canesm:', ret_model
	
# CNRM-CM5	
	model_mask=choose_mask(natural_files,'z38g','z3dz')+choose_mask(natural_files,'z8nk','z9fb')
	print sum(model_mask)
	ret_model=ret_time(natural_data[model_mask,:,:],obs)
	print '2014, CNRM-CM5:', ret_model

	
# CSIRO-Mk3-6-0
	model_mask=choose_mask(natural_files,'z3e0','z3jk')+choose_mask(natural_files,'z9fc','za73')
	print sum(model_mask)
	ret_model=ret_time(natural_data[model_mask,:,:],obs)
	print '2014 CSIRO-Mk3-6-0 return time:',ret_model
	
# GFDL-CM3
	model_mask=choose_mask(natural_files,'z2m8','z2rr')+choose_mask(natural_files,'z5kg','z6c7')
	print sum(model_mask)
	ret_model=ret_time(natural_data[model_mask,:,:],obs)
	print '2014 GFDL-CM3 return time:',ret_model

# GISS-E2-H
	model_mask=choose_mask(natural_files,'z2rs','z2xb')+choose_mask(natural_files,'z6c8','z73z')
	print sum(model_mask)
	ret_model=ret_time(natural_data[model_mask,:,:],obs)
	print '2014 GISS-E2-H return time:',ret_model
	
# GISS-E2-R
	model_mask=choose_mask(natural_files,'z3jk','z3p3')+choose_mask(natural_files,'za74','zayv')
	print sum(model_mask)
	ret_model=ret_time(natural_data[model_mask,:,:],obs)
	print '2014 GISS-E2-R return time:',ret_model
	
# HadGEM2-ES
	model_mask=choose_mask(natural_files,'z2xc','z32v')+choose_mask(natural_files,'z740','z7vr')
	print sum(model_mask)
	ret_model=ret_time(natural_data[model_mask,:,:],obs)
	print '2014 HadGEM2-ES return time:',ret_model
	
# HadGEM2-ES
	model_mask=choose_mask(natural_files,'z2xc','z32v')+choose_mask(natural_files,'z740','z7vr')
	print sum(model_mask)
	ret_model=ret_time(natural_data[model_mask,:,:],obs)
	print '2014 HadGEM2-ES return time:',ret_model
	
# IPSL-CM5A-LR
	model_mask=choose_mask(natural_files,'z3p4','z3un')+choose_mask(natural_files,'zayw','zbqn')
	print sum(model_mask)
	ret_model=ret_time(natural_data[model_mask,:,:],obs)
	print '2014 IPSL-CM5A-LR return time:',ret_model
	
# IPSL-CM5A-MR
	model_mask=choose_mask(natural_files,'z3uo','z407')+choose_mask(natural_files,'zbqo','zcif')
	print sum(model_mask)
	ret_model=ret_time(natural_data[model_mask,:,:],obs)
	print '2014 IPSL-CM5A-MR return time:',ret_model
	
# MIROC-ESM
	model_mask=choose_mask(natural_files,'z408','z45r')+choose_mask(natural_files,'zcig','zda7')
	print sum(model_mask)
	ret_model=ret_time(natural_data[model_mask,:,:],obs)
	print '2014 MIROC-ESM return time:',ret_model
	
# MMM
	model_mask=choose_mask(natural_files,'z45s','z4bb')+choose_mask(natural_files,'zda8','ze1z')
	print sum(model_mask)
	ret_model=ret_time(natural_data[model_mask,:,:],obs)
	print '2014 MMM return time:',ret_model
	
# Gets return time for slice of the grid. 
# Roughly corresponding to:
# (ymin,xmin) being the north western corner
# (ymax,xmax) being the south easter norner
# These directions are not exact on the rotated grid 
def print_region_returntimes(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,xmin,xmax,ymin,ymax):
	print_returntimes(historical_data[:,ymin:ymax,xmin:xmax],natural_data[:,ymin:ymax,xmin:xmax],clim_hist_data[:,ymin:ymax,xmin:xmax],clim_nat_data[:,ymin:ymax,xmin:xmax],obs[ymin:ymax,xmin:xmax])


def get_region_returntimes(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,xmin,xmax,ymin,ymax):
	return get_returntimes(historical_data[:,ymin:ymax,xmin:xmax],natural_data[:,ymin:ymax,xmin:xmax],clim_hist_data[:,ymin:ymax,xmin:xmax],clim_nat_data[:,ymin:ymax,xmin:xmax],obs[ymin:ymax,xmin:xmax])


def get_returntimes(historical_data,natural_data,clim_hist_data,clim_nat_data,obs):
	ret_hist=ret_time(historical_data,obs)
	ret_nat=ret_time(natural_data,obs)
	ret_clim_hist=ret_time(clim_hist_data,obs)
	ret_clim_nat=ret_time(clim_nat_data,obs)
	lp=land_points(obs)
	return ret_hist,ret_nat,ret_clim_hist,ret_clim_nat,lp
	
def print_returntimes(historical_data,natural_data,clim_hist_data,clim_nat_data,obs):

	ret_hist=ret_time(historical_data,obs)
	print '2014 All Forcings return time:',ret_hist

	ret_nat=ret_time(natural_data,obs)
	print '2014 Natural return time:',ret_nat
	
	# Print return time separated by model
	print_model_returntimes(natural_data,obs)

######## Climatologies

	ret_clim_hist=ret_time(clim_hist_data,obs)
	print 'Climatology All Forcings return time:',ret_clim_hist	

	ret_clim_nat=ret_time(clim_nat_data,obs)
	print 'Climatology Natural return time:',ret_clim_nat

##############################################################

# Searches for closest point to p_lat,p_lon and returns grid reference
def find_closest(p_lat,p_lon,lat_coord,lon_coord):
	ny,nx=lat_coord.shape
	min_dist=100000000
	min_point=(-1,-1)
	for j in range(ny):
		for i in range(nx):
			dist= (lat_coord[j,i]-p_lat)**2+((lon_coord[j,i]-p_lon)%360)**2
#			print lat_coord[j,i],p_lat,lon_coord[j,i],p_lon,dist
			if dist<min_dist:
				min_point=(j,i)
				min_dist=dist
	return min_point

#############################################################

def plot_spacial(historical_data,natural_data,obs,lat_coord,lon_coord,folder):	
	# Spatially resolved count of ensembles hotter than obs:	
	condition=historical_data>obs
	count_hist=condition.sum(0)/(historical_data.shape[0]*1.0)
	
	condition=natural_data>obs
	count_nat=condition.sum(0)/(natural_data.shape[0]*1.0)

	# Change in likelihood
	change=(count_hist-count_nat)/(count_nat+1e-5)
	
	# Anomaly between obs and ensemble mean of historical simulations
	anom=obs-(historical_data.mean(0))

	# Anomaly between obs and ensemble mean of nat simulations
	anom_nat=obs-(natural_data.mean(0))

	# FAR
        far=1.0-((count_nat)/(count_hist))
#	far=1.0-((count_nat)/(count_hist+1e-15))
#	far=far*(np.logical_or(count_hist!=0.0,count_nat!=0.0)*1.0) # Set to zero where count_hist and count_nat==0
#	far=far*((count_hist!=0.0)*1.0)
	#sanity check
#	print 'sanity check of FAR, should be 0', np.logical_and(count_hist==0.0 and count_nat!=0.0).sum()
################################################
# Plot stuff

	if not os.path.exists(folder):
		os.mkdir(folder)
	folder=folder+'/' # Make sure folder has a slash


#First get rid of annoying Mac hidden files
	delfiles=glob.glob(folder+'._*.png')
	for f in delfiles:
		os.remove(f)
	ny,nx=lat_coord.shape
		
	# Stereographic projection. lon_0 and lat_0 are the central points
#	m = Basemap(projection='stere',lon_0=lon_coord[ny/2,nx/2]+7.0,lat_0=lat_coord[ny/2,nx/2]-5.0,resolution='c',llcrnrlon=lon_coord[-1,0],llcrnrlat=lat_coord[-1,0],urcrnrlon=lon_coord[0,-1],urcrnrlat=lat_coord[0,-1])
	m = Basemap(projection='stere',lon_0=lon_coord[ny/2,nx/2]+7.0,lat_0=lat_coord[ny/2,nx/2],resolution='c',llcrnrlon=lon_coord[-1,0],llcrnrlat=lat_coord[-1,0]+2.0,urcrnrlon=lon_coord[0,-1]-7.0,urcrnrlat=lat_coord[0,-1])
	x,y=m(lon_coord[:],lat_coord[:])
	plt.set_cmap('coolwarm')		
	circles=[30,40,50,60,70]
	meridians=[-40,-30,-20,-10,0,10,20,30,40,50,60]
	
	plt.figure(1)
	plt.clf()
	m.contourf(x,y,(count_hist),np.arange(0,1.05,.05),extend='both')
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	plt.colorbar()
	plt.title('Probability of temperature being above 2014 temp = '+str(count_hist.mean()))
	plt.savefig(folder+'prob_hist.png')
		
	plt.figure(2)
	plt.clf()
	m.contourf(x,y,(count_nat),np.arange(0,0.35,.01),extend='both')
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	plt.colorbar()
	plt.title('Probability of temperature being above 2014 temp \nwith natural forcings only = '+str(count_nat.mean()))
	plt.savefig(folder+'prob_nat.png')
	
	plt.figure(3)
	plt.clf()
	plt.title('Observational dataset of Temperature for 2014 = '+str(obs.mean()))
	m.contourf(x,y,obs,20)
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	plt.colorbar()
	plt.savefig(folder+'obs.png')
	
	plt.figure(4)
	plt.clf()
	plt.title('2014 Anomaly (Observations - Ensemble mean) = '+str(anom.mean()))
	limit=max([anom.max(),abs(anom.min())])
	contours=np.linspace(-limit,limit,20)
	m.contourf(x,y,anom,contours)
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	plt.colorbar()
	plt.savefig(folder+'anom_hist.png')
	
	plt.figure(8)
	plt.clf()
	plt.title('2014 Anomaly (Observations - Nat Ensemble mean) = '+str(anom_nat.mean()))
	limit=max([anom_nat.max(),abs(anom_nat.min())])
	contours=np.linspace(-limit,limit,20)
	m.contourf(x,y,anom_nat,contours)
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	plt.colorbar()
	plt.savefig(folder+'anom_nat.png')	
	
	plt.figure(5)
	plt.clf()
	plt.title('Climatalogical bias used to bias correct model = '+str(bias.mean()))
	limit=max([bias.max(),abs(bias.min())])
	contours=np.linspace(-limit,limit,20)
	m.contourf(x,y,bias,contours)
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	plt.colorbar()
	plt.savefig(folder+'bias.png')
	
	plt.figure(6)
	plt.clf()
	plt.title('FAR = '+str(far.mean()))
	m.contourf(x,y,far,np.arange(.40,1.05,.05),extend='both')
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	plt.colorbar()
	plt.savefig(folder+'far.png')	
	
	plt.figure(7)
	plt.clf()
	plt.title('Change in likelihood of temp greater than 2014 ='+str(change.mean()))
	m.contourf(x,y,change,np.arange(0,22,2),extend='both')
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	plt.colorbar()
	plt.savefig(folder+'change_risk.png')	
	
#	plt.figure(1)
#	plt.clf()
#	m.contourf(x,y,historical_data.min(0))
#	m.drawcoastlines()
#	m.drawcountries()
#	plt.colorbar()
#	plt.title('Ensemble min of temperature')
#	plt.savefig(folder+'hist_min.png')
#	
#	plt.figure(1)
#	plt.clf()
#	m.contourf(x,y,historical_data.max(0))
#	m.drawcoastlines()
#	m.drawcountries()
#	plt.colorbar()
#	plt.title('Ensemble min of temperature')
#	plt.savefig(folder+'hist_max.png')
#	
#	plt.figure(1)
#	plt.clf()
#	m.contourf(x,y,natural_data.min(0))
#	m.drawcoastlines()
#	m.drawcountries()
#	plt.colorbar()
#	plt.title('Ensemble min of temperature')
#	plt.savefig(folder+'nat_min.png')
#	
#	plt.figure(1)
#	plt.clf()
#	m.contourf(x,y,natural_data.max(0))
#	m.drawcoastlines()
#	m.drawcountries()
#	plt.colorbar()
#	plt.title('Ensemble min of temperature')
#	plt.savefig(folder+'nat_max.png')
		
        plt.figure(1)
        plt.clf()
        m.contourf(x,y,historical_data.mean(0)-natural_data.mean(0))
        m.drawcoastlines()
        m.drawcountries()
        plt.colorbar()
        plt.title('Temperature: Historical - Nat')
        plt.savefig(folder+'hist-nat.png')	
	
	
def plot_region(data,lat_coord,lon_coord,lat_s,lat_e,lon_s,lon_e,region_name):
	plt.clf()
	circles=[30,40,50,60,70]
	meridians=[-40,-30,-20,-10,0,10,20,30,40,50,60]
	ny,nx=lat_coord.shape
	m = Basemap(projection='stere',lon_0=lon_coord[ny/2,nx/2]+7.0,lat_0=lat_coord[ny/2,nx/2],resolution='c',llcrnrlon=lon_coord[-1,0],llcrnrlat=lat_coord[-1,0]+2.0,urcrnrlon=lon_coord[0,-1]-7.0,urcrnrlat=lat_coord[0,-1])
	x,y=m(lon_coord[:],lat_coord[:])
	plt.title('Plotting data for region: '+region_name)
	try:
		m.contourf(x[lat_s:lat_e,lon_s:lon_e],y[lat_s:lat_e,lon_s:lon_e],data[lat_s:lat_e,lon_s:lon_e],20)
		plt.colorbar()
	except:
		m.plot(x[lat_s:lat_e,lon_s:lon_e],y[lat_s:lat_e,lon_s:lon_e],'r.') # Just put a dot at the x y point
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	try:
		os.remove('check_regions/._'+region_name+'.png')
	except:
		pass
	plt.savefig('check_regions/'+region_name+'.png')
		
#########################################################################
# Main function
#########################################################################

if __name__=='__main__':


###############  Model climatological bias

	bias_files='/home/cenv0437/scratch/data_from_ouce/EU_???_temp-cruts_0.44_mean.nc'
	bias = ave_bias(bias_files)[:,:]
	
	
################  Observations
	
	# Obs using higher resolution cru.ts climatology
	obs_p5deg='/home/cenv0437/scratch/data_from_ouce/CRU_TS_dec13-nov14_crut4anomalies.nc'

	# Regrid obs to rotated regional grid
	rot_template='/home/cenv0437/scratch/data_from_ouce/hadam3p_eu_z4ao_2013_1_009238311_0_tasmean.nc'
	
############ Get rotated grid info
	
	f_rot=netcdf_file(rot_template,'r')
	lat_coord=f_rot.variables['global_latitude0'][4:-7,4:-4]
	lon_coord=f_rot.variables['global_longitude0'][4:-7,4:-4]
	f_rot.close()
	
############# Remap obs to rotated grid	
	
	obs=remap_p5deg_to_rotated(obs_p5deg,rot_template)[:,:]


	
	# Print debug info of lat and lon indices
	lon_indices,lat_indices=np.meshgrid(np.arange(lat_coord.shape[1]),np.arange(lat_coord.shape[0]))
	print lat_coord.shape
	print lat_indices.shape, lat_indices.max()
	print lon_indices.shape, lon_indices.max()
	
	plot_region(lat_indices,lat_coord,lon_coord,0,-1,0,-1,'Lat_indices')
	plot_region(lon_indices,lat_coord,lon_coord,0,-1,0,-1,'Lon_indices')
	plot_region(lat_coord,lat_coord,lon_coord,0,-1,0,-1,'Lat_coords')
	plot_region((lon_coord-180)%360-180,lat_coord,lon_coord,0,-1,0,-1,'Lon_coords')
	reg_karoly=np.logical_or(np.logical_or(lat_coord<30 ,lat_coord>75) , np.logical_and(lon_coord >45 , lon_coord<348))
#	reg_karoly=np.logical_or(lat_coord<30 ,lat_coord>75)
	plot_region(np.ma.masked_where(reg_karoly,lat_coord),lat_coord,lon_coord,0,-1,0,-1,'EU Karoly Actual')
	reg_geertjan=np.logical_or(np.logical_or(lat_coord<30 ,lat_coord>76) , np.logical_and(lon_coord >45 , lon_coord<335))
	plot_region(np.ma.masked_where(reg_geertjan,lat_coord),lat_coord,lon_coord,0,-1,0,-1,'EU Geert Jan Actual')

#################  Model Data:

	read_data=False
	if read_data or not os.path.exists('historical_data.pkl'):

		infiles=glob.glob('/home/cenv0437/scratch/batch_100/hadam3p_eu_*_tasmean.nc')
		# Have to choose ensemble members by umid
		historical_files=choose_tasknames(infiles,'z200','z2gn')+choose_tasknames(infiles,'z4c0','z4sn')
		natural_files=[x for x in infiles if x not in historical_files]

		clim_hist_files=glob.glob('/home/cenv0437/scratch/batch_43/hadam3p_eu_*_tasmean.nc')
		clim_nat_files=glob.glob('/home/cenv0437/scratch/batch_45/hadam3p_eu_*_tasmean.nc')
	
#################################################
# Load Data	

	# Load historical data into single array and bias correct
		print '2014 historical files:',len(historical_files)
		hist_raw=np.ma.array(load_ensemble(historical_files))[:,:,:]
		# Apply mask (only need to do this when not bias correcting)
		for i in range(len(historical_files)):
				hist_raw[i,:]=np.ma.masked_where(bias.mask,hist_raw[i,:])
		historical_data=hist_raw-bias
		print 'loaded all forcings 2014 data...'
	
	
	# Load natural simulation data into single array and bias correct
		print '2014 natural files:',len(natural_files)
		natural_raw=np.ma.array(load_ensemble(natural_files))[:,:,:]
		# Apply mask (only need to do this when not bias correcting)
		for i in range(len(natural_files)):
				natural_raw[i,:]=np.ma.masked_where(bias.mask,natural_raw[i,:])
		natural_data=natural_raw-bias
		print 'loaded natural 2014 data...'

	######## Climatologies

	# Ensemble of historical simulations between 2000-2012
		print '1985-2011 all forcings files:',len(clim_hist_files)
		clim_hist_data=np.ma.array(load_ensemble(clim_hist_files)[:,:,:])-bias
		print 'loaded all forcings 1985-2011 data...'

	# Ensemble of natural simulations between 2000-2012
		print '1985-2011 natural files:',len(clim_nat_files)
		clim_nat_data=np.ma.array(load_ensemble(clim_nat_files)[:,:,:])-bias
		print 'loaded natural 1985-2011 data...'
		
		# Write data to pickle
		fdata=open('historical_data.pkl','wb')
		pickle.dump(historical_data,fdata,-1)
		fdata.close()
		fdata=open('natural_data.pkl','wb')
		pickle.dump(natural_data,fdata,-1)
		fdata.close()
		fdata=open('clim_hist_data.pkl','wb')
		pickle.dump(clim_hist_data,fdata,-1)
		fdata.close()
		fdata=open('clim_nat_data.pkl','wb')
		pickle.dump(clim_nat_data,fdata,-1)
		fdata.close()

	else: #load from pickle
		fdata=open('historical_data.pkl','rb')
		historical_data=pickle.load(fdata)
		fdata.close()
		fdata=open('natural_data.pkl','rb')
		natural_data=pickle.load(fdata)
		fdata.close()
		fdata=open('clim_hist_data.pkl','rb')
		clim_hist_data=pickle.load(fdata)
		fdata.close()
		fdata=open('clim_nat_data.pkl','rb')
		clim_nat_data=pickle.load(fdata)
		fdata.close()
		print 'loaded data from pkl files'
		
		
#		for i in range(clim_nat_data.shape[0]):
#				clim_nat_data[i,:]=np.ma.masked_where(bias.mask,clim_nat_data[i,:])
#				if np.any(clim_nat_data[i,:]>350.0) or np.any(clim_nat_data[i,:]<170.):
#					clim_nat_data[i,:]=np.ma.masked

#		for i in range(clim_hist_data.shape[0]):
#				clim_hist_data[i,:]=np.ma.masked_where(bias.mask,clim_hist_data[i,:])
#		clim_hist_data=np.ma.masked_invalid(clim_hist_data)
#		clim_hist_data=np.ma.masked_invalid(clim_nat_data)
#		print 'remasked clim data'
		

########################################################
# Print Diagnostics

	plt.figure(1)
	
	hist=[]
	nat=[]
	clim_hist=[]
	clim_nat=[]
	lp=[]
	
	#All of Europe
	vals=get_returntimes(historical_data,natural_data,clim_hist_data,clim_nat_data,obs)
	print vals
	hist.append(vals[0])
	nat.append(vals[1])
	clim_hist.append(vals[2])
	clim_nat.append(vals[3])
	lp.append(vals[4])
	
	ensemble_ave=historical_data.mean(0)

	
	
	
	#Specific regions:
	print 'calcluating region for England...'
	lat_s,lon_s=find_closest(59,-11,lat_coord,lon_coord)
	lat_e,lon_e=find_closest(50,1.2,lat_coord,lon_coord)
	print 'veryfy start:',lat_coord[lat_s,lon_s],lon_coord[lat_s,lon_s]
	print 'veryfy end:',lat_coord[lat_e,lon_e],lon_coord[lat_e,lon_e]
	plot_region(ensemble_ave,lat_coord,lon_coord,lat_s,lat_e,lon_s,lon_e,'England')
	vals=get_region_returntimes(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,lon_s,lon_e,lat_s,lat_e)
	print 'return times',vals
	print 'FAR:',1-vals[0]/vals[1],1-vals[2]/vals[3]
	
	print 'calculating region for Germany...'
	lat_s,lon_s=find_closest(54.7,6.0,lat_coord,lon_coord)
	lat_e,lon_e=find_closest(48,14.,lat_coord,lon_coord)
#	lat_s,lon_s=find_closest(53.5,6.0,lat_coord,lon_coord)
#	lat_e,lon_e=find_closest(46,19,lat_coord,lon_coord)
	print 'boxes:',lat_s,lat_e,lon_s,lon_e
	print 'veryfy start:',lat_coord[lat_s,lon_s],lon_coord[lat_s,lon_s]
	print 'veryfy end:',lat_coord[lat_e,lon_e],lon_coord[lat_e,lon_e]
	plot_region(ensemble_ave,lat_coord,lon_coord,lat_s,lat_e,lon_s,lon_e,'Germany')
	vals=get_region_returntimes(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,lon_s,lon_e,lat_s,lat_e)
        print 'return times',vals
        print 'FAR:',1-vals[0]/vals[1],1-vals[2]/vals[3]
	
	print 'Mediterranean...'
	lat_s,lon_s=find_closest(45.7,0,lat_coord,lon_coord)
	lat_e,lon_e=find_closest(34.8,28,lat_coord,lon_coord)
	print 'veryfy start:',lat_coord[lat_s,lon_s],lon_coord[lat_s,lon_s]
	print 'veryfy end:',lat_coord[lat_e,lon_e],lon_coord[lat_e,lon_e]
	plot_region(ensemble_ave,lat_coord,lon_coord,lat_s,lat_e,lon_s,lon_e,'Mediterranean')
	vals=get_region_returntimes(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,lon_s,lon_e,lat_s,lat_e)
        print 'return times',vals
        print 'FAR:',1-vals[0]/vals[1],1-vals[2]/vals[3]
	
	print 'Scandanavia...'
	lat_s,lon_s=find_closest(71,1.2,lat_coord,lon_coord)
	lat_e,lon_e=find_closest(58.7,30.5,lat_coord,lon_coord)
	print 'veryfy start:',lat_coord[lat_s,lon_s],lon_coord[lat_s,lon_s]
	print 'veryfy end:',lat_coord[lat_e,lon_e],lon_coord[lat_e,lon_e]
	plot_region(ensemble_ave,lat_coord,lon_coord,lat_s,lat_e,lon_s,lon_e,'Scandanavia')
	vals=get_region_returntimes(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,lon_s,lon_e,lat_s,lat_e)
        print 'return times',vals
        print 'FAR:',1-vals[0]/vals[1],1-vals[2]/vals[3]
	
	print 'Central europe'
	lat_s,lon_s=find_closest(54.3,2,lat_coord,lon_coord)
	lat_e,lon_e=find_closest(45.,35.0,lat_coord,lon_coord)
	print 'veryfy start:',lat_coord[lat_s,lon_s],lon_coord[lat_s,lon_s]
	print 'veryfy end:',lat_coord[lat_e,lon_e],lon_coord[lat_e,lon_e]
	plot_region(ensemble_ave,lat_coord,lon_coord,lat_s,lat_e,lon_s,lon_e,'Central Europe')
	vals=get_region_returntimes(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,lon_s,lon_e,lat_s,lat_e)
        print 'return times',vals
        print 'FAR:',1-vals[0]/vals[1],1-vals[2]/vals[3]
	
	print 'Greece:'
	lat_s,lon_s=find_closest(40.8,19.4,lat_coord,lon_coord)
	lat_e,lon_e=find_closest(34.4,27.5,lat_coord,lon_coord)
	print 'veryfy start:',lat_coord[lat_s,lon_s],lon_coord[lat_s,lon_s]
	print 'veryfy end:',lat_coord[lat_e,lon_e],lon_coord[lat_e,lon_e]
	plot_region(ensemble_ave,lat_coord,lon_coord,lat_s,lat_e,lon_s,lon_e,'Greece')
	vals=get_region_returntimes(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,lon_s,lon_e,lat_s,lat_e)
        print 'return times',vals
        print 'FAR:',1-vals[0]/vals[1],1-vals[2]/vals[3]
	
	print 'Italy:'	
	lat_s,lon_s=find_closest(45.8,7.3,lat_coord,lon_coord)
	lat_e,lon_e=find_closest(36.4,17.4,lat_coord,lon_coord)
	print 'veryfy start:',lat_coord[lat_s,lon_s],lon_coord[lat_s,lon_s]
	print 'veryfy end:',lat_coord[lat_e,lon_e],lon_coord[lat_e,lon_e]
	plot_region(ensemble_ave,lat_coord,lon_coord,lat_s,lat_e,lon_s,lon_e,'Italy')
	vals=get_region_returntimes(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,lon_s,lon_e,lat_s,lat_e)
        print 'return times',vals
        print 'FAR:',1-vals[0]/vals[1],1-vals[2]/vals[3]
	
	print 'Oxford:'	
	lat_s,lon_s=find_closest(51.7,-1.2,lat_coord,lon_coord)
	print 'veryfy start:',lat_coord[lat_s,lon_s],lon_coord[lat_s,lon_s]
	plot_region(ensemble_ave,lat_coord,lon_coord,lat_s,lat_s+1,lon_s,lon_s+1,'Oxford')
	vals=get_region_returntimes(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,lon_s,lon_s+1,lat_s,lat_s+1)
        print 'return times',vals
        print 'FAR:',1-vals[0]/vals[1],1-vals[2]/vals[3]
	
	
#TODO: Fix these so they match a bit better!
	print 'EU Geert Jan:'
	lat_s,lon_s=find_closest(76,-25,lat_coord,lon_coord)
	lat_e,lon_e=find_closest(30,45,lat_coord,lon_coord)
	plot_region(ensemble_ave,lat_coord,lon_coord,lat_s,lat_e,lon_s,lon_e,'EU Geert Jan')
	print 'veryfy start:',lat_coord[lat_s,lon_s],lon_coord[lat_s,lon_s]
	print 'veryfy end:',lat_coord[lat_e,lon_e],lon_coord[lat_e,lon_e]
#	plot_region(ensemble_ave,lat_coord,lon_coord,lat_s,lat_e,lon_s,lon_e,'EU Geert Jan')
	vals=get_region_returntimes(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,lon_s,lon_e,lat_s,lat_e)
	print vals	
	
	print 'EU Karoly:'	
	lat_s,lon_s=find_closest(75,-12,lat_coord,lon_coord)
	lat_e,lon_e=find_closest(30,45,lat_coord,lon_coord)
	print 'veryfy start:',lat_coord[lat_s,lon_s],lon_coord[lat_s,lon_s]
	print 'veryfy end:',lat_coord[lat_e,lon_e],lon_coord[lat_e,lon_e]
#	reg_karoly=np.nan*(lat_coord<30 or lat_coord>75 or lon_coord >45 and lon_coord<348)
#	plot_region(ensemble_ave*reg_karoly,lat_coord,lon_coord,0,-1,0,-1,'EU Karoly')
	vals=get_region_returntimes(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,lon_s,lon_e,lat_s,lat_e)
	print vals	
	
	print 'lat_coords,lon_coords:',lat_coord[0,0],lon_coord[0,0],lat_coord[-1,-1],lon_coord[-1,-1]		
	
	
	# Some random regions
#	for i in range(50):
#		try:
#			vals=get_region_returntimes(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,50-i,51+i,0,-1)
#			vals=get_region_returntimes(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,0,-1,50-i,51+i)
#			vals=get_region_returntimes(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,50-i,51+i,50-i,51+i)
#			vals=get_region_returntimes(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,70,71+i,39-i,40+i)
#		except: 
#			continue
#		if vals[4]>0:
#			hist.append(vals[0])
#			nat.append(vals[1])
#			clim_hist.append(vals[2])
#			clim_nat.append(vals[3])
#			lp.append(vals[4])
			
	# Some other random regions
#	for i in range(50):
#		try:
#			vals=get_region_returntimes(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,50-1,51+i,0,-1)
#		except: 
#			continue
#		if vals[4]>0:
#			hist.append(vals[0])
#			nat.append(vals[1])
#			clim_hist.append(vals[2])
#			clim_nat.append(vals[3])
#			lp.append(vals[4])
		
######################################################
# Plot stuff

#First get rid of annoying Mac hidden files
	delfiles=glob.glob('._*.png')
	for f in delfiles:
		os.remove(f)
		
	plt.figure(1)
	plt.plot(lp,hist,'.',label='hist')
	plt.plot(lp,nat,'.',label='nat')
	plt.plot(lp,clim_hist,'.',label='clim_hist')
	plt.plot(lp,clim_nat,'.',label='clim_nat')
	plt.title('return time vs size of region')
	plt.legend()
	plt.ylim([0,700])
	plt.savefig('figure1.png')
		

	plot_spacial(historical_data,natural_data,obs,lat_coord,lon_coord,'spacial_figs_2014')
	plot_spacial(clim_hist_data,clim_nat_data,obs,lat_coord,lon_coord,'spacial_figs_clim')

	#plt.show()
