#!/usr/local/bin/python2.7
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
def create_netcdf_stripped(template,data,outname,lat_nwtart=4,lat_send=7,lon_nwtart=4,lon_send=4):
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
			leng=int(template.dimensions[dim])-lat_nwtart-lat_send
			outfile.createDimension(dim,leng)
			outfile.createVariable(dim,'f',(dim,))
			outfile.variables[dim][:]=template.variables[dim][lat_nwtart:-lat_send]
			outfile.variables[dim].__setattr__('axis','Y')
		elif dim[:3]=='lon' : 
			leng=int(template.dimensions[dim])-lon_nwtart-lon_send
			outfile.createDimension(dim,leng)
			outfile.createVariable(dim,'f',(dim,))
			outfile.variables[dim][:]=template.variables[dim][lon_nwtart:-lon_send]
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
	outfile.variables['global_latitude0'][:]=template.variables['global_latitude0'][lat_nwtart:-lat_send,lon_nwtart:-lon_send]
	for att in template.variables['global_latitude0']._attributes:
			outfile.variables['global_latitude0'].__setattr__(att,template.variables['global_latitude0'].__getattribute__(att))
			
	# Global longitude
	outfile.createVariable('global_longitude0','f',template.variables['global_longitude0'].dimensions)
	outfile.variables['global_longitude0'][:]=template.variables['global_longitude0'][lat_nwtart:-lat_send,lon_nwtart:-lon_send]
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

def plot_spacial(historical_data,natural_data,obs,lat_coord,lon_coord,subfig):	
	folder='spacial_figs'
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
	plt.subplot(subfig)
	m.contourf(x,y,(count_hist),np.arange(0,1.05,.05),extend='both')
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	plt.suptitle('Probability of temperature being above 2014 temp = '+str(count_hist.mean()))
        if subfig%2==1:
                plt.title('a)')
        else:
                plt.title('b)')
                plt.colorbar()
		plt.savefig(folder+'prob_hist.png')
		
	plt.figure(2)
	plt.subplot(subfig)
	m.contourf(x,y,(count_nat),np.arange(0,0.35,.01),extend='both')
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	plt.suptitle('Probability of temperature being above 2014 temp \nwith natural forcings only = '+str(count_nat.mean()))
        if subfig%2==1:
                plt.title('a)')
        else:
                plt.title('b)')
		plt.colorbar()
		plt.savefig(folder+'prob_nat.png')
	
	plt.figure(3)
	plt.subplot(subfig)
	plt.suptitle('Observational dataset of Temperature for 2014 = '+str(obs.mean()))
	m.contourf(x,y,obs,20)
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
        if subfig%2==1:
                plt.title('a)')
        else:
                plt.title('b)')
		plt.colorbar()
		plt.savefig(folder+'obs.png')
	
	plt.figure(4)
	plt.subplot(subfig)
	plt.suptitle('2014 Anomaly (Observations - Ensemble mean) = '+str(anom.mean()))
	limit=max([anom.max(),abs(anom.min())])
	contours=np.linspace(-limit,limit,20)
	m.contourf(x,y,anom,contours)
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
        if subfig%2==1:
                plt.title('a)')
        else:
                plt.title('b)')
		plt.colorbar()
		plt.savefig(folder+'anom_hist.png')
	
	plt.figure(8)
	plt.subplot(subfig)
	plt.suptitle('2014 Anomaly (Observations - Nat Ensemble mean) = '+str(anom_nat.mean()))
	limit=max([anom_nat.max(),abs(anom_nat.min())])
	contours=np.linspace(-limit,limit,20)
	m.contourf(x,y,anom_nat,contours)
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
        if subfig%2==1:
                plt.title('a)')
        else:
                plt.title('b)')
		plt.colorbar()
		plt.savefig(folder+'anom_nat.png')	
	
	plt.figure(5)
	plt.subplot(subfig)
	plt.suptitle('Climatalogical bias used to bias correct model = '+str(bias.mean()))
	limit=max([bias.max(),abs(bias.min())])
	contours=np.linspace(-limit,limit,20)
	m.contourf(x,y,bias,contours)
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
        if subfig%2==1:
                plt.title('a)')
        else:
                plt.title('b)')
		plt.colorbar()
		plt.savefig(folder+'bias.png')
	
	plt.figure(6)
	ga=plt.gca()
	plt.subplot(subfig)
#	plt.title('FAR = '+str(far.mean()))
	c=m.contourf(x,y,far,np.arange(.40,1.05,.05),extend='both')
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	if subfig%2==1:
		plt.title('a)')
	else:
		plt.title('b)')
		# Make an axis for the colorbar on the right side
		cax = plt.gcf().add_axes([0.9, 0.1, 0.03, 0.8])
		plt.colorbar(c, cax=cax)
		plt.savefig(folder+'far.png')	
	
	plt.figure(7)
	plt.subplot(subfig)
	plt.suptitle('Change in likelihood of temp greater than 2014 ='+str(change.mean()))
	c=m.contourf(x,y,change,np.arange(0,22,2),extend='both')
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
        if subfig%2==1:
                plt.title('a)')
        else:
                plt.title('b)')
                cax = plt.gcf().add_axes([0.9, 0.1, 0.03, 0.8])
		plt.colorbar(c,cax=cax)
		plt.tight_layout()
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
		
        plt.figure(9)
        plt.subplot(subfig)
        m.contourf(x,y,historical_data.mean(0)-natural_data.mean(0))
        m.drawcoastlines()
        m.drawcountries()
        if subfig%2==1:
                plt.title('a)')
        else:
                plt.title('b)')
		plt.colorbar()
		plt.tight_layout()
        	plt.savefig(folder+'hist-nat.png')	
	plt.suptitle('Temperature: Historical - Nat')
	
def plot_region(data,lat_coord,lon_coord,lat_nw,lat_se,lon_nw,lon_se,region_name):
	plt.clf()
	circles=[30,40,50,60,70]
	meridians=[-40,-30,-20,-10,0,10,20,30,40,50,60]
	ny,nx=lat_coord.shape
	m = Basemap(projection='stere',lon_0=lon_coord[ny/2,nx/2]+7.0,lat_0=lat_coord[ny/2,nx/2],resolution='c',llcrnrlon=lon_coord[-1,0],llcrnrlat=lat_coord[-1,0]+2.0,urcrnrlon=lon_coord[0,-1]-7.0,urcrnrlat=lat_coord[0,-1])
	x,y=m(lon_coord[:],lat_coord[:])
	plt.title('Plotting data for region: '+region_name)
	try:
		m.contourf(x[lat_nw:lat_se,lon_nw:lon_se],y[lat_nw:lat_se,lon_nw:lon_se],data[lat_nw:lat_se,lon_nw:lon_se],20)
		plt.colorbar()
	except:
		m.plot(x[lat_nw:lat_se,lon_nw:lon_se],y[lat_nw:lat_se,lon_nw:lon_se],'r.') # Just put a dot at the x y point
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	try:
		os.remove('check_regions/._'+region_name+'.png')
	except:
		pass
	plt.savefig('check_regions/'+region_name+'.png')


def print_region_vals(region,lat_nw,lon_nw,lat_se,lon_se,historical_data,natural_data,clim_hist_data,clim_nat_data,lat_coord,lon_coord,table):
	#Specific regions:
	print 'calcluating region: '+region
	y_nw,x_nw=find_closest(lat_nw,lon_nw,lat_coord,lon_coord)
	y_se,x_se=find_closest(lat_se,lon_se,lat_coord,lon_coord)
	# Debug info
	print 'veryfy start:',lat_coord[y_nw,x_nw],lon_coord[y_nw,x_nw]
	print 'veryfy end:',lat_coord[y_se,x_se],lon_coord[y_se,x_se]
	plot_region(historical_data[0,:],lat_coord,lon_coord,y_nw,y_se,x_nw,x_se,region)
	# Calculate the return times
	vals=get_region_returntimes(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,x_nw,x_se,y_nw,y_se)
	print 'return times',vals
	print 'FAR:',1-vals[0]/vals[1],1-vals[2]/vals[3]
#        table.write(region+' & '+str(1-vals[0]/vals[1])+' & '+str(1-vals[2]/vals[3])+'\\\\\n')
	table.write(region+' & '+str(vals[0])+' & '+str(vals[1])+' & '+str(vals[2])+' & '+str(vals[3])+'\\\\\n')
	
def print_point_vals(region,lat_nw,lon_nw,historical_data,natural_data,clim_hist_data,clim_nat_data,lat_coord,lon_coord,table):
	#Specific regions:
	print 'calcluating region: '+region
	lat_nw,lon_nw=find_closest(lat_nw,lon_nw,lat_coord,lon_coord)
	lat_se=lat_nw+1
	lon_se=lon_nw+1
	# Debug info
	print 'veryfy start:',lat_coord[lat_nw,lon_nw],lon_coord[lat_nw,lon_nw]
	print 'veryfy end:',lat_coord[lat_se,lon_se],lon_coord[lat_se,lon_se]
	plot_region(historical_data[0,:],lat_coord,lon_coord,lat_nw,lat_se,lon_nw,lon_se,region)
	# Calculate the return times	
	vals=get_region_returntimes(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,lon_nw,lon_se,lat_nw,lat_se)
	print 'return times',vals
	print 'FAR:',1-vals[0]/vals[1],1-vals[2]/vals[3]
#        table.write(region+' & '+str(1-vals[0]/vals[1])+' & '+str(1-vals[2]/vals[3])+'\\\\\n')
	table.write(region+' & '+str(vals[0])+' & '+str(vals[1])+' & '+str(vals[2])+' & '+str(vals[3])+'\\\\\n')
		
		
		
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
	

# Latex table of return times
	np.set_printoptions(precision=3) #reduce number of sig figures printed
	table=open('table.txt','w')
	table.write('\\begin{tabular}{l c c c c}\n')
#	table.write('Region & FAR 2014 & FAR Clim \\\\\n')
	table.write('Region & All2014 & Nat2014 &AllClim & NatClim \\\\\n')
	table.write('\\hline\n')
	
	#All of Europe
	vals=get_returntimes(historical_data,natural_data,clim_hist_data,clim_nat_data,obs)
	print vals
	print 'FAR:',1-vals[0]/vals[1],1-vals[2]/vals[3]
	region= 'Europe'
#	table.write(region+' & '+str(1-vals[0]/vals[1])+' & '+str(1-vals[2]/vals[3])+'\\\\\n')
	table.write(region+' & '+str(vals[0])+' & '+str(vals[1])+' & '+str(vals[2])+' & '+str(vals[3])+'\\\\\n')
	
	#Specific regions:
	
	print_region_vals('England',59,-11,50,1.2,historical_data,natural_data,clim_hist_data,clim_nat_data,lat_coord,lon_coord,table)
	print_region_vals('Germany',54.7,6.0,48,14.,historical_data,natural_data,clim_hist_data,clim_nat_data,lat_coord,lon_coord,table)
	print_region_vals('Mediterranean',43.7,-5.,34.8,28,historical_data,natural_data,clim_hist_data,clim_nat_data,lat_coord,lon_coord,table)
	print_region_vals('Scandanavia',71,1.2,58.7,30.5,historical_data,natural_data,clim_hist_data,clim_nat_data,lat_coord,lon_coord,table)
	print_region_vals('Central EU',54.3,2,45.,35.0,historical_data,natural_data,clim_hist_data,clim_nat_data,lat_coord,lon_coord,table)
	print_region_vals('Greece',40.8,19.4,34.4,27.5,historical_data,natural_data,clim_hist_data,clim_nat_data,lat_coord,lon_coord,table)
	print_region_vals('Italy',45.8,7.3,36.4,17.4,historical_data,natural_data,clim_hist_data,clim_nat_data,lat_coord,lon_coord,table)
	print_point_vals('Oxford',51.7,-1.2,historical_data,natural_data,clim_hist_data,clim_nat_data,lat_coord,lon_coord,table)

# Finish up table
	table.write('\\end{tabular}\n')
	table.close()


#TODO: Fix these so they match a bit better!
#	print 'EU Geert Jan:'
#	lat_nw,lon_nw=find_closest(76,-25,lat_coord,lon_coord)
#	lat_se,lon_se=find_closest(30,45,lat_coord,lon_coord)
#	plot_region(ensemble_ave,lat_coord,lon_coord,lat_nw,lat_se,lon_nw,lon_se,'EU Geert Jan')
#	print 'veryfy start:',lat_coord[lat_nw,lon_nw],lon_coord[lat_nw,lon_nw]
#	print 'veryfy end:',lat_coord[lat_se,lon_se],lon_coord[lat_se,lon_se]
#	plot_region(ensemble_ave,lat_coord,lon_coord,lat_nw,lat_se,lon_nw,lon_se,'EU Geert Jan')
#	vals=get_region_returntimes(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,lon_nw,lon_se,lat_nw,lat_se)
#	print vals	
	
#	print 'EU Karoly:'	
#	lat_nw,lon_nw=find_closest(75,-12,lat_coord,lon_coord)
#	lat_se,lon_se=find_closest(30,45,lat_coord,lon_coord)
#	print 'veryfy start:',lat_coord[lat_nw,lon_nw],lon_coord[lat_nw,lon_nw]
#	print 'veryfy end:',lat_coord[lat_se,lon_se],lon_coord[lat_se,lon_se]
#	reg_karoly=np.nan*(lat_coord<30 or lat_coord>75 or lon_coord >45 and lon_coord<348)
#	plot_region(ensemble_ave*reg_karoly,lat_coord,lon_coord,0,-1,0,-1,'EU Karoly')
#	vals=get_region_returntimes(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,lon_nw,lon_se,lat_nw,lat_se)
#	print vals	
	
#	print 'lat_coords,lon_coords:',lat_coord[0,0],lon_coord[0,0],lat_coord[-1,-1],lon_coord[-1,-1]		
	
########################################
#  Plot the spacial plots

	plot_spacial(historical_data,natural_data,obs,lat_coord,lon_coord,121)
	plot_spacial(clim_hist_data,clim_nat_data,obs,lat_coord,lon_coord,122)

