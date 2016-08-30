#!/usr/local/bin/python2.7
# NOTE python2.7 must be used for this script

import glob,os,sys
from sort_umids import sort_tasknames,choose_mask,choose_tasknames
import numpy as np
from netcdf_file import netcdf_file
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pickle
from mpl_toolkits.basemap import Basemap,cm
import random
from scipy import stats

pkl_dir='/gpfs/projects/cpdn/scratch/cenv0437/pickle_data/'

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
	i=0
	for f in filenames:
		try:
			tmp=netcdf_file(f,'r').variables['tas'][0,0,4:-7,4:-4] # Cut off border points
			if tmp.max()>350.0 or tmp.min()<170. or not np.all(np.isfinite(tmp)):
				print 'error: wierd vals',f
				continue
			else:
				data[i,:]=tmp
				i=i+1
		except:
			print 'Error, cannot load files',f
			continue
	return data[:i,:]
	
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
	p5deg_nc=netcdf_file(p5deg_file,'r').variables['tmp']
	p5deg_data=p5deg_nc[:12,:].mean(0)
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

# Works with E-OBS data
def remap_p5deg_to_rotated2(p5deg_data,rot_template,fillValue):
        f_rot=netcdf_file(rot_template,'r')
        g_lat=f_rot.variables['global_latitude0'][:] #[4:-7,4:-4]
        g_lon=f_rot.variables['global_longitude0'][:] #[4:-7,4:-4]
        # g_lon2 is in range -180,180
        g_lon2=g_lon.copy()
        g_lon2[g_lon>180.]=g_lon[g_lon>180.]-360.
        print 'lon2:',g_lon2.min(),g_lon2.max()
        print 'local coords',g_lon.shape
        print 'orig coords',p5deg_data.shape
        dataout=np.zeros([108,114])
        for i in range(114):
                for j in range(108):
                	# Mask data outside e-obs region
                	if g_lat[j,i] < 25.25 or g_lat[j,i] > 75.25 or g_lon2[j,i] < -40.25 or g_lon2[j,i] > 75.25:
                		dataout[j,i]=fillValue
                	else:
				# Calculate nearest neighbor indices
				index_lat=int(((g_lat[j,i]-25.25)/0.5)+0.5)
				index_lon=int(((g_lon2[j,i]+40.25)/0.5)+0.5)
				if not p5deg_data.mask[index_lat,index_lon]:
					dataout[j,i]=p5deg_data[index_lat,index_lon]
				else:
					dataout[j,i]=fillValue
#				print g_lat[j,i],g_lon2[j,i],'to',index_lat*.5+25.25,index_lon*.5-40.25
        dataout=np.ma.masked_values(dataout,fillValue)
        return dataout
	
def remap_2deg_to_rotated(p5deg_file,rot_template):
	f_rot=netcdf_file(rot_template,'r')
	g_lat=f_rot.variables['global_latitude0'][4:-7,4:-4]
	g_lon=f_rot.variables['global_longitude0'][4:-7,4:-4]
	p5deg_nc=netcdf_file(p5deg_file,'r').variables['temperature_anomaly']
	p5deg_data=p5deg_nc[:12,:].mean(0)
	print p5deg_data.shape
	dataout=np.zeros([108,114])
	for i in range(114):
		for j in range(108):
			# Calculate nearest neighbor indices 
			index_lat=int(((g_lat[j,i]+89.)/2.)+0.5) 
			index_lon=int((g_lon[j,i]/2.)+0.5)%180
			dataout[j,i]=p5deg_data[index_lat,index_lon]
	dataout=np.ma.masked_values(dataout,p5deg_nc._FillValue)
	return dataout	
	
def remap_1deg_to_rotated(data,rot_template,fillValue):
	f_rot=netcdf_file(rot_template,'r')
	g_lat=f_rot.variables['global_latitude0']#[4:-7,4:-4]
	g_lon=f_rot.variables['global_longitude0']#[4:-7,4:-4]

	print data.shape
	dataout=np.zeros([108,114])
	for i in range(114):
		for j in range(108):
			# Calculate nearest neighbor indices 
			index_lat=int(((g_lat[j,i]+89.5))+1.) 
			index_lon=int((g_lon[j,i])+179.5)%360
			dataout[j,i]=data[index_lat,index_lon]
	dataout=np.ma.masked_values(dataout,fillValue)
	dataout=np.ma.masked_invalid(dataout)
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

#################

# Calculate the return time of the mean of the data being greater than an observational dataset
def ret_time(data,obs):

	ensembles=data.shape[0]*1.0 # number of ensembles. 
	try:
		return ensembles/(data.mean(1).mean(1)>obs.mean(0).mean(0)).sum(0)
	except:
		return np.inf

#################

# Calculate the return time of the mean of the data being greater than an observational dataset
def ret_time_sampled(data,obs):

	ensembles=data.shape[0]*1.0 # number of ensembles. 
	samples=200 # Number of subsamples to try
	samplenumber=ensembles*2/3 # number of ensemble members in each sample
	samplecounts=np.zeros([samples]) # array for counts of each sample
	samplechoice=np.concatenate((np.ones([samples]),np.zeros([ensembles-samples])))==1.

	counts=data.mean(1).mean(1)>obs.mean(0).mean(0)
	# Now subsample counts to get an idea of the uncertainty:
	for i in range(samples):
		np.random.shuffle(samplechoice) # random choice of ensemble members
		samplecounts[i]=counts[samplechoice].sum(0)
		
	samplecounts=samples/samplecounts # convert to ret time
	samplecounts=np.sort(samplecounts)
	try:
		val=ensembles/counts.sum(0)
	except:
		val=1.e20
	min=samplecounts[int(10*samples/100)] # 10th percentile
	max=samplecounts[-int(10*samples/100)] # 90th percentile
	
#	if val==np.inf: val=1.e20
	if max==np.inf or max ==np.nan: max=1.e20
	if min==np.inf or max ==np.nan: min=1.e20
	return min,val,max
#	mean=samplecounts.mean()
#	std=samplecounts.std()
#	return mean-2.*std,ensembles/counts.sum(0),mean+2.*std

#################

# Calculate the return time of the mean of the data being greater than an observational dataset
def ret_time_bootstrap(data,obs):
	
	ensembles=data.shape[0] # number of ensembles. 
	samples=200 # Number of subsamples to try
	samplenumber=ensembles # number of ensemble members in each sample
	samplecounts=np.zeros([samples]) # array for counts of each sample
	sample_data=np.zeros([ensembles]) 

	event=data.mean(1).mean(1)>obs.mean(0).mean(0) # True if event else false
	# Now subsample counts to get an idea of the uncertainty:
	for i in range(samples):
		sample_data=np.zeros([ensembles]) 
		for y in range(ensembles):
			x = random.uniform(0, ensembles)
			sample_data[y] = event[x]
		samplecounts[i]=sample_data.sum(0) # Number of events in sample
		
	samplecounts=ensembles/samplecounts # convert to ret time
	samplecounts=np.sort(samplecounts)
	try:
		val=ensembles*1./event.sum(0)
	except:
		val=1.e20
	min=samplecounts[int(10*samples/100)] # 10th percentile
	max=samplecounts[-int(10*samples/100)] # 90th percentile
	
#	if val==np.inf: val=1.e20
	if max==np.inf or max==np.nan: max=1.e20
	if min==np.inf or max==np.nan: min=1.e20
	return min,val,max
#	mean=samplecounts.mean()
#	std=samplecounts.std()
#	return mean-2.*std,ensembles/counts.sum(0),mean+2.*std

#################

# Calculate the return time of the mean of the data being greater than an observational dataset
def ret_time_bootstrap_stats(data,obs):

	ensembles=data.shape[0] # number of ensembles. 
	samples=1000 # Number of subsamples to try
	samplenumber=ensembles # number of ensemble members in each sample
	samplecounts=np.zeros([samples]) # array for counts of each sample
	sample_data=np.zeros([ensembles])
	if len(data.shape)>1: # Spatial region
		event=data.mean(1).mean(1)>obs.mean(0).mean(0) # True if event else false
	else: # Single data point
		event=data>obs
	# Now subsample counts to get an idea of the uncertainty:
	for i in range(samples):
		picks=np.random.random_integers(0,ensembles-1,ensembles)
		samplecounts[i]=event[picks].sum(0) # Number of events in sample
		
	samplecounts=ensembles/samplecounts # convert to ret time
	samplecounts=np.sort(samplecounts)
	try:
		val=(ensembles*1.)/event.sum(0)
	except:
		val=1.e20
		
	min=stats.scoreatpercentile(samplecounts, 5)
	max=stats.scoreatpercentile(samplecounts, 95)

	if max==np.inf or max==np.nan: max=1.e20
	if min==np.inf or min==np.nan: min=1.e20
	return min,val,max

###############################################################################

# Calculate the return time of the mean of the data being greater than an observational dataset
def ret_time_bootstrap_stats_2d(data,obs):

	ensembles=data.shape[0] # number of ensembles. 
	samples=1000 # Number of subsamples to try
	samplenumber=ensembles # number of ensemble members in each sample
	samplecounts=np.zeros([samples]) # array for counts of each sample
	sample_data=np.zeros([ensembles])
	if len(data.shape)>1: # Spatial region
		event=data.mean(1).mean(1)>obs.mean(0).mean(0) # True if event else false
	else: # Single data point
		event=data>obs
	# Now subsample counts to get an idea of the uncertainty:
	for i in range(samples):
		sample_data=np.zeros([ensembles]) 
		for y in range(ensembles):
			x = random.uniform(0, ensembles)
			sample_data[y] = event[x]
		samplecounts[i]=sample_data.sum(0) # Number of events in sample
		
	samplecounts=ensembles/samplecounts # convert to ret time
	samplecounts=np.sort(samplecounts)
	try:
		val=(ensembles*1.)/event.sum(0)
	except:
		val=1.e20
		
	min=stats.scoreatpercentile(samplecounts, 5)
	max=stats.scoreatpercentile(samplecounts, 95)

	if max==np.inf or max==np.nan: max=1.e20
	if min==np.inf or min==np.nan: min=1.e20
	return min,val,max

###############################################################################

# Calculate the return time of the mean of the data being greater than an observational dataset
def ret_time_bootstrap_stats_spatial(data,obs,samples=1000):

	ensembles=data.shape[0] # number of ensembles. 
	#samples=1000 # Number of subsamples to try
	samplenumber=ensembles # number of ensemble members in each sample
	samplecounts=np.zeros([samples,data.shape[1],data.shape[2]]) # array for counts of each sample
	event=data>obs # True if event else false
	
	# Now subsample counts to get an idea of the uncertainty:
	for i in range(samples):
		sample_data=np.zeros(data.shape) 
		for y in range(ensembles):
			x = random.uniform(0, ensembles)
			sample_data[y] = event[x]
		samplecounts[i,:]=sample_data.sum(0) # Number of events in sample
		
	samplecounts=ensembles/samplecounts # convert to ret time
	samplecounts=np.sort(samplecounts,axis=0)
	try:
		val=(ensembles*1.)/event.sum(0)
	except:
		val=1.e20
		
	min=stats.scoreatpercentile(samplecounts, 5,axis=0)
	max=stats.scoreatpercentile(samplecounts, 95,axis=0)

	max[max==np.inf]=1.e20
	val[val==np.inf]=1.e20
	min[min==np.inf]=1.e20
	max=np.ma.masked_where(obs.mask,max)
	min=np.ma.masked_where(obs.mask,min)
	val=np.ma.masked_where(obs.mask,val)
	return min,val,max
	
#############################################################

def ratio_bootstrap_stats_spatial(hist,nat,obs,samples=1000):

	ensembles_hist=hist.shape[0] # number of hist ensembles.
	ensembles_nat=nat.shape[0]

	counts_hist=np.zeros([samples,hist.shape[1],hist.shape[2]]) # array for counts of each sample
	counts_nat=np.zeros([samples,hist.shape[1],hist.shape[2]]) # array for counts of each sample
	event_hist=hist>obs # True if event else false
	event_nat=nat>obs
	
	# Now resample to get an idea of the uncertainty:
	for i in range(samples):
	
#		sample_hist=np.zeros(hist.shape)
#		for y in range(ensembles_hist):
#			x = random.uniform(0, ensembles_hist)
#			sample_hist[y] = event_hist[x]

#		sample_nat=np.zeros(nat.shape) 
#		for y in range(ensembles_nat):
#			x=random.uniform(0, ensembles_nat)
#			sample_nat[y] = event_nat[x]

		sample_hist=event_hist[np.random.randint(0,high=ensembles_hist,size=ensembles_hist),:]
		sample_nat=event_nat[np.random.randint(0,high=ensembles_nat,size=ensembles_nat),:]

		# Number of events in sample
		counts_hist[i,:]=sample_hist.sum(0) 
		counts_nat[i,:]=sample_nat.sum(0)
		
	# Normalise counts into a return time
	counts_hist=ensembles_hist/counts_hist
	counts_nat=ensembles_nat/counts_nat
	
	# Calculate ratio
	ratios = counts_nat/counts_hist
	# Sort ratios
	ratios = np.sort(ratios,axis=0)
	# Calculate percentiles
	min=stats.scoreatpercentile(ratios, 5,axis=0)
	max=stats.scoreatpercentile(ratios, 95,axis=0)

	# Mask out ocean and set infinity to 1e20
	max[max==np.inf]=1.e20
	min[min==np.inf]=1.e20
	max=np.ma.masked_where(obs.mask,max)
	min=np.ma.masked_where(obs.mask,min)
	return min,max


###############################################################################

def calc_return_time_confidences_2d(ey_data, direction="ascending", c=[0.05, 0.95], bsn=1e5):
	# c = confidence intervals (percentiles) to calculate
	# bsn = boot strap number, number of times to resample the distribution
	
	samples=ey_data.shape[0]*ey_data.shape[1]
	# create the store
	sample_store = numpy.zeros((bsn,samples), 'f')
	# do the resampling
	for s in range(0, int(bsn)):
		t_data = numpy.zeros(ey_data.shape, 'f')
		for y in range(0, ey_data.shape[0]):
			x = random.uniform(0, ey_data.shape[0])
			t_data[y] = ey_data[x,:]
		t_data=t_data.flatten()
		t_data.sort()
		# reverse if necessary
		if direction == "descending":
			t_data = t_data[::-1]
		sample_store[s,:] = t_data
	# now for each confidence interval find the  value at the percentile
	conf_inter = numpy.zeros((len(c), ey_data.shape[0]*ey_data.shape[1]), 'f')
	for c0 in range(0, len(c)):
		for y in range(0, samples):
			data_slice = sample_store[:,y]
			conf_inter[c0,y] = stats.scoreatpercentile(data_slice, c[c0]*100)
	return conf_inter

################################################

def get_model_fnames(natural_files):
	natural_files_dict={}
	natural_files=np.array(natural_files)
	
	model = 'CCSM4'	
	model_mask = choose_mask(natural_files,'z2go','z2m7')+choose_mask(natural_files,'z4so','z5kf')
	print model , sum(model_mask)
	natural_files_dict[model]=natural_files[model_mask]
	
	model = 'CanESM'
	model_mask=choose_mask(natural_files,'z32w','z38f')+choose_mask(natural_files,'z7vs','z8nj')
	print model , sum(model_mask)
	natural_files_dict[model]=natural_files[model_mask]
	
	model = 'CNRM-CM5'
	model_mask=choose_mask(natural_files,'z38g','z3dz')+choose_mask(natural_files,'z8nk','z9fb')
	print model , sum(model_mask)
	natural_files_dict[model]=natural_files[model_mask]

	model = 'CSIRO-Mk3-6-0'
	model_mask=choose_mask(natural_files,'z3e0','z3jk')+choose_mask(natural_files,'z9fc','za73')
	print model , sum(model_mask)
	natural_files_dict[model]=natural_files[model_mask]
	
	model = 'GFDL-CM3'
	model_mask=choose_mask(natural_files,'z2m8','z2rr')+choose_mask(natural_files,'z5kg','z6c7')
	print model , sum(model_mask)
	natural_files_dict[model]=natural_files[model_mask]

	model = 'GISS-E2-H'
	model_mask=choose_mask(natural_files,'z2rs','z2xb')+choose_mask(natural_files,'z6c8','z73z')
	print model , sum(model_mask)
	natural_files_dict[model]=natural_files[model_mask]
	
	model = 'GISS-E2-R'
	model_mask=choose_mask(natural_files,'z3jk','z3p3')+choose_mask(natural_files,'za74','zayv')
	print model , sum(model_mask)
	natural_files_dict[model]=natural_files[model_mask]
	
	model = 'HadGEM2-ES'
	model_mask=choose_mask(natural_files,'z2xc','z32v')+choose_mask(natural_files,'z740','z7vr')
	print model , sum(model_mask)
	natural_files_dict[model]=natural_files[model_mask]
	
	model = 'IPSL-CM5A-LR'
	model_mask=choose_mask(natural_files,'z3p4','z3un')+choose_mask(natural_files,'zayw','zbqn')
	print model , sum(model_mask)
	natural_files_dict[model]=natural_files[model_mask]
	
	model = 'IPSL-CM5A-MR'
	model_mask=choose_mask(natural_files,'z3uo','z407')+choose_mask(natural_files,'zbqo','zcif')
	print model , sum(model_mask)
	natural_files_dict[model]=natural_files[model_mask]
	
	model = 'MIROC-ESM'
	model_mask=choose_mask(natural_files,'z408','z45r')+choose_mask(natural_files,'zcig','zda7')
	print model , sum(model_mask)
	natural_files_dict[model]=natural_files[model_mask]
	
	model = 'MMM'
	model_mask=choose_mask(natural_files,'z45s','z4bb')+choose_mask(natural_files,'zda8','ze1z')
	print model , sum(model_mask)
	natural_files_dict[model]=natural_files[model_mask]
	
	return natural_files_dict


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

def get_region_returntimes2(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,xmin,xmax,ymin,ymax):
	a,b,c,d,e=get_returntimes(historical_data[:,ymin:ymax,xmin:xmax],natural_data[:,ymin:ymax,xmin:xmax],clim_hist_data[:,ymin:ymax,xmin:xmax],clim_nat_data[:,ymin:ymax,xmin:xmax],obs[ymin:ymax,xmin:xmax])
	f=historical_data[:,ymin:ymax,xmin:xmax].mean(1).mean(1).std(0)
	g=natural_data[:,ymin:ymax,xmin:xmax].mean(1).mean(1).std(0)
	return a,b,c,d,e,f,g

def get_region_returntimes_sampled(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,xmin,xmax,ymin,ymax):
	return get_returntimes_sampled(historical_data[:,ymin:ymax,xmin:xmax],natural_data[:,ymin:ymax,xmin:xmax],clim_hist_data[:,ymin:ymax,xmin:xmax],clim_nat_data[:,ymin:ymax,xmin:xmax],obs[ymin:ymax,xmin:xmax])

def get_region_returntimes_bootstrap(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,xmin,xmax,ymin,ymax):
	return get_returntimes_bootstrap(historical_data[:,ymin:ymax,xmin:xmax],natural_data[:,ymin:ymax,xmin:xmax],clim_hist_data[:,ymin:ymax,xmin:xmax],clim_nat_data[:,ymin:ymax,xmin:xmax],obs[ymin:ymax,xmin:xmax])

def get_region_returntimes_bootstrap_stats(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,xmin,xmax,ymin,ymax):
	return get_returntimes_bootstrap_stats(historical_data[:,ymin:ymax,xmin:xmax],natural_data[:,ymin:ymax,xmin:xmax],clim_hist_data[:,ymin:ymax,xmin:xmax],clim_nat_data[:,ymin:ymax,xmin:xmax],obs[ymin:ymax,xmin:xmax])

def get_returntimes_sampled(historical_data,natural_data,clim_hist_data,clim_nat_data,obs):
	ret_hist=ret_time_sampled(historical_data,obs)
	ret_nat=ret_time_sampled(natural_data,obs)
	ret_clim_hist=ret_time_sampled(clim_hist_data,obs)
	ret_clim_nat=ret_time_sampled(clim_nat_data,obs)
	lp=land_points(obs)
	return ret_hist,ret_nat,ret_clim_hist,ret_clim_nat,lp
	
def get_returntimes_bootstrap(historical_data,natural_data,clim_hist_data,clim_nat_data,obs):
	ret_hist=ret_time_bootstrap(historical_data,obs)
	ret_nat=ret_time_bootstrap(natural_data,obs)
	ret_clim_hist=ret_time_bootstrap(clim_hist_data,obs)
	ret_clim_nat=ret_time_bootstrap(clim_nat_data,obs)
	lp=land_points(obs)
	return ret_hist,ret_nat,ret_clim_hist,ret_clim_nat,lp
	
def get_returntimes_bootstrap_stats(historical_data,natural_data,clim_hist_data,clim_nat_data,obs):
	ret_hist=ret_time_bootstrap_stats(historical_data,obs)
	ret_nat=ret_time_bootstrap_stats(natural_data,obs)
	ret_clim_hist=ret_time_bootstrap_stats(clim_hist_data,obs)
	ret_clim_nat=ret_time_bootstrap_stats(clim_nat_data,obs)
	lp=land_points(obs)
	return ret_hist,ret_nat,ret_clim_hist,ret_clim_nat,lp

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
	
##############################################################

# Searches for closest point to p_lat,p_lon and returns grid reference
def find_closest_1d(p_lat,p_lon,lat_coord,lon_coord):
	ny=lat_coord.shape[0]
	nx=lon_coord.shape[0]
	min_dist=100000000
	min_point=(-1,-1)
	for j in range(ny):
		for i in range(nx):
			dist= (lat_coord[j]-p_lat)**2+(np.abs(lon_coord[i]-p_lon)%360)**2
#			print lat_coord[j,i],p_lat,lon_coord[j,i],p_lon,dist
			if dist<min_dist:
				min_point=(j,i)
				min_dist=dist
	print min_point,min_dist
	return min_point
	
# Searches for closest point to p_lat,p_lon and returns grid reference
def find_closest_1d_v2(p_lat,p_lon,lat_coord,lon_coord):
	ny=lat_coord.shape[0]
	nx=lon_coord.shape[0]
	min_dist=100000000
	minx=-1
	miny=-1
	for j in range(ny):
		dist=(lat_coord[j]-p_lat)**2
		if dist<min_dist:
			miny=j
			min_dist=dist
	min_dist=100000000
	print 'plon',p_lon
	for i in range(nx):
		dist= (np.abs(lon_coord[i]-p_lon)%360)**2
		if dist<min_dist:
				minx=i
				min_dist=dist
	return miny,minx

#############################################################

def plot_spacial(historical_data,natural_data,obs,lat_coord,lon_coord,subfig,obsname):	
	folder='spacial_figs_'+obsname
	# Spatially resolved count of ensembles hotter than obs:	
	condition=historical_data>obs
	count_hist=condition.sum(0)/(historical_data.shape[0]*1.0)
	
	condition=natural_data>obs
	count_nat=condition.sum(0)/(natural_data.shape[0]*1.0)

	# Change in likelihood
	change=(count_hist-count_nat)/(count_nat+1e-5)
	
	# Ratio 
	ratio = count_hist / (count_nat +1e-16)
	
	# Anomaly between obs and ensemble mean of historical simulations
	anom=obs-(historical_data.mean(0))

	# Anomaly between obs and ensemble mean of nat simulations
	anom_nat=obs-(natural_data.mean(0))

	# FAR
	far=1.0-((count_nat)/(count_hist)) # Areas where count_hist =0 -> NaN

#	far=1.0-((count_nat)/(count_hist+1e-15))
##	far=far*(np.logical_or(count_hist!=0.0,count_nat!=0.0)*1.0) # Set to zero where count_hist and count_nat==0
#	far=far*((count_hist!=0.0)*1.0)   # Set to zero where count_hist=0
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
	print 'lon0=',lon_coord[ny/2,nx/2]+7.0
	print 'lat0=',lat_coord[ny/2,nx/2]
	print 'llcrnrlon=',lon_coord[-1,0]
	print 'llcrnrlat=',lat_coord[-1,0]+2.0
	print 'urcrnrlon=',lon_coord[0,-1]-7.0
	print 'urcrnrlat=',lat_coord[0,-1]
	
	x,y=m(lon_coord[:],lat_coord[:])
	plt.set_cmap('coolwarm')		
	circles=[30,40,50,60,70]
	meridians=[-40,-30,-20,-10,0,10,20,30,40,50,60]
	
	plt.figure(1)
	plt.subplot(subfig)
	m.pcolor(x,y,(count_hist),vmin=0,vmax=1.)#np.arange(0,1.05,.05),extend='both')
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	plt.suptitle('Probability of temperature being above 2014 temp = '+'{:.2f}'.format(count_hist.mean()))
	if subfig%2==1:
		plt.title('a)')
	else:
		plt.title('b)')
		plt.colorbar(extend='max')
		plt.savefig(folder+'prob_hist.png')
		
	plt.figure(2)
	plt.subplot(subfig)
	m.pcolor(x,y,(count_nat),vmin=0,vmax=0.3)#,np.arange(0,0.35,.01),extend='both')
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	plt.suptitle('Probability of temperature being above 2014 temp \nwith natural forcings only = '+'{:.2f}'.format(count_nat.mean()))
	if subfig%2==1:
		plt.title('a)')
	else:
		plt.title('b)')
		plt.colorbar(extend='max')
		plt.savefig(folder+'prob_nat.png')
	
	plt.figure(3)
	plt.subplot(subfig)
	plt.suptitle('Observational dataset of Temperature for 2014 = '+'{:.2f}'.format(obs.mean()))
	m.pcolor(x,y,obs)
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
	plt.suptitle('2014 Anomaly (Observations - Ensemble mean) = '+'{:.2f}'.format(anom.mean()))
	limit=max([anom.max(),abs(anom.min())])
	contours=np.linspace(-limit,limit,20)
	m.pcolor(x,y,anom,vmin=-3,vmax=3)#,contours)
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	if subfig%2==1:
		plt.title('a)')
	else:
		plt.title('b)')
		plt.colorbar(extend='both')
		plt.savefig(folder+'anom_hist.png')
	
	plt.figure(8)
	plt.subplot(subfig)
	plt.suptitle('2014 Anomaly (Observations - Nat Ensemble mean) = '+'{:.2f}'.format(anom_nat.mean()))
	limit=max([anom_nat.max(),abs(anom_nat.min())])
	contours=np.linspace(-limit,limit,20)
	m.pcolor(x,y,anom_nat,vmin=-3,vmax=3)#,contours)
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	if subfig%2==1:
		plt.title('a)')
	else:
		plt.title('b)')
		plt.colorbar(extend='both')
		plt.savefig(folder+'anom_nat.png')	
	
	plt.figure(5)
	plt.subplot(subfig)
	plt.suptitle('Climatalogical bias used to bias correct model = '+'{:.2f}'.format(bias.mean()))
	limit=max([bias.max(),abs(bias.min())])
	contours=np.linspace(-limit,limit,20)
	m.pcolor(x,y,bias,vmin=-limit/2.,vmax=limit/2.)#,contours)
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	if subfig%2==1:
		plt.title('a)')
	else:
		plt.title('b)')
		plt.colorbar(extend='both')
		plt.savefig(folder+'bias.png')
	
	plt.figure(6)
	plt.subplot(subfig)
	far_axes.append(plt.gca())
#	plt.title('FAR = '+'{:.2f}'.format(far.mean()))
	c=m.pcolor(x,y,far,vmin=0.6,vmax=1.)#,np.arange(.40,1.05,.05),extend='both')
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	if subfig%2==1:
		plt.title('a)')
	else:
		plt.title('b)')
		# Make an axis for the colorbar on the right side
#		cax = plt.gcf().add_axes([0.9, 0.1, 0.03, 0.8])
#		plt.colorbar(c, cax=cax)
		plt.colorbar(c,orientation='horizontal',ax=far_axes,aspect=40,pad=0.05)
		plt.savefig(folder+'far.png')	
	
	plt.figure(7)
	plt.subplot(subfig)
	plt.suptitle('Ratio of likelihoods of temp greater than 2014 ='+'{:.2f}'.format(ratio.mean()))
#	change[change==0.0]=1.e-16
	print 'ratio',ratio.min(),ratio.max()
	c=m.pcolor(x,y,change,norm=LogNorm(vmin=1., vmax=1.e3))
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
		plt.savefig(folder+'ratio_risk.png')	
		
	plt.figure(13)
	plt.subplot(subfig)
	plt.suptitle('Relative trend in temperature')
	hist_mean=historical_data.mean(0)
	nat_mean=natural_data.mean(0)
	c=m.pcolor(x,y,(hist_mean-nat_mean)/(hist_mean.mean()-nat_mean.mean()) ,vmin=0,vmax=2)
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	if subfig%2==1:
		plt.title('a)')
	else:
		plt.title('b)')
		cax = plt.gcf().add_axes([0.9, 0.1, 0.03, 0.8])
		plt.colorbar(c,cax=cax,extend='max')
		plt.tight_layout()
		plt.savefig(folder+'trend.png')	
		
	plt.figure(17)
	plt.subplot(subfig)
	hist_mean=historical_data.mean(0)
	nat_mean=natural_data.mean(0)
	plt.suptitle('Absolute trend in temperature: '+str((hist_mean-nat_mean).mean()))
	c=m.pcolor(x,y,(hist_mean-nat_mean) ,vmin=0,vmax=2)
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	if subfig%2==1:
		plt.title('a)')
	else:
		plt.title('b)')
		cax = plt.gcf().add_axes([0.9, 0.1, 0.03, 0.8])
		plt.colorbar(c,cax=cax,extend='max')
		plt.tight_layout()
		plt.savefig(folder+'trend_abs.png')	
		
	plt.figure(14)
	plt.subplot(subfig)
	plt.suptitle('Standard deviation of temperature (hist)')
	c=m.pcolor(x,y,(historical_data.std(0)) ,vmin=0,vmax=2)
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	if subfig%2==1:
		plt.title('a)')
	else:
		plt.title('b)')
		cax = plt.gcf().add_axes([0.9, 0.1, 0.03, 0.8])
		plt.colorbar(c,cax=cax,extend='max')
		plt.tight_layout()
		plt.savefig(folder+'stdev_hist.png')
	
	plt.figure(15)
	plt.subplot(subfig)
	plt.suptitle('Standard deviation of temperature (nat)')
	c=m.pcolor(x,y,(natural_data.std(0)) ,vmin=0,vmax=2)
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	if subfig%2==1:
		plt.title('a)')
	else:
		plt.title('b)')
		cax = plt.gcf().add_axes([0.9, 0.1, 0.03, 0.8])
		plt.colorbar(c,cax=cax,extend='max')
		plt.tight_layout()
		plt.savefig(folder+'stdev_nat.png')


	plt.figure(16)
	plt.subplot(subfig)
	plt.suptitle('Standard deviation of temperature (all)')
	all_detrend=np.ma.concatenate((natural_data,historical_data-(hist_mean-nat_mean)),axis=0)
	c=m.pcolor(x,y,all_detrend.std(0),vmin=0,vmax=2)
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	if subfig%2==1:
		plt.title('a)')
	else:
		plt.title('b)')
		cax = plt.gcf().add_axes([0.9, 0.1, 0.03, 0.8])
		plt.colorbar(c,cax=cax,extend='max')
		plt.tight_layout()
		plt.savefig(folder+'stdev_all.png')


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
		plt.suptitle('Temperature: Historical - Nat')
		plt.savefig(folder+'hist-nat.png')	

	

	
def plot_region(data,lat_coord,lon_coord,lat_nw,lat_se,lon_nw,lon_se,region_name):
	plt.clf()
	circles=[30,40,50,60,70]
	meridians=[-40,-30,-20,-10,0,10,20,30,40,50,60]
	ny,nx=lat_coord.shape
	m = Basemap(projection='stere',lon_0=lon_coord[ny/2,nx/2]+7.0,lat_0=lat_coord[ny/2,nx/2],resolution='c',llcrnrlon=lon_coord[-1,0],llcrnrlat=lat_coord[-1,0]+2.0,urcrnrlon=lon_coord[0,-1]-7.0,urcrnrlat=lat_coord[0,-1])
	x,y=m(lon_coord[:],lat_coord[:])
	try:
		m.pcolor(x[lat_nw:lat_se,lon_nw:lon_se],y[lat_nw:lat_se,lon_nw:lon_se],data[lat_nw:lat_se,lon_nw:lon_se])
		plt.colorbar()
	except:
		plt.clf()
		m.plot(x[lat_nw:lat_se,lon_nw:lon_se],y[lat_nw:lat_se,lon_nw:lon_se],'r.') # Just put a dot at the x y point
	plt.title('Plotting data for region: '+region_name)
	m.drawcoastlines()
	m.drawcountries()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	try:
		os.remove('check_regions/._'+region_name+'.png')
	except:
		pass
	plt.savefig('check_regions/'+region_name+'.png')


def print_region_vals(region,lat_nw,lon_nw,lat_se,lon_se,historical_data,natural_data,clim_hist_data,clim_nat_data,lat_coord,lon_coord):
	#Specific regions:
	print 'calcluating region: '+region
	y_nw,x_nw=find_closest(lat_nw,lon_nw,lat_coord,lon_coord)
	y_se,x_se=find_closest(lat_se,lon_se,lat_coord,lon_coord)
	# Debug info
	print 'region:',y_nw,y_se,x_nw,x_se 
	print 'veryfy start:',lat_coord[y_nw,x_nw],lon_coord[y_nw,x_nw]
	print 'veryfy end:',lat_coord[y_se,x_se],lon_coord[y_se,x_se]
	plot_region(historical_data[0,:],lat_coord,lon_coord,y_nw,y_se,x_nw,x_se,region)
	# Calculate the return times
	vals=get_region_returntimes(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,x_nw,x_se,y_nw,y_se)
	print 'return times',vals
	print 'FAR:',1-vals[0]/vals[1],1-vals[2]/vals[3]
#	table.write(region+' & '+'{:.2f}'.format(1-vals[0]/vals[1])+' & '+'{:.2f}'.format(1-vals[2]/vals[3])+'\\\\\n')
	table.write(region+' & '+'{:.2f}'.format(vals[0])+' & '+'{:.2f}'.format(vals[1])+' & '+'{:.3f}'.format(1-vals[0]/vals[1])+' & '+'{:.2f}'.format(vals[2])+' & '+'{:.2f}'.format(vals[3])+' & '+'{:.3f}'.format(1-vals[2]/vals[3])+'\\\\\n')
	
def print_region_vals2(region,lat_nw,lon_nw,lat_se,lon_se,historical_data,natural_data,clim_hist_data,clim_nat_data,lat_coord,lon_coord):
	#Specific regions:
	print 'calcluating region: '+region
	y_nw,x_nw=find_closest(lat_nw,lon_nw,lat_coord,lon_coord)
	y_se,x_se=find_closest(lat_se,lon_se,lat_coord,lon_coord)
	# Debug info
	print 'veryfy start:',lat_coord[y_nw,x_nw],lon_coord[y_nw,x_nw]
	print 'veryfy end:',lat_coord[y_se,x_se],lon_coord[y_se,x_se]
	plot_region(historical_data[0,:],lat_coord,lon_coord,y_nw,y_se,x_nw,x_se,region)
	
	# Print mean and standard dev over region:
	obsmean=obs[y_nw:y_se,x_nw:x_se].mean()
	print 'obs:',obsmean
	print 'historical2014: mean, anomaly,stdev',historical_data[:,y_nw:y_se,x_nw:x_se].mean(),historical_data[:,y_nw:y_se,x_nw:x_se].mean(1).mean(1).std(0)
	print 'historicalClim: mean, anomaly,stdev',clim_hist_data[:,y_nw:y_se,x_nw:x_se].mean(),clim_hist_data[:,y_nw:y_se,x_nw:x_se].mean(1).mean(1).std(0)
	print 'natural2014: mean, anomaly,stdev',natural_data[:,y_nw:y_se,x_nw:x_se].mean(),natural_data[:,y_nw:y_se,x_nw:x_se].mean(1).mean(1).std(0)
	print 'historicalClim: mean, anomaly,stdev',clim_nat_data[:,y_nw:y_se,x_nw:x_se].mean(),clim_nat_data[:,y_nw:y_se,x_nw:x_se].mean(1).mean(1).std(0)
	# Calculate the return times
	vals=get_region_returntimes_bootstrap_stats(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,x_nw,x_se,y_nw,y_se)
	print 'return times 2014: actual,natural',vals[0],vals[1]
	print 'return times clim: actual,natural',vals[2],vals[3]
	
	plt.figure()
	plt.hist(clim_hist_data[:,y_nw:y_se,x_nw:x_se].mean(1).mean(1),200,normed=1,facecolor='blue',alpha=0.4,label='actual')
	plt.hist(clim_nat_data[:,y_nw:y_se,x_nw:x_se].mean(1).mean(1),200,normed=1,facecolor='green',alpha=0.4,label='natural')
	plt.axvline(x=obsmean,color='r')
	plt.legend()
	plt.savefig('histo_clim.png')
	
	
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
#	table.write(region+' & '+'{:.2f}'.format(1-vals[0]/vals[1])+' & '+'{:.2f}'.format(1-vals[2]/vals[3])+'\\\\\n')
	table.write(region+' & '+'{:.2f}'.format(vals[0])+' & '+'{:.2f}'.format(vals[1])+' & '+'{:.3f}'.format(1-vals[0]/vals[1])+' & '+'{:.2f}'.format(vals[2])+' & '+'{:.2f}'.format(vals[3])+' & '+'{:.3f}'.format(1-vals[2]/vals[3])+'\\\\\n')	
		
		
#########################################################################
# Main function
#########################################################################

if __name__=='__main__':


###############  Model climatological bias

	bias_files='/home/cenv0437/scratch/data_from_ouce/EU_???_temp-cruts_0.44_mean.nc'
	bias = ave_bias(bias_files)[:,:]
	print 'bias mask',bias.mask.sum()
	
	
################  Observations
	
	# Obs using higher resolution cru.ts climatology
#	obs_p5deg='/home/cenv0437/scratch/data_from_ouce/CRU_TS_dec13-nov14_crut4anomalies.nc'


	obs_p5deg= '/home/cenv0437/scratch/data_from_ouce/CRU_TS_Absolute_plus_BEST_anomaly_201312_201412.nc'
	obsname='BEST'
	
#	obs_p5deg= '/home/cenv0437/scratch/data_from_ouce/CRU_TS_Absolute_plus_CRU4_anomaly_201312_201412.nc'
#	obsname='CRU4'
	
#	obs_p5deg= '/home/cenv0437/scratch/data_from_ouce/CRU_TS_Absolute_plus_GISS_anomaly_201312_201412.nc'
#	obsname='GISS'	

	# Regrid obs to rotated regional grid
	rot_template='/home/cenv0437/scratch/data_from_ouce/hadam3p_eu_z4ao_2013_1_009238311_0_tasmean.nc'
	rot_template='/home/cenv0437/scratch/data_from_ouce/hadam3p_eu_template.nc'
	
############ Get rotated grid info
	
	f_rot=netcdf_file(rot_template,'r')
	lat_coord=f_rot.variables['global_latitude0'][4:-7,4:-4]
	lon_coord=f_rot.variables['global_longitude0'][4:-7,4:-4]
	f_rot.close()
	
############# Remap obs to rotated grid	
	
	obs=remap_p5deg_to_rotated(obs_p5deg,rot_template)[:,:]
	# Mask out area masked in the bias correction file
	obs=np.ma.masked_where(bias.mask,obs)
	print 'obs mask',obs.mask.sum()
	


	
	# Print debug info of lat and lon indices
	lon_indices,lat_indices=np.meshgrid(np.arange(lat_coord.shape[1]),np.arange(lat_coord.shape[0]))
	print lat_coord.shape
	print lat_indices.shape, lat_indices.max()
	print lon_indices.shape, lon_indices.max()
	
#	plot_region(lat_indices,lat_coord,lon_coord,0,-1,0,-1,'Lat_indices')
#	plot_region(lon_indices,lat_coord,lon_coord,0,-1,0,-1,'Lon_indices')
#	plot_region(lat_coord,lat_coord,lon_coord,0,-1,0,-1,'Lat_coords')
#	plot_region((lon_coord-180)%360-180,lat_coord,lon_coord,0,-1,0,-1,'Lon_coords')
	
	reg_central_eng=~np.logical_and(np.logical_and(lat_coord>51 ,lat_coord<54) , lon_coord>357)
#	plot_region(np.ma.masked_where(reg_central_eng,lon_coord),lat_coord,lon_coord,0,-1,0,-1,'Central England')
	
	lsm = netcdf_file('/home/cenv0437/scratch/data_from_ouce/land_mask_eu50.nc').variables['lsm'][0,0,4:-7,4:-4]+1
	reg_north_eng=~np.logical_and(np.logical_and(lat_coord>54 ,lat_coord<57) , np.logical_or(lon_coord>354,lon_coord<2))
	reg_north_eng=np.logical_or(reg_north_eng,lsm)
#	plot_region(np.ma.masked_where(reg_north_eng,lon_coord),lat_coord,lon_coord,0,-1,0,-1,'Northern England')
	
	# Netherlands point
	buffer=0
	lat_i,lon_i=find_closest(52,5,lat_coord,lon_coord)
	print lsm[lat_i-buffer:lat_i+1+buffer,lon_i-buffer:lon_i+1+buffer]
	plot_region(lsm,lat_coord,lon_coord,lat_i-buffer,lat_i+1+buffer,lon_i-buffer,lon_i+1+buffer,'52N5E')
	
	reg_karoly=np.logical_or(np.logical_or(lat_coord<30 ,lat_coord>75) , np.logical_and(lon_coord >45 , lon_coord<348))
#	reg_karoly=np.logical_or(lat_coord<30 ,lat_coord>75)
	plot_region(np.ma.masked_where(reg_karoly,lat_coord),lat_coord,lon_coord,0,-1,0,-1,'EU Karoly Actual')
	reg_geertjan=np.logical_or(np.logical_or(lat_coord<30 ,lat_coord>76) , np.logical_and(lon_coord >45 , lon_coord<335))
	plot_region(np.ma.masked_where(reg_geertjan,lat_coord),lat_coord,lon_coord,0,-1,0,-1,'EU Geert Jan Actual')

#################  Model Data:

	read_data=False
	if read_data or not os.path.exists(pkl_dir+'historical_data.pkl'):

		infiles=glob.glob('/home/cenv0437/scratch/batch_100/hadam3p_eu_*_tasmean.nc')
		# Have to choose ensemble members by umid
		historical_files=choose_tasknames(infiles,'z200','z2gn')+choose_tasknames(infiles,'z4c0','z4sn')
		natural_files=[x for x in infiles if x not in historical_files]

		print 'len historical batch100',len(historical_files)
		historical_files=historical_files+glob.glob('/home/cenv0437/scratch/batch_166/hadam3p_eu_????_*_tasmean.nc') #addition from batch 166


		clim_hist_files=glob.glob('/home/cenv0437/scratch/batch_43/hadam3p_eu_[op]???_20*_tasmean.nc')
		clim_nat_files=glob.glob('/home/cenv0437/scratch/batch_45/hadam3p_eu_q???_20*_tasmean.nc')
	
#################################################
# Load Data	

	# Load historical data into single array and bias correct
		print '2014 historical files:',len(historical_files)
		hist_raw=np.ma.array(load_ensemble(historical_files))[:,:,:]
		# Apply mask from obs
		for i in range(len(historical_files)):
				hist_raw[i,:]=np.ma.masked_where(obs.mask,hist_raw[i,:])
		historical_data=hist_raw-bias
		print 'loaded all forcings 2014 data...'
	
	
	# Load natural simulation data into single array and bias correct
		print '2014 natural files:',len(natural_files)
		natural_raw=np.ma.array(load_ensemble(natural_files))[:,:,:]
		# Apply mask from obs
		for i in range(len(natural_files)):
				natural_raw[i,:]=np.ma.masked_where(obs.mask,natural_raw[i,:])
		natural_data=natural_raw-bias
		print 'loaded natural 2014 data...'

	######## Climatologies

	# Ensemble of historical simulations between 2000-2012
		print '1985-2011 all forcings files:',len(clim_hist_files)
		clim_hist_data=np.ma.array(load_ensemble(clim_hist_files)[:,:,:])-bias
		for i in range(len(clim_hist_data)):
			clim_hist_data[i,:]=np.ma.masked_where(obs.mask,clim_hist_data[i,:])
		print 'loaded all forcings 1985-2011 data...'

	# Ensemble of natural simulations between 2000-2012
		print '1985-2011 natural files:',len(clim_nat_files)
		clim_nat_data=np.ma.array(load_ensemble(clim_nat_files)[:,:,:])-bias
		# Mask out areas masked in obs
		for i in range(len(clim_nat_data)):
			clim_nat_data[i,:]=np.ma.masked_where(obs.mask,clim_nat_data[i,:])
		print 'loaded natural 1985-2011 data...'
		
		# Write data to pickle
		fdata=open(pkl_dir+'historical_data.pkl','wb')
		pickle.dump(historical_data,fdata,-1)
		fdata.close()
		fdata=open(pkl_dir+'natural_data.pkl','wb')
		pickle.dump(natural_data,fdata,-1)
		fdata.close()
		fdata=open(pkl_dir+'clim_hist_data.pkl','wb')
		pickle.dump(clim_hist_data,fdata,-1)
		fdata.close()
		fdata=open(pkl_dir+'clim_nat_data.pkl','wb')
		pickle.dump(clim_nat_data,fdata,-1)
		fdata.close()

	else: #load from pickle
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
		print 'data masks',historical_data.mask[0,:].sum()
		print 'data masks',natural_data.mask[0,:].sum()
		print 'data masks',clim_hist_data.mask[0,:].sum()
		print 'data masks',clim_nat_data.mask[0,:].sum()
			
		
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
#	np.set_printoptions(precision=3) #reduce number of sig figures printed
#	table=open('table.txt','w')
#	table.write('\\begin{tabular}{l c c c c c c}\n')
#	table.write('Region & FAR 2014 & FAR Clim \\\\\n')
#	table.write('Region & All2014 & Nat2014 & FAR2014 & AllClim & NatClim & FARClim \\\\\n')
#	table.write('\\hline\n')
	
	#All of Europe
#	print "Vals"
#	print get_returntimes(historical_data,natural_data,clim_hist_data,clim_nat_data,obs)
#	print "sampled"
#	print get_returntimes_sampled(historical_data,natural_data,clim_hist_data,clim_nat_data,obs)
#	print "bootstrap"
#	print get_returntimes_bootstrap(historical_data,natural_data,clim_hist_data,clim_nat_data,obs)
#	print "bootstrap_stats"
	print get_returntimes_bootstrap_stats(historical_data,natural_data,clim_hist_data,clim_nat_data,obs)
	
	# France:
#	print_region_vals2('France',49.6,-5.3,43.6,7.5,historical_data,natural_data,clim_hist_data,clim_nat_data,lat_coord,lon_coord)
	
	
	#Heat wave region
#	print_region_vals('HeatWave',53.2,-5.5,43.8,15.4,historical_data,natural_data,clim_hist_data,clim_nat_data,lat_coord,lon_coord)
#	print_region_vals('HeatWave_larger',55.,-5.,42.,20.,historical_data,natural_data,clim_hist_data,clim_nat_data,lat_coord,lon_coord)
	
	
#	print vals
#	print 'FAR:',1-vals[0]/vals[1],1-vals[2]/vals[3]
#	region= 'Europe'
#	table.write(region+' & '+'{:.2f}'.format(1-vals[0]/vals[1])+' & '+'{:.2f}'.format(1-vals[2]/vals[3])+'\\\\\n')
#	table.write(region+' & '+'{:.2f}'.format(vals[0])+' & '+'{:.2f}'.format(vals[1])+' & '+'{:.3f}'.format(1-vals[0]/vals[1])+' & '+'{:.2f}'.format(vals[2])+' & '+'{:.2f}'.format(vals[3])+' & '+'{:.3f}'.format(1-vals[2]/vals[3])+'\\\\\n')	
	#Specific regions:
	
	#print_region_vals('England',59,-11,50,1.2,historical_data,natural_data,clim_hist_data,clim_nat_data,lat_coord,lon_coord,table)
	#print_region_vals('Germany',54.7,6.0,48,14.,historical_data,natural_data,clim_hist_data,clim_nat_data,lat_coord,lon_coord,table)
	#print_region_vals('Mediterranean',43.7,-5.,34.8,28,historical_data,natural_data,clim_hist_data,clim_nat_data,lat_coord,lon_coord,table)
	#print_region_vals('Scandanavia',71,1.2,58.7,30.5,historical_data,natural_data,clim_hist_data,clim_nat_data,lat_coord,lon_coord,table)
#	print_region_vals('Central EU',54.3,2,45.,35.0,historical_data,natural_data,clim_hist_data,clim_nat_data,lat_coord,lon_coord,table)
#	print_region_vals('Greece',40.8,19.4,34.4,27.5,historical_data,natural_data,clim_hist_data,clim_nat_data,lat_coord,lon_coord,table)
	#print_region_vals('Italy',45.8,7.3,36.4,17.4,historical_data,natural_data,clim_hist_data,clim_nat_data,lat_coord,lon_coord,table)
	#print_point_vals('Oxford',51.7,-1.2,historical_data,natural_data,clim_hist_data,clim_nat_data,lat_coord,lon_coord,table)

# Finish up table
#	table.write('\\end{tabular}\n')
#	table.close()


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

	far_axes=[]
	plot_spacial(historical_data,natural_data,obs,lat_coord,lon_coord,121,obsname)
	plot_spacial(clim_hist_data,clim_nat_data,obs,lat_coord,lon_coord,122,obsname)

