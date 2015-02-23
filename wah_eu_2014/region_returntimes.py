import glob,os
from sort_umids import sort_tasknames,choose_mask,choose_tasknames
import numpy as np
from netcdf_file import netcdf_file
import matplotlib.pyplot as plt
import pickle

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
	data=np.zeros([len(filenames),108,114])
	for i,f in enumerate(filenames):
		tmp=netcdf_file(f,'r').variables['tas'][0,0,4:-7,4:-4] # Cut off border points
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
	return ensembles/(data.mean(1).mean(1)>obs.mean(0).mean(0)).sum(0)

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


#########################################################################
# Main function
#########################################################################

if __name__=='__main__':


###############  Model climatological bias

	bias_files='/ouce-home/staff/cenv0270/CPDN/Europe_2014/observations/CRUT4/EU_???_temp-cruts_0.44_mean.nc'
	bias = ave_bias(bias_files)[:,:]
	
	
################  Observations
	
	# CRUTEM observational set
	obs_N96='/ouce-home/staff/cenv0270/CPDN/Europe_2014/observations/final/CRU_T2m_dec-nov_N96_interp_lsm.nc'

	# Obs using GFS analysis
#	obs_p5deg='/ouce-home/staff/cenv0270/CPDN/Europe_2014/observations/HALF_final/GFS_dec13-nov14_biascorrected.nc'

	# Obs using higher resolution cru.ts climatology
	obs_p5deg='/ouce-home/staff/cenv0270/CPDN/Europe_2014/observations/HALF_final/CRU_TS_dec13-nov14_crut4anomalies.nc'

	# Regrid obs to rotated regional grid
	rot_template='/ouce-home/staff/cenv0437/data/batch_100/hadam3p_eu_z4ao_2013_1_009238311_0_tasmean.nc'
#	obs=remap_N96_to_rotated(obs_N96,rot_template)
	obs=remap_p5deg_to_rotated(obs_p5deg,rot_template)[:,:]


#################  Model Data:

	infiles=glob.glob('/ouce-home/staff/cenv0437/data/batch_100/hadam3p_eu_*_tasmean.nc')
	# Have to choose ensemble members by umid
	historical_files=choose_tasknames(infiles,'z200','z2gn')+choose_tasknames(infiles,'z4c0','z4sn')
	natural_files=[x for x in infiles if x not in historical_files]

	clim_hist_files=glob.glob('/ouce-home/staff/cenv0437/data/batch_43/hadam3p_eu_*_tasmean.nc')
	clim_nat_files=glob.glob('/ouce-home/staff/cenv0437/data/batch_45/hadam3p_eu_*_tasmean.nc')

#################################################
# Load Data	

	read_data=False
	if read_data:
		
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
		clim_hist_data=load_ensemble(clim_hist_files)[:,:,:]-bias
		print 'loaded all forcings 1985-2011 data...'

	# Ensemble of natural simulations between 2000-2012
		print '1985-2011 natural files:',len(clim_nat_files)
		clim_nat_data=load_ensemble(clim_nat_files)[:,:,:]-bias
		print 'loaded natural 1985-2011 data...'
		
		# Write data to pickle
		fdata=open('data.pkl','wb')
		pickle.dump(historical_data,fdata,-1)
		pickle.dump(natural_data,fdata,-1)
		pickle.dump(clim_hist_data,fdata,-1)
		pickle.dump(clim_nat_data,fdata,-1)
		fdata.close()

	else: #load from pickle
		fdata=open('data.pkl','rb')
		historical_data=pickle.load(fdata)
		natural_data=pickle.load(fdata)
		clim_hist_data=pickle.load(fdata)
		clim_nat_data=pickle.load(fdata)
		fdata.close()	
		print 'loaded data from data.pkl'
		

########################################################
# Print Diagnostics

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
	
	# Some random regions
	for i in range(50):
		try:
#			vals=get_region_returntimes(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,50-i,51+i,0,-1)
#			vals=get_region_returntimes(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,0,-1,50-i,51+i)
			vals=get_region_returntimes(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,50-i,51+i,50-i,51+i)
#			vals=get_region_returntimes(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,70,71+i,39-i,40+i)
		except: 
			continue
		if vals[4]>0:
			hist.append(vals[0])
			nat.append(vals[1])
			clim_hist.append(vals[2])
			clim_nat.append(vals[3])
			lp.append(vals[4])
			
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
		

	#plt.show()
