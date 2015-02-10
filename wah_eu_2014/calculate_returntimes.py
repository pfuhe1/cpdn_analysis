import glob,os
from sort_umids import sort_tasknames
import numpy as np
from netcdf_file import netcdf_file
import matplotlib.pyplot as plt

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


if __name__=='__main__':
	infiles='/ouce-home/staff/cenv0437/data/batch_100/hadam3p_eu_*_tasmean.nc'
	bias_files='/ouce-home/staff/cenv0270/CPDN/Europe_2014/observations/CRUT4/EU_???_temp-cruts_0.44_mean.nc'
	bias = ave_bias(bias_files)[:,:-10]
	(historical_files,natural_files)=sort_tasknames(infiles)

	# Load historical data into single array and bias correct
	print len(historical_files)
	hist_raw=np.ma.array(load_ensemble(historical_files))[:,:,:-10]
	# Apply mask (only need to do this when not bias correcting)
	for i in range(len(historical_files)):
	        hist_raw[i,:]=np.ma.masked_where(bias.mask,hist_raw[i,:])
	historical_data=hist_raw-bias
	
	# Load natural simulation data into single array and bias correct
	print len(natural_files)
	natural_raw=np.ma.array(load_ensemble(natural_files))[:,:,:-10]
	# Apply mask (only need to do this when not bias correcting)
	for i in range(len(natural_files)):
	        natural_raw[i,:]=np.ma.masked_where(bias.mask,natural_raw[i,:])
	natural_data=natural_raw-bias
	
	# CRUTEM observational set
	obs_N96='/ouce-home/staff/cenv0270/CPDN/Europe_2014/observations/final/CRU_T2m_dec-nov_N96_interp_lsm.nc'

	# Obs using GFS analysis
#	obs_p5deg='/ouce-home/staff/cenv0270/CPDN/Europe_2014/observations/HALF_final/GFS_dec13-nov14_biascorrected.nc'

	# Obs using higher resolution cru.ts climatology
	obs_p5deg='/ouce-home/staff/cenv0270/CPDN/Europe_2014/observations/HALF_final/CRU_TS_dec13-nov14_crut4anomalies.nc'

	# Regrid obs to rotated regional grid
	rot_template='/ouce-home/staff/cenv0437/data/batch_100/hadam3p_eu_z4ao_2013_1_009238311_0_tasmean.nc'
#	obs=remap_N96_to_rotated(obs_N96,rot_template)
	obs=remap_p5deg_to_rotated(obs_p5deg,rot_template)[:,:-10]

	oldfiles=glob.glob('/ouce-home/staff/cenv0437/data/batch_43/hadam3p_eu_*_tasmean.nc')
	print len(oldfiles)
        old_raw=load_ensemble(oldfiles)[:,:,:-10]
        old_data=old_raw-bias

#################################################
# Calculate diagnostics	
	
	# Count number of model runs which are above the obs
	condition=historical_data>obs
	count_hist=condition.sum(0)/(len(historical_files)*1.0)
	
	condition=natural_data>obs
	count_nat=condition.sum(0)/(len(natural_files)*1.0)
	
	chist2=(historical_data.mean(1).mean(1)>obs.mean(0).mean(0)).sum(0)/(len(historical_files)*1.0)
	cnat2=(natural_data.mean(1).mean(1)>obs.mean(0).mean(0)).sum(0)/(len(natural_files)*1.0)

	chist3=(hist_raw.mean(1).mean(1)>obs.mean(0).mean(0)).sum(0)/(len(historical_files)*1.0)
	cnat3=(natural_raw.mean(1).mean(1)>obs.mean(0).mean(0)).sum(0)/(len(natural_files)*1.0)

# Old ensemble of historical simulations between 2000-2008
	cold=(old_data.mean(1).mean(1)>obs.mean(0).mean(0)).sum(0)/(len(oldfiles)*1.0)

# Region to look at germany:
	hmean=hist_raw[:,50:60,10:50].mean(1).mean(1)
	nmean=natural_raw[:,50:60,10:50].mean(1).mean(1)	
	obsmean=obs[50:60,10:50].mean(0).mean(0)
	chist4=(hmean>obsmean).sum(0)/(len(historical_files)*1.0)
	cnat4=(nmean>obsmean).sum(0)/(len(natural_files)*1.0)

	print 'return times: all forcings=',1/chist2,'natural forcings=',1/cnat2,'2000-2008=',1/cold
	print 'far',1-(cnat2/chist2)
	print 'change in risk', (chist2-cnat2)/cnat2
	
	print 'Without bias correction:'
	print 'return times: all forcings=',1/chist3,'natural forcings=',1/cnat3
	print 'far',1-(cnat3/chist3)
	print 'change in risk', (chist3-cnat3)/cnat3
	
	print 'Germany:'
	print 'return times: all forcings=',1/chist4,'natural forcings=',1/cnat4
	print 'far',1-(cnat4/chist4)
	print 'change in risk', (chist4-cnat4)/cnat4
	
	# Non vectorised way of calculating count 
#	count = np.zeros([119,122])
#	for i in range(122):
#		for j in range(119):
#			count[j,i] = (historical_data[:,j,i]>obs).sum()/(len(historical_files)*1.0)
	
	# Change in likelihood
	change=(count_hist-count_nat)/(count_nat+1e-5)
	
	# Anomaly between obs and ensemble mean of historical simulations
	anom=obs-(historical_data.mean(0))

	# Anomaly between obs and ensemble mean of nat simulations
	anom_nat=obs-(natural_data.mean(0))

	# FAR
	far=1.0-((count_nat)/(count_hist+1e-5))
	far=far*((count_hist!=0.0)*1.0) # Set to zero where count_hist==0

####################################################
	# Write out data to netcdf
#	create_netcdf_stripped(netcdf_file(rot_template),obs,'CRU_T2m_dec-nov_regional.nc')
#	create_netcdf_stripped(netcdf_file(rot_template),hist_raw,'hist_raw.nc')
#	create_netcdf_stripped(netcdf_file(rot_template),historical_data,'hist_corrected.nc')
#	create_netcdf_stripped(netcdf_file(rot_template),natural_raw,'natural_raw.nc')
#	create_netcdf_stripped(netcdf_file(rot_template),natural_data,'natural_corrected.nc')
#	create_netcdf_stripped(netcdf_file(rot_template),anom,'anom_corrected.nc')
#	create_netcdf_stripped(netcdf_file(rot_template),count_hist,'prob_historical_corrected.nc')
#	create_netcdf_stripped(netcdf_file(rot_template),count_nat,'prob_natural_corrected.nc')
#	create_netcdf_stripped(netcdf_file(rot_template),far,'far_corrected.nc')

######################################################
# Plot stuff

#First get rid of annoying Mac hidden files
	delfiles=glob.glob('._*.png')
	for f in delfiles:
		os.remove(f)
		
	plt.figure(1)
	plt.contourf((count_hist[::-1,:]),np.arange(0,1.05,.05),extend='both')
	plt.colorbar()
	plt.title('Probability of temperature being above 2014 temp = '+str(count_hist.mean()))
	plt.savefig('figure1.png')
		
	plt.figure(2)
	plt.contourf((count_nat[::-1,:]),np.arange(0,0.25,.01),extend='both')
	plt.colorbar()
	plt.title('Probability of temperature being above 2014 temp \nwith natural forcings only = '+str(count_nat.mean()))
	plt.savefig('figure2.png')
	
	plt.figure(3)
	plt.title('Observational dataset of Temperature for 2014 = '+str(obs.mean()))
	plt.contourf(obs[::-1,:],20)
	plt.colorbar()
	plt.savefig('figure3.png')
	
	plt.figure(4)
	plt.title('2014 Anomaly (Observations - Ensemble mean) = '+str(anom.mean()))
	plt.contourf(anom[::-1,:],20)
	plt.colorbar()
	plt.savefig('figure4.png')
	
	plt.figure(5)
	plt.title('Climatalogical bias used to bias correct model = '+str(bias.mean()))
	plt.contourf(bias[::-1,:],20)
	plt.colorbar()
	plt.savefig('figure5.png')
	
	plt.figure(6)
	plt.title('FAR = '+str(far.mean()))
	plt.contourf(far[::-1,:],np.arange(0,1.05,.05),extend='both')
	plt.colorbar()
	plt.savefig('figure6.png')	
	
	plt.figure(7)
	plt.title('Change in likelihood of temp greater than 2014 ='+str(change.mean()))
	plt.contourf(change[::-1,:],np.arange(0,22,2),extend='both')
	plt.colorbar()
	plt.savefig('figure7.png')	
	
	plt.figure(8)
	plt.title('2014 Anomaly (Observations - Nat Ensemble mean) = '+str(anom_nat.mean()))
	plt.contourf(anom_nat[::-1,:],20)
	plt.colorbar()
	plt.savefig('figure8.png')	
		
	#plt.show()
