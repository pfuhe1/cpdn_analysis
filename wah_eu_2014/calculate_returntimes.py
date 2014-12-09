import glob,os
from sort_umids import sort_tasknames
import numpy as np
from netcdf_file import netcdf_file
import matplotlib.pyplot as plt
from extract_zips import create_netcdf
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
#			print j,i
			index_lat=int(((g_lat[j,i]+90)/1.25)+0.5)%145
			index_lon=int((g_lon[j,i]/1.875)+0.5)%192
#			print index_lat,index_lon
			dataout[j,i]=N96_data[index_lat,index_lon]
	dataout=np.ma.masked_values(dataout,0.0)
	return dataout			


# Create netcdf file, removing the 'sponge layer' at the edge of the domain
def create_netcdf_stripped(template,data,outname,lat_start=4,lat_end=7,lon_start=4,lon_end=4):
	# create outfile object
	outfile=netcdf_file(outname,'w')
	
	# Create dimensions copied from template file
	temp=template.variables['tas']
	for dim in temp.dimensions:
		if dim=='time1': 
			leng=0
			outfile.createDimension(dim,leng)
			outfile.createVariable(dim,'f',(dim,))
			outfile.variables[dim][:]=template.variables[dim][:]
		elif dim[:3]=='lat': 
			leng=int(template.dimensions[dim])-lat_start-lat_end
			outfile.createDimension(dim,leng)
			outfile.createVariable(dim,'f',(dim,))
			outfile.variables[dim][:]=template.variables[dim][lat_start:-lat_end]
		elif dim[:3]=='lon' : 
			leng=int(template.dimensions[dim])-lon_start-lon_end
			outfile.createDimension(dim,leng)
			outfile.createVariable(dim,'f',(dim,))
			outfile.variables[dim][:]=template.variables[dim][lon_start:-lon_end]
		else: 
			leng=int(template.dimensions[dim])
			outfile.createDimension(dim,leng)
			outfile.createVariable(dim,'f',(dim,))
			outfile.variables[dim][:]=template.variables[dim][:]
			
		for att in template.variables[dim]._attributes:
			outfile.variables[dim].__setattr__(att,template.variables[dim].__getattribute__(att))
	
	# Rotated pole variable
	outfile.createVariable('rotated_pole0','c',())
	for att in template.variables['rotated_pole0']._attributes:
			outfile.variables['rotated_pole0'].__setattr__(att,template.variables['rotated_pole0'].__getattribute__(att))
	
	# Global latitude
	outfile.createVariable('global_latitude0','f',template.variables['global_latitude0'].dimensions)
	outfile.variables['global_latitude0'][:]=template.variables['global_latitude0'][lat_start:-lat_end,lon_start-lon_end]
	for att in template.variables['global_latitude0']._attributes:
			outfile.variables['global_latitude0'].__setattr__(att,template.variables['global_latitude0'].__getattribute__(att))
			
	# Global longitude
	outfile.createVariable('global_longitude0','f',template.variables['global_longitude0'].dimensions)
	outfile.variables['global_longitude0'][:]=template.variables['global_longitude0'][lat_start:-lat_end,lon_start-lon_end]
	for att in template.variables['global_longitude0']._attributes:
			outfile.variables['global_longitude0'].__setattr__(att,template.variables['global_longitude0'].__getattribute__(att))
	
	# Create data variable (named tas)
	outfile.createVariable('tas','f',temp.dimensions)
	outfile.variables['tas'][:]=data
	for att in temp._attributes:
		outfile.variables['tas'].__setattr__(att,temp.__getattribute__(att))
	
	outfile.flush()
	outfile.close()


if __name__=='__main__':
	infiles='/ouce-home/staff/cenv0437/data/batch_100/hadam3p_eu_*_tasmean.nc'
	bias_files='/ouce-home/staff/cenv0270/CPDN/Europe_2014/observations/CRUT4/EU_???_temp-cruts_0.44_mean.nc'
	bias = ave_bias(bias_files)
	(historical_files,natural_files)=sort_tasknames(infiles)
	
	# Load historical data into single array and bias correct
	print len(historical_files)
	hist_raw=load_ensemble(historical_files)
	historical_data=hist_raw-bias
	
	# Load natural simulation data into single array and bias correct
	print len(natural_files)
	natural_data=load_ensemble(natural_files)-bias
	
	
	# TODO Get real OBS!!!
	obs_N96='/ouce-home/staff/cenv0270/CPDN/Europe_2014/observations/final/CRU_T2m_dec-nov_N96_interp_lsm.nc'
	rot_template='/ouce-home/staff/cenv0437/data/batch_100/hadam3p_eu_z4ao_2013_1_009238311_0_tasmean.nc'
	obs=remap_N96_to_rotated(obs_N96,rot_template)
	
	# Write out calculated obs on rotated grid
	create_netcdf_stripped(netcdf_file(rot_template),obs,'CRU_T2m_dec-nov_regional.nc')
	
	# Count number of model runs which are above the obs
	condition=historical_data>obs
	count_hist=condition.sum(0)/(len(historical_files)*1.0)
	
	condition=natural_data>obs
	count_nat=condition.sum(0)/(len(natural_files)*1.0)

# Non vectorised way of calculating count 
#	count = np.zeros([119,122])
#	for i in range(122):
#		for j in range(119):
#			count[j,i] = (historical_data[:,j,i]>obs).sum()/(len(historical_files)*1.0)

# Plot stuff
	plt.figure(1)
#	plt.contourf((count_hist[::-1,:]))
	plt.contourf((historical_data.mean(0)-obs)[::-1,:],20)
	plt.colorbar()
	
	plt.figure(2)
#	plt.contourf((count_nat[::-1,:]))
	plt.contourf(bias[::-1,:],20)
	plt.colorbar()
	
	plt.figure(3)
	plt.contourf(obs[::-1,:],20)
	plt.colorbar()
	plt.show()	
	
#	far=1.0-(count_nat/count_hist)
#	plt.figure(3)
#	plt.contourf(far[::-1,:])
#	plt.colorbar()
#	plt.show()
	