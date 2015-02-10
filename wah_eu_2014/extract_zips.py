import os,shutil,zipfile,glob
import numpy as np
#from scipy.io.netcdf import netcdf_file
from netcdf_file import netcdf_file



# Function to calculate yearly mean from monthly means
def yearly_mean(ncfiles):
	# Open first nc file to get dimensions
	tmp=netcdf_file(ncfiles[0],'r')
	meantemp=np.zeros(tmp.variables['field16'].shape)
	tmp.close()

	i=0
	# Loop over list of months 
	for filename in ncfiles:
		ncfile=netcdf_file(filename,'r')
		data=ncfile.variables['field16'][:]
		
		# Data consistency checking (for Temperature only!!!)
		if data.min()<170.0 or data.max()>350.0:
			raise Exception('Data outside reasonable bounds')
		meantemp=meantemp+data
		i=i+1
	if i!=12: # Need yearly mean
		raise Exception('Less than a years data is present')
	else:
		return meantemp/12
 
#function to return zip name 
def nc_filename(umid,year,zip_num,pp_code):
	mon_map={1:'dec',2:'jan',3:'feb',4:'mar',5:'apr',6:'may',7:'jun',8:'jul',9:'aug',10:'sep',11:'oct',12:'nov'}
	if zip_num>1: year=year+1
	decade=(int(year)-1800)/10
	year_code=str(year%10)
	if decade<10:   dec_code=str(decade) #decade is 0 to 9
	else:           dec_code=chr(decade-10+ord('a')) #decade needs to be converted to char
	
	return umid+pp_code+dec_code+year_code+mon_map[zip_num]+'.nc'
 
 # Extract one ncfile from each monthly zip and return list
def extract_zips(taskdir,outpath):
#	file_map={1:'l3dec',2:'l4jan',3:'l4feb',4:'l4mar',5:'l4apr',6:'l4may',7:'l4jun',8:'l4jul',9:'l4aug',10:'l4sep',11:'l4oct',12:'l4nov'}
#	mon_map={'jan':1,'feb':2,'mar':3,'apr':4,'may':5,'jun':6,'jul':7,'aug':8,'sep':9,'oct':10,'nov':11,'dec':12}	

	# Make tmp folder
	tmpdir=os.path.join(outpath,'tmp')
	if not os.path.exists(tmpdir):
		os.makedirs(tmpdir)		
	taskname = taskdir.split('/')[-1]
	[umid,year]=taskname.split('_')[2:4]
	extracted=[]
	for i in range(1,13):
		zipname=os.path.join(taskdir,taskname+'_'+str(i)+'.zip')
		if os.path.exists(zipname):
			zip=zipfile.ZipFile(zipname,'r')
			# Just extract the monthly means from regional model
#			ncfile=umid+'ga.pe'+file_map[i]+'.nc'
			ncfile=nc_filename(umid,int(year),i,'ga.pe')
			zip.extract(ncfile,tmpdir) 
			extracted.append(os.path.join(tmpdir,ncfile))
		# Only return if all of the zips are present
		else: raise Exception('Not all zips are present')
	return extracted
	
def create_netcdf(template,data,outname):
	# create outfile object
	outfile=netcdf_file(outname,'w')
	
	# Create dimensions copied from template file
	temp=template.variables['field16']
	for dim in temp.dimensions:
		if dim=='time1': leng=0
		else: leng=int(template.dimensions[dim])
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
	outfile.variables['global_latitude0'][:]=template.variables['global_latitude0'][:]
	for att in template.variables['global_latitude0']._attributes:
			outfile.variables['global_latitude0'].__setattr__(att,template.variables['global_latitude0'].__getattribute__(att))
			
	# Global longitude
	outfile.createVariable('global_longitude0','f',template.variables['global_longitude0'].dimensions)
	outfile.variables['global_longitude0'][:]=template.variables['global_longitude0'][:]
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
	incoming='/gpfs/projects/cpdn/storage/boinc/upload/hadam3p_eu/batch100/hadam3p_eu_*'
	outpath='/gpfs/projects/cpdn/scratch/cenv0437/batch_100'
#	incoming='/gpfs/projects/cpdn/storage/boinc/upload/hadam3p_eu/batch43/hadam3p_eu_????_2*'
#	outpath='/gpfs/projects/cpdn/scratch/cenv0437/batch_43'
	# Make sure output exists
	if not os.path.exists(outpath):
		os.makedirs(outpath)
		
	tasks=glob.glob(incoming)
	for taskdir in tasks:
		taskname = taskdir.split('/')[-1]
		outfile=os.path.join(outpath,taskname+'_tasmean.nc')
		print outfile
		
		if not os.path.exists(outfile):
			try:
				ncfiles=extract_zips(taskdir,outpath)
				meantemp=yearly_mean(ncfiles)
				tmp=netcdf_file(ncfiles[0],'r') # use as a template for output file
				create_netcdf(tmp,meantemp,outfile)
				# Remove temporary ncfiles
				for f in ncfiles:
					os.remove(f)
			# Copy data to arcus
#				os.system('scp '+outfile+' cenv0437@arcus.oerc.ox.ac.uk:data/batch_100')
			except Exception as e:
				print 'Error for',taskname,e
		else:
			print taskname, 'already converted to yearly mean'
