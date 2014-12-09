import os,shutil,zipfile,glob
import numpy as np
#from scipy.io.netcdf import netcdf_file
from netcdf_file import netcdf_file
file_map={1:'l3dec',2:'l4jan',3:'l4feb',4:'l4mar',5:'l4apr',6:'l4may',7:'l4jun',8:'l4jul',9:'l4aug',10:'l4sep',11:'l4oct',12:'l4nov'}

def yearly_mean(ncfiles,outfile):
	print ncfiles[0]
	# Open first nc file to get dimensions
	tmp=netcdf_file(ncfiles[0],'r')
	nt,nl,nlat,nlon=tmp.variables['field16'].shape
	meantemp=np.zeros([1,1,nlat,nlon])
	i=0
	# Loop over list of tasks and copy to scratch space
	for filename in ncfiles:
		ncfile=netcdf_file(filename,'r')
		data=ncfile.variables['field16'][:]
		if data.min()<170.0 or data.max()>350.0:
			return False
		meantemp=meantemp+data
		i=i+1
	if i!=12: # Need yearly mean
		return False
	else:
		meantemp=meantemp/12
		create_netcdf(tmp,meantemp,outfile)
 
def extract_zips(taskdir,outpath):
	# Make tmp folder
	tmpdir=os.path.join(outpath,'tmp')
	if not os.path.exists(tmpdir):
		os.makedirs(tmpdir)		
	taskname = taskdir.split('/')[-1]
	umid=taskname.split('_')[2]
	extracted=[]
	for i in range(1,13):
		zipname=os.path.join(taskdir,taskname+'_'+str(i)+'.zip')
		if os.path.exists(zipname):
			zip=zipfile.ZipFile(zipname,'r')
			# Just extract the monthly means from regional model
			ncfile=umid+'ga.pe'+file_map[i]+'.nc'
			zip.extract(ncfile,tmpdir) 
			print zipname,i
			extracted.append(os.path.join(tmpdir,ncfile))
		# Only return if all of the zips are present
		else: return []
	return extracted
	
def create_netcdf(template,data,outname):
	# create outfile object
	outfile=netcdf_file(outname,'w')
	
	# Create dimimport os,shutil,zipfile,glob
import numpy as np
#from scipy.io.netcdf import netcdf_file
from netcdf_file import netcdf_file
file_map={1:'l3dec',2:'l4jan',3:'l4feb',4:'l4mar',5:'l4apr',6:'l4may',7:'l4jun',8:'l4jul',9:'l4aug',10:'l4sep',11:'l4oct',12:'l4nov'}

def yearly_mean(ncfiles,outfile):
	# Open first nc file to get dimensions
	tmp=netcdf_file(ncfiles[0],'r')
	nt,nl,nlat,nlon=tmp.variables['field16'].shape
	meantemp=np.zeros([1,1,nlat,nlon])
	i=0
	# Loop over list of tasks and copy to scratch space
	for filename in ncfiles:
		ncfile=netcdf_file(filename,'r')
		data=ncfile.variables['field16'][:]
		if data.min()<170.0 or data.max()>350.0:
			return False
		meantemp=meantemp+data
		i=i+1
	if i!=12: # Need yearly mean
		return False
	else:
		meantemp=meantemp/12
		create_netcdf(tmp,meantemp,outfile)
		print outfile
 
def extract_zips(taskdir,outpath):
	# Make tmp folder
	tmpdir=os.path.join(outpath,'tmp')
	if not os.path.exists(tmpdir):
		os.makedirs(tmpdir)		
	taskname = taskdir.split('/')[-1]
	umid=taskname.split('_')[2]
	extracted=[]
	for i in range(1,13):
		zipname=os.path.join(taskdir,taskname+'_'+str(i)+'.zip')
		ncfile=umid+'ga.pe'+file_map[i]+'.nc'
		if os.path.exists(zipname):
			zip=zipfile.ZipFile(zipname,'r')
			# Just extract the monthly means from regional model
			zip.extract(ncfile,tmpdir)
			extracted.append(os.path.join(tmpdir,ncfile))
		# Only return if all of the zips are present
		else: return []
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
	# Make sure output exists
	if not os.path.exists(outpath):
		os.makedirs(outpath)
		
	tasks=glob.glob(incoming)
	for taskdir in tasks:
		taskname = taskdir.split('/')[-1]
		ncfiles=extract_zips(taskdir,outpath)
		outfile=os.path.join(outpath,taskname+'_tasmean.nc')
		yearly_mean(ncfiles,outfile)