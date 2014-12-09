import glob,os
#import numpy as np
#from scipy.io.netcdf import netcdf_file
#from netcdf_file import netcdf_file


def umid_to_int(umid_string):
# Test cases
#	list= ['0001','000a','000z','0010','00a0']
#	for umid in list:
#		print 'umid',umid, umid_to_int(umid)
# Expected output:
# umid 0001 1
# umid 000a 10
# umid 000z 35
# umid 0010 36
# umid 00a0 360
	id_int=0
	for i in range(4):
		c=umid_string[-i-1]
		if c.isdigit():
			id_int=id_int+int(c)*(36**i)
		else:
			id_int=id_int+(ord(c)-ord('a')+10)*36**i
	return id_int

def sort_tasknames(infiles):
	natural_files=[]
	historical_files=[]
	for f in glob.glob(infiles):
		umid=os.path.basename(f).split('_')[2]
		if umid_to_int(umid)>=umid_to_int('z2go'):
			natural_files.append(f)
		else:
			historical_files.append(f)
	return (historical_files,natural_files)

if __name__=='__main__':
	infiles='/home/ouce-cpdn/cenv0437/data/batch_100/hadam3p_eu_*_tasmean.nc'
	(historical_files,natural_files)=sort_tasknames(infiles)
	
	# Print out lists		
	print 'historical files:',len(historical_files)
	for f in historical_files:
		print f
		
	print 'natural files:',len(natural_files)
	for f in natural_files:
		print f