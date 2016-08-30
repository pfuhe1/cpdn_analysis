import glob,os
import numpy as np
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
	return int(umid_string,36) # Use inbuilt python function

def sort_tasknames(infiles):
	#If necessary glob infiles from a string into a list
	if isinstance(infiles,basestring):infiles=glob.glob(infiles)
	natural_files=[]
	historical_files=[]
	for f in infiles:
		umid=os.path.basename(f).split('_')[2]
		if umid_to_int(umid)>=umid_to_int('z2go'):
			natural_files.append(f)
		else:
			historical_files.append(f)
	return (historical_files,natural_files)

def choose_tasknames(infiles,startid,endid):
        choosefiles=[]
        #If necessary glob infiles from a string into a list
        if isinstance(infiles,basestring):infiles=glob.glob(infiles)
        for f in infiles:
                umid=os.path.basename(f).split('_')[2]
                if umid_to_int(startid)<=umid_to_int(umid) and umid_to_int(umid)<=umid_to_int(endid):
                        choosefiles.append(f)
        return choosefiles
        
        
def choose_mask(infiles,startid,endid):
	#If necessary glob infiles from a string into a list
	if isinstance(infiles,basestring):infiles=glob.glob(infiles)
	choosefiles=np.zeros(len(infiles),dtype=bool)
	for i,f in enumerate(infiles):
		umid=os.path.basename(f).split('_')[2]
		if umid_to_int(startid)<=umid_to_int(umid) and umid_to_int(umid)<=umid_to_int(endid):
			choosefiles[i]=True
        return choosefiles
        
        
############
#The following functions assume that the umid is the second element rather than the third


def choose_tasknames1(infiles,startid,endid):
        choosefiles=[]
        #If necessary glob infiles from a string into a list
        if isinstance(infiles,basestring):infiles=glob.glob(infiles)
        for f in infiles:
                umid=os.path.basename(f).split('_')[1]
                if umid_to_int(startid)<=umid_to_int(umid) and umid_to_int(umid)<=umid_to_int(endid):
                        choosefiles.append(f)
        return choosefiles
        
        
def choose_mask1(infiles,startid,endid):
	#If necessary glob infiles from a string into a list
	if isinstance(infiles,basestring):infiles=glob.glob(infiles)
	choosefiles=np.zeros(len(infiles),dtype=bool)
	for i,f in enumerate(infiles):
		umid=os.path.basename(f).split('_')[1]
		if umid_to_int(startid)<=umid_to_int(umid) and umid_to_int(umid)<=umid_to_int(endid):
			choosefiles[i]=True
        return choosefiles
        

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
