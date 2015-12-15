#!/bin/env python2.7
import glob
import zipfile,gzip
import shutil,sys,os

sys.path.append('../cpdn_ancil_maker')
from checkdate_ancil_dump import checkdate

if __name__=="__main__":
	in_dir='/gpfs/projects/cpdn/storage/boinc/upload'
	in_umids='/home/cenv0437/cpdn_analysis/eof_analysis/eof_plots_hist/umids_chosen.txt'
	outdir='../../scratch/restarts_221/hist/'
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	for umid in open(in_umids):
		umid=umid.strip()
		tmp='/gpfs/projects/cpdn/storage/boinc/upload/hadam3p_eu/batch221/hadam3p_eu_'+umid+'_*/hadam3p_eu_'+umid+'*_12.zip'
		zipname=glob.glob(tmp)[0]
		print zipname
		f=zipfile.ZipFile(zipname)
		f.extractall(outdir)
		f.close()
		f_regional=outdir+'restart_eu_'+umid+'_2015_11_01.gz'
		f_atmos=outdir+'restart_atmos_'+umid+'_2015_11_01.gz'
		with open(outdir+'region_restart.day', 'rb') as f_in, gzip.open(f_regional, 'wb') as f_out:
			shutil.copyfileobj(f_in, f_out)
			if checkdate(f_regional):
				print 'created',f_regional
			else:
				print 'error, wrong date for',f_regional
		with open(outdir+'atmos_restart.day', 'rb') as f_in, gzip.open(f_atmos, 'wb') as f_out:
			shutil.copyfileobj(f_in, f_out)
			if checkdate(f_atmos):
				print 'created',f_atmos
			else:
				print 'error, wrong date for',f_atmos
		
		os.remove(outdir+'region_restart.day')
		os.remove(outdir+'atmos_restart.day')