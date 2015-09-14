import glob
import zipfile,gzip
import shutil,sys,os

sys.path.append('../cpdn_ancil_maker')
from checkdate_ancil_dump import checkdate

if __name__=="__main__":
	in_dir='/gpfs/projects/cpdn/storage/boinc/upload'
	in_umids='/home/cenv0437/cpdn_analysis/eof_analysis/eof_plots_nat/umids_chosen.txt'
	outdir='../../scratch/restarts_181/nat/'
	for umid in open(in_umids):
		umid=umid.strip()
		tmp='/gpfs/projects/cpdn/storage/boinc/upload/hadam3p_pnw/batch181/sub-batch*/hadam3p_pnw_'+umid+'_*/hadam3p_pnw_'+umid+'*_19.zip'
		zipname=glob.glob(tmp)[0]
		print zipname
		f=zipfile.ZipFile(zipname)
		f.extractall(outdir)
		f.close()
		f_regional=outdir+'restart_pnw_'+umid+'_2015_06_01.gz'
		f_atmos=outdir+'restart_atmos_'+umid+'_2015_06_01.gz'
		with open('region_restart.day', 'rb') as f_in, gzip.open(f_regional, 'wb') as f_out:
			shutil.copyfileobj(f_in, f_out)
			if checkdate(f_regional):
				print 'created',f_regional
			else:
				print 'error, wrong date for',f_regional
		with open('atmos_restart.day', 'rb') as f_in, gzip.open(f_atmos, 'wb') as f_out:
			shutil.copyfileobj(f_in, f_out)
			if checkdate(f_atmos):
				print 'created',f_atmos
			else:
				print 'error, wrong date for',f_atmos
		