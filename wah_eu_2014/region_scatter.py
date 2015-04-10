#!/usr/local/bin/python2.7
# NOTE python2.7 must be used for this script

import glob,os
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle
from mpl_toolkits.basemap import Basemap,cm
from region_returntimes import get_region_returntimes2,remap_p5deg_to_rotated
from scipy.io.netcdf import netcdf_file

pkl_dir='/gpfs/projects/cpdn/scratch/cenv0437/pickle_data/'

if __name__=='__main__':

# Load model data from pickle files
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

	obs_p5deg='/home/cenv0437/scratch/data_from_ouce/CRU_TS_Absolute_plus_BEST_anomaly_201312_201412.nc'
	obsname='BEST'
	#obs_p5deg='/home/cenv0437/scratch/data_from_ouce/CRU_TS_Absolute_plus_CRU4_anomaly_201312_201412.nc'
#	obsname='CRU4'

	#obs_p5deg='/home/cenv0437/scratch/data_from_ouce/CRU_TS_Absolute_plus_GISS_anomaly_201312_201412.nc'
#	obsname='GISS'




	calculate_regions=False

	if calculate_regions:
	# Cacluate return times for many regions
		print historical_data.shape	
		nens,nlat,nlon=historical_data.shape
		hist=np.zeros([30,nlat/5+1,nlon/5+1])
		nat=np.zeros([30,nlat/5+1,nlon/5+1])
		clim_hist=np.zeros([30,nlat/5+1,nlon/5+1])
		clim_nat=np.zeros([30,nlat/5+1,nlon/5+1])
		lp=np.zeros([30,nlat/5+1,nlon/5+1])	
		hist_var=np.zeros([30,nlat/5+1,nlon/5+1]) #variance of hist ensemble 
		nat_var=np.zeros([30,nlat/5+1,nlon/5+1]) # variance of nat ensemble
		hist_cor=np.zeros([30,nlat/5+1,nlon/5+1])
		nat_cor=np.zeros([30,nlat/5+1,nlon/5+1])
		
		#obs data
		# Obs using higher resolution cru.ts climatology
#		obs_p5deg='/home/cenv0437/scratch/data_from_ouce/CRU_TS_dec13-nov14_crut4anomalies.nc'


		# Regrid obs to rotated regional grid
		rot_template='/home/cenv0437/scratch/data_from_ouce/hadam3p_eu_z4ao_2013_1_009238311_0_tasmean.nc'
		f_rot=netcdf_file(rot_template,'r')
		lat_coord=f_rot.variables['global_latitude0'][4:-7,4:-4]
		lon_coord=f_rot.variables['global_longitude0'][4:-7,4:-4]
		f_rot.close()	
		obs=remap_p5deg_to_rotated(obs_p5deg,rot_template)[:,:]
		
		# Some random regions
		for n in range(30): #width/height of domain
			print n
			for i in range((nlon-n)/5):
				for j in range((nlat-n)/5):
					try:
						vals=get_region_returntimes2(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,i*5,i*5+(n*3+2),j*5,j*5+(n*3+2))
						vin=historical_data[j*5:j*5+(n*3+2),i*5:i*5+(n*3+2)]
						v1=historical_data[:,j*4,i*5:i*5+(n*3+2)]
						v2=historical_data[:,j*5+(n*3+3),i*5:i*5+(n*3+2)]
						v3=historical_data[:,j*5:j*5+(n*3+2),i*4]
						v4=historical_data[:,j*5:j*5+(n*3+2),i*5+(n*3+3)]
					except Exception,e:
						print e
						continue
					if vals[4]>0:
						hist[n,j,i]=vals[0]
						nat[n,j,i]=vals[1]
						clim_hist[n,j,i]=vals[2]
						clim_nat[n,j,i]=vals[3]
						lp[n,j,i]=vals[4]
						hist_var[n,j,i]=vals[5]
						nat_var[n,j,i]=vals[6]
						vout=np.ma.concatenate((v1,v2,v3,v4),axis=1) # chose points around boundary of region
						hist_cor[n,j,i]=np.corrcoef(v1.mean(1),vout.mean(1))[0,1]
						
	#Write out these values to a pickle file
		data2=open(pkl_dir+'data2'+obsname+'.pkl','wb')
		pickle.dump(hist,data2,-1)
		pickle.dump(nat,data2,-1)
		pickle.dump(clim_hist,data2,-1)
		pickle.dump(clim_nat,data2,-1)
		pickle.dump(lp,data2,-1)
		pickle.dump(hist_var,data2,-1)
		pickle.dump(nat_var,data2,-1)
		pickle.dump(hist_cor,data2,-1)
		data2.close()
	else:
	# read existing values from pickle	
		data2=open(pkl_dir+'data2'+obsname+'.pkl','rb')
		hist=pickle.load(data2)
		nat=pickle.load(data2)
		clim_hist=pickle.load(data2)
		clim_nat=pickle.load(data2)
		lp=pickle.load(data2)
		hist_var=pickle.load(data2)
		nat_var=pickle.load(data2)
		hist_cor=pickle.load(data2)
		data2.close()

###############################################
# Plotting	

#First get rid of annoying Mac hidden files
	delfiles=glob.glob('scatter_figs/._*.png')
	for f in delfiles:
		os.remove(f)
		
# Some data manipulation...
	lp=lp.flatten()
	hist=hist.flatten()
	nat=nat.flatten()
	clim_hist=clim_hist.flatten()
	clim_nat=clim_nat.flatten()
	

	
	condition=np.logical_or(np.logical_not(np.isfinite(hist)),hist>1000)
	hist=np.ma.masked_where(condition,hist)
	

#	condition=np.logical_or(np.logical_not(np.isfinite(nat)),nat>4000)
	print (~np.isfinite(nat)).sum()
	nat=np.ma.masked_where(nat==np.NaN,nat) # mask out any possible nan's
	print (~np.isfinite(nat)).sum()
	nat[~np.isfinite(nat)]=natural_data.shape[0]	
	print (~np.isfinite(nat)).sum()
	print nat.min(),nat.max()
	condition=np.logical_or(np.logical_not(np.isfinite(clim_hist)),clim_hist>4000)
	clim_hist=np.ma.masked_where(condition,clim_hist)
	
#	condition=np.logical_or(np.logical_not(np.isfinite(clim_nat)),clim_nat>4000)

	print (~np.isfinite(clim_nat)).sum()
	clim_nat=np.ma.masked_where(clim_nat==np.NaN,clim_nat) # mask out any possible nan's
	print (~np.isfinite(clim_nat)).sum()
	clim_nat[~np.isfinite(clim_nat)]=clim_nat_data.shape[0]
	print (~np.isfinite(clim_nat)).sum()
	print clim_nat.min(),clim_nat.max()


	print '###########################################################\n'

	plt.figure(10)
	plt.plot(lp,nat.flatten(),'.',label='nat')
#	plt.plot(lp,clim_hist.flatten(),'g.',label='clim_hist')
#	plt.plot(lp,hist.flatten(),'.b',label='hist')
	plt.plot(lp,clim_nat.flatten(),'.',label='clim_nat')
	plt.title('return time vs size of region')
	plt.legend()
	plt.ylim([0,3000])
	plt.xlim([0,5000])
	plt.xlabel('Size of region (grid-points on land)')
	plt.ylabel('Return time (years)')
	plt.savefig('figure_nat.png')
	
	plt.figure(2)
	plt.plot(lp,(1-clim_hist/nat).flatten(),'.c',label='Clim hist,nat2014 vs 2014 obs')
	plt.plot(lp,(1-clim_hist/clim_nat).flatten(),'.g',label='clim hist, clim nat, vs 2014 obs')
	plt.plot(lp,(1-hist/nat).flatten(),'.b',label='2014')
#	plt.plot(curve_fit(np.exp2,lp,(1-hist/nat)),'-',label='2014 fit')
	plt.title('FAR vs size of region')
	plt.xlabel('Size of region (grid-points on land)')
	plt.ylabel('FAR')
	plt.ylim([0.6,1])
	plt.xlim([0,5000])
	plt.legend(loc='best')
	plt.savefig('scatter_figs/figure_far_'+obsname+'.png')

	plt.clf()

# Histogram bin sizes
	bin_width=40
	centre_size=np.arange(5+bin_width/2,5000,bin_width)
	hist_bins=np.arange(0.6,1.02,.02)
	
# Calculate histogram:
	histo2d=np.zeros([len(hist_bins)-1,len(centre_size)])
	for i,size in enumerate(centre_size):
		condition=np.logical_or(lp<size-bin_width/2,lp>=size+bin_width/2)
		tmp=np.ma.masked_where(condition,1-hist/nat,copy=True)
		histo,bins=np.histogram(tmp[~tmp.mask],hist_bins)
		centre_far=bins[:-1]+(bins[1:]-bins[:-1])/2.
		histo2d[:,i]=histo/((~condition).sum()*1.0) # Normalise histogram for each region size bin
	plt.contourf(centre_size,centre_far,histo2d,100)#,np.arange(0,1.01,.01))
	plt.title('Proportion of subregions: FAR for 2014 vs '+obsname)
	plt.xlabel('Size of region (grid-points on land)')
	plt.ylabel('FAR')
	plt.colorbar()
	plt.savefig('scatter_figs/figure_farhist_'+obsname+'1.png')

	plt.clf()
# Calculate histogram:

	histo2d=np.zeros([len(hist_bins)-1,len(centre_size)])
	for i,size in enumerate(centre_size):
		condition=np.logical_or(lp<size-bin_width/2,lp>=size+bin_width/2)
		tmp=np.ma.masked_where(condition,1-clim_hist/clim_nat,copy=True)
		histo,bins=np.histogram(tmp[~tmp.mask],hist_bins)
		centre_far=bins[:-1]+(bins[1:]-bins[:-1])/2.
		histo2d[:,i]=histo/((~condition).sum()*1.0) # Normalise histogram for each region size bin
	plt.contourf(centre_size,centre_far,histo2d,100)#,np.arange(0,1.01,.01))
	plt.title('Proportion of subregions: FAR for Climatology vs '+obsname)
	plt.xlabel('Size of region (grid-points on land)')
	plt.ylabel('FAR')
	plt.colorbar()
	plt.savefig('scatter_figs/figure_farhist_'+obsname+'2.png')
	
	plt.clf()
# Calculate histogram:

	histo2d=np.zeros([len(hist_bins)-1,len(centre_size)])
	for i,size in enumerate(centre_size):
		condition=np.logical_or(lp<size-bin_width/2,lp>=size+bin_width/2)
		tmp=np.ma.masked_where(condition,1-clim_hist/nat,copy=True)
		histo,bins=np.histogram(tmp[~tmp.mask],hist_bins)
		centre_far=bins[:-1]+(bins[1:]-bins[:-1])/2.
		histo2d[:,i]=histo/((~condition).sum()*1.0) # Normalise histogram for each region size bin
	plt.contourf(centre_size,centre_far,histo2d,100)#,np.arange(0,1.01,.01))
	plt.title('Proportion of subregions: FAR for Climhist/Nat2014 vs '+obsname)
	plt.xlabel('Size of region (grid-points on land)')
	plt.ylabel('FAR')
	plt.colorbar()
	plt.savefig('scatter_figs/figure_farhist_'+obsname+'3.png')



	
	plt.clf()
#	plt.scatter(lp,(1-hist/nat).flatten(),c=np.log(nat.max()/nat),label='2014')
	plt.scatter(lp,(1-hist/nat).flatten(),c=hist_cor)
	plt.colorbar()
#	plt.plot(curve_fit(np.exp2,lp,(1-hist/nat)),'-',label='2014 fit')
	plt.title('FAR vs size of region')
	plt.xlabel('Size of region (grid-points on land)')
	plt.ylabel('FAR ')
	plt.ylim([0.6,1])
	plt.xlim([0,5000])
#	plt.legend(loc='best')
	plt.savefig('scatter_figs/figure_far_'+obsname+'1.png')
	

	
	
	plt.clf()
	plt.scatter(lp,(1-clim_hist/clim_nat).flatten(),c=np.log(clim_nat.max()/clim_nat),label='clim hist, clim nat, vs 2014 obs')
#	plt.scatter(lp,(1-clim_hist/clim_nat).flatten(),c=nat_var)
	plt.colorbar()
#	plt.plot(curve_fit(np.exp2,lp,(1-hist/nat)),'-',label='2014 fit')
	plt.title('FAR vs size of region: clim hist v clim nat against '+obsname)
	plt.xlabel('Size of region (grid-points on land)')
	plt.ylabel('FAR ')
	plt.ylim([.6,1])
	plt.xlim([0,5000])
#	plt.legend(loc='best')
	plt.savefig('scatter_figs/figure_far_'+obsname+'2.png')
	
	plt.clf()
	plt.scatter(lp,(1-clim_hist/nat).flatten(),c=np.log(nat.max()/nat),label='Clim hist,nat2014 vs 2014 obs')
	plt.colorbar()
#	plt.plot(curve_fit(np.exp2,lp,(1-hist/nat)),'-',label='2014 fit')
	plt.title('FAR vs size of region: Clim hist,nat2014 vs against '+obsname)
	plt.xlabel('Size of region (grid-points on land)')
	plt.ylabel('FAR ')
	plt.ylim([0.6,1])
	plt.xlim([0,5000])
#	plt.legend(loc='best')
	plt.savefig('scatter_figs/figure_far_'+obsname+'3.png')

	plt.figure(4)
	plt.semilogy(lp,1-(1-clim_hist/nat).flatten(),'.c',label='Clim hist,nat2014 vs 2014 obs')
	plt.semilogy(lp,1-(1-clim_hist/clim_nat).flatten(),'.g',label='clim hist, clim nat, vs 2014 obs')
	plt.semilogy(lp,1-(1-hist/nat).flatten(),'.b',label='2014')
#	plt.plot(curve_fit(np.exp2,lp,(1-hist/nat)),'-',label='2014 fit')
	plt.title('FAR vs size of region')
	plt.xlabel('Size of region (grid-points on land)')
	plt.ylabel('1-FAR ')
	plt.ylim([0,1])
	plt.xlim([0,5000])
	plt.legend(loc='best')
	plt.savefig('scatter_figs/figure_farlog_'+obsname+'.png')
	
	plt.clf()
	ax = plt.gca()
	ax.scatter(lp,1-(1-hist/nat).flatten(),c=np.log(nat.max()/nat),label='2014')
	ax.set_yscale('log')
#	plt.plot(curve_fit(np.exp2,lp,(1-hist/nat)),'-',label='2014 fit')
	plt.title('FAR vs size of region')
	plt.xlabel('Size of region (grid-points on land)')
	plt.ylabel('1-FAR ')
	ax.set_ylim(0.0001,1)
	ax.set_xlim(0,5000)
	plt.legend(loc='best')
	plt.savefig('scatter_figs/figure_farlog_'+obsname+'1.png')
	
	plt.clf()
	ax = plt.gca()
	plt.scatter(lp,1-(1-clim_hist/clim_nat).flatten(),c=np.log(clim_nat.max()/clim_nat),label='clim hist, clim nat, vs 2014 obs')
	ax.set_yscale('log')
#	plt.plot(curve_fit(np.exp2,lp,(1-hist/nat)),'-',label='2014 fit')
	plt.title('FAR vs size of region')
	plt.xlabel('Size of region (grid-points on land)')
	plt.ylabel('1-FAR ')
	plt.ylim([0.0001,1])
	plt.xlim([0,5000])
	plt.legend(loc='best')
	plt.savefig('scatter_figs/figure_farlog_'+obsname+'2.png')
	
	plt.clf()
	ax = plt.gca()
	plt.scatter(lp,1-(1-clim_hist/nat).flatten(),c=np.log(nat.max()/nat),label='Clim hist,nat2014 vs 2014 obs')
	ax.set_yscale('log')
#	plt.plot(curve_fit(np.exp2,lp,(1-hist/nat)),'-',label='2014 fit')
	plt.title('FAR vs size of region')
	plt.xlabel('Size of region (grid-points on land)')
	plt.ylabel('1-FAR ')
	plt.ylim([0.0001,1])
	plt.xlim([0,5000])
	plt.legend(loc='best')
	plt.savefig('scatter_figs/figure_farlog_'+obsname+'3.png')
	
	plt.clf()
#	plt.scatter(lp,(1-hist/nat).flatten(),c=np.log(nat.max()/nat),label='2014')
	plt.scatter(lp,1/nat,c=nat_var)
	plt.colorbar()
#	plt.plot(curve_fit(np.exp2,lp,(1-hist/nat)),'-',label='2014 fit')
	plt.title('P0 vs size of region')
	plt.xlabel('Size of region (grid-points on land)')
	plt.ylabel('FAR ')
#	plt.ylim([0.6,1])
	plt.xlim([0,5000])
#	plt.legend(loc='best')
	plt.savefig('scatter_figs/figure_p0_'+obsname+'.png')
	

	
	
	plt.clf()
	plt.scatter(lp,1/hist.flatten(),c=hist_var,label='clim hist, clim nat, vs 2014 obs')
#	plt.scatter(lp,(1-clim_hist/clim_nat).flatten(),c=nat_var)
	plt.colorbar()
#	plt.plot(curve_fit(np.exp2,lp,(1-hist/nat)),'-',label='2014 fit')
	plt.title('P1 '+obsname)
	plt.xlabel('Size of region (grid-points on land)')
	plt.ylabel('return time ')
#	plt.ylim([.6,1])
	plt.xlim([0,5000])
#	plt.legend(loc='best')
	plt.savefig('scatter_figs/figure_p1_'+obsname+'.png')
	
	plt.clf()
	
	plt.figure(3)
	hist=np.ma.masked_where(np.logical_or(lp<300,lp>800),hist)
	print 'non masked values',hist.shape[0] -hist.mask.sum()
	histo,bins=np.histogram((1-hist/nat).flatten(),[0,.5,.6,.7,.8,.9,.91,.92,.93,.94,.95,.96,.97,.98,.99,1.0,1.1])
	print histo,bins
	centre=bins[:-1]+(bins[1:]-bins[:-1])/2.
	plt.plot(centre[:-1],histo[:-1])
	plt.title("Histogram of FAR values")
#	plt.savefig('histogram.png')
		
	plt.show()
