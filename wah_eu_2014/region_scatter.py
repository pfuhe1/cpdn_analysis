#!/usr/local/bin/python2.7
# NOTE python2.7 must be used for this script

import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle
from mpl_toolkits.basemap import Basemap,cm
from region_returntimes import get_region_returntimes



if __name__=='__main__':

# Load model data from pickle files
	fdata=open('historical_data.pkl','rb')
	historical_data=pickle.load(fdata)
	fdata.close()
	fdata=open('natural_data.pkl','rb')
	natural_data=pickle.load(fdata)
	fdata.close()
	fdata=open('clim_hist_data.pkl','rb')
	clim_hist_data=pickle.load(fdata)
	fdata.close()
	fdata=open('clim_nat_data.pkl','rb')
	clim_nat_data=pickle.load(fdata)
	fdata.close()
	print 'loaded data from pkl files'

	calculate_regions=False

	if calculate_regions:
	# Cacluate return times for many regions	
		nlat,nlon=lon_coord.shape
		hist=np.zeros([30,nlat/5+1,nlon/5+1])
		nat=np.zeros([30,nlat/5+1,nlon/5+1])
		clim_hist=np.zeros([30,nlat/5+1,nlon/5+1])
		clim_nat=np.zeros([30,nlat/5+1,nlon/5+1])
		lp=np.zeros([30,nlat/5+1,nlon/5+1])	
		# Some random regions
		for n in range(30): #width/height of domain
			print n
			for i in range((nlon-n)/5):
				for j in range((nlat-n)/5):
					try:
						vals=get_region_returntimes(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,i*5,i*5+(n*3+2),j*5,j*5+(n*3+2))
					except:
						continue
					if vals[4]>0:
						hist[n,j,i]=vals[0]
						nat[n,j,i]=vals[1]
						clim_hist[n,j,i]=vals[2]
						clim_nat[n,j,i]=vals[3]
						lp[n,j,i]=vals[4]
						
	#Write out these values to a pickle file
		data2=open('data2.pkl','wb')
		pickle.dump(hist,data2,-1)
		pickle.dump(nat,data2,-1)
		pickle.dump(clim_hist,data2,-1)
		pickle.dump(clim_nat,data2,-1)
		pickle.dump(lp,data2,-1)
		data2.close()
	else:
	# read existing values from pickle	
		data2=open('data2.pkl','rb')
		hist=pickle.load(data2)
		nat=pickle.load(data2)
		clim_hist=pickle.load(data2)
		clim_nat=pickle.load(data2)
		lp=pickle.load(data2)
		data2.close()

###############################################
# Plotting	

#First get rid of annoying Mac hidden files
	delfiles=glob.glob('._*.png')
	for f in delfiles:
		os.remove(f)
		
# Some data manipulation...
	lp=lp.flatten()
	hist=hist.flatten()
	nat=nat.flatten()
	clim_hist=clim_hist.flatten()
	clim_nat=clim_nat.flatten()
	hist=np.ma.masked_where(np.logical_not(np.isfinite(hist)),hist)
	nat=np.ma.masked_where(np.logical_not(np.isfinite(nat)),nat)
	

	plt.figure(10)
	plt.plot(lp,nat.flatten(),'.',label='nat')
#	plt.plot(lp,clim_hist.flatten(),'g.',label='clim_hist')
#	plt.plot(lp,hist.flatten(),'.b',label='hist')
	plt.plot(lp,clim_nat.flatten(),'.',label='clim_nat')
	plt.title('return time vs size of region')
	plt.legend()
	plt.ylim([0,3000])
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
	plt.ylim([0,1])
	plt.legend(loc='best')
	plt.savefig('figure_far.png')
	
	plt.figure(3)
	hist=np.ma.masked_where(np.logical_or(lp<300,lp>800),hist)
	print 'non masked values',hist.shape[0] -hist.mask.sum()
	histo,bins=np.histogram((1-hist/nat).flatten(),[0,.5,.6,.7,.8,.9,.91,.92,.93,.94,.95,.96,.97,.98,.99,1.0,1.1])
	print histo,bins
	centre=bins[:-1]+(bins[1:]-bins[:-1])/2.
	plt.plot(centre[:-1],histo[:-1])
	plt.title("Histogram of FAR values")
	plt.savefig('histogram.png')
		
	#plt.show()
