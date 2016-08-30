#!/usr/local/bin/python2.7
# NOTE python2.7 must be used for this script

import glob,os
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle
from mpl_toolkits.basemap import Basemap,cm
from region_returntimes import get_region_returntimes_sampled, get_region_returntimes_bootstrap_stats,remap_p5deg_to_rotated, find_closest
from scipy.io.netcdf import netcdf_file
import matplotlib.patches as mpatches

pkl_dir='/gpfs/projects/cpdn/scratch/cenv0437/pickle_data/'

class Arraylist:
	def __init__(self):
		self.lort1900=np.zeros(40)
		self.rt1900=np.zeros(40)
		self.lort2014=np.zeros(40)
		self.rt2014=np.zeros(40)
		self.loratio=np.zeros(40)
		self.ratio=np.zeros(40)

def parse(var):
	if var=='infinity':
		return 1.e20
	else:
		return float(var)


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

#	obs_p5deg= '/home/cenv0437/scratch/data_from_ouce/CRU_TS_Absolute_plus_BEST_anomaly_201312_201412.nc'
#	obsname='BEST'

#	obs_p5deg= '/home/cenv0437/scratch/data_from_ouce/CRU_TS_Absolute_plus_CRU4_anomaly_201312_201412.nc'
#	obsname='CRU4'

	obs_p5deg= '/home/cenv0437/scratch/data_from_ouce/CRU_TS_Absolute_plus_GISS_anomaly_201312_201412.nc'
	obsname='GISS'

#	obs_p5deg= '/home/cenv0437/scratch/data_from_ouce/CRU_TS_Absolute_plus_NCEP_anomaly_201312_201412_all.nc'
#	obsname='NCEP'

#	obs_p5deg= '/home/cenv0437/scratch/data_from_ouce/CRU_TS_Absolute_plus_CFSR_anomaly_201312_201412_all.nc'
#	obsname='CFSR'


	# Some specific regions:
	regions=[(51.2, 13.8),(62.561671, 13.958159),(51.749181, -1.281326),(55.685268, 37.274331),(43.959139, 20.910196),(40.402960, -3.732137)]
	rnames=['Germany','Sweden','England','Russia','Serbia','Spain']
	rcolors=['k','g','b','r','c','m']

# Cacluate return times for many regions

	print historical_data.shape	
	nens,nlat,nlon=historical_data.shape
	hist=np.zeros([len(regions),50])
	hist_up=np.zeros([len(regions),50])
	hist_low=np.zeros([len(regions),50])
	nat=np.zeros([len(regions),50])
	nat_up=np.zeros([len(regions),50])
	nat_low=np.zeros([len(regions),50])
	clim_hist_low=np.zeros([len(regions),50])
	clim_hist=np.zeros([len(regions),50])
	clim_hist_up=np.zeros([len(regions),50])
	clim_nat_low=np.zeros([len(regions),50])
	clim_nat=np.zeros([len(regions),50])
	clim_nat_up=np.zeros([len(regions),50])
	lp=np.zeros([len(regions),50])
	hist_var=np.zeros([len(regions),50]) #variance of hist ensemble 
	nat_var=np.zeros([len(regions),50]) # variance of nat ensemble
	

	# Regrid obs to rotated regional grid
	rot_template='/home/cenv0437/scratch/data_from_ouce/hadam3p_eu_z4ao_2013_1_009238311_0_tasmean.nc'
	f_rot=netcdf_file(rot_template,'r')
	lat_coord=f_rot.variables['global_latitude0'][4:-7,4:-4]
	lon_coord=f_rot.variables['global_longitude0'][4:-7,4:-4]
	f_rot.close()	
	obs=remap_p5deg_to_rotated(obs_p5deg,rot_template)[:,:]
		

	for reg,rval in enumerate(regions):
		j0,i0=find_closest(rval[0],rval[1],lat_coord,lon_coord)
		print rnames[reg],':',lat_coord[j0,i0],lon_coord[j0,i0]
		for n in range(50): #width/height of domain
			print n,
			try:
					vals=get_region_returntimes_bootstrap_stats(historical_data,natural_data,clim_hist_data,clim_nat_data,obs,i0-n,i0+n+1,j0-n,j0+n+1)
			except Exception,e:
				print e
				break
			if vals[4]>0: # If there are any land points to analyse
				hist_low[reg,n]=vals[0][0]
				hist[reg,n]=vals[0][1]
				hist_up[reg,n]=vals[0][2]
				nat_low[reg,n]=vals[1][0]
				nat[reg,n]=vals[1][1]
				nat_up[reg,n]=vals[1][2]
				clim_hist_low[reg,n]=vals[2][0]
				clim_hist[reg,n]=vals[2][1]
				clim_hist_up[reg,n]=vals[2][2]
				clim_nat_low[reg,n]=vals[3][0]
				clim_nat[reg,n]=vals[3][1]
				clim_nat_up[reg,n]=vals[3][2]
			print ''

###############################################
# Plotting	

#First get rid of annoying Mac hidden files
	delfiles=glob.glob('scatter_figs/._*.png')
	for f in delfiles:
		os.remove(f)
		
# Some data manipulation...
#	lp=lp.flatten()
#	hist=hist.flatten()
#	nat=nat.flatten()
#	clim_hist=clim_hist.flatten()
#	clim_nat=clim_nat.flatten()
	

	
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

# Empirical data:

	# Open file:
	fnames=glob.glob('../../scratch/data_from_ouce/empirical_fig2/*')
	regiondata={}
	
	for f in fnames:
		region=os.path.basename(f)[:-4]
		print region
		regiondata[region]=Arraylist() #Initialise arrays
		for line in open(f):
			if line.strip()[0]=='#':
				continue #skip commented lines
			try:
				vars=line.split()
				n=int(vars[0])
				regiondata[region].lort1900[n]=parse(vars[1])
				regiondata[region].rt1900[n]=parse(vars[2])
				regiondata[region].lort2014[n]=parse(vars[3])
				regiondata[region].rt2014[n]=parse(vars[4])
				regiondata[region].loratio[n]=parse(vars[5])
				regiondata[region].ratio[n]=parse(vars[6])
			except Exception as e:
				print 'Error, failure to parse line:',line
				print e

#Approx width of box at 50 degrees North	
	xvals=(2.*np.arange(40)+1.)*40075/360.*np.cos(0.87) 
# Colors for the plots
	dict_colors={'Germany':'k', 'Sweden':'g', 'England':'b', 'Russia':'r', 'Serbia':'c', 'Spain':'m'}
	
	print '###########################################################\n'
	
# Plotting stuff

#First get rid of annoying Mac hidden files
	delfiles=glob.glob('region_figs/._*.png')
	for f in delfiles:
		os.remove(f)
		

	
plt.figure(1)
plt.clf()
#	plt.subplot(131)
#	plt.plot([],[],'k--',label='Region:')
#	plt.figure(2)
#	plt.subplot(132)
#	plt.plot([],[],'k--',label='Region:')
#	plt.figure(3)
#	plt.subplot(133)
#	plt.plot([],[],'k--',label='Region:')

box_width=np.arange(50,5000,100)
reg_lines=[]
hatch_hist='oo'
hatch_nat='OO'
hatch_2014='//'
hatch_clim='\\\\'
hatch_emp='..'

# Plot the data for return times / probabilities for each region
for i,reg in enumerate(rnames):
	print reg
	
	plt.subplot(321)
	plt.semilogy(box_width,1/hist[i,:], rcolors[i]+'.-')#,label=reg+'-hist')
	plt.fill_between( box_width,1/hist_up[i,:],1/hist_low[i,:],color=rcolors[i],alpha=0.1,hatch=hatch_hist)
	plt.semilogy(box_width,1/nat[i,:], rcolors[i]+'-')#,label=reg)#+'-nat')
	plt.fill_between( box_width,1/nat_up[i,:],1/nat_low[i,:],color=rcolors[i],alpha=0.1,hatch=hatch_nat)
	
	plt.subplot(323)
	plt.semilogy(box_width,1/clim_hist[i,:], rcolors[i]+'.-')#,label=reg+'-hist')
	plt.fill_between( box_width,1/clim_hist_up[i,:],1/clim_hist_low[i,:], color=rcolors[i],alpha=0.1,hatch=hatch_hist)
	plt.semilogy(box_width,1/clim_nat[i,:], rcolors[i]+'-')#,label=reg)#+'-nat')
	plt.fill_between( box_width,1/clim_nat_up[i,:],1/clim_nat_low[i,:], color=rcolors[i],alpha=0.1,hatch=hatch_nat)
	
	plt.subplot(325)
	plt.plot(box_width,(1-hist[i,:]/nat[i,:]), rcolors[i]+'.-')#,label=reg+'-hist')
	plt.plot(box_width,(1-clim_hist[i,:]/clim_nat[i,:]), rcolors[i]+'-')#,label=reg)#+'-nat')
	plt.fill_between( box_width,1-hist_up[i,:]/nat_low[i,:],1-hist_low[i,:]/nat_up[i,:], color=rcolors[i],alpha=0.1,hatch=hatch_2014)
	plt.fill_between( box_width,1-clim_hist_up[i,:]/clim_nat_low[i,:],1-clim_hist_low[i,:]/clim_nat_up[i,:], color=rcolors[i],alpha=0.1,hatch=hatch_clim)
	
	plt.subplot(322)
	reg_lines.append(plt.plot([],[],rcolors[i],label=reg)[0])
	
	plt.subplot(324)
#		plt.semilogy(xvals,1/regiondata[reg].lort2014, '.-'+dict_colors[reg])#,linewidth=0.5,markersize=2.5)
	plt.semilogy(xvals,1/regiondata[reg].rt2014, '.-'+dict_colors[reg],linewidth=1)
#		plt.semilogy(xvals,1/regiondata[reg].lort1900, '-'+dict_colors[reg])#,linewidth=0.5)
	plt.semilogy(xvals,1/regiondata[reg].rt1900, '-'+dict_colors[reg],linewidth=1)
	plt.fill_between( xvals,1/regiondata[reg].rt1900,1/regiondata[reg].lort1900, alpha=0.1,color=dict_colors[reg],hatch=hatch_nat)
	plt.fill_between( xvals,1/regiondata[reg].rt2014,1/regiondata[reg].lort2014, alpha=0.1,color=dict_colors[reg],hatch=hatch_hist)
	
	plt.subplot(326)
	plt.plot(xvals,1-1/regiondata[reg].ratio,'-'+dict_colors[reg])
#		plt.plot(xvals,1-1/regiondata[reg].loratio,dict_colors[reg])
	plt.fill_between( xvals,1-1/regiondata[reg].loratio,1-1/regiondata[reg].ratio,alpha=0.1,color=dict_colors[reg],hatch=hatch_emp)

###################################
# Plot formatting

plt.subplot(321)
#	plt.plot([],[],'k--',label='Ensemble:')
plt.plot([],[],'k.-',label='hist2014')
plt.plot([],[],'k-',label='nat2014')
#	plt.title('Probablility using 2014 ensembles')
plt.title('a)')
#	plt.legend(loc='best')
plt.ylim([0.001,1])
plt.xlim([0,4000])
plt.xlabel('Scale of region (km)')
plt.ylabel('Probability of event')
ax=plt.gca()
ticks=ax.get_xticks()
ax.set_xticks(ticks[::2])

plt.subplot(323)
#	plt.plot([],[],'k--',label='Ensemble:')
plt.plot([],[],'k.-',label='histClim')
plt.plot([],[],'k-',label='natClim')
plt.title('b)')
#	plt.title('Probablility using Clim ensembles')
#	plt.legend(loc='best')
plt.ylim([0.001,1])
plt.xlim([0,4000])
plt.xlabel('Scale of region (km)')
plt.ylabel('Probability of event')
ax=plt.gca()
ticks=ax.get_xticks()
ax.set_xticks(ticks[::2])

plt.subplot(325)
#	plt.plot([],[],'k--',label='Ensemble:')
plt.plot([],[],'k.-',label='2014')
plt.plot([],[],'k-',label='Clim')
plt.title('d)')
#	plt.title('FAR')
#	plt.legend(loc='best')
plt.ylim([0.7,1])
plt.xlim([0,4000])
plt.xlabel('Scale of region (km)')
plt.ylabel('FAR')
ax=plt.gca()
ticks=ax.get_xticks()
ax.set_xticks(ticks[::2])


plt.subplot(324)
plt.title('c)')
#plt.plot([],[],'k.-',label='lort2014')
#plt.plot([],[],'k.-.',label='rt2014')
#plt.plot([],[],'k-',label='lort1901')
#plt.plot([],[],'k',label='rt1901')
plt.ylim([0.001,1])
plt.xlim([0,4000])
plt.xlabel('Scale of region (km)')
plt.ylabel('Probability of event')
ax=plt.gca()
ticks=ax.get_xticks()
ax.set_xticks(ticks[::2])

plt.subplot(326)
plt.title('e)')
plt.plot([],[],'k.-',label='loratio')
plt.plot([],[],'k-.',label='ratio')
plt.ylim([0.7,1])
plt.xlim([0,4000])
plt.xlabel('Scale of region (km)')
plt.ylabel('FAR')
ax=plt.gca()
ticks=ax.get_xticks()
ax.set_xticks(ticks[::2])

plt.subplot(322)
#plt.title('b)')
# Put in legend for colors
plt.axis('off')
dotline=plt.plot([],[],'k.-')[0]
line=plt.plot([],[],'k-')[0]
dashdot=plt.plot([],[],'k-.')[0]

leg1= plt.legend(handles=reg_lines,loc='right',fontsize='small')
plt.gca().add_artist(leg1)

patches=[]
patchlabels=['Hist Prob','Nat Prob','2014 FAR','Clim FAR','Empirical FAR']
patches.append((dotline,mpatches.Patch(color='k',hatch=hatch_hist,alpha=0.3)))
patches.append((line,mpatches.Patch(color='k',hatch=hatch_nat,alpha=0.3)))
patches.append((dotline,mpatches.Patch(color='k',hatch=hatch_2014,alpha=0.3)))
patches.append((line,mpatches.Patch(color='k',hatch=hatch_clim,alpha=0.3)))
patches.append((line,mpatches.Patch(color='k',hatch=hatch_emp,alpha=0.3)))

plt.legend(patches,patchlabels,loc='center left',fontsize='small')

plt.tight_layout()
plt.savefig('combined_figs/figure_combined_'+obsname+'.png')


#####################################
# Far plot
#	plt.figure(3)
	#plt.plot(lp.flatten()**0.5*50.,(1-hist.flatten()/nat.flatten()),'.',label='2014')
#	plt.plot(lp,clim_hist.flatten(),'g.',label='clim_hist')
#	plt.plot(lp,hist.flatten(),'.b',label='hist')
	#plt.plot(lp.flatten()**0.5*50.,(1-clim_hist.flatten()/clim_nat.flatten()),'.',label='clim')
#	plt.title('relative risk vs size of region')
#	plt.legend()
#	plt.ylim([0,3000])
#	plt.xlim([0,5000])
#	plt.xlabel('Scale of region (km)')
#	plt.ylabel('RR')
#	plt.savefig('figure_far_specific.png')

	
