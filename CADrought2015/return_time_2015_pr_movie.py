#############################################################################
# Program : return_time_2014_bootstrap.py
# Author  : Sarah Sparrow (assisted by nathalie Schaller)
# Date	: 05/03/2014
# Purpose : Plot return time periods for 2014 all forcing and natural data
# Updates : 02/04/2014 Updated from return_time_2014.py to include bootstrap 
#		   error bars
#		   25/04/2014 Updated to include CMIP5 model bootstraps
#############################################################################

import sys
import os
import numpy as np
import math
import datetime
import fnmatch
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import glob
from scipy.io.netcdf import netcdf_file
import matplotlib.animation as animat
import random

EU = os.path.expanduser
sys.path.append(EU("~massey/Coding/cpdn_analysis/"))


from return_time_plot import *

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset


def return_time_pr(i,fig,hist_data,nat_data,out_dir=False):

	# Number of ensemble members
	nens_hist=hist_data.shape[0]
	nens_nat=nat_data.shape[0]

	# Reduce number of ensemble members
	nens_hist=min(nens_hist,i)
	nens_nat=min(nens_nat,i)
	
	hist_data=hist_data[:nens_hist]
	nat_data=nat_data[:nens_nat]
	
	# Number of simulations used here
	nsims=nens_hist +nens_nat
	

	#----------------------------------------------------------------------
	# Set the plot fint size
	#---------------------------------------------------------------------
	font = {'family' : 'sans-serif',
			'size'   : 20}

	matplotlib.rc('font', **font)

	#--------------------------------------------------------------------------
	# Define model names for files and type
	#--------------------------------------------------------------------------
	types=["all_forcings","natural"]
	models=["CanESM2", "CCSM4", "CNRM-CM5", "CSIRO-Mk3-6-0", 
			"GFDL-CM3", "GISS-E2-H", "GISS-E2-R", "HadGEM2-ES", "IPSL-CM5A-LR", "IPSL-CM5A-MR", "MIROC-ESM"]
	
#	season_col_allf="RoyalBlue" 
#	allf_mec="MediumBlue"
	season_col_allf='orange'
	allf_mec='darkorange'
	
	season_col_nat="MediumSeaGreen"
	nat_mec="DarkGreen"
	
	season_col_clim='RoyalBlue'
	clim_mec='MediumBlue'

	plt.clf()
	fig.set_size_inches(8,8)
	ax = fig.add_subplot(1,1,1)
	fig.subplots_adjust(bottom=0.15)
	ax.set_ylabel("Threshold of total precipitaion (mm)",fontsize=16)
	ax.set_xlabel("Chance of having precipitation less than threshold",fontsize=16)
	plt.setp(ax.get_xticklabels(),fontsize=12)
	plt.setp(ax.get_yticklabels(),fontsize=12)


	# Hist data
	y_data_all, x_data_all = calc_return_times(hist_data,
											   direction="ascending", 
											   period=1)
#	print 'histy',y_data_all.min(),y_data_all.max()
#	print 'histx',x_data_all.min(),x_data_all.max()
	conf_all = calc_return_time_confidences(hist_data,
											direction="ascending",bsn=1e3)  
	l1=ax.semilogx(x_data_all,y_data_all, marker='o',markersize=7,
				   linestyle='None',mec=allf_mec,mfc=season_col_allf,
				   color=season_col_allf,fillstyle='full',
				   label="Actual",zorder=5)
	conf_all_5=conf_all[0,:].squeeze()
	conf_all_95=conf_all[1,:].squeeze()
	cl1=ax.fill_between(x_data_all,conf_all_5,conf_all_95,facecolor=allf_mec,edgecolor=allf_mec,alpha=0.3,linewidth=1.5,zorder=4)

	# Nat data
	y_data_all, x_data_all = calc_return_times(nat_data,
											   direction="ascending", 
											   period=1)
	conf_all = calc_return_time_confidences(nat_data,
											direction="ascending",bsn=1e3)  
	l1=ax.semilogx(x_data_all,y_data_all, marker='o',markersize=7,
				   linestyle='None',mec=nat_mec,mfc=season_col_nat,
				   color=season_col_nat,fillstyle='full',
				   label="Natural",zorder=5)
	conf_all_5=conf_all[0,:].squeeze()
	conf_all_95=conf_all[1,:].squeeze()
	cl1=ax.fill_between(x_data_all,conf_all_5,conf_all_95,facecolor=nat_mec,edgecolor=nat_mec,alpha=0.3,linewidth=1.5,zorder=4)
	
	
	ax.set_title("California Nov2014-Apr2015 Precipitation\n"+str(nsims).zfill(5)+" Simulations")

	ax.set_ylim(0,1000) #0,450
#	labels=['1','10','100','1000']
#	ax.set_yticklabels(labels)

	ax.set_xlim(1,1e3) 
	labels=['','','1/10','1/100','1/1000']
	ax.set_xticklabels(labels)

	plt.legend(loc='upper right')
	if out_dir:
		if not os.path.exists(out_dir):
        		os.mkdir(out_dir)
		fname_out=os.path.join(out_dir,"return_time_pr_"+str(nsims).zfill(5)+".png")
		fig.savefig(fname_out,dpi=28.75*2)



def calif_mean(pnw_data):
	calif_indices=[6318,6319,6320,6321,6322,6323,6324,6325,6326,6327,6328,6329,6330,6331
,6428,6429,6430,6431,6432,6433,6434,6435,6436,6437,6438,6439,6440,6441
,6538,6539,6540,6541,6542,6543,6544,6545,6546,6547,6548,6549,6550,6551
,6648,6649,6650,6651,6652,6653,6654,6655,6656,6657,6658,6659,6660,6661
,6758,6759,6760,6761,6762,6763,6764,6765,6766,6767,6768,6769,6770,6771
,6867,6868,6869,6870,6871,6872,6873,6874,6875,6876,6877,6878,6879,6880
,6881,6977,6978,6979,6980,6981,6982,6983,6984,6985,6986,6987,6988,6989
,6990,6991,7087,7088,7089,7090,7091,7092,7093,7094,7095,7096,7097,7098
,7099,7100,7101,7198,7199,7200,7201,7202,7203,7204,7205,7206,7207,7208
,7209,7210,7211,7309,7310,7311,7312,7313,7314,7315,7316,7317,7318,7319
,7320,7321,7419,7420,7421,7422,7423,7424,7425,7426,7427,7428,7429,7430
,7431,7529,7530,7531,7532,7533,7534,7535,7536,7537,7538,7539,7540,7541
,7639,7640,7641,7642,7643,7644,7645,7646,7647,7648,7649,7650,7651,7749
,7750,7751,7752,7753,7754,7755,7756,7757,7758,7759,7760,7761,7860,7861
,7862,7863,7864,7865,7866,7867,7868,7869,7870,7871,7872,7971,7972,7973
,7974,7975,7976,7977,7978,7979,7980,7981,7982,7983,7984,8081,8082,8083
,8084,8085,8086,8087,8088,8089,8090,8091,8092,8093,8094,8095,8192,8194
,8195,8196,8197,8198,8199,8200,8201,8202,8203,8204,8205,8206,8304,8305
,8306,8307,8308,8309,8310,8311,8312,8313,8314,8315,8316,8317,8413,8414
,8415,8416,8417,8418,8419,8420,8421,8422,8423,8424,8425,8426,8427,8428
,8523,8524,8525,8526,8527,8528,8529,8530,8531,8532,8533,8534,8535,8536
,8537,8538,8539,8540,8634,8635,8636,8637,8638,8639,8640,8641,8642,8643
,8644,8645,8646,8647,8648,8649,8650,8651,8744,8745,8746,8747,8748,8749
,8750,8751,8752,8753,8754,8755,8756,8757,8758,8759,8760,8761,8762,8856
,8857,8858,8859,8860,8861,8862,8863,8864,8865,8866,8867,8868,8869,8870
,8871,8872,8873,8965,8966,8967,8968,8969,8970,8971,8972,8973,8974,8975
,8976,8977,8978,8979,8980,8981,8982,8983,8984,9075,9076,9077,9078,9079
,9080,9081,9082,9083,9084,9085,9086,9087,9088,9089,9090,9091,9092,9093
,9094,9095,9186,9187,9188,9189,9190,9191,9192,9193,9194,9195,9196,9197
,9198,9199,9200,9201,9202,9203,9204,9205,9206,9207,9297,9298,9299,9300
,9301,9302,9303,9304,9305,9306,9307,9308,9309,9310,9311,9312,9313,9314
,9315,9316,9317,9318,9408,9409,9410,9411,9412,9413,9414,9415,9416,9417
,9418,9419,9420,9421,9422,9423,9424,9425,9426,9427,9428,9429,9519,9520
,9521,9522,9523,9524,9525,9526,9527,9528,9529,9530,9531,9532,9533,9534
,9535,9536,9537,9538,9539,9540,9629,9630,9631,9632,9633,9634,9635,9636
,9637,9638,9639,9640,9641,9642,9643,9644,9645,9646,9647,9648,9649,9650
,9651,9740,9741,9742,9743,9744,9745,9746,9747,9748,9749,9750,9751,9752
,9753,9754,9755,9756,9757,9758,9759,9760,9761,9762,9850,9851,9852,9853
,9854,9855,9856,9857,9858,9859,9860,9861,9862,9863,9864,9865,9866,9867
,9868,9869,9870,9871,9872,9873,9960,9961,9962,9963,9964,9965,9966,9967
,9968,9969,9970,9971,9972,9973,9974,9975,9976,9977,9978,9979,9980,9981
,9982,9983,9984,10075,10076,10077,10078,10079,10080,10081,10082,10083,10084,10085
,10086,10087,10088,10089,10090,10091,10092,10093,10186,10187,10188,10189,10190,10191
,10192,10193,10194,10195,10196,10197,10198,10199,10200,10201,10202,10299,10300,10301
,10302,10303,10304,10305,10306,10307,10308,10309,10310,10311,10312,10410,10411,10412
,10413,10414,10415,10416,10417,10418,10419,10420,10421,10422,10522,10523,10524,10525
,10526,10527,10528,10529,10530,10531,10532,10633,10634,10635,10636,10637,10638,10639
,10640,10641,10642,10643,10743,10744,10745,10746,10747,10748,10749,10750,10751,10752
,10753,10853,10854,10855,10856,10857,10858,10859]
	calif_mask=np.ones(pnw_data.shape).flatten() # Start with mask everywhere
	for pt in calif_indices:
		calif_mask[pt]=0 #unmask points in california
	calif_mask=np.reshape(calif_mask,pnw_data.shape)
	calif_data=np.ma.masked_where(calif_mask,pnw_data)
	return calif_data.mean()

#Main controling function
def main():

	# Actual
	fnames_pr_hist=glob.glob('/home/cenv0437/scratch/hadam3p_pnw/batch181/sub-batch1/ga.pe/field90/field90_????_2014-11.nc')
	random.shuffle(fnames_pr_hist)
	pr_hist=np.zeros([len(fnames_pr_hist)])
	for i,fname in enumerate(fnames_pr_hist):
		tmpsum=0
		monthfiles=glob.glob(fname[:-10]+'*') # TODO, could loop over specific months instead
		if len(monthfiles)==6:		
			for fname2 in monthfiles:
				tmp=netcdf_file(fname,'r').variables['field90'][:,0,:].mean(0) # # Ave precip rate over month 
				tmpsum=tmpsum+calif_mean(tmp)
			pr_hist[i]=tmpsum*2592000.
		else:
			pr_hist[i]=np.nan
	
	
	# Natural Forcings
	nat_subbatches=range(3,15)
	fnames_pr_nat=[]
	for subbatch in nat_subbatches:
		fnames_pr_nat=fnames_pr_nat+glob.glob('/home/cenv0437/scratch/hadam3p_pnw/batch181/sub-batch'+str(subbatch)+'/ga.pe/field90/field90_????_2014-11.nc')
	random.shuffle(fnames_pr_nat)
	pr_nat=np.zeros([len(fnames_pr_nat)])
	for i,fname in enumerate(fnames_pr_nat):
		tmpsum=0
		monthfiles=glob.glob(fname[:-10]+'*') # TODO, could loop over specific months instead
		if len(monthfiles)==6:
			for fname2 in monthfiles:
				tmp=netcdf_file(fname,'r').variables['field90'][:,0,:].mean(0) # Ave precip rate over month 
				tmpsum=tmpsum+calif_mean(tmp)
			pr_nat[i]=tmpsum*2592000. # Convert from rate per second to total over month
		else:
#			raise Exception('wrong number of files: '+str(len(monthfiles)))
			pr_nat[i]=np.nan


	make_video=True
	
	if make_video:
		# Create list of frames:
		imax=100 #Maximum number of ensembles members to include	
#		imax=max(pr_hist.shape[0],pr_nat.shape[0])
#		imax=min(pr_hist.shape[0],pr_nat.shape[0])
		
		ilist=[imax] # First frame has the max number of ensembles
		alpha=0.05
		for i in range(3,imax):
			idx=i+int(math.exp((i)*alpha))
			ilist.append(idx)
			if idx>=imax:
				break

		print 'Making video with ',str(len(ilist)),'frames'
		fig=plt.figure()			
		anim = animat.FuncAnimation(fig, return_time_pr,
									   fargs=(fig,pr_hist,pr_nat),
									   frames=ilist, interval=20, blit=True)
	 
	 	fout='/home/cenv0437/cpdn_analysis/CADrought2015/web_content/pr_movie.mp4'
		print 'created animation...',
		anim.save(fout,writer='ffmpeg',fps=20,extra_args=['-vcodec', 'libx264'])
		print "saved as",fout
	
	else:
		return_time_pr(imax,plt.figure(),pr_hist,pr_nat,'.')
		
	print 'Finished!'
	

#Washerboard function that allows main() to run on running this file
if __name__=="__main__":
  main()