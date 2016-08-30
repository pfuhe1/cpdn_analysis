#!/usr/local/bin/python2.7
#############################################################################
# Program : return_time_2014_bootstrap.py
# Author  : Peter Uhe, based on script by Sarah Sparrow (assisted by nathalie Schaller)
# Date	: 03/07/2015
# Purpose : Plot return time periods for 2015 all forcing, climatology and natural data
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
import time
import pickle

from return_time_plot import *

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

pkl_dir='/home/cenv0437/cpdn_analysis/CADrought2015'

def return_time_swe(i,fig,datad,mecs,cols,region,fname_out=False):
	print i
	
	# Number of simulations used here
	nsims=0

	#----------------------------------------------------------------------
	# Set the plot fint size
	#---------------------------------------------------------------------
	font = {'family' : 'sans-serif','size'   : 20}
	matplotlib.rc('font', **font)


	### Set up plot
	plt.clf()
	fig.set_size_inches(8,8)
	ax = fig.add_subplot(1,1,1)

	for label,data in datad.iteritems():
		# Get number of ensembles members and reduce down to i if necessary
		nens=min(data.shape[0],i)
		data=data[:nens]
		nsims=nsims+nens
		
		# Return Times
		y_data, x_data = calc_return_times(data,  direction="ascending",  period=1)
		l1=ax.semilogx(x_data,y_data, marker='o',markersize=5,
				   linestyle='None',mec=mecs[label],mfc=cols[label],
				   color=cols[label],fillstyle='full',
				   label=label,zorder=5)
		# Confidence Interval	   
		conf = calc_return_time_confidences(data,direction="ascending",bsn=1e3)  
		conf_5=conf[0,:].squeeze()
		conf_95=conf[1,:].squeeze()
		cl1=ax.fill_between(x_data,conf_5,conf_95,facecolor=mecs[label],edgecolor=mecs[label],alpha=0.3,linewidth=1.5,zorder=4)
	
	
	# Plotting Stuff
	fig.subplots_adjust(bottom=0.15)
	ax.set_ylabel("Threshold of snow amount (mm)",fontsize=16)
	ax.set_xlabel("Chance of having snow less than threshold",fontsize=16)
	plt.setp(ax.get_xticklabels(),fontsize=12)
	plt.setp(ax.get_yticklabels(),fontsize=12)
	ax.set_title(region+" April 1 Snow Amount\n"+str(nsims).zfill(5)+" Simulations")

	ylims={'California':(0,50),'Oregon':(0,30),'Washington':(0,100)}
	ax.set_ylim(ylims[region][0],ylims[region][1]) 
	ax.set_xlim(1,1e3) 
	labels=['','','1/10','1/100','1/1000']
	ax.set_xticklabels(labels)
	plt.legend(loc='upper right',numpoints=1)
	
	# Save figure if fname_out is defined
	if fname_out:
		fig.savefig(fname_out,dpi=28.75*2)

def region_mean(pnw_data,region):
	indices={}
	indices['California']=[
6318,6319,6320,6321,6322,6323,6324,6325,6326,6327,6328,6329,6330,6331
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
	indices['Washington']=[
2803,2804,2805,2806,2807,2808,2809,2810,2811,2812,2813,2814,2815,2816,2817
,2818,2819,2914,2915,2916,2917,2918,2919,2920,2921,2922,2923,2924,2925,2926
,2927,2928,2929,3023,3024,3025,3026,3027,3028,3029,3030,3031,3032,3033,3034
,3035,3036,3037,3038,3039,3127,3128,3134,3135,3136,3137,3138,3139,3140,3141
,3142,3143,3144,3145,3146,3147,3148,3149,3237,3238,3239,3240,3241,3242,3245
,3246,3247,3248,3249,3250,3251,3252,3253,3254,3255,3256,3257,3258,3259,3348
,3349,3350,3351,3352,3353,3354,3355,3356,3357,3358,3359,3360,3361,3362,3363
,3364,3365,3366,3367,3368,3369,3370,3458,3459,3460,3461,3462,3463,3464,3465
,3466,3467,3468,3469,3470,3471,3472,3473,3474,3475,3476,3477,3478,3479,3480
,3568,3569,3570,3571,3572,3573,3574,3575,3576,3577,3578,3579,3580,3581,3582
,3583,3584,3585,3586,3587,3588,3589,3590,3678,3679,3680,3681,3682,3683,3684
,3685,3686,3687,3688,3689,3690,3691,3692,3693,3694,3695,3696,3697,3698,3699
,3700,3789,3790,3791,3792,3793,3794,3795,3796,3797,3798,3799,3800,3801,3802
,3803,3804,3805,3806,3807,3808,3809,3810,3899,3900,3901,3902,3903,3904,3905
,3906,3907,3908,3909,3910,3911,3912,3913,3914,3915,3916,3917,3918,3919,3920
,4009,4010,4011,4012,4013,4014,4015,4016,4017,4018,4019,4020,4021,4022,4023
,4024,4025,4026,4027,4028,4029,4030,4121,4122,4123,4124,4125,4126,4127,4128
,4129,4130,4131,4132,4133,4134,4135,4136,4137,4138,4139,4140,4232,4233,4234
,4235,4236,4237,4238,4239,4240,4241,4242,4243,4244,4343,4344,4345,4346,4347
,4348,4349,4350,4454]
	indices['Oregon']=[
4120,4229,4230,4231,4245,4246,4247,4248,4249,4250,4339,4340,4341,4342,4351
,4352,4353,4354,4355,4356,4357,4358,4359,4360,4361,4362,4449,4450,4451,4452
,4453,4455,4456,4457,4458,4459,4460,4461,4462,4463,4464,4465,4466,4467,4468
,4469,4470,4471,4472,4559,4560,4561,4562,4563,4564,4565,4566,4567,4568,4569
,4570,4571,4572,4573,4574,4575,4576,4577,4578,4579,4580,4581,4669,4670,4671
,4672,4673,4674,4675,4676,4677,4678,4679,4680,4681,4682,4683,4684,4685,4686
,4687,4688,4689,4690,4691,4779,4780,4781,4782,4783,4784,4785,4786,4787,4788
,4789,4790,4791,4792,4793,4794,4795,4796,4797,4798,4799,4800,4801,4889,4890
,4891,4892,4893,4894,4895,4896,4897,4898,4899,4900,4901,4902,4903,4904,4905
,4906,4907,4908,4909,4910,4998,4999,5000,5001,5002,5003,5004,5005,5006,5007
,5008,5009,5010,5011,5012,5013,5014,5015,5016,5017,5018,5019,5020,5108,5109
,5110,5111,5112,5113,5114,5115,5116,5117,5118,5119,5120,5121,5122,5123,5124
,5125,5126,5127,5128,5129,5130,5131,5218,5219,5220,5221,5222,5223,5224,5225
,5226,5227,5228,5229,5230,5231,5232,5233,5234,5235,5236,5237,5238,5239,5240
,5241,5328,5329,5330,5331,5332,5333,5334,5335,5336,5337,5338,5339,5340,5341
,5342,5343,5344,5345,5346,5347,5348,5349,5350,5351,5438,5439,5440,5441,5442
,5443,5444,5445,5446,5447,5448,5449,5450,5451,5452,5453,5454,5455,5456,5457
,5458,5459,5460,5461,5547,5548,5549,5550,5551,5552,5553,5554,5555,5556,5557
,5558,5559,5560,5561,5562,5563,5564,5565,5566,5567,5568,5569,5570,5571,5657
,5658,5659,5660,5661,5662,5663,5664,5665,5666,5667,5668,5669,5670,5671,5672
,5673,5674,5675,5676,5677,5678,5679,5680,5681,5767,5768,5769,5770,5771,5772
,5773,5774,5775,5776,5777,5778,5779,5780,5781,5782,5783,5784,5785,5786,5787
,5788,5789,5790,5791,5877,5878,5879,5880,5881,5882,5883,5884,5885,5886,5887
,5888,5889,5890,5891,5892,5893,5894,5895,5896,5897,5898,5899,5900,5901,5987
,5988,5989,5990,5991,5992,5993,5994,5995,5996,5997,5998,5999,6000,6001,6002
,6003,6004,6005,6006,6007,6008,6009,6010,6011,6097,6098,6099,6100,6101,6102
,6103,6104,6105,6106,6107,6108,6109,6110,6111,6112,6113,6114,6115,6116,6117
,6118,6119,6120,6121,6208,6209,6210,6211,6212,6213,6214,6215,6216,6217,6218
,6219,6220,6221,6222,6223,6224,6225,6226,6227,6228,6229]

	region_mask=np.ones(pnw_data.shape).flatten() # Start with mask everywhere	
	for pt in indices[region]:
		region_mask[pt]=0 #unmask points in california
	region_mask=np.reshape(region_mask,pnw_data.shape)
	region_data=np.ma.masked_where(region_mask,pnw_data)
	region_data=np.ma.masked_values(region_data,-1.073742e+09) # mask missing values
	return region_data.mean()


#Main controling function
def main():
	
	months=['2014-12','2015-01','2015-02']
	regions=['California','Oregon','Washington']
	
	f_pickle=pkl_dir+'/swe_results.pkl'
	if os.path.exists(f_pickle):
		fdata=open(f_pickle,'rb')
		data_dict_hist=pickle.load(fdata)
		data_dict_clim=pickle.load(fdata)
		data_dict_nat=pickle.load(fdata)
		fdata.close()
	else:
		data_dict_hist={}
		data_dict_clim={}
		data_dict_nat={}
		
	f_pickle2=pkl_dir+'/swe_results2.pkl'
	if os.path.exists(f_pickle2):
		fdata=open(f_pickle2,'rb')
		data_dict_nat2=pickle.load(fdata)
	else:
		data_dict_nat2={}


	# Actual
	fnames_data_hist=glob.glob('/home/cenv0437/scratch/hadam3p_pnw/batch181/sub-batch15/ga.pd/field93/field93_????_2015-04.nc')
	fnames_data_hist=fnames_data_hist+glob.glob('/home/cenv0437/scratch/hadam3p_pnw/batch181/sub-batch1/ga.pd/field93/field93_????_2015-04.nc')
#	random.shuffle(fnames_data_hist)
	print len(fnames_data_hist)
	for i,fname in enumerate(fnames_data_hist):
		umid=os.path.basename(fname).split('_')[1]
		if umid not in data_dict_hist.keys():
			tmplist=np.zeros([len(regions)])
			try:	
				print os.path.basename(fname)
				tmp=netcdf_file(fname,'r').variables['field93'][0,0,:] # Take first of april only
				for j,region in enumerate(regions):
					tmplist[j]=region_mean(tmp,region)
				data_dict_hist[umid]=tmplist
			except Exception,e:
				print "Error loading files",e
				continue

				
			
	# Climatology
	fnames_data_clim=glob.glob('/home/cenv0437/scratch/hadam3p_pnw/batch181/sub-batch29/ga.pd/field93/field93_????_2015-04.nc')
	random.shuffle(fnames_data_clim)
	print len(fnames_data_clim)
	data_clim=np.zeros([len(regions),len(fnames_data_clim)])
	for i,fname in enumerate(fnames_data_clim):
		umid=os.path.basename(fname).split('_')[1]
		if umid not in data_dict_clim.keys():
			tmplist=np.zeros([len(regions)])
			try:	
				print os.path.basename(fname)
				tmp=netcdf_file(fname,'r').variables['field93'][0,0,:] # Take first of april only
				for j,region in enumerate(regions):
					tmplist[j]=region_mean(tmp,region)
				data_dict_clim[umid]=tmplist
			except Exception,e:
				print "Error loading files",e
				continue


	# Natural Forcings
	nat_subbatches=range(3,15)+range(17,29)
	fnames_data_nat=[]
	for subbatch in nat_subbatches:
		fnames_data_nat=fnames_data_nat+glob.glob('/home/cenv0437/scratch/hadam3p_pnw/batch181/sub-batch'+str(subbatch)+'/ga.pd/field93/field93_????_2015-04.nc')
	random.shuffle(fnames_data_nat)
	print len(fnames_data_nat)

	for i,fname in enumerate(fnames_data_nat):
		umid=os.path.basename(fname).split('_')[1]
		if umid not in data_dict_nat.keys():
			tmplist=np.zeros([len(regions)])
			try:	
				print os.path.basename(fname)
				tmp=netcdf_file(fname,'r').variables['field93'][0,0,:] # Take first of april only
				for j,region in enumerate(regions):
					tmplist[j]=region_mean(tmp,region)
				data_dict_nat[umid]=tmplist
			except Exception,e:
				print "Error loading files",e
				continue

	# Natural Forcings
	nat_subbatches2=range(30,42)
	fnames_data_nat2=[]
	for subbatch in nat_subbatches2:
		fnames_data_nat2=fnames_data_nat2+glob.glob('/home/cenv0437/scratch/hadam3p_pnw/batch181/sub-batch'+str(subbatch)+'/ga.pd/field93/field93_????_2015-04.nc')
	random.shuffle(fnames_data_nat2)
	print len(fnames_data_nat2)

	for i,fname in enumerate(fnames_data_nat2):
		umid=os.path.basename(fname).split('_')[1]
		if umid not in data_dict_nat2.keys():
			tmplist=np.zeros([len(regions)])
			try:	
				print os.path.basename(fname)
				tmp=netcdf_file(fname,'r').variables['field93'][0,0,:] # Take first of april only
				for j,region in enumerate(regions):
					tmplist[j]=region_mean(tmp,region)
				data_dict_nat2[umid]=tmplist
			except Exception,e:
				print "Error loading files",e
				continue

		
	# Save data for reuse later
	fdata=open(f_pickle,'wb')
	pickle.dump(data_dict_hist,fdata,-1)
	pickle.dump(data_dict_clim,fdata,-1)			
	pickle.dump(data_dict_nat,fdata,-1)
	fdata.close()
	
	fdata=open(f_pickle2,'wb')
	pickle.dump(data_dict_nat2,fdata,-1)
	fdata.close()
	
	mecs={}
	cols={}
	
	# Define color scheme:
	cols['Actual']='orange'
	mecs['Actual']='darkorange'
	
	cols['Natural']="MediumSeaGreen"
	mecs['Natural']="DarkGreen"
	
	cols['Climatology']='RoyalBlue'
	mecs['Climatology']='MediumBlue'
	
	cols['NaturalFix']='Purple'
	mecs['NaturalFix']='DarkMagenta'
	
	
	# Set Maximum number of ensembles members to include
#		imax=100	
	imax=max([len(data_dict_hist.keys()),len(data_dict_clim.keys()),len(data_dict_nat2.keys())])
	nens_total=sum([len(data_dict_hist.keys()),len(data_dict_clim.keys()),len(data_dict_nat2.keys())])
	
	# Create list of frames:
	ilist=[imax] # First frame has the max number of ensembles
	alpha=0.07
	for i in range(3,imax):
		idx=i+int(math.exp((i)*alpha))
		ilist.append(idx)
		if idx>=imax:
			break
		
	for i,region in enumerate(regions):
		datad={}
		datad['Actual']=np.array(data_dict_hist.values())[:,i]
		datad['Climatology']=np.array(data_dict_clim.values())[:,i]
		datad['Natural']=np.array(data_dict_nat2.values())[:,i]
		#datad['NaturalFix']=np.array(data_dict_nat2.values())[:,i]	

		fout='/home/cenv0437/cpdn_analysis/CADrought2015/archive_fix/swe_'+region+'_n'+str(nens_total)
		webfile='/home/cenv0437/cpdn_analysis/CADrought2015/web_content/swe_'+region

		# Plot histograms of data
#		plt.figure()
#		for label,data in datad.iteritems():
#			plt.hist(data,50,histtype='stepfilled',normed=1,facecolor=cols[label],edgecolor=cols[label],alpha=0.2,label=label)
#			plt.axvline(x=data.mean(),color=cols[label])
#		plt.legend(loc='upper left')
#		plt.savefig('archive_extra/swe_hist_'+region+'.png')
		
		# Make still image
		f_ext='.png'
		return_time_swe(imax,plt.figure(),datad,mecs,cols,region,fname_out=fout+f_ext)
		print "saved as",fout+f_ext
		# Link file in web_content folder
		if os.path.exists(webfile+f_ext):
			os.remove(webfile+f_ext)
		os.symlink(fout+f_ext,webfile+f_ext)

#		# Make video
		f_ext='.mp4'
		print 'Making video with ',len(ilist),'frames'
		fig=plt.figure()
		anim = animat.FuncAnimation(fig, return_time_swe,
			fargs=(fig,datad,mecs,cols,region),
			frames=ilist, interval=20, blit=True)
		print 'created animation...',
		anim.save(fout+f_ext,writer='ffmpeg',fps=20,extra_args=['-vcodec', 'libx264'])
		print "saved as",fout+f_ext
		# Link file in web_content folder
		if os.path.exists(webfile+f_ext):
			os.remove(webfile+f_ext)
		os.symlink(fout+f_ext,webfile+f_ext)
		
	
	print 'Finished!'

#Washerboard function that allows main() to run on running this file
if __name__=="__main__":
  main()
