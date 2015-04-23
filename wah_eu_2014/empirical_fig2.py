import numpy as np
import glob,os
from netcdf_file import netcdf_file
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle
from mpl_toolkits.basemap import Basemap,cm	

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
		return np.inf
	else:
		return float(var)

if __name__=='__main__':


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
	
	for reg,data in regiondata.iteritems():

		plt.figure(1)
		plt.plot(xvals,1-1/data.ratio,':'+dict_colors[reg])
		plt.plot(xvals,1-1/data.loratio,dict_colors[reg])

		
		plt.figure(2)
		plt.semilogy(xvals,1/data.lort2014,'.-'+dict_colors[reg])#,linewidth=0.5,markersize=2.5)
		plt.semilogy(xvals,1/data.rt2014,'-.'+dict_colors[reg],linewidth=1)
	

		plt.semilogy(xvals,1/data.lort1900,'-'+dict_colors[reg])#,linewidth=0.5)
		plt.semilogy(xvals,1/data.rt1900,':'+dict_colors[reg],linewidth=1)

		plt.figure(3)
		plt.semilogy(xvals,1/data.rt1900,'-'+dict_colors[reg],linewidth=1)
		plt.semilogy(xvals,1/data.rt2014,'.-'+dict_colors[reg],linewidth=1)
	
	plt.figure(1)
	plt.savefig('tmp_figs/emipircal_fig21.png')
	plt.figure(2)	
	plt.ylim([0.001,1])
	plt.savefig('tmp_figs/emipircal_fig22.png')
	plt.figure(3)
	plt.ylim([0.001,1])
	plt.savefig('tmp_figs/emipircal_fig23.png')
