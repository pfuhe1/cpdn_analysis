# Eof analysis of soil moisture from dumps
#soil levels
#10,30,60,200  (thickness, cm)
#0,10,40,100,300 (boundary cm)

import sys,os
import numpy as np
import glob
from scipy.io import netcdf_file
#eofs installed to /gpfs/home/cenv0437/.local/lib/python2.7/site-packages/
sys.path.append('/gpfs/home/cenv0437/.local/lib/python2.7/site-packages/')
from eofs.standard import Eof
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap,cm

# Load ensemble of files into a single array
def load_regional(filenames):
	data=np.ma.zeros([len(filenames),105,98])
	i=0
	for f in filenames:
		try:
			print i
			tmp=netcdf_file(f,'r').variables['sm'][:]
			if not np.all(np.isfinite(tmp)):
				print 'error: wierd vals',f
				continue
			else:
				# Integral over top 1m of soil
				data[i,:]=tmp[0,0,6:-6,6:-6]*.1+tmp[0,1,6:-6,6:-6]*.3+tmp[0,2,6:-6,6:-6]*.6 
				i=i+1
		except Exception,e:
			print 'Error, cannot load files',f,e
			continue
	return data[:i,:]
	
	
# Load ensemble of files into a single array
def load_global(filenames):
	data=np.ma.zeros([len(filenames),145,192])
	i=0
	for f in filenames:
		try:
			tmp=netcdf_file(f,'r').variables['sm']
			if not np.all(np.isfinite(tmp)):
				print 'error: wierd vals',f
				continue
			else:
				# Integral over top 1m of soil
				data[i,:]=tmp[0,0,6:-6,6:-6]*.1+tmp[0,0,6:-6,6:-6]*.3+tmp[0,0,6:-6,6:-6]*.6 
				i=i+1
		except:
			print 'Error, cannot load files',f
			continue
	return data[:i,:]
	
def plot_region_pnw(data,lat_coord,lon_coord,lat_nw,lat_se,lon_nw,lon_se,title,varfrac):
	plt.clf()
	circles=[30,40,50,60,70]
	meridians=[-130,-120,-110]
	ny,nx=lat_coord.shape
	m = Basemap(projection='stere',lon_0=lon_coord[ny/2,nx/2],lat_0=lat_coord[ny/2,nx/2],resolution='l',llcrnrlon=lon_coord[-1,0],llcrnrlat=lat_coord[-1,0],urcrnrlon=lon_coord[0,-1],urcrnrlat=lat_coord[0,-1])
	x,y=m(lon_coord[:],lat_coord[:])
	
	plt.title('Plotting data for: '+title+'\nvariance:'+str(varfrac))
	try:
		max=np.absolute(data[lat_nw:lat_se,lon_nw:lon_se]).max()
		m.pcolor(x[lat_nw:lat_se,lon_nw:lon_se],y[lat_nw:lat_se,lon_nw:lon_se],data[lat_nw:lat_se,lon_nw:lon_se],vmin=-max,vmax=max)
		plt.colorbar()
	except:
		raise
		m.plot(x[lat_nw:lat_se,lon_nw:lon_se],y[lat_nw:lat_se,lon_nw:lon_se],'r.') # Just put a dot at the x y point
	m.drawcoastlines()
	m.drawcountries()
	y=m.drawstates()
	m.drawparallels(circles)
	m.drawmeridians(meridians)
	try:
		os.remove(output_dir+'/._'+title+'.png')
	except:
		pass
	plt.savefig(output_dir+'/'+title+'.png')
	
if __name__=='__main__':
		
		output_dir='eof_plots_hist'
		if not os.path.exists(output_dir):
			os.mkdir(output_dir)
		filenames=glob.glob('/gpfs/projects/cpdn/scratch/cenv0437/dump_sm/hadam3p_pnw/batch181/sub-batch15/region_*.nc')
#		for batchid in range(30,42):
#			filenames.extend(glob.glob('/gpfs/projects/cpdn/scratch/cenv0437/dump_sm/hadam3p_pnw/batch181/sub-batch'+str(batchid)+'/region_*.nc'))
		fexample=netcdf_file('/home/cenv0437/scratch/tmp2/pfzkga.pel3dec.nc','r')
		lat_coord=fexample.variables['global_latitude0'][6:-6,6:-6]
		lon_coord=fexample.variables['global_longitude0'][6:-6,6:-6]
	
		#np.random.shuffle(filenames)
		data=load_regional(filenames)
		data=np.ma.masked_values(data,2.e+20)
		print "data loaded",data.shape
	
		# Set up info
		plt.set_cmap('RdBu')
		neofs=5
		nens=data.shape[0]
		nwanted=57
	
		print "\nWhole ensemble..."
		solver=Eof(data)
		print 'set up EOF solver'
		pcs=solver.pcs(npcs=neofs,pcscaling=1)
		eofs=solver.eofs(neofs=neofs)
		varfrac=solver.varianceFraction(neigs=neofs)
		print 'calculated EOFs'
		print 'printing EOFs'
		for i in range(neofs):
			print 'EOF',i
			plt.clf()
			plot_region_pnw(eofs[i,:],lat_coord,lon_coord,0,-1,0,-1,'EOF'+str(i),varfrac[i])
		print "plotting histograms of PCs"
		for i in range(3):
			plt.clf()
			plt.hist(pcs[:,i],200,range=(-4,4),normed=1,alpha=0.4,label='pc'+str(i))
			plt.ylim([0,.6])
			plt.savefig(output_dir+'/histogram_pc'+str(i)+'.png')
		print "plotting mean and stdev of ensemble"
		plot_region_pnw(data[:].mean(0),lat_coord,lon_coord,0,-1,0,-1,'mean',data.mean())
		plot_region_pnw(data[:].std(0),lat_coord,lon_coord,0,-1,0,-1,'stdev',data.std())

		print "\nSubsampled ensemble..."
		radius2D=(pcs[:,0]**2+pcs[:,1]**2)**.5#+pcs[:,2]**2)**.5
		radius3D=(pcs[:,0]**2+pcs[:,1]**2+pcs[:,2]**2)**.5
		radius=radius3D # Use first 2 EOFs to choose subsample

		# TODO: change from using radius to probability distribution
		#maskranges=[0.0,0.5,1.,1.5,2.,2.5]
		#maskranges=np.arange(0,2.8,.4)
		#np.append(maskranges,[10.])
		maskranges=[0,.5,.8,1.,1.2,1.4,1.6,1.8,2.0,2.3,2.6,10.]
		colors=['r','orange','brown','purple','green','r','orange','brown','purple','green','r','orange','brown','purple','green','r','orange','brown','purple','green']

		print "selecting sub sample and ploting as PC scatterplot"
		selected=0
		masksubset=np.zeros([nens],dtype='bool')
		plt.clf()
		plt.scatter(pcs[:,0],pcs[:,1])
		for i in range(len(maskranges)-1):
			mask=np.logical_and(radius>maskranges[i],radius<=maskranges[i+1])
			nvals=pcs[mask,0].shape[0]
		#	nselect=min(10,nvals)
			nselect= (nvals* nwanted)/nens
			print 'range',maskranges[i],maskranges[i+1]
			print 'selected',nvals,nselect
			selected=selected+ nselect
			select=np.zeros([nvals],dtype='bool')				
			select[:nselect]=1
			np.random.shuffle(select)
			k=0
			for j in range(nens):
				if mask[j]:
					masksubset[j]=select[k]
					k=k+1
			plt.scatter(pcs[mask,0][select],pcs[mask,1][select],color=colors[i])

		#plt.scatter(pcs[masksubset,0],pcs[masksubset,1],color='r')
		print 'number of ensembles selected:',selected
		plt.xlabel('pc0')
		plt.ylabel('pc1')
		plt.title('scatter of pc0 vs pc1')
		plt.savefig(output_dir+'/scatter.png')
		# Calculate EOFs
		solver2=Eof(data[masksubset,:])
		pcs2=solver2.pcs(npcs=neofs,pcscaling=1)
		eofs2=solver2.eofs(neofs=neofs)
		varfrac2=solver2.varianceFraction(neigs=neofs)
		print 'calculated EOFs!'
		print "plotting mean and stdev of subsample..."
		plot_region_pnw(data[masksubset,:].mean(0),lat_coord,lon_coord,0,-1,0,-1,'mean_subset3D',data[masksubset,:].mean())
		plot_region_pnw(data[masksubset,:].std(0),lat_coord,lon_coord,0,-1,0,-1,'stdev_subset3D',data[masksubset,:].std())
		plot_region_pnw(data.mean(0)-data[masksubset,:].mean(0),lat_coord,lon_coord,0,-1,0,-1,'mean_diff3D',-1)
		plot_region_pnw(data.std(0)-data[masksubset,:].std(0),lat_coord,lon_coord,0,-1,0,-1,'stdev_diff3D',-1)
		print 'plotting EOFs for subsample...'
		for i in range(neofs):
			print 'EOF',i
			plt.clf()
			plot_region_pnw(eofs2[i,:],lat_coord,lon_coord,0,-1,0,-1,'EOF'+str(i)+'_subset3D',varfrac2[i])
	
		print "Files chosen:"
		fout=open(output_dir+'/umids_chosen.txt','w')
		for i in np.array(filenames)[masksubset]:
			fout.write(os.path.basename(i).split('_')[1]+'\n') # Print Umid for each dump selected
		fout.close()

		# Choose a random subset as a comparison to our other method
		#
# 		print "\nRadom subsampling..."
# 		maskrandom=np.zeros([nens],dtype='bool')
# 		maskrandom[:51]=1
# 		np.random.shuffle(maskrandom)
# 		# Calculate EOFs
# 		solver3=Eof(data[maskrandom,:])
# 		pcs3=solver3.pcs(npcs=neofs,pcscaling=1)
# 		eofs3=solver3.eofs(neofs=neofs)
# 		varfrac3=solver3.varianceFraction(neigs=neofs)
# 		print 'calculated EOFs'
# 		print "plotting mean and stdev"
# 		plot_region_pnw(data[masksubset,:].mean(0),lat_coord,lon_coord,0,-1,0,-1,'mean_subsetrand',data[masksubset,:].mean())
# 		plot_region_pnw(data[masksubset,:].std(0),lat_coord,lon_coord,0,-1,0,-1,'stdev_subsetrand',data[masksubset,:].std())
# 		plot_region_pnw(data.mean(0)-data[masksubset,:].mean(0),lat_coord,lon_coord,0,-1,0,-1,'mean_diffrand',-1)
# 		plot_region_pnw(data.std(0)-data[masksubset,:].std(0),lat_coord,lon_coord,0,-1,0,-1,'stdev_diffrand',-1)
# 		print "plotting EOFs"
# 		for i in range(neofs):
# 			print 'EOF',i
# 			plt.clf()
# 			plot_region_pnw(eofs3[i,:],lat_coord,lon_coord,0,-1,0,-1,'EOF'+str(i)+'_subsetrand',varfrac2[i])
# 		print 'plotting as PC scatter plot'
# 		plt.clf()
# 		plt.scatter(pcs[:,0],pcs[:,1])
# 		plt.scatter(pcs[maskrandom,0],pcs[masksubset,1],color='r')
# 		plt.xlabel('pc0')
# 		plt.ylabel('pc1')
# 		plt.title('scatter of pc0 vs pc1')
# 		plt.savefig(output_dir+'/scatter_rand.png')