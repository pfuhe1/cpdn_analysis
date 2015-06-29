#!/usr/bin/env python

###############################################################################
# Program : wah_extract.py
# Author  : Neil Massey
# Date	  : 04/06/14
# Purpose : download w@h output files from a list of urls and extract the 
#           requested fields into separate netCDF files
# Requires: scipy.io.netcdf
###############################################################################

import sys, getopt, os, ast
import urllib, tempfile, zipfile, shutil
import numpy, numpy.ma
import math
import gzip
#from scipy.io.netcdf import netcdf_file
from netcdf_file import * 		# uncomment this line if you don't have scipy
from conv_rot_grid import rot2glob, glob2rot

###############################################################################

def strip(s):
	return s.strip()

###############################################################################

def read_urls(urls_file):
	fh = gzip.open(urls_file, 'r')
	urls = fh.readlines()
	fh.close()
	urls = map(strip, urls)
	return urls
	
def read_urls_gpfs(urls_file,gpfsdir):
	fh = gzip.open(urls_file, 'r')
	urls=[]
	for url in gzip.open(urls_file, 'r'):
		if url.strip():
			urls.append(gpfsdir+url.strip().split('results')[1])
	print urls
	return urls

###############################################################################

def get_umid_boinc_year(url):
	# get the umid, boinc id and year from url
	# get the last file
	lf = url.split("/")[-1][:-4]
	boinc = url.split("/")[-2]
	umid = lf.split("_")[2]
	year = lf.split("_")[3]
	return umid, boinc, year

###############################################################################

def get_coord_string(coord_list):
	coord_string = ""
	for p in range(0,4):
		V = str(coord_list[p])
		if (p == 1 or p == 3):
			if V[0] == "-":
				coord_string += "S" + V[1:] + "_"
			else:
				coord_string += "N" + V[0:] + "_"
		else:
			if V[0] == "-":
				coord_string += "W" + V[1:] + "_"
			else:
				coord_string += "E" + V[0:] + "_"
	return coord_string[:-1]

###############################################################################

def get_output_field_name(field):
	fname = field[1]
	if field[2] != []:
		coord_string = get_coord_string(field[2])
		fname += "_" + coord_string
	if field[3] != 'all':
		fname += "_" + field[3]
	return fname

###############################################################################

def make_directories(output_dir, year, boinc, field_list, n_valid):
	# make the output top directory
	c_path = output_dir
	download=False
	if not os.path.exists(c_path):
		os.mkdir(c_path)
	# make the year directory
	c_path += "/" + year
	if not os.path.exists(c_path):
		os.mkdir(c_path)
	# make the boinc / umid directory
	c_path += "/" + boinc
	if not os.path.exists(c_path):
		os.mkdir(c_path)
	# finally - make the field directories
	for field in field_list:
		# make the pattern directory (i.e. pattern for global / regional / regional daily etc.)
		f_path = c_path + "/" + field[0]
		if not os.path.exists(f_path):
			os.mkdir(f_path)
		# if there is not subsetting or processing of the field then the name is 
		# just the field name
		f_path = f_path + "/" + get_output_field_name(field)
			
		if not os.path.exists(f_path):
			os.mkdir(f_path)
		if len(os.listdir(f_path)) != n_valid:
			download=True

	return c_path, download
	
def make_directories2(output_dir, field_list, n_valid):
	# make the output top directory
	c_path = output_dir
	download=False
	if not os.path.exists(c_path):
		os.makedirs(c_path)
	# finally - make the field directories
	for field in field_list:
		# make the pattern directory (i.e. pattern for global / regional / regional daily etc.)
		f_path = c_path + "/" + field[0]
		if not os.path.exists(f_path):
			os.mkdir(f_path)
		# if there is not subsetting or processing of the field then the name is 
		# just the field name
		f_path = f_path + "/" + get_output_field_name(field)
			
		if not os.path.exists(f_path):
			os.mkdir(f_path)
		if len(os.listdir(f_path)) != n_valid:
			download=True

	return c_path, download

###############################################################################

def get_idx(point, array):
	D = array[1] - array[0]
	S = array[0]
	I = int(0.5+(point - S) / D)
	return I

###############################################################################

def subset_dimensions(in_dimensions, field, plon, plat):
	fs = field[2]
	# fs has elements [lon_l, lat_l, lon_r, lat_r]
	out_dims = []
	subset_dims = []
	# check whether we need to translate from global to rotated grid
	if plon == 0.0 and plat == 90.0:
		if fs != []:
			rlon_s = fs[0]
			rlon_e = fs[2]
			rlat_s = fs[1]
			rlat_e = fs[3]
	else:
		if fs != []:
			rlon_s, rlat_s = glob2rot(fs[0], fs[1], plon, plat)
			rlon_e, rlat_e = glob2rot(fs[2], fs[3], plon, plat)
	# get the lon / lat indices
	remap_data = False
	for d in in_dimensions:
		if "longitude" in d[0]:
			if fs == []:
				lon_idx_s = 0
				lon_idx_e = d[1].shape[0]
			else:
				lon_idx_s = get_idx(rlon_s, d[1])
				lon_idx_e = get_idx(rlon_e, d[1]) + 1
			if lon_idx_s < 0:
				remap_data = True
			if lon_idx_e > d[1].shape[0]:
				remap_data = True
			# if the data needs remapping
			if remap_data:
				d_len = d[1].shape[0]				# get the length
				d_len_d2 = d_len / 2				# get the length div 2
				lon_idx_s += d_len_d2				# remap the start index
				lon_idx_e += d_len_d2				# remap the end index
				d_new = numpy.zeros([d_len], 'f')	# create a new dimension variable
				d_new[0:d_len_d2] = d[1][d_len_d2:d_len]	# copy the right hand half to the left hand
				d_new[d_len_d2:d_len] = d[1][0:d_len_d2]	# copy the left hand half to the right hand
				d[1] = d_new						# reassign
			# if processing data then do not append dimension to output dimensions
			if field[3] == "all":
				out_dims.append([d[0], d[1][lon_idx_s:lon_idx_e]])
			subset_dims.append([d[0], d[1][lon_idx_s:lon_idx_e]])

		elif "latitude" in d[0]:
			if fs == []:
				lat_idx_s = 0
				lat_idx_e = d[1].shape[0]
			else:
				lat_idx_s = get_idx(rlat_s, d[1])
				lat_idx_e = get_idx(rlat_e, d[1]) +1
			if lat_idx_e < lat_idx_s:
				temp = lat_idx_e
				lat_idx_e = lat_idx_s
				lat_idx_s = temp
			if lat_idx_s < 0:
				lat_idx_s = 0
			if lat_idx_e > d[1].shape[0]:
				lat_idx_e = d[1].shape[0]
			# if processing data then do not append dimension
			if field[3] == "all":
				out_dims.append([d[0], d[1][lat_idx_s:lat_idx_e]])
			subset_dims.append([d[0], d[1][lat_idx_s:lat_idx_e]])
		else:
			out_dims.append([d[0], d[1]])
			subset_dims.append([d[0], d[1]])
	# if processing data then need a "pt" dimension
	if field[3] != "all":
		out_dims.append(["pt", numpy.array([1.0])])
	return out_dims, subset_dims, [lon_idx_s, lat_idx_s, lon_idx_e, lat_idx_e], remap_data

###############################################################################

def calc_grid_box_area(lon, lat, lon_d, lat_d, plon, plat):
	# calculate the area of a grid box.  This is only an approximation as it
	# assumes the grid box is square, rather than the polygonal shape the
	# rotated boxes actually are.
	lon1_rad, lat1_rad = rot2glob(lon-lon_d,lat+lat_d, plon, plat)
	lon2_rad, lat2_rad = rot2glob(lon+lon_d,lat-lat_d, plon, plat)
		
   	area = math.fabs(lon2_rad-lon1_rad) *\
		   math.fabs(math.sin(lat2_rad) - math.sin(lat1_rad))
	return area

###############################################################################

def aamean(var_ma_data, plon, plat, subset_dims):
	# calculate the weights
	if plon == 0.0 and plat == 90.0:
		# just a regular grid so the weights are simply the cosine of the latitude
		for dim in subset_dims:
			if "latitude" in dim[0]:
				weights = numpy.cos(numpy.radians(dim[1]))
		lon_means = numpy.ma.mean(var_ma_data, axis=3)
		means = numpy.ma.average(lon_means, weights=weights, axis=2)
	else:
		# identify the latitude and longitude dimensions and get the coordinates
		for dim in subset_dims:
			if "latitude" in dim[0]:
				lats = dim[1]
			if "longitude" in dim[0]:
				lons = dim[1]
		# calculate the weights
		weights = numpy.zeros([lats.shape[0], lons.shape[0]], 'f')
		lon_d = 0.5 * (lons[1] - lons[0])
		lat_d = 0.5 * (lats[1] - lats[0])
		for j in range(0, lats.shape[0]):
			for i in range(0, lons.shape[0]):
				# weights are area grid box
				weights[j,i] = calc_grid_box_area(lons[i], lats[j], lon_d, lat_d, plon, plat)
		# multiply through each timestep / level by the weights
		weighted_data = var_ma_data[:,:] * weights
		# calculate the means - sums first
		means = numpy.ma.sum(numpy.ma.sum(weighted_data, axis=3), axis=2)
		# calculate the means
		for t in range(0, var_ma_data.shape[0]):
			for z in range(0, var_ma_data.shape[1]):
				# only want to use the sum of the weights of the non missing values
				idx = ~var_ma_data[t,z].mask
				# now calculate the sum of the weights
				sum_weights = numpy.sum(weights[idx])		
				means[t,z] =  means[t,z] / sum_weights
	return means

###############################################################################

def process_data(var_in_data, process, mv, plon, plat, subset_dims, valid_min, valid_max):
	# mask array before processing
	# mask is either the missing value or where the data is outside the
	# valid range
	var_ma_data = numpy.ma.masked_where((var_in_data == mv) | \
										(var_in_data < valid_min) | \
										(var_in_data > valid_max), var_in_data)
	# do the various processes based on what has been requested
	if process == "min":
		out_data = numpy.ma.min(numpy.ma.min(var_ma_data, axis=2), axis=2)
	elif process == "max":
		out_data = numpy.ma.max(numpy.ma.max(var_ma_data, axis=2), axis=2)
	elif process == "mean":
		out_data = aamean(var_ma_data, plon, plat, subset_dims)
	elif process == "sum":
		out_data = numpy.ma.sum(numpy.ma.sum(var_ma_data, axis=2), axis=2)
	if len(out_data.shape) != len(var_in_data.shape)-1:
		out_data = out_data.reshape(out_data.shape[0], out_data.shape[1], 1)
	return out_data

###############################################################################

def get_missing_value(attrs):
	if "missing_value" in attrs.keys():
		mv = attrs["missing_value"]
	if "_FillValue" in attrs.keys():
		mv = attrs["_FillValue"]
	return mv

###############################################################################

def get_rotated_pole(attrs, nc_in_file):
	# get the rotated pole longitude / latitude (for calculating weights)
	if "grid_mapping" in attrs.keys():
		grid_map_name = attrs["grid_mapping"]
		grid_map_var = nc_in_file.variables[grid_map_name]	
		plon = grid_map_var._attributes["grid_north_pole_longitude"]
		plat = grid_map_var._attributes["grid_north_pole_latitude"]
	else:
		plon = 0.0
		plat = 90.0
	return plon, plat

###############################################################################

def extract_netcdf_var(nc_in_file, nc_out_file, field):
	nc_in_var = nc_in_file.variables[field[1]]
	in_dimensions = []
	process = field[3]
	v_min = field[4]
	v_max = field[5]
	# now copy the dimensions from input netcdf
	for d in nc_in_var.dimensions:
		# get the input dimension and the data
		dim_in_var = nc_in_file.variables[d]
		dim_in_data = dim_in_var[:]
		in_dimensions.append([d, dim_in_data])

	# get the rotated pole definition	
	plon, plat = get_rotated_pole(nc_in_var._attributes, nc_in_file)
	# subset the dimensions to create the out_dimensions
	out_dims, subset_dims, lon_lat_idxs, remap_data = subset_dimensions(in_dimensions, field, plon, plat)
	# if the longitude and latitude indexes are < 0 then we need to remap the data so that 0deg is
	# in the middle of the field, not at the beginning
	if remap_data:
		in_data = nc_in_var[:,:,lon_lat_idxs[1]:lon_lat_idxs[3],:]	# get the input data - can subset latitude early
		new_data = numpy.zeros(in_data.shape, 'f') # create a new store
		d_len = in_data.shape[3]					# get the longitude length
		d_len_d2 = d_len / 2						# lon length div 2
		new_data[:,:,:,0:d_len_d2] = in_data[:,:,:,d_len_d2:d_len]	# copy right hand half to left hand
		new_data[:,:,:,d_len_d2:d_len] = in_data[:,:,:,0:d_len_d2]  # copy left hand half to right hand
		var_out_data = new_data[:,:,:,lon_lat_idxs[0]:lon_lat_idxs[2]] # get the subset data
	else:		
		var_out_data = nc_in_var[:,:,lon_lat_idxs[1]:lon_lat_idxs[3], lon_lat_idxs[0]:lon_lat_idxs[2]]
	
	# if the data is going to be processed then do the processing here
	if process != "all":
		# get the missing value first
		mv = get_missing_value(nc_in_var._attributes)
		# get the rotated pole
		var_out_data = process_data(var_out_data, process, mv, plon, plat, subset_dims, v_min, v_max)
		
	for d in out_dims:
		# create the output dimension and variable
		nc_out_file.createDimension(d[0], d[1].shape[0])
		dim_out_var = nc_out_file.createVariable(d[0], d[1].dtype, (d[0],))
		# assign the output variable data and attributes from the input
		if d[0] in nc_in_file.variables.keys():
			dim_in_var = nc_in_file.variables[d[0]]
			dim_out_var._attributes = dim_in_var._attributes
		elif d[0] == "pt":
			# if it's the "pt" dimension then create an attribute indicating the domain of the
			# mean-ed / max-ed / min-ed variable
			dom_str = ""
			if field[2] == []:
				dom_str = "global  "
			else:
				for i in range(0, 4):
					dom_str += str(field[2][i]) + ", "
			dim_out_var._attributes["domain"] = dom_str[:-2]
		dim_out_var[:] = d[1][:]
		
	# create the variable
	out_dim_names = [d[0] for d in out_dims]
	nc_out_var = nc_out_file.createVariable(field[1], var_out_data.dtype, out_dim_names)
	# assign the attributes
	nc_out_var._attributes = nc_in_var._attributes
	# remove the grid mapping and coordinates from the dictionary if they exist and process is not all
	if process != "all":
		if "grid_mapping" in nc_out_var._attributes:
			del nc_out_var._attributes["grid_mapping"]
		if "coordinates" in nc_out_var._attributes:
			del nc_out_var._attributes["coordinates"]
		if "cell_method" in nc_out_var._attributes:
			nc_out_var._attributes["cell_method"] += ", area: " + process + " "
	# assign the data
	nc_out_var[:] = var_out_data
	# check for rotated pole and copy variable if it exists
	if "grid_mapping" in nc_out_var._attributes and len(out_dims) == 4:
		grid_map_name = nc_out_var._attributes["grid_mapping"]
		grid_map_var = nc_in_file.variables[grid_map_name]
		grid_map_out_var = nc_out_file.createVariable(grid_map_name, 'c', ())
		grid_map_out_var._attributes = grid_map_var._attributes
		# get the global longitude / global latitude vars
		coords = (nc_out_var._attributes["coordinates"]).split(" ");
		global_lon_var = nc_in_file.variables[coords[0]]
		global_lat_var = nc_in_file.variables[coords[1]]
		global_lon_data = global_lon_var[lon_lat_idxs[1]:lon_lat_idxs[3], lon_lat_idxs[0]:lon_lat_idxs[2]]
		global_lat_data = global_lat_var[lon_lat_idxs[1]:lon_lat_idxs[3], lon_lat_idxs[0]:lon_lat_idxs[2]]
		# create the global latitude / global longitude variables
		out_global_lon_var = nc_out_file.createVariable(coords[0], global_lon_data.dtype, (out_dims[2][0], out_dims[3][0]))
		out_global_lon_var[:] = global_lon_data
		out_global_lat_var = nc_out_file.createVariable(coords[1], global_lat_data.dtype, (out_dims[2][0], out_dims[3][0]))
		out_global_lat_var[:] = global_lat_data
		out_global_lon_var._attributes = global_lon_var._attributes
		out_global_lat_var._attributes = global_lat_var._attributes
		
###############################################################################

def extract_netcdf(zf, zf_file, base_path, field, temp_dir):
	# extract the file to a tempfile
	tmp_file_path = temp_dir + "/" + zf_file
	zf.extract(zf_file, temp_dir)
	# open as netCDF to a temporary file
	nc_file = netcdf_file(tmp_file_path,'r')
	# create the output netCDF file
	o_field_name = get_output_field_name(field)
	out_name = base_path + "/" + field[0] + "/" + o_field_name + "/" + zf_file[:-3] + "_" + o_field_name + ".nc"
	# check whether it exists
	if os.path.exists(out_name):
		return
	# get the variable from the input
	if not field[1] in nc_file.variables.keys():
		print "Could not extract field: " + field[1] + " from file: " + zf_file
		return
	out_ncf = netcdf_file(out_name, "w")
	extract_netcdf_var(nc_file, out_ncf, field)
	out_ncf.close()
	nc_file.close()
	# delete the temp directory and its contents
	os.remove(tmp_file_path)

###############################################################################

def extract(url, field_list, output_dir, temp_dir, n_valid,gpfs):
	# fetch the url to a temp (zip) file using urllib
	try:
		# get the umid, boinc id and year from the url
		umid, boinc, year = get_umid_boinc_year(url)
		# make directories to hold the output
		if gpfs:
			base_path, download = make_directories2(output_dir, field_list, n_valid)
		else:
			base_path, download = make_directories(output_dir, year, boinc, field_list, n_valid)
		# check whether this has already been processed
		if not download:
			return
		zf_fh = tempfile.NamedTemporaryFile(mode='w', delete=False,dir=temp_dir)
		if gpfs:
			zf = zipfile.ZipFile(url,'r')
			#urlret=shutil.copyfile(url,zf_fh.name)
		else:
		# open the zip
			urlret = urllib.urlretrieve(url, zf_fh.name)
			zf = zipfile.ZipFile(zf_fh.name,'r')
		# list the zip contents
		zf_list = zf.namelist()
		# check to see which files match the pattern
		for field in field_list:
			pattern = field[0]
			for zf_file in zf_list:
				if pattern in zf_file:
					extract_netcdf(zf, zf_file, base_path, field, temp_dir)
		if not gpfs:
			os.remove(zf_fh.name)
	except Exception,e:
		print "Could not extract url: " + url
		print e

###############################################################################

def usage():
	print "Usage : wah_extract.py --urls|-u  --fields|-f --n_valid|-n --output_dir|-o"
	print "      : fields has the format : [pattern,field_name,[region],process,valid_min,valid_max]"
	print "      : where [region] = [lon_l,lat_l,lon_r,lat_r]"
	print "      :        process = min|max|mean|sum|all"

###############################################################################

if __name__ == "__main__":
	urls_file = ""
	fields = ""
	output_dir = ""
	n_valid = 0
	gpfs=True # Use this if the gpfs is mounted
	gpfs_dir='/gpfs/projects/cpdn/storage/boinc/upload'
	opts, args = getopt.getopt(sys.argv[1:], 'u:p:f:o:n:', 
							   ['urls=', 'fields=', 'output_dir=', 'n_valid='])

	for opt, val in opts:
		if opt in ['--urls', '-u']:
			urls_file = val
		if opt in ['--fields', '-f']:
			fields = val
		if opt in ['--n_valid', '-n']:
			n_valid = int(val)
		if opt in ['--output_dir', '-o']:
			output_base = val
				
	if urls_file == "" or fields == "" or output_base == "" or n_valid == 0:
		usage()
		exit()
		
	# split the field list up
	field_list = ast.literal_eval(fields)
	for field in field_list:
		if len(field) != 6:
			print "Fields argument not formatted correctly"
			print field
			usage()
			exit()
	
	# get the urls from the file
	#if gpfs:
	#	urls=read_urls_gpfs(urls_file,gpfs_dir)
	#else:
	urls = read_urls(urls_file)
	
	# create a temporary directory - do we have permission?
	temp_dir = tempfile.mkdtemp(dir='/gpfs/projects/cpdn/scratch/cenv0437/')
	
	# download each url
	for u in urls:
		print u
		if gpfs:
			try:
				u=u.split('results')[1]
				gpfs_dir='/gpfs/projects/cpdn/storage/boinc/upload'
			
				output_dir=output_base+os.path.dirname(u) #Add subdirectories for batch and workunit
				u=gpfs_dir+u # Actual location of the file
			except:
				continue
		extract(u, field_list, output_dir, temp_dir, n_valid,gpfs)
		
	# remove the temporary directory
	shutil.rmtree(temp_dir)
