# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 10:39:03 2023

@author: 11913
"""

import numpy as np
import math as mt
import subprocess,os
#import radio correction as rrd
#importing correlation correction
import string 
import xml.etree.ElementTree as ET
import pdb
from scipy import signal

import remove_corr_bias as rcb
import scipy.io as sio
import h5py

import subprocess







#croppying the coherence;

def Cropping_coh (coh_2d, lat_max, lat_min, lon_max, lon_min):
    [M,N] = np.shape(coh_2d)
    lat_axis = np.linspace(lat_max, lat_min, M)
    lon_axis = np.linspace(lon_min, lon_max, N)
    [LON, LAT] = np.meshgrid(lon_axis, lat_axis)
    
    LON[np.isnan(coh_2d) or coh_2d ==0] = np.nan
    LAT[np.isnan(coh_2d) or coh_2d ==0] = np.nan
    
    
    LAT_1d = LAT.flatten()
    LON_1d = LON.flatten()
    lat_valid_max = np.nanmax(LAT_1d)
    lat_valid_min = np.nanmin(LAT_1d)
    lon_valid_max = np.nanmax(LON_1d)
    lon_valid_min = np.nanmin(LON_1d)
    
    line_min = np.where(lat_axis == lat_valid_max)[0]
    line_max = np.where(lat_axis == lat_valid_min)[0]
    pixel_min = np.where(lon_axis == lon_valid_min)[0]
    pixel_max = np.where(lon_axis == lon_valid_max)[0]
    
    lat_ax_valid = lat_axis[line_min: line_max+1]
    lon_ax_valid = lon_axis[pixel_min: pixel_max+ 1]
    
    coh_valid = coh_2d[line_min:line_max+1, pixel_min: pixel_max+ 1]
    
    
    return coh_valid, lat_ax_valid, lon_ax_valid


def Extract_PreProc(subdir_path, outpath='cor_input.mat'):
    
    ResFileDir = subdir_path
    
    
    
    
    if 'forFSH' not in subdir_path:
      ResFileDir = os.path.join(subdir_path,'forFSH')
      
    if len(os.listdir(ResFileDir))!=7 and len(os.listdir(ResFileDir))!=8:
        print('this folder is invalid:' + ResFileDir)
        return
      
    xmlfile_2_read = [f for f in os.listdir(ResFileDir) if f.endswith('Proc.xml')][0]
    
    
    path_2_xml =  os.path.join(subdir_path, xmlfile_2_read)
    xml_tree = ET.parse(path_2_xml)
    root = xml_tree.getroot()
    root_tag = root.tag
    
    
    #reading the relevant parameters;
    try:
      range_pixel_res = float(root.findall("./master/instrument/range_pixel_size")[0].text)
    except:
      range_pixel_res = float(root.findall("./reference/instrument/range_pixel_size")[0].text)
    try:
      llambda = float(root.findall("./master/instrument/radar_wavelength")[0].text)
    except:
      llambda = float(root.findall("./reference/instrument/radar_wavelength")[0].text)
    try:
      first_range = float(root.findall("./runTopo/inputs/range_first_sample")[0].text)
    except:
      first_range = float(root.findall("./reference/frame/starting_range")[0].text)
    try:
      num_range_bin = int(root.findall("./frame/inputs/width")[0].text)
    except:
      num_range_bin = int(root.findall("./reference/frame/number_of_samples")[0].text)
    center_range = first_range + (num_range_bin/2-1)*range_pixel_res
    try:
      incid_angle = float(root.findall("./master/instrument/incidence_angle")[0].text)
    except:
      incid_angle = float(root.findall("./reference/instrument/incidence_angle")[0].text)
      baseline_top = float(root.findall("./baseline/perp_baseline_top")[0].text)
      baseline_bottom = float(root.findall("./baseline/perp_baseline_bottom")[0].text)
      baseline = (baseline_bottom+baseline_top)/2
    try:
      sensor = root.findall("./master/platform/mission")[0].text
    except:
      sensor = root.findall("./reference/platform/mission")[0].text
  
    if sensor == 'ALOS':
      numLooks = 20
    elif sensor == 'ALOS2':
      numLooks = 30
    else:
      raise Exception("invalid sensor: supported sensors include ALOS and ALOS-2 only")
      
      
    #reading the geocoding information
    xmlfile_2_read = [f for f in os.listdir(ResFileDir) if f.endswith('topophase.cor.geo.xml')][0]
    path_2_xml = os.path.join(ResFileDir, xmlfile_2_read)
    
    
    tree = ET.parse(path_2_xml)
    root = tree.getroot()
    delta_array = np.array([])
    start_array = np.array([])
    size_array = np.array([], dtype=np.int32)
    
    for size in root.iter('property'):
       if size.get('name') == 'size':
           size_array = np.append(size_array, int(size.find('value').text))
    for delta_val in root.iter('property'):
       if delta_val.get('name') == 'delta':
           delta_array = np.append(delta_array, float(delta_val.find('value').text))
    for start_val in root.iter('property'):
       if start_val.get('name') == 'startingvalue':
           start_array = np.append(start_array, float(start_val.find('value').text))
    end_array = start_array + size_array * delta_array
    north = max(start_array[1],end_array[1])
    south = min(start_array[1],end_array[1])
    east = max(start_array[0],end_array[0])
    west = min(start_array[0],end_array[0])
    coords = [north, south, west, east]
    
    #cropping the image here
    #....
    geo_nwidth = size_array[0]
    geo_nlines = size_array[1]
    corner_lat = north
    corner_lon = west
    step_lat = delta_array[1]
    step_lon = delta_array[0]
    
    
    
    #reading the resampleOnlyImage files
    #xmlfilet = [f for f in os.listdir(ResFileDir) if f.endswith('resampOnlyImage.amp.geo.xml')][0]
    #xmlfile = os.path.join(ResFileDir, xmlfilet)
    #tree = ET.parse(xmlfile)
    #root = tree.getroot()
    #delta_array = np.array([])# Simon Kraatz, UMass Amherst
# April 28, 2020
    #start_array = np.array([])
    #size_array = np.array([], dtype=np.int32)
    
    #reading the geolocated cor and amplitude files;
    finnn = os.path.join(ResFileDir, 'topophase.cor.geo')
    fid_cor = open(finnn, "rb")
    #might get failure
    #fid_cor = open(directory + "int_"+date1+"_"+date2+"/topophase.cor.geo", "rb")
    cor_file = np.fromfile(fid_cor, dtype=np.dtype('<f'))
##    pdb.set_trace()
    corr = cor_file.reshape(2*geo_nwidth, -1, order='F')
    corr = corr[:,0:geo_nlines]
    corr_mag = corr[geo_nwidth:2*geo_nwidth,:]
    
    
    #reading amplitude files
    finnn = os.path.join(ResFileDir, 'resampOnlyImage.amp.geo')
    fid_amp = open(finnn, "rb")
    #fid_amp = open(directory + "int_"+date1+"_"+date2+"/resampOnlyImage.amp.geo", "rb")
    amp_file = np.fromfile(fid_amp, dtype=np.dtype('<f'))
    inty = amp_file.reshape(2*geo_nwidth, -1, order='F')
    inty = inty[:,0:geo_nlines]
    inty1 = inty[::2,:]
    inty2 = inty[1::2,:]
    
    #handling the outliers
    inty1 = np.power(inty1,2)                   # Hardcoded based on 2 range looks and 10 azimuth looks
    inty2 = np.power(inty2,2)

    inty1[inty1 <= 0] = np.NaN
    inty2[inty2 <= 0] = np.NaN
    corr_mag[corr_mag <= 0] = np.NaN
    
    #define the noise level
    if sensor == 'ALOS':
        if root_tag == 'insarProc':
            ####### ALOS thermal noise level (insarApp)
            N1 = 55.5**2
            N2 = 55.5**2
        elif root_tag == 'stripmapProc':
            ####### ALOS thermal noise level (stripmapApp)
            N1 = (55.5/81)**2
            N2 = (55.5/81)**2
                    
        else:
            raise Exception("invalid *Proc.xml file!!!")
    elif sensor == 'ALOS2':
        if root_tag == 'insarProc':
            ####### ALOS-2 thermal noise level (insarApp)
            N1 = 25848**2
            N2 = 25848**2
        elif root_tag == 'stripmapProc':
            ####### ALOS-2 thermal noise level (stripmapApp)
            N1 = 18114**2
            N2 = 18114**2
        else:
            raise Exception("invalid *Proc.xml file!!!")
    else:
        raise Exception("invalid sensor: supported sensors include ALOS and ALOS-2 only")
            
        
    #thermal correction:
    S1 = inty1 - N1
    g_th_1 = np.zeros(S1.shape)
    g_th_1[S1>N1] = np.sqrt(S1[S1>N1] / (S1[S1>N1] + N1))
    g_th_1[np.isnan(S1)] = np.NaN
    g_th_1[S1 <= N1] = np.NaN

    S2 = inty2-N2
    g_th_2 = np.zeros(S2.shape)
    g_th_2[S2>N2] = np.sqrt(S2[S2>N2] / (S2[S2>N2] + N2))
    g_th_2[np.isnan(S2)] = np.NaN
    g_th_2[S2 <= N2] = np.NaN

    g_th = g_th_1 * g_th_2

    corr_mag[corr_mag<0] = 0
    corr_mag[corr_mag>1] = 1
    #corr_mag = rcb.remove_corr_bias(corr_mag,numLooks)
    corr_mag[corr_mag<0] = 0

    corr_vs = corr_mag / g_th
    
    
    pi = mt.pi
    #correcting the geocoding shifts
    #....
    
    
    
    #correcting the geometric deccorelation by a simple spatially invariant equations
    gamma_factor = 1 - (2 * mt.fabs(baseline) * mt.cos(incid_angle / 180 * pi) * range_pixel_res / mt.sin(incid_angle / 180 * pi) / llambda / center_range)
    corr_vs = corr_vs/gamma_factor
    corr_vs[corr_vs>1] = 1
    corr_vs[corr_vs<0] = 0
    
    flag_grad_cal = 0
    
    #It might be useful for ALOS2
    if flag_grad_cal == 1:
        x = np.linspace(north, south, geo_nlines)
        y = np.linspace(west, east, geo_nwidth)
        [X,Y] = np.meshgrid(x,y)
        A = np.vstack([X[~np.isnan(corr_vs)], Y[~np.isnan(corr_vs)], np.ones(np.size(corr_vs[~np.isnan(corr_vs)]))]).T
        coeff = np.linalg.lstsq(A, corr_vs[~np.isnan(corr_vs)])[0]
        corr_vs = corr_vs - X*coeff[0] - Y*coeff[1]
        corr_vs[corr_vs>1] = 1
        corr_vs[corr_vs<0] = 0
    
    
    kz = -2 * pi * 2 / llambda / center_range / mt.sin(incid_angle/180*pi) * baseline
    kz = mt.fabs(kz)
    
    file_name = os.path.join(ResFileDir, outpath)
    if '.mat' not in file_name:
        file_name = file_name + '.mat'
    sio.savemat(file_name, {'corr_vs':corr_vs,'kz':kz,'coords':coords})
    
    return  1 
    
    
    
    # sd
    


def preproc_batch (workdir):
    files = os.listdir(workdir)
    output_name = 'cor_input.mat'
    for file in files:
        subdir_loop = os.path.join(workdir, file)
        if os.path.isdir(subdir_loop) and file.find('_')!=-1:
            exe_path = os.path.join(subdir_loop,'isce_proc_dir')
            if os.path.exists(exe_path):
                files_sub = os.listdir(exe_path)
                for file_sub in files_sub:
                   judge_dir = os.path.join(exe_path, file_sub,'forFSH')
                   if os.path.isdir(judge_dir):
                       subdir_loop = os.path.join(exe_path, file_sub,'forFSH')
                       #print(subdir_loop)
                       Extract_PreProc(subdir_loop, output_name)


def write_to_txt(content, file_path):
    file = open(file_path, 'a')
    file.write(content)
    file.close()
        
def runCmd(cmd):
    out = subprocess.getoutput(cmd)
    return out


def gdal_mosaic_each_dir(workdir,backscatter_source,year, alos_or_alos2 ):
    #flag alos_or_alos2 to indicate the backscatter generated based on alos or alos2: 1 for alos2, 0 for alos
    matfile = os.path.join(workdir,'Height_inversion.mat')
    if os.path.exists(os.path.join(workdir,'Height_inversion_rh100.mat')):
        matfile = os.path.join(workdir,'Height_inversion_rh100.mat')
    if os.path.exists(os.path.join(workdir,'Height_inversion_rh98.mat')):
        matfile = os.path.join(workdir,'Height_inversion_rh98.mat')
        
    if os.path.exists(matfile):
        try:
            conts = h5py.File(matfile)
            #conts = sio.loadmat(matfile)
            coords = conts['GEDI_coords']
            latmax = coords[0]
            latmin = coords[1]
            lonmin = coords[2]
            lonmax = coords[3]    
        except:
            conts = sio.loadmat(matfile)
            coords = conts['GEDI_coords'][0]
            latmax = coords[0]
            latmin = coords[1]
            lonmin = coords[2]
            lonmax = coords[3] 
    else:
        print('This directory is not valid! Please check it!')
        return None
            
        #coords = conts['GEDI_coords'][0]
        #latmax = coords[0]
        #latmin = coords[1]
        #lonmin = coords[2]
        #lonmax = coords[3] 
        
    lat_prefix = 'N'
    lon_prefix = 'W'
        
    y_idx_max = int(np.ceil(latmax))
    y_idx_min = int(np.ceil(latmin))
    x_idx_max = int(np.floor(lonmax))
    x_idx_min = int(np.floor(lonmin))
        
        
    backscatter_output = 'backscatter_mosaic.tif'
        
    path_to_output = os.path.join(workdir,backscatter_output)
    #if os.path.exists(path_to_output):
        #return 0
        #os.remove(path_to_output)
        
    tiff_lists = 'tiff_lists.txt';
    path_to_lists =os.path.join(workdir,tiff_lists)
    if os.path.exists(path_to_lists):
        os.remove(path_to_lists)
                
        
    numfiles = (y_idx_max-y_idx_min+1) *(x_idx_max - x_idx_min+1) 
    for i in range(y_idx_min,y_idx_max+1):
        lat_prefix = 'N'
        if i<0:
            lat_prefix = 'S'
        y_idx_loop = np.abs(i)
            
        for j in range(x_idx_min,x_idx_max+1):
            lon_prefix = 'E'
            if j<0:
                lon_prefix = 'W'
                    
            x_idx_loop = np.abs(j)
            #alos
            if alos_or_alos2 == 0:
                file_name = lat_prefix + '{0}'.format(y_idx_loop) + lon_prefix + '{0}'.format(x_idx_loop)+'_' + '20{0}_sl_HV_F__DAR.tif'.format(year)
            #or alos2
            else:
                file_name = lat_prefix + '{0}'.format(y_idx_loop) + lon_prefix + '{0}'.format(x_idx_loop)+'_' + '{0}_sl_HV_F02DAR.tif'.format(year)
                
            if x_idx_loop<100:
                #alos
                if alos_or_alos2 == 0:
                    file_name = lat_prefix + '{0}'.format(y_idx_loop) + lon_prefix + '0{0}'.format(x_idx_loop)+'_' + '20{0}_sl_HV_F__DAR.tif'.format(year)
                #alos2
                else:
                    file_name = lat_prefix + '{0}'.format(y_idx_loop) + lon_prefix + '0{0}'.format(x_idx_loop)+'_' + '{0}_sl_HV_F02DAR.tif'.format(year)
                #write_to_txt
            path_name = os.path.join(backscatter_source,file_name)
            if os.path.exists(path_name):
                write_to_txt(path_name + '\n',path_to_lists)
            else:
                print('Require access to ' + path_name  +' in path:' + workdir)
                return None
        
        
    cmd = 'gdal_merge.py -o ' +  path_to_output  +' --optfile '+ path_to_lists
    runCmd(cmd)
    #os.system(cmd)
        
        
                
            
        
               
                       
def backscatter_mosaic_v1(workdir,backscatter_source,year,alos_or_alos2 = 1):
    
    
    if (alos_or_alos2!=0 and alos_or_alos2!=1):
        print('The flag alos_or_alos2 should be either 0 or 1')
        return None
    
    files = os.listdir(workdir)
    for file in files:
        subdir_loop = os.path.join(workdir, file)
        if os.path.isdir(subdir_loop) and file.find('_')!=-1:
            gdal_mosaic_each_dir(subdir_loop,backscatter_source,year,alos_or_alos2)
            
            
def backscatter_mosaic_v2(workdir,backscatter_source,year, alos_or_alos2 = 1):
    
    if (alos_or_alos2!=0 and alos_or_alos2!=1):
        print('The flag alos_or_alos2 should be either 0 or 1')
        return None
    
    files = os.listdir(workdir)
    for file in files:
        subdir_loop = os.path.join(workdir, file)
        if os.path.isdir(subdir_loop) and file.find('_')!=-1:
            exe_path = os.path.join(subdir_loop,'isce_proc_dir')
        if os.path.exists(exe_path):
            files_sub = os.listdir(exe_path)
            for file_sub in files_sub:
                judge_dir = os.path.join(exe_path, file_sub,'forFSH')
                if os.path.isdir(judge_dir):
                    subdir_loop = os.path.join(exe_path, file_sub,'forFSH')
                    print('processing dir:' + subdir_loop)
                    gdal_mosaic_each_dir(subdir_loop,backscatter_source, year,alos_or_alos2)
                   
                   

    
if __name__ == '__main__':
    #workdir = '/media/yanglei/B87CC1AA7CC163AA/alos_hainan_1.1/468_370/isce_proc_dir/20070709_20070824/forFSH'
    #workdir = '/media/yanglei/B87CC1AA7CC163AA/alos_hainan_1.1/468_360/isce_proc_dir/20070709_20070824'
    #workdir = '/media/yanglei/B87CC1AA7CC163AA/alos_hainan_1.1/468_360/isce_proc_dir/20070824_20071009/forFSH'
    #workdir = '/media/yanglei/B87CC1AA7CC163AA/alos_hainan_1.1/468_360/isce_proc_dir/20071009_20071124/forFSH'
    #workdir = '/media/yanglei/B87CC1AA7CC163AA/alos_hainan_1.1/469_370/isce_proc_dir/20070610_20070726/forFSH'
    #workdir = '/media/yanglei/C656BE2356BE13E1/amazon_one_pair_test/alos_pair_2/815_930/forFSH'
    #workdir = '/media/yanglei/C656BE2356BE13E1/amazon_one_pair_test/alos_pair_2/630_815/forFSH'
    #workdir = '/media/yanglei/B87CC1AA7CC163AA/alos_hainan_1.1/469_360/isce_proc_dir/20070726_20070910/forFSH'
    #workdir = '/media/yanglei/B87CC1AA7CC163AA/370_v3/forFSH'
    
    #generate the cor_input.mat for the further processing;
    #Extract_PreProc('/media/yanghai/disk5/alos_2_china_us/neimeng_liaoning_heilongjiang/alos_pair3_431_830/pair_2/forFSH', 'cor_input.mat')
    
    #Or all the preprocessed InSAR data are placed in a root like
    #workdir = '/media/yanghai/disk3/dalian_reproc/proc'
    #The preproc function will go to each directory of InSAR pair in each path_frame directory to generate the cor_input.mat 
    #preproc_batch(workdir) 
      
      
    
    # Another function is provided here for mosaicing the ALOS backscatter tiles for each InSAR pair
    #processing one single InSAR pair 
    #the path to one single InSAR directory 
    #input path = /media/yanghai/disk3/china_ne_unzip/416_930/isce_proc_dir/20090819_20091004/forFSH
      
    #gdal_mosaic_each_dir('/media/yanghai/disk3/china_ne_unzip/416_930/isce_proc_dir/20090819_20091004/forFSH','/media/yanghai/disk4/china_ne_19_backscatter',19, 1 )
    
    workdir = '/media/yanghai/disk3/dalian_reproc/proc'
    #preproc_batch(workdir)
    #preproc_batch(workdir)
     
    #backscatter_debug
    workdir = '/media/yanghai/disk2/saihanba_reproc'
    #backscatter_source = '/media/yanghai/disk5/howland_backscatter'
    
    #backscatter_mosaic_v1(workdir,backscatter_source,'09')
    
    
    workdir = '/media/yanghai/disk3/426_2009'
    workdir = '/media/yanghai/disk3/china_ne_unzip'
    workdir = '/media/yanghai/disk2/saihanba_reproc'
    backscatter_source = '/media/yanghai/disk2/saihanba_backscatter_21'
    
    #gdal_mosaic_each_dir(workdir,backscatter_source,19,1)
    #backscatter_mosaic_v2(workdir,backscatter_source,21,1)
    
    workdir = '/media/yanghai/disk4/wnmf_reproc/126_870/isce_proc_dir/20100814_20100929/forFSH'
    backscatter_source = '/media/yanghai/disk4/wnmf_reproc/backscatter_10'
    #gdal_mosaic_each_dir(workdir,backscatter_source,10,0)
    
    
    
    
    workdir = '/media/yanghai/disk2/NewHampshire/forFSH'
    #workdir = '/media/yanghai/disk2/ma_vermont'
    output_name = 'cor_input.mat'
    #Extract_PreProc(workdir, output_name)
    #workdir = '/media/yanghai/disk5/me_07'
    
    
    backscatter_source = '/media/yanghai/disk4/china_ne_19_backscatter'
    #backscatter_source = '/media/yanghai/disk4/china_ne_21_backscatter'
    workdir = '/media/yanghai/disk3/china_ne_unzip'
    workdir = '/media/yanghai/disk3/reproc_v3/proc'
    workdir = '/media/yanghai/disk3/backscatter_reproc_v2'
    workdir = '/media/yanghai/disk3/dalian_reproc/proc'
    backscatter_mosaic_v2(workdir,backscatter_source,19,1)
    #backscatter_mosaic_v2(workdir,backscatter_source,21,1)
    #preproc_batch(workdir)

            
    
