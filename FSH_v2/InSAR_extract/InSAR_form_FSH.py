#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 09:54:48 2023

@author: yanghai
"""
import numpy as np 
import scipy.io as sio
#import commands
import subprocess
import os
import time
import argparse
import pdb
import os
import isce
import shelve
import string
import sys
import shutil


def DeleteFiles(path):
    dirs = []
    files = os.listdir(path)
    for file in files:
        subdir_loop = os.path.join(path, file)
        if os.path.isdir(subdir_loop):
            #files_inner = os.listdir(subdir_loop)
            #dest_folder = os.path.join(subdir_loop, 'forFSH')
            
            
           files_2 = os.listdir(subdir_loop)
            
           for file_2 in files_2:
               file_path_2 = os.path.join(subdir_loop, file_2)
               if(file_2 != 'forFSH'):
                   if(os.path.isdir(file_path_2)):
                       shutil.rmtree(file_path_2, ignore_errors=True)
                   else:
                       os.remove(file_path_2)

def runCmd(cmd):
    out = subprocess.getoutput(cmd)
    return out

def Extract_files(path):
    dirs = []
    files = os.listdir(path)
    for file in files:
        subdir_loop = os.path.join(path, file)
        if os.path.isdir(subdir_loop):
            #files_inner = os.listdir(subdir_loop)
            dest_folder = os.path.join(subdir_loop, 'forFSH')
            
            
            if (os.path.exists(dest_folder)):
                continue
        
            if (not(os.path.exists(dest_folder))):
                os.mkdir(dest_folder)
                
            
            if os.path.exists(os.path.join(subdir_loop, 'stripmapProc.xml')):
               shutil.copy (os.path.join(subdir_loop, 'stripmapProc.xml'),dest_folder)
            else:
                continue
            
            
            
            subdir_intf = os.path.join(subdir_loop, 'interferogram')
            
            #resampOnlyImg
            path_tmp = os.path.join(subdir_intf, 'resampOnlyImage.amp.geo')
            if os.path.exists(path_tmp):
               shutil.copy (path_tmp,dest_folder)
            else:
                continue
            
            path_tmp = os.path.join(subdir_intf, 'resampOnlyImage.amp.geo.vrt')
            if os.path.exists(path_tmp):
               shutil.copy (path_tmp,dest_folder)
            else:
                continue
            
            path_tmp = os.path.join(subdir_intf, 'resampOnlyImage.amp.geo.xml')
            if os.path.exists(path_tmp):
               shutil.copy (path_tmp,dest_folder)
            else:
                continue
            

            #topophase.cor.geo
            path_tmp = os.path.join(subdir_intf, 'topophase.cor.geo')
            if os.path.exists(path_tmp):
               shutil.copy (path_tmp,dest_folder)
            else:
                continue
            
            path_tmp = os.path.join(subdir_intf, 'topophase.cor.geo.vrt')
            if os.path.exists(path_tmp):
               shutil.copy (path_tmp,dest_folder)
            else:
                continue
            
            path_tmp = os.path.join(subdir_intf, 'topophase.cor.geo.xml')
            if os.path.exists(path_tmp):
               shutil.copy (path_tmp,dest_folder)
            else:
                continue


            
            
            
       
            
    #return dirs
    
    
def extract_files_single(path):
    files = os.listdir(path)
    subdir_loop = path
    
    if os.path.isdir(subdir_loop):
            #files_inner = os.listdir(subdir_loop)
            dest_folder = os.path.join(subdir_loop, 'forFSH')
            
            if (not(os.path.exists(dest_folder))):
                os.mkdir(dest_folder)
                    
            shutil.copy (os.path.join(subdir_loop, 'stripmapProc.xml'),dest_folder)
            
            
            subdir_intf = os.path.join(subdir_loop, 'interferogram')
            shutil.copy (os.path.join(subdir_intf, 'resampOnlyImage.amp.geo'),dest_folder)
            shutil.copy (os.path.join(subdir_intf, 'resampOnlyImage.amp.geo.vrt'),dest_folder)
            shutil.copy (os.path.join(subdir_intf, 'resampOnlyImage.amp.geo.xml'),dest_folder)
            
            shutil.copy (os.path.join(subdir_intf, 'topophase.cor.geo'),dest_folder)
            shutil.copy (os.path.join(subdir_intf, 'topophase.cor.geo.vrt'),dest_folder)
            shutil.copy (os.path.join(subdir_intf, 'topophase.cor.geo.xml'),dest_folder)



def DeleteFiles_single(path):
    dirs = []
    subdir_loop = path
    if os.path.isdir(subdir_loop):
            files_2 = os.listdir(subdir_loop)
            
            for file_2 in files_2:
               file_path_2 = os.path.join(subdir_loop, file_2)
               if(file_2 != 'forFSH'):
                   if(os.path.isdir(file_path_2)):
                       shutil.rmtree(file_path_2, ignore_errors=True)
                   else:
                       os.remove(file_path_2)           
            


        
                       

def runISCE_all_subdir(path,alos_or_alos1_1=0):
    dirs = []
    files = os.listdir(path)
    for file in files:
        subdir_loop = os.path.join(path, file)
        print(subdir_loop)
        if os.path.isdir(subdir_loop):
            if not(os.path.isdir(os.path.join(subdir_loop,'interferogram')) or os.path.isdir(os.path.join(subdir_loop,'forFSH'))):
                print(path)
                if alos_or_alos1_1 ==0:
                    cmd ='single_scene_stripmapApp_alos_slc.py -f {0}'.format(os.path.join(path, file))
                else:
                    cmd ='single_scene_stripmapApp.py -f {0}'.format(os.path.join(path, file))
                runCmd(cmd)
                #subprocess.call(["single_scene_stripmapApp_alos_slc.py","-f",os.path.join(path, file)],shell=True)
                
#                extract_files_single(os.path.join(path, file))
#                DeleteFiles_single(os.path.join(path, file))
            

def runISCE_batch(rootpath,alos_or_alos1_1=0):
    files = os.listdir(rootpath)
    for file in files:
        subdir_loop = os.path.join(rootpath, file)
        if os.path.isdir(subdir_loop) and file.find('_')!=-1:
            exe_path = os.path.join(subdir_loop,'isce_proc_dir')
            if os.path.exists(exe_path):
                runISCE_all_subdir(exe_path,alos_or_alos1_1)
    


def runExtractFile_batch (rootpath):
    files = os.listdir(rootpath)
    for file in files:
        subdir_loop = os.path.join(rootpath, file)
        if os.path.isdir(subdir_loop) and file.find('_')!=-1:
            exe_path = os.path.join(subdir_loop,'isce_proc_dir')
            if os.path.exists(exe_path):
                Extract_files(exe_path)
    
    

                       
                       
def runDeleteFiles_batch (rootpath):
    files = os.listdir(rootpath)
    for file in files:
        subdir_loop = os.path.join(rootpath, file)
        if os.path.isdir(subdir_loop) and file.find('_')!=-1:
            exe_path = os.path.join(subdir_loop,'isce_proc_dir')
            if os.path.exists(exe_path):
                DeleteFiles(exe_path)
                       

def InSAR_Form(path, flag_intf_proc = 0 ):
    
    files = os.listdir(path)
    print(path) 
    path2file = []
    file_name = [] 
    date_str = []
    year_num = []
    month_num = []
    
    kkk = 0
    alos1_or_alos1_1 = 0
    for file in files:
        subdir = os.path.join(path,file)
        if os.path.isdir(subdir) and file.find('ALPSRP')!=-1:
          os.chdir(subdir)   
          kkk = kkk+1
           #entry(kkk). name = file
           #code = file[6:15]
          try:
              filepath = str.split(runCmd('find `pwd` -name "summary.txt"'))[0]
              alos1_or_alos1_1 = 0#alos1.1
          except:
              filepath = str.split(runCmd('find `pwd` -name "workreport"'))[0]
              alos1_or_alos1_1 = 1 #alos1.0
          #filepath = str.split(runCmd('find `pwd` -name "workreport"'))
         
         
        
          if alos1_or_alos1_1 ==0:
              line = runCmd('fgrep Lbi_ObservationDate '+filepath)
          else:
              line = runCmd('fgrep Img_SceneCenterDateTime '+filepath)

          #name = path[0:-11]
          path2file.append(subdir)
          file_name.append(file)
          if alos1_or_alos1_1 ==1:
              date = str.split(line,'=')[1][2:10] #alos l0
          else:
              date = str.split(line,'=')[1][1:9] #alos l1
          date_str.append(date )
          
          year_num.append(np.float64(date[0:4]))
          month_num.append(np.float64(date[-4:-2]))
          os.chdir(path) 
          
        
        
    n = np.size(file_name)
    
    
    name_dir = 'isce_proc_dir'
    path2dir = os.path.join(path,name_dir)
    if not(os.path.exists(path2dir)):
        os.mkdir(path2dir)
    else:
        return 0
    
    path2_proc_dir = os.path.join(path,name_dir)    
    for i in range(n):
        mon_1 = month_num[i] + year_num[i]*12
        for j in range(i+1,n):
            mon_2 = month_num[j] + year_num[j]*12
            
            
            if abs(mon_2 - mon_1)<3.5:
                
            
                subsubdir_name = date_str[i] + '_' + date_str[j]
            
                path2intf = os.path.join(path2_proc_dir, subsubdir_name)
            
                if not(os.path.exists(path2intf)):
                    os.mkdir(path2intf)
                    
                subdir_i = os.path.join(path2intf, file_name[i])
                 
                if not(os.path.exists(subdir_i)):
                    #os.mkdir(subdir_i)
                    shutil.copytree(path2file[i],subdir_i)
            
                subdir_j = os.path.join(path2intf, file_name[j])
                 
                if not(os.path.exists(subdir_j)):
                    #os.mkdir(subdir_j)
                    shutil.copytree(path2file[j],subdir_j)
                
                if (flag_intf_proc == 1):
                    cmd ='single_scene_stripmapApp_alos_slc.py -f {0}'.format(path2intf)
                    runCmd(cmd)
                    
                    
    return alos1_or_alos1_1
                    
def InSAR_form_batch(rootpath, flag_intf_proc=0):
    files = os.listdir(rootpath)
    
    for file in files:
        subdir = os.path.join(rootpath,file)
        os.chdir(rootpath)
        if os.path.isdir(subdir) and file.find('_')!=-1:
            InSAR_Form(subdir,flag_intf_proc)
        

    return 0
                    
                    
                    
    
if __name__ == '__main__':
    #forming the InSAR pairs.....
    #exe_path = '/media/yanghai/disk4/426_2009_proc'
    #exe_path ='/media/yanghai/disk4/431_2009_proc'
    #exe_path = '/media/yanghai/disk4/wnmf_reproc'
    exe_path = '/media/yanghai/disk4/423_reproc'
    exe_path = '/media/yanghai/disk4/418_920_reproc/'
    exe_path = '/media/yanghai/disk3/reproc_massive_v2'
    exe_path = '/media/yanghai/disk3/dalian_reproc/proc'
    #exe_path = '/media/yanghai/disk3/china_ne_unzip/443_940'
    #InSAR_Form(exe_path)
    #runISCE_all_subdir(os.path.join(exe_path,'isce_proc_dir'))
    #runDeleteFiles_batch(exe_path)
    code = InSAR_form_batch(exe_path)
    
    #executing ISCE processing
    runISCE_batch(exe_path,code)
    
    exe_path  = '/media/yanghai/disk3/china_ne_unzip'
    #extract_files_single(work_dir)
    #DeleteFiles_single(work_dir)
    #extracting the files 
    #runExtractFile_batch(exe_path)
    #runDeleteFiles_batch(exe_path)
 
    
    
    #exe_path = '/media/yanglei/B87CC1AA7CC163AA/alos_hainan_1.1/468_360'
    #InSAR_Form(exe_path)
    #exe_path = '/media/yanglei/B87CC1AA7CC163AA/alos_hainan_1.1/468_360/isce_proc_dir'
    #runISCE_all_subdir(exe_path)           
            