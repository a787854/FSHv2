# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 23:24:26 2024

A script to show how to use the relevant function to extract the needed data.

@author: yanghai
"""


import InSAR_form_FSH as Iff
import FSH_prepro as Fp

if __name__ == '__main__':


    #specify your working directory containing all the classified ALOS InSAR pair foler
    exe_path = '/media/yanghai/Proc_dir'
    
    #forming all the InSAR pairs within three months
    code = Iff.InSAR_form_batch(exe_path)
    
    #run ISCE software to process each subdir
    Iff.runISCE_batch(exe_path,code)
    
    #extract the information after the processing, you can also use the preproc function to extract the needed files in specific folder
    Fp.preproc_batch(exe_path)    
    
    
    #if you have all the backscatter files (should be .tif) in one folder, you can mosaic the backscatter files to cover each InSAR scene
    backscatter_source = '/media/yanghai/china_ne_21_backscatter'
    alos2_flag = 1;
    alos_flag = 0;
    # 21 is the last two letters of 2021
    Fp.backscatter_mosaic_v2(exe_path,backscatter_source,21,alos2_flag)
    
    
    

    
    
    
        

