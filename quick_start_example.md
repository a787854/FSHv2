   ## Download the example data including ALOS files, extracted GEDI data, backscatter product, and forest/Non-forest masking file from the [baiduNetdisk (code: 7gy9)](https://pan.baidu.com/s/1i8TLB8vmYJE_6xVxpUv77Q) or the [google drive link] 
  (https://drive.google.com/drive/folders/1L2rMJK9SsR7YuBTz3BUw0nmz_yhkCVK5?usp=drive_link;

   ## Open the “insar_proc_data” folder, and unzip all the ALOS files into the folder.

   ## Invoke the preproc_batch.py script to perform following operations:
   * Use the ClassfiyALOS function to classify the ALOS (Level-1 or Level-1.1) data according to its frame and orbit numbers, and to move the data into the corresponding “frame_orbit” folder;
   
   * Use the InSAR_form_batch function to form the InSAR pairs with a temporal baseline less than three months in each folder. The formed InSAR pairs are stored in the folders under the “isce_proc” folder;

   * Use the runISCEbatch function to perform the interferometric processing for each pair. The generated InSAR results stored in the forFSH folder under each InSAR folder;

   * Use the coh_extract_batch function to extract the coherence magnitude map and relevant parameters;

   * If you have all the needed backscatter products in one folder (like the backscatter folder in the example folder), please use backscatter_mosaic_v2 function to create backscatter_mosaic.tif file covering each InSAR scene in each directory;

  ## The following operations are executed in the Matlab environment:
  * Use the gedi_matching_batch.m script to the sparse rh95, rh98, or rh100 GEDI samples with the InSAR coherence magnitude map. Please note you should specify the path to GEDI mat file (e.g., fsh_demo\gedi_data\gedi_howland.mat in the example data);
  * Use the fsh_inversion_batch.m function to obtain the InSAR coherence-based forest height estimates, please specify the path to the folder including all the preprocessed data (e.g., fsh_demo\insar_proc_data);
  * Use the backscatter_inv_batch.m function to replace the FSH estimates of short vegetation with backscatter estimates (e.g., fsh_demo\insar_proc_data);
  * Invoke the mosaic_scenes.m function to performing mosaicking inside the folder.
