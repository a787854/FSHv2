%define the bbox of the scene
clear,clc

%reading the mask file 
% [MASK,R] = geotiffread('Maine_NLCD2011_nonwildland.tif');
% [N,M] = size(MASK);
% %construct the grid in an agreement with classify map;
% lat_ref_axis = linspace(R.LatitudeLimits(2),R.LatitudeLimits(1),N);
% lon_ref_axis = linspace(R.LongitudeLimits(1),R.LongitudeLimits(2),M);


%china_ne_bbox

workdir = '/media/yanghai/disk3/china_ne_unzip';
cd(workdir);

extra_radius = 0.01;
lat_min = 38.6 - extra_radius;
lat_max = 51.6 + extra_radius;
lon_min = 115 - extra_radius;
lon_max = 134.93 + extra_radius;
    
delta_lat = -0.000277777*3;
delta_lon = -delta_lat;
%lat_max = 47.47;
%lat_min = 40.4;
%lon_max = -66.94;
%lon_min = -74.5;

%orginal lat/lon grids in masking images
%idx_lat = and(lat_ref_axis>=lat_min,lat_ref_axis<=lat_max);
%lat_ax_mosaic = lat_ref_axis(idx_lat);
%idx_lon =  and(lon_ref_axis>=lon_min,lon_ref_axis<=lon_max);
%lon_ax_mosaic = lon_ref_axis(idx_lon);

lat_ax_mosaic = lat_max:delta_lat:lat_min;
lon_ax_mosaic = lon_min:delta_lon:lon_max;

load mask_china_ne.mat
[N,M] = size(mask_map);
lat_mask_axis = linspace(R.LatitudeLimits(2),R.LatitudeLimits(1),N);
lon_mask_axis = linspace(R.LongitudeLimits(1),R.LongitudeLimits(2),M);
mask_map_interp = interp2(lon_mask_axis,lat_mask_axis(:),single(mask_map),lon_ax_mosaic,lat_ax_mosaic(:));

mask_map_interp = mask_map_interp>0.5;

%lat_ax_mosaic = lat_max: lat_delta:lat_min;
%lon_ax_mosaic = lon_min:lon_delta:lon_max;
height_mosaic = zeros(length(lat_ax_mosaic),length(lon_ax_mosaic))+ nan;
error_mosaic =  zeros(length(lat_ax_mosaic),length(lon_ax_mosaic)) + nan;

row_idx = 1:length(lat_ax_mosaic);
col_idx = 1:length(lon_ax_mosaic);

%clear up_mat middle_mat left_mat right_mat down_mat


info = dir;
diary on

for i = 3:length (info)
    
    if(~exist(info(i).name,'dir'))
        continue
    end
    
    isce_proc_dir = [workdir '/' info(i).name '/isce_proc_dir' ];
    if (not(exist(isce_proc_dir,'dir')))
        continue
    end
    
    mat_name_loop = find_valid_mat(isce_proc_dir);
    
    if(isempty(mat_name_loop))
        continue;
    end
    
    
    %mat_name_loop = [info(i).name '/Height_inv_with_backscatter.mat' ];
    
    load(mat_name_loop);
    H_inverted(or(H_inverted<0,H_inverted>60))  = nan;
    [N,M] = size(H_inverted);
    lat_axis_loop = linspace(GEDI_coords(1),GEDI_coords(2),N);
    lon_axis_loop = linspace(GEDI_coords(3),GEDI_coords(4),M);
    
    %determine the row indexH_inverted_comb
    idx_lat_loop = and(lat_ax_mosaic>=GEDI_coords(2),lat_ax_mosaic<=GEDI_coords(1));
    lat_interp_loop = lat_ax_mosaic(idx_lat_loop);
    row_idx_loop  = row_idx(idx_lat_loop);
    if(isempty(row_idx_loop))
        continue;
    end
    y0 = row_idx_loop(1);
    yL = length(row_idx_loop);
    
    %determine the col index
    idx_lon_loop = and(lon_ax_mosaic>=GEDI_coords(3),lon_ax_mosaic<=GEDI_coords(4));
    lon_interp_loop = lon_ax_mosaic(idx_lon_loop);
    col_idx_loop = col_idx(idx_lon_loop);
    if(isempty(col_idx_loop))
        continue;
    end
    x0 = col_idx_loop(1);
    xL = length(col_idx_loop);
    
    %
    height_interp = interp2(lon_axis_loop,lat_axis_loop',H_inverted,lon_interp_loop,lat_interp_loop');
    error_interp  = interp2(lon_axis_loop,lat_axis_loop',Resi_op,lon_interp_loop,lat_interp_loop');
    
    
    for kk= 1:yL
        for jj = 1:xL
            if(isnan(height_interp(kk,jj)))
                continue
            end
            
            input_h = height_interp(kk,jj);
            output_h = height_mosaic(kk-1+y0, jj -1 + x0);
            
            input_err = error_interp(kk,jj);
            output_err = error_mosaic(kk-1+y0, jj -1 + x0);
            
            if(isnan(output_h))
                height_mosaic(kk-1+y0, jj -1 + x0) = input_h;
            else
                if(input_err<output_err)
                    height_mosaic(kk-1+y0, jj -1 + x0) = input_h; % (input_h+output_h)/2.0;
                end
            end
            
        end
    end
    
    
    
    
    
    
end

diary off
%Mask_mosaic = MASK(idx_lat,idx_lon);

%height_mosaic(MASK==1) = nan;

figure,imagesc(height_mosaic,[0 35]),colormap('jet')
%figure,imagesc(height_mosaic,'AlphaData',~isnan(height_mosaic),[0 35]),colormap('jet')

height_write = height_mosaic;

height_write(mask_map_interp==0)  = nan;
%height_write(isnan(height_write)) = 255;
figure,imagesc(height_write,[0 35]),colormap('jet')
R_new = georasterref('RasterSize', size(height_mosaic), ...
      'RasterInterpretation', 'cells', 'ColumnsStartFrom', 'north', ...
      'LatitudeLimits', [min(lat_ax_mosaic) max(lat_ax_mosaic)], 'LongitudeLimits', [min(lon_ax_mosaic) max(lon_ax_mosaic) ])
  
  geotiffwrite('China_NE_90m.tif',height_write,R_new,'TiffType','bigtiff','CoordRefSysCode', 4326);
%figure,imagesc(lon_vec,(lat_vec),B,[0 40]),axis xy

function [slope_val] = read_slope_txt(file_name)
    fid = fopen(file_name);
    oneLine = fgetl(fid);
    tmp = split(oneLine,':');
    slope_val = tmp{2};
    fclose (fid);
end

function [slope_val,fitting_error] = read_inv_txt(file_name)
    fid = fopen(file_name);
    oneLine = fgetl(fid);
    tmp = split(oneLine,':');
    slope_val = str2num(tmp{2});
    oneLine = fgetl(fid);
    tmp = split(oneLine,':');
    fitting_error = str2num(tmp{2});
    fclose (fid);
end


function [str_to_mat_name_loop]  = find_valid_mat(sub_dir)
    info = dir(sub_dir);
    
    
    
    str_to_mat_name_loop = [];
    
    slope_save = 99999;

    for i = 3:length (info)
        
        if(~exist([sub_dir '/' info(i).name],'dir'))
            continue
        end
        
        
        dir_loop = [sub_dir '/' info(i).name '/forFSH'];
        
        if (exist(dir_loop, 'dir'))
            height_mat = [dir_loop '/Height_inv_with_backscatter_rh98.mat' ];
            
            if(not(exist(height_mat,'file')))
                continue;
            end
            
            
            slope_val_loop = nan;
            slope_file = [dir_loop '/slope.txt' ];
            if ((exist(slope_file,'file')))
                
                %height_mat = [info(i).name '/Height_inv_with_backscatter_rh98.mat' ];
                
                slope_val_loop = str2num(read_slope_txt(slope_file));
            end
            
            slope_file2 = [dir_loop '/inv_ac_rh98.txt' ];
            if ((exist(slope_file2,'file')))
                
                %height_mat = [info(i).name '/Height_inv_with_backscatter_rh98.mat' ];
                
                %slope_val_loop = str2num(read_slope_txt(slope_file2));
                [slope_val_loop,fitting_error] = read_inv_txt(slope_file2);
%                 if(fitting_error>6)
%                     slope_val_loop = nan;
%                 end
                
                
            end

            if(isnan(slope_val_loop))
                warning(['Invalid Pair. Please check:' dir_loop]);
                continue;
            end
            disp(dir_loop);
            if (slope_val_loop<slope_save)
                slope_save = slope_val_loop;
                str_to_mat_name_loop = height_mat;
            end
            
            
        end
            
        
        
        
        
        
    end
    
    if(slope_save==99999)
         warning(['This InSAR scene is not valid:' sub_dir]);
    end
    
    

end