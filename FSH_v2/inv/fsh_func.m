function [] = fsh_func(rh_mode)
%clear,clc
tic
%load LIDAR_SAR_sel_98.mat
%load LIDAR_SAR_98_new.mat
%load LIDAR_SAR_rh100.mat
%load LIDAR_SAR_rh98.mat
%%load LIDAR_SAR_rh98.mat
%load LIDAR_SAR_rh98_full.mat 

%load LIDAR_SAR_rh98_fullpow_night.mat 
%load LIDAR_SAR_rh98_6alg.mat
%load LIDAR_SAR_rh100.mat

%1:rh100, 2:rh98, 3:rh95

if nargin==0
    rh_mode = 1;
end



switch rh_mode
    case 1
        if (not(exist('LIDAR_SAR_rh100.mat','file')))
            now_dir = pwd;
            warning(['No available LIDAR_SAR_rh100.mat. Invaild InSAR pair:' now_dir ]);
            return;
        end
        
        load LIDAR_SAR_rh100.mat;
        fsh_inv_filename = 'fsh_inv_rh100.tif';
        coh_inv_filename = 'coh_inv_rh100.tif';
        mat_file_name = 'Height_inversion_rh100.mat';
        txt_output = 'inv_ac_rh100.txt';
    case 2
        if (not(exist('LIDAR_SAR_rh98.mat','file')))
            now_dir = pwd;
            warning(['No available LIDAR_SAR_rh98.mat. Invaild InSAR pair:' now_dir ]);
            return;
        end
        
        load LIDAR_SAR_rh98.mat;
        fsh_inv_filename = 'fsh_inv_rh98.tif';
        coh_inv_filename = 'coh_inv_rh98.tif';
        mat_file_name = 'Height_inversion_rh98.mat';
        txt_output = 'inv_ac_rh98.txt';
    case 3
        
        if (not(exist('LIDAR_SAR_rh95.mat','file')))
            now_dir = pwd;
            warning(['No available LIDAR_SAR_rh95.mat. Invaild InSAR pair:' now_dir ]);
            return;
        end
        
        load LIDAR_SAR_rh95.mat;
        fsh_inv_filename = 'fsh_inv_rh95.tif';
        coh_inv_filename = 'coh_inv_rh95.tif';
        mat_file_name = 'Height_inversion_rh95.mat';
        txt_output = 'inv_ac_rh95.txt';
    otherwise
        error('Not supported input for rh_mode!');
        
end


%load LIDAR_SAR_rh98.mat;
%load LIDAR_SAR_rh98_fullpow_new.mat
%load LIDAR_SAR_rh98_6alg_v3.mat
rh98 = GEDI;
 %load LIDAR_SAR_no_sel.mat
 
 rh100  = GEDI ;
 %GEDI(GEDI>50) = nan;

%sfigure,imagesc(rh100 - rh98, 'AlphaData',~isnan(rh100),[-3,3]),colormap('turbo')
figure,imagesc(GEDI, 'AlphaData',~isnan(GEDI),[0,40]),colormap('jet')

% if(sum(isnan(GEDI(:)))== size(GEDI,1)*size(GEDI,2))
%     return
% end
% 
% if (exist('fsh_inv.tif','file'))
%     return
% end
hv_coh = hv_coh;

lat_axis = linspace(GEDI_coords(1), GEDI_coords(2) ,  size(hv_coh,1));
lon_axis = linspace(GEDI_coords(3), GEDI_coords(4) ,  size(hv_coh,2));

[LAT,LON] = ndgrid(lat_axis, lon_axis);


if 0

    [FNF,R] = readgeoraster('FNF_mosaic_v2.tif');
    [N,M] = size(FNF);
% %construct the grid in an agreement with classify map;
    lat_FNF_axis = linspace(R.LatitudeLimits(2),R.LatitudeLimits(1),N);
    lon_FNF_axis = linspace(R.LongitudeLimits(1),R.LongitudeLimits(2),M);
    mask = double(FNF<2);
    mask_interp = interp2(lon_FNF_axis,lat_FNF_axis(:),mask,lon_axis,lat_axis(:));
    mask_interp = mask_interp >=0.3;
    
end






[nRow, nCol] = size(hv_coh);

if (0)
    valididx=and(and(~isnan(hv_coh),hv_coh>=0.1),GEDI>13);
    coh_plot = hv_coh(valididx);
    GEDI_plot = GEDI(valididx);
    DensScat(coh_plot,GEDI_plot);xlim([0 1]),ylim([10 35]),xlabel('coh'),ylabel('GEDI RH98');,ylim([ 13 40]), xlim([0.1 1])
    coef = polyfit(coh_plot,GEDI_plot,1)
    
    x = [0.1 1];
    y = x*coef(1) + coef(2);
    hold on, plot(x,y,'r','LineWidth',1.5)
    
    Hm_trunc = coh_plot;
    Pm_trunc = GEDI_plot;
    save slope_assess.mat Hm_trunc Pm_trunc coef
end

hv_coh(hv_coh <=0.10) = nan;
valididx=and(~isnan(hv_coh),hv_coh>=0.1);
GEDI_global = GEDI(valididx);
coh_global = hv_coh(valididx);

if rh_mode == 1
    idx = and(GEDI_global>=12,GEDI_global<=50);
    if sum(idx(:)>0)<4e4
        idx = and(GEDI_global>=5,GEDI_global<=50);
    end
elseif rh_mode ==2
    idx = and(GEDI_global>=13,GEDI_global<=50);
    if sum(idx(:)>0)<5e4
        idx = and(GEDI_global>=7,GEDI_global<=50);
    end
elseif rh_mode ==3
    idx = and(GEDI_global>=2,GEDI_global<=50);
end

if exist('mask_interp','var')
    idx(mask_interp==0) = 0;
end


GEDI_global = GEDI_global(idx);
coh_global = coh_global(idx);


LUT = load('arc_sinc_data.mat');
 yy_array_1 = linspace(LUT.Y(1),LUT.Y(end), length(LUT.Y));
 xx_array_1 = interp1(LUT.Y,LUT.X,yy_array_1 );
 LUT_yy_array = yy_array_1;
 LUT_xx_array = xx_array_1;
%[sol] = search_2D_kb(double(coh_global), double(GEDI_global), 0.4,1.3,3,17,LUT.Y,LUT.X);
%




idx = and(~isnan(GEDI_global), and(GEDI_global<60, coh_global>0.05));
nvalid_gedi = sum(idx(:));

nvalid_gedi = sum(~isnan(GEDI(:)));

if(nvalid_gedi<10000)
    warning('This pair does not have enough gedi samples, we will not process it any more!');
    close all
    return
    
end



Smin_global = 0.5;
Smax_global = 1.5;
Cmin_global = 7;
Cmax_global = 25;
nS = 64;
nC = 64;
height_upper = 55;

[S_global,C_global] = global_cuda_interface(coh_global, GEDI_global, Smin_global, Smax_global,nS, Cmin_global,Cmax_global,nC,LUT_xx_array, LUT_yy_array,height_upper);


coef = polyfit(coh_global,GEDI_global,1);
if((coef(1))>0)
    DensScat(coh_global,GEDI_global);xlim([0 1]),ylim([10 35]),xlabel('coh'),ylabel('gedi height');
    warning('This pair show less sensitivity!!!');
end
DensScat(coh_global,GEDI_global);xlim([0 1]),ylim([10 35]),xlabel('coh'),ylabel('gedi height');
% 
if(or(S_global < 0.6 ,S_global >1.7))
    S_global= 0.9;
    C_global = 9;
end

if(or(C_global == Cmin_global,S_global >19))
    S_global= 0.9;
    C_global = 9;
end





sol = [S_global,C_global];
%sol = [   0.885714285714286   8.952755905511811];


   
 SS_op =  sol(1);
 CC_op = sol(2);
  
  % 16.4187
%0.4746   16.4187
input_coh = hv_coh;
input_coh(input_coh<0.03) = nan;
win = 5;
H_op=arc_sinc(input_coh./SS_op).*CC_op;H_op=nanmedfilt2(H_op,[win win]);
%H_op(GEDI_full<13) = nan;
figure,imagesc(H_op,'AlphaData',~isnan(H_op),[0 35]),colormap('jet'),title(['S-global =' num2str(S_global) ' , C-gloabl = ' num2str(C_global)]);

if 0
    fid = fopen('slope.txt','w');
    fprintf(fid,'slope:%f\n',coef(1));
    %fprintf(fid,'Mean_fitting_error:%f\n',Mean_fit_error);
    fclose(fid);
    close all
    return 
end


S_global = sol(1);
C_global = sol(2);
S_mat = zeros(size(hv_coh)) + S_global;
C_mat = zeros(size(hv_coh)) + C_global;
S_mat(isnan(hv_coh)) = nan; 
C_mat(isnan(hv_coh)) = nan;

residual = zeros(size(hv_coh));
gedi_canopy_height = GEDI;
gedi_canopy_height(gedi_canopy_height>60) = nan;
gamma = hv_coh;

win_nx = 32; win_ny = 32;
vec = 1:win_nx+1;
vec =abs( vec- win_nx/2-1);
distance =sqrt( vec.^2 + vec.'.^2);
x_central = win_nx/2 + 1;
y_central = win_ny/2 + 1;


min_val = min(distance(:));
max_val = max(distance(:));
%30 as 
weight = zeros(size(distance))+0.25;

distance_thres_1 = min_val + 0.3* (max_val - min_val);
weight(distance<=distance_thres_1) = 1;

distance_thres_2 =  min_val + 0.6* (max_val - min_val);
weight(and(distance>=distance_thres_1,distance<=distance_thres_2)) = 0.5;

weight(y_central,x_central) = 11;

ny = size(gedi_canopy_height,1);
nx = size(gedi_canopy_height,2);
y_seq = 1:ny;
x_seq = 1:nx;
max_valid = -999;


S_search = [S_global - 0.5 1.3];
if S_search(1)<0.55
    S_search(1) = 0.55;
end
C_search = [C_global - 4 C_global + 7 ];
S_min = S_search(1);
S_max = S_search(end);
C_min = C_search(1);
C_max = C_search(end);


%gedi_canopy_height(gedi_canopy_height) = nan;

 
 idx = and(~isnan(gedi_canopy_height), and(gedi_canopy_height<60, gamma>0.05));
nvalid_gedi = sum(idx(:));

if(nvalid_gedi<10000)
    warning('This pair does not have enough gedi samples, we will not process it any more!');
end
 
 nS = 64;
 nC = 64;
 
 
 
  [h_resi,S_par, C_par, lines_central,pixels_central] = local_op_cuda_interface(gamma,gedi_canopy_height,  nvalid_gedi, S_min, S_max,...
    nS, C_min, C_max, nC , LUT_xx_array, LUT_yy_array, win_nx, win_ny, weight, 60, 0, 0.05, 11);
 

[Y,X] = ndgrid(y_seq,x_seq);

Xsub = double(pixels_central);
Ysub = double(lines_central);

method='natural';

F = scatteredInterpolant(Ysub,Xsub,S_par,method);
SS_op=F(Y,X); SS_ = nanmedfilt2(SS_op,[win win]);

F = scatteredInterpolant(Ysub,Xsub,C_par,method);
CC_op=F(Y,X); CC_ = nanmedfilt2(CC_op,[win win]);

F = scatteredInterpolant(Ysub,Xsub,h_resi,method);
Resi_op=F(Y,X);


mask = isnan(input_coh);
SS_op(mask)=nan;CC_op(mask)=nan;
input_coh = hv_coh;
input_coh(input_coh<0.02) = nan;
H_op=arc_sinc(input_coh./SS_op).*CC_op;H_op=nanmedfilt2(H_op,[3 3]);
figure,imagesc(H_op,'AlphaData',~isnan(H_op),[0,35]),colormap('jet')
H_inverted = H_op;
H_op(or(H_op<0,H_op>60)) = nan;
hv_coh(or(H_op<0,H_op>60)) = nan;



R_new = georasterref('RasterSize', size(H_op), ...
    'RasterInterpretation', 'cells', 'ColumnsStartFrom', 'north', ...
    'LatitudeLimits', [min(lat_axis) max(lat_axis)], 'LongitudeLimits', [min(lon_axis) max(lon_axis) ]);

EPSG = 4326;



geotiffwrite(fsh_inv_filename,H_op,R_new,'CoordRefSysCode', EPSG);
geotiffwrite(coh_inv_filename,hv_coh,R_new,'CoordRefSysCode', EPSG);
%S_discrete = S_par;

% CC_ = nanmedfilt2(CC_op,[5 5]);
% SS_ = nanmedfilt2(SS_op,[5 5]);
% H_op=arc_sinc(input_coh./SS_op).*CC_;H_op=nanmedfilt2(H_op,[5 5]);
% figure,imagesc(H_op,'AlphaData',~isnan(H_op),[0,35]),colormap('jet'),title('nan median filtering on S,C par')
% win1 = 16;
% S_filter2 = nanconv(SS_op,ones(win1)./win1^2);
% C_filter2 = nanconv(CC_op,ones(win1)./win1^2);
% H_op=arc_sinc(input_coh./S_filter2).*C_filter2;H_op=nanmedfilt2(H_op,[5 5]);
% figure,imagesc(H_op,'AlphaData',~isnan(H_op),[0,35]),colormap('jet'),title('nanconv (16) on S,C par')
Resi_op(isnan(hv_coh)) = nan;
Resi_op(not(and(Resi_op>0,Resi_op<150))) = nan;
Mean_fit_error = mean(Resi_op(:),'omitnan');



save(mat_file_name,'H_inverted','GEDI','GEDI_coords','SS_op','CC_op','Resi_op','Mean_fit_error');



fid = fopen(txt_output,'w');
fprintf(fid,'slope:%f\n',coef(1));
fprintf(fid,'Mean_fitting_error:%f\n',Mean_fit_error);
fclose(fid);

close all
return

toc


C_discrete = C_par;
S_full = SS_op;
C_full = CC_op;
GEDI_discrete = GEDI;
Line_discrete = Ysub;
Pixel_discrete = Xsub;

save test_rh98.mat input_coh S_discrete C_discrete S_full C_full GEDI_full GEDI_discrete Line_discrete Pixel_discrete


%write a geotiff files

% S_discrete = S_par;
% C_discrete = C_par;
% L_discrete = lines_central;
% P_discrete = pixels_central;

%save Height_inversion_pycuda_check.mat S_discrete C_discrete L_discrete P_discrete lat_axis lon_axis


%save Height_inversion_rh98.mat H_inverted GEDI GEDI_coords SS_op CC_op Resi_op

end

