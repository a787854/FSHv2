function [] = backscatter_inv_func(rh_mode)
%generate gamma_naought


if nargin ==0
    rh_mode = 1
end



[DN,R] = readgeoraster('backscatter_mosaic.tif');
DN(DN==0) = nan;

gamma0 = 20*log10(double(DN)) - 83;
gamma0 =10.^(gamma0./10);
gamma0(gamma0==0) = nan;
%gamma0=nanmedfilt2(gamma0,[3 3]);
kernel = ones(3,3)./ 9;
gamma0 = conv2(gamma0,kernel,'same');
DN = conv2(DN,kernel,'same');

%geotiffwrite('amazon_gamma0.tif',gamma0,R);


[N,M] = size(DN);
% %construct the grid in an agreement with classify map;
 lat_ref_gamma = linspace(R.LatitudeLimits(2),R.LatitudeLimits(1),N);
 lon_ref_gamma = linspace(R.LongitudeLimits(1),R.LongitudeLimits(2),M);
 
 replace_value = 0;
 
 switch rh_mode
     case 1
         if not(exist('fsh_inv_rh100.tif','file'))
            [FSH,R]  = readgeoraster('fsh_inv.tif');
         else
            [FSH,R]  = readgeoraster('fsh_inv_rh100.tif');
         end
          replace_value = 17;
          load LIDAR_SAR_rh100.mat
     case 2
         [FSH,R]  = readgeoraster('fsh_inv_rh98.tif');
          load LIDAR_SAR_rh98.mat
          replace_value = 17;
     case 3
         [FSH,R]  = readgeoraster('fsh_inv_rh95.tif');
          load LIDAR_SAR_rh95.mat
          replace_value = 13;
     otherwise
          [FSH,R]  = readgeoraster('fsh_inv.tif');
          load LIDAR_SAR_rh100.mat
 end
 
 %[FSH,R]  = readgeoraster('fsh_inv_rh100.tif');
 FSH(FSH==nan) = 0;
  [N,M] = size(FSH);
% %construct the grid in an agreement with classify map;
 lat_FSH_axis = linspace(R.LatitudeLimits(2),R.LatitudeLimits(1),N);
 lon_FSH_axis = linspace(R.LongitudeLimits(1),R.LongitudeLimits(2),M);
 
 
 back_interp = interp2(lon_ref_gamma, lat_ref_gamma(:), gamma0,lon_FSH_axis,lat_FSH_axis(:),'linear');
 DN_interp = interp2(lon_ref_gamma, lat_ref_gamma(:), double(DN),lon_FSH_axis,lat_FSH_axis(:),'linear');
 
 

 
 idx = and (~isnan(back_interp),and(GEDI>=0.5,GEDI<=20));
 
A  = 0.10;
B = 0.0622;
C = 1.0143;
%fixing C and 
%A = 0.0759;
%B =  0.1566;
%C =  1.0143;
 sol(1) = A;
 sol(2) = B;


A_search = [A-A/3,  A+A/3];
B_search = [B-B/3,  B+B/3];
%[sol] = kb_fitting_exp(double(GEDI(idx)),double(back_interp(idx)),A_search,B_search);
sol = backscatter_search_2D(double(back_interp(idx)), double(GEDI(idx)), A_search(1), A_search(2), B_search(1),B_search(2),1);
A = sol(1);
B = sol(2);
 
h_inv = exp_height_inv(double(back_interp),sol(1),sol(2),C);
 
 
 

win = 3;
h_inv=nanmedfilt2(h_inv,[win win]);


idx = and (~isnan(back_interp),GEDI<=13);




coef = polyfit(double(back_interp(idx)), double(GEDI(idx)),1);
h_inv2 = back_interp.*coef(1) +coef(2);


coef2 = polyfit(double(DN_interp(idx)), double(GEDI(idx)),1);
h_inv3 = DN_interp.*coef2(1) +coef2(2);

% figure,
% ax(1) = subplot(1,3,1),imagesc(h_inv,[0,30]),colormap('jet');
% ax(2) = subplot(1,3,2),imagesc(h_inv2,[0,30]),colormap('jet'); 
% ax(3) = subplot(1,3,3),imagesc(h_inv3,[0,30]),colormap('jet');
% 
% linkaxes(ax);





if exist('inv_ac_rh98.txt','file')
    
    [slope_val,fitting_error] = read_inv_txt('inv_ac_rh98.txt');
    if(fitting_error>5)
        replace_value = 70;
    end
end

inv__ = h_inv;

idx = and(inv__<=10,FSH<=60);

FSH(idx) = inv__(idx);
%FSH = h_inv;
FSH(FSH<0)= nan;
FSH(FSH>50)= nan;
switch rh_mode
    case 1
        geotiffwrite('rh100_with_backscatter.tif',FSH,R);
    case 2
        geotiffwrite('rh98_with_backscatter.tif',FSH,R);
    case 3
        geotiffwrite('rh95_with_backscatter.tif',FSH,R);
end
%geotiffwrite('fsh100_with_backscatter.tif',FSH,R);
 
% if not(exist('Height_inversion_rh100.mat','file'))
%    load Height_inversion.mat;
% else
%    load Height_inversion_rh100.mat;
% end
% 
% H_inverted = FSH;
%save Height_inv_with_backscatter.mat H_inverted CC_op SS_op Resi_op GEDI GEDI_coords 


switch rh_mode
    case 1
        %geotiffwrite('rh100_with_backscatter.tif',FSH,R);
        load Height_inversion_rh100.mat;
        H_inverted = FSH;
        save Height_inv_with_backscatter_rh100.mat H_inverted CC_op SS_op Resi_op GEDI GEDI_coords
    case 2
        load Height_inversion_rh98.mat;
        H_inverted = FSH;
        save Height_inv_with_backscatter_rh98.mat H_inverted CC_op SS_op Resi_op GEDI GEDI_coords
    case 3
        load Height_inversion_rh95.mat;
        H_inverted = FSH;
        save Height_inv_with_backscatter_rh95.mat H_inverted CC_op SS_op Resi_op GEDI GEDI_coords
end
close all
 
end


function [slope_val,fitting_error] = read_inv_txt(file_name)
    fid = fopen(file_name);
    oneLine = fgetl(fid);
    tmp = split(oneLine,':');
    slope_val = str2double(tmp{2});
    oneLine = fgetl(fid);
    tmp = split(oneLine,':');
    fitting_error = str2double(tmp{2});
    fclose (fid);
end

