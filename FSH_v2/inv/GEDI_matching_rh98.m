function [file_name_save] = GEDI_matching_rh98(workdir,path_to_gedi_mat, gedi_int_flag, GEDI_bigger)

path_record = pwd;
file_name_save = 'LIDAR_SAR_rh98.mat';
cd (workdir);

if(not(exist('cor_input.mat','file')))
    error('We cannot find cor_input.mat file in the current directory!');
end


if not(exist('gedi_int_flag','var'))
    gedi_int_flag = 0;
end

load ('cor_input.mat');
corr_vs = corr_vs';
lat_axis = linspace(coords(1),coords(2),size(corr_vs,1));
lon_axis = linspace(coords(3),coords(4),size(corr_vs,2));

[LAT,LON] = ndgrid(lat_axis,lon_axis);

%crop the image
mask_ = isnan(corr_vs);
LAT(mask_) = nan;
LON(mask_) = nan; 


%starting row
[max_lat,id1] = max(LAT(:),[],'omitnan');
[max_lat_row,col] = ind2sub(size(LAT),id1);

%ending row 
[min_lat,id1] = min(LAT(:),[],'omitnan');
[min_lat_row,col] = ind2sub(size(LAT),id1);

%staring col 
[min_lon,id1] = min(LON(:),[],'omitnan');
[row,min_lon_col] = ind2sub(size(LAT),id1);

%ending col 
[max_lon,id1] = max(LON(:),[],'omitnan');
[row,max_lon_col] = ind2sub(size(LAT),id1);

gamma_valid = corr_vs(max_lat_row:min_lat_row,min_lon_col:max_lon_col);
lat_vec_valid = lat_axis(max_lat_row:min_lat_row);
lon_vec_valid = lon_axis(min_lon_col:max_lon_col);

if(length(lat_vec_valid)<2)
    warning(['Invliad InSAR pair: ' workdir '\n'])
    return;
end

if(length(lon_vec_valid)<2)
    warning(['Invliad InSAR pair: ' workdir '\n'])
    return;
end

N = lat_vec_valid(1);
S = lat_vec_valid(end);
W = lon_vec_valid(1);
E = lon_vec_valid(end);
%loading the gedi points

[LAT_IN, LON_IN] = ndgrid(lat_vec_valid, lon_vec_valid);


if(not(exist(path_to_gedi_mat,'file')))
    error('Sorry! We cannot find the gedi file!');
end



if not(exist('GEDI_bigger','var'))
    GEDI_bigger = load (path_to_gedi_mat);
end
idx = and(and(GEDI_bigger.lat>=coords(2),GEDI_bigger.lat<=coords(1)),and(GEDI_bigger.lon>=coords(3),GEDI_bigger.lon<=coords(4)));




%RH100_in = GEDI_bigger.rh100(idx);
RH100_in = GEDI_bigger.rh98(idx);
lat_in = GEDI_bigger.lat(idx);
lon_in = GEDI_bigger.lon(idx);

if (isempty(RH100_in) || isempty(lat_in) || isempty(lon_in))
    warning('No available GEDI samples! Please check the GEDI file!!!');
    return
end

%nearst neiborhood is firstly carried out to determine which pixel is to be interpolated 
%N = coords(1);
%W = coords(3);
resY = lat_vec_valid(1) - lat_vec_valid(2);
resX = lon_vec_valid(2) - lon_vec_valid(1);
dimY = length(lat_vec_valid);
dimX = length(lon_vec_valid);
B = zeros(size(gamma_valid));

for ii=1:length(lat_in)
    lon_loop=lon_in(ii);
    lat_loop=lat_in(ii);
    rh100_loop=RH100_in(ii);
    idY=floor((N-lat_loop)/resY)+1;
    idX=floor((lon_loop-W)/resX)+1;
    if (idY>0)&&(idY<dimY+1)&&(idX>0)&&(idX<dimX+1)
        if ~isnan(B(idY,idX))
            B(idY,idX)=(B(idY,idX)+rh100_loop)/2;
        else
            B(idY,idX)=rh100_loop;
        end
    end
%     if mod(ii,1000)==0
%         ii,num
%     end
end

B(B==0) = nan;
idx = ~isnan(B);

lat_valid = LAT_IN(idx);
lon_valid = LON_IN(idx);



grid_rh100 = griddata(lon_in,lat_in,RH100_in,lon_valid, lat_valid,'natural');




if gedi_int_flag
    GEDI_full = griddata(lon_in, lat_in, RH100_in, LON_IN,LAT_IN,'linear');
    R_new = georasterref('RasterSize', size(GEDI_full), ...
    'RasterInterpretation', 'cells', 'ColumnsStartFrom', 'north', ...
    'LatitudeLimits', [min(LAT_IN(:)) max(LAT_IN(:))], 'LongitudeLimits', [min(LON_IN(:)) max(LON_IN(:))]);

    EPSG = 4326;



    geotiffwrite('GEDI_rh98_full.tif',GEDI_full,R_new,'CoordRefSysCode', EPSG);
%figure,imagesc(GEDI_full, 'alphaData',~isnan(gamma_valid),[0 35]),colormap('jet')

end

B_new = zeros(size(B))+nan;
B_new(idx) =  grid_rh100;
B_new_1 = B_new;
%figure, imagesc(lon_axis, lat_axis, B_new,[0 50]),colormap(turbo)
B = B_new;
GEDI = B_new;
GEDI_coords = [N S W E];
hv_coh = gamma_valid;

save(file_name_save,'GEDI','GEDI_coords','hv_coh');

cd(path_record);
end

