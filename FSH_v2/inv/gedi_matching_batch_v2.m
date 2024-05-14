


%for processing each dir: a compact version

%rootdir = '/media/yanghai/disk4/wnmf_reproc/126_870/isce_proc_dir';
%path_to_gedi_mat = '/media/yanghai/disk2/china_gedi_folder/gedi_china_ne.mat';
path_to_gedi_mat = 'new_england_v2.mat';
info_root = dir(rootdir);

for i = 1:length(info_root)
    if(~contains(info_root(i).name,'_'))
        continue;
    end
    
    folder_ = [rootdir '/' info_root(i).name '/forFSH'  ];
    
    if(exist(folder_,'dir'))
        
        
        cor_input_file = [folder_ '/cor_input.mat' ];
        
        
        if (exist(cor_input_file,'file'))
            
            [file_name_save] = GEDI_matching_rh98(folder_,path_to_gedi_mat)
        else 
            warning(['This directory does not have a cor_input.mat file! pelase check it:\n ',folder_ ]);
        
        end
        
   
        
        
        
        
    end
    
    
    
end
