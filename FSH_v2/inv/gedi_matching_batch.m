%gedi points matching 


rootdir = '/media/yanghai/disk5/me_10/';
path_to_gedi_mat = 'maine_2state_v2.mat';
info_root = dir(rootdir);

for i = 1:length(info_root)
    if(~contains(info_root(i).name,'_'))
        continue;
    end
    
    folder_ = [rootdir info_root(i).name '/isce_proc_dir' ];
    
    if(exist(folder_,'dir'))
        info_subdir1 = dir(folder_);
        
        
        for j = 1:length(info_subdir1)
            
            if(~contains(info_subdir1(j).name,'_'))
                continue;
            end
            
            folder_2 = [folder_ '/' info_subdir1(j).name '/forFSH'  ];
            
            
            if(exist(folder_2,'dir'))
                
                cor_input_file = [folder_2 '/cor_input.mat' ];
                
                if (exist(cor_input_file,'file'))
                    
                    [file_name_save] = GEDI_matching(folder_2,path_to_gedi_mat)
                end
                
                
            end
            
            
        end
        
    end
    
    
        
end


