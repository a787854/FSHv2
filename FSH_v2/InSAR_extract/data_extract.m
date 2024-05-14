%put the files into corresponding folders
%this_dir = pwd;
%cmd = 'find . -name "ALPSRP*" -print';
%[s,w] = unix(cmd);

% go through the list and extract all entries
%nentries = length(find(fix(w)==10));   % this is the number of carriage returns in the returned string from the unix command
entry = dir;

for ientry= 3:length(entry)
  %[entry(ientry).name,w] = strtok(w);
 
  % get the code that relates to the row/path/pass numbers
  idex = findstr(entry(ientry).name,'ALPSRP');
  if(isempty(idex))
    continue
  end
  code = entry(ientry).name((idex+6):(idex+14));
  entry(ientry).row = str2num(code(6:9));   % the row number
  entry(ientry).path_plus = rem(str2num(code(1:5)),671);   % get the number related to the path number
  entry(ientry).pass = floor(str2num(code(1:5))/671);
  
  
  % the reminder can be found from the asf archive. for 
  if(entry(ientry).path_plus == 555)
    path = 470;
  elseif(entry(ientry).path_plus  == 621)   % 
    path = 469;
  elseif(entry(ientry).path_plus  == 373)   % resevoir
    path = 468;
  elseif(entry(ientry).path_plus == 125)    % manaus area
    path = 467;
  elseif (entry(ientry).path_plus  == 548)
    path = 466;
  elseif (entry(ientry).path_plus == 300)
    path = 465;
  elseif (entry(ientry).path_plus == 446)
    path = 471;
  elseif (entry(ientry).path_plus == 555)
    path = 117;
  elseif (entry(ientry).path_plus == 132)
    path = 118;
  elseif (entry(ientry).path_plus == 380)
    path = 119;
  elseif (entry(ientry).path_plus == 628)
    path = 120;  
  elseif (entry(ientry).path_plus == 205)
    path = 121;
  elseif (entry(ientry).path_plus == 453)
    path = 122;
  elseif (entry(ientry).path_plus == 30)
    path = 123;  
  elseif (entry(ientry).path_plus == 278)
    path = 124;
  elseif (entry(ientry).path_plus == 526)
    path = 125;
  elseif (entry(ientry).path_plus == 103)
    path = 126;  
    
   
  else
      path = 0;
  end

  
  dir_name = [num2str(path) '_' num2str(entry(ientry).row)]
  if (not(exist(dir_name,'dir')))
    mkdir(dir_name)
  end
  
  if not(exist([dir_name, entry(ientry).name ],'file'))
    movefile(entry(ientry).name,dir_name)
  end
  
% %   entry(ientry).dir = [this_dir entry(ientry).name(2:(idex-5))];
%   entry(ientry).dir = [this_dir entry(ientry).name(2:(idex(2)-5))];
%   %cmd = ['fgrep "CenterDateTime" ' entry(ientry).dir 'summary.txt'];  % get the observation date
%   cmd = ['fgrep "CenterDateTime" ' entry(ientry).dir 'workreport'];  % get the observation date
%   [s,w2] = unix(cmd);
%   entry(ientry).year = str2num(w2(26:29));
%   entry(ientry).month = str2num(w2(30:31));
%   entry(ientry).day = str2num(w2(32:33));
% %   entry(ientry).date_str = w2(26:33);
%   entry(ientry).date_str = w2(39:46);
% %   cmd = ['fgrep "ProcessLevel" ' entry(ientry).dir 'summary.txt'];  % get the observation date
% %   [s,w3] = unix(cmd);
% %   entry(ientry).proc = w3(19:21);
end