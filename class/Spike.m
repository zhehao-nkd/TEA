classdef Spike < handle
    
    properties
        
        T
        time
        pc
        channel
        unit
        waveform
        
    end
    
    
    methods
        
        function s = Spike(T)
            
            s.T = T;
            
            s.time = s.T.timestamp; % convert the timestamps into a matrix
                %data.waveform = T{:,4:end}; % number of columns may change
                
                vname = s.T.Properties.VariableNames;% name of vareiables
                idxwaveform = find (~cellfun(@isempty, regexp(vname, '^Var\d+$')));
                s.waveform = table2array(s.T(:,idxwaveform));
                
                idxpc = find (~cellfun(@isempty, regexp(vname, '^pc\d+$'))); % idx of PCs
                s.pc = table2array(s.T(:,idxpc));
                
                idxchannel = ~cellfun(@isempty, regexp(vname, 'channelname'));
                s.channel = unique(table2array(s.T(:,idxchannel)));
                
                idxunit = ~cellfun(@isempty, regexp(vname, 'unit'));
                s.unit = unique(table2array(s.T(:,idxunit)));
            
            
        end
        
        
      
        
        
        function result = getResult(s)
            
            result = s.result;
            
        end
        
        function time = getTime(s, whichunit,range)
            
            if ~exist('whichunit','var')
                time = {s.result.time}.';
            elseif exist('range','var')
                time = s.result(whichunit).time;
                time = time(time>range(1)& time < range(2));
            else
                 time = s.result(whichunit).time;
            end      
                
            
            
        end
        
        function properties = getsingle(s,which)
            properties.T = s.separatedT{which};
            properties.result = s.result{which};
        end
        
        
            
        
    end
    
    methods(Static)
        
        function separatedT = split(path_txt)
            
            rawT = readtable(path_txt);
            
            % sorted plus unsorted
            
            % Read first line as the title
            %?????????
            fid = fopen(path_txt);
            frewind(fid);
            first_line = deblankl(fgetl(fid)); % deblankl for removal of blanks, as VariableNames can't include blanks
            titlecell = split(first_line,',');
            lentitle = length(titlecell);
            rawT.Properties.VariableNames(1:lentitle) = titlecell;
            
            
            T = rawT(startsWith(rawT.channelname,'SPKC'),:); % Remove non-SPKC channels
            if isempty(T)
                 T = rawT(startsWith(rawT.channelname,'SPK'),:);
                 warning('Dangerous !!! No SPKC detected, using SPK instead.'); % to fit jelena's data
            end
            channel_name = unique(T.('channelname')); % Get SPKC channel names
            
            separatedT = {};
            k = 1;
            for n = 1:length(channel_name)
                
                temp = T(startsWith(T.channelname,channel_name{n}),:);
                
                unit_name = unique(temp.('unit'));
                
                for m = 1:length(unit_name)
                    
                    separatedT{k}= temp((temp.unit == unit_name(m)),:);
                    k = k + 1;
                end
              %disp('MAGA!')  
                
            end
            
            separatedT = Spike.rm0(separatedT);
            
        end
        
        
        function separatedT = rm0(separatedT) % remove unit 0, which are unsorted
            
            not0 = [];
            k = 0;
            for idx = 1: length(separatedT)
            
               if unique(separatedT{idx}.unit)~=0
                   k = k + 1;
                   not0(k) = idx;
               end
            end
            
            separatedT = separatedT(not0);  % 
                   
        end
            
    end
    
    
end
