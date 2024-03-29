classdef Spike < handle
    
    properties    
        T
        time
        pc
        channel
        unit
        waveform % Here returns the waveforms
        peak
        valley
        ephys_fs
        
    end
    
    
    methods
        
        function s = Spike(T)
            
            s.T = T; % T is the imcoming table
            
            s.time = s.T.timestamp; % Convert the timestamps into a matrix
            %data.waveform = T{:,4:end}; % number of columns may change
            
            vname = s.T.Properties.VariableNames;% name of vareiables
            % Extract waveforms
            idxwaveform = find (~cellfun(@isempty, regexp(vname, '^Var\d+$')));
            tabel_contain_Var = s.T(:,idxwaveform);
            tabel_contain_Var(:,all(ismissing(tabel_contain_Var)))=[];
            s.waveform = table2array(tabel_contain_Var);
            
            % Extract PC values
            idxpc = find (~cellfun(@isempty, regexp(vname, '^pc\d+$'))); % idx of PCs
            s.pc = table2array(s.T(:,idxpc));
            
            % Extract peak
            idxpeak = find (~cellfun(@isempty, regexp(vname, 'peak'))); % idx of peaks
            s.peak = table2array(s.T(:,idxpeak));
            
             % Extract valley
            idxvalley = find (~cellfun(@isempty, regexp(vname, 'valley'))); % idx of valleys
            s.valley = table2array(s.T(:,idxvalley));
            
            % Extract channelname
            idxchannel = ~cellfun(@isempty, regexp(vname, 'channelname'));
            s.channel = unique(table2array(s.T(:,idxchannel)));
            
            % Extract unit name
            idxunit = ~cellfun(@isempty, regexp(vname, 'unit'));
            s.unit = unique(table2array(s.T(:,idxunit)));
            
            
        end
          
       
        function set_ephys_fs(s,ephys_fs)
            s.ephys_fs = ephys_fs;
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
        
        function mean_waveform(s) % calculate the mean waveform
            
            this_waveform = s.waveform;
            figure('Color','w');
            hold on
            plot(this_waveform.','Color',[.5,.5,.5]);
            plot(max(this_waveform),'--','Color','blue');
            plot(min(this_waveform),'--','Color','blue');
            plot(mean(this_waveform),'Color','red');
            set(gcf,'defaultTextInterpreter','none');
            title(sprintf('%s_%u',s.channel{1},s.unit));
            hold off
            
        end
        
        function meandiff = cal_wavelength(s)
            waveforms = s.waveform.';
           
           [valley,Ivalley] =  min (waveforms);
           
           for k = 1: length(Ivalley)
               
           [peak(k),Ipeak(k)] =  max (waveforms(Ivalley(k):end,k));
           Ipeak(k) = Ipeak(k) + Ivalley(k);
           
           end
          
           diffs = (Ipeak - Ivalley)/s.ephys_fs*1000; % ephys_fs is the sa,pling rate of ephys signal
           
           meandiff = mean(diffs); % 单位是 ms
           
            
            
        end
        
            
        
    end
    
    methods(Static)
        
        function which(path_txt) % to find out how many channel-unit exist in the path_txt file
            separatedT = Spike.split(path_txt);
            
            for k = 1: length(separatedT)
                temp = unique(separatedT{1}.channelname);
                channel = temp{1};
                temp = unique(separatedT{1}.unit);
                unit = temp;
                channel_unit = sprintf('%s-%u',channel,unit);
                disp(channel_unit);
                newline;
            end
        end
        
        
        function separatedT = split(path_txt_or_pl2) % split exported data to single channel-unit table, unit0 (unsorted) are removed
            dbstop if error
            [~,~,ext] = fileparts(path_txt_or_pl2);
            switch ext
                case '.txt'
                    separatedT = fromtxt(path_txt_or_pl2);
                case '.pl2'
                    separatedT = frompl2(path_txt_or_pl2);
            end

            function separatedT = frompl2(path_txt_or_pl2)
                pl2info = PL2GetFileIndex(path_txt_or_pl2);
                separatedT = {};
                count = 0;
                cated_spchannels = vertcat(pl2info.SpikeChannels{:});
                cated_spchannels = cated_spchannels(find(~cellfun(@isempty, regexp(cellstr({cated_spchannels.Name}.'),'SPK_SPKC\d+'))));
                for k = 1:length(cated_spchannels)
                    numunits = min(setdiff(find(cated_spchannels(k).UnitCounts == 0),[1]))-1; % can't be the first unit (unit 1), because unsorted migght be zero
                    chaname = regexp(cated_spchannels(k).Name,'SPKC\d+','match');
                    chaname = chaname{1};
                    unitsnames = 1:1:numunits-1; % unsorted (==0) not included
                    
                    collectunit = {};
                    wvswvs = {};
                    tsts = {};
                    for kk = 1:length(unitsnames)
                         wvs= PL2Waves(path_txt_or_pl2, cated_spchannels(k).Name,unitsnames(kk));
                         wvswvs{kk} = wvs.Waves;
                         tsts{kk} = wvs.Ts;
                         collectunit{kk} = repmat(unitsnames(kk), size(wvs.Waves,1),1);
                    end

                   unit = vertcat(collectunit{:});
                   if isempty(unit)
                       continue
                   end
                   timestamp = vertcat(tsts{:});
                   allwvs = vertcat(wvswvs{:});
                   [~,score,~] = pca(allwvs);

                   % calculate features
                   pc1 = score(:,1);
                   pc2 = score(:,2);
                   pc3 = score(:,3);
                   peak = max(allwvs,[],2);
                   valley = min(allwvs,[],2);

                   channelname = cellstr(repmat(convertCharsToStrings(chaname ),length(unit),1));
                   channeltable = table(channelname,unit,timestamp,pc1,pc2,pc3,peak,valley);

                   wvstable = array2table(allwvs);
                   wvstable.Properties.VariableNames = ...% rename variable names
                       cellfun(@(x)strrep(x,'allwvs','Var'),wvstable.Properties.VariableNames,'Uni',0);

                   channeltable = horzcat(channeltable ,wvstable);
                   uniqueunits = unique([channeltable.unit]);

                   for kkk = 1:length(uniqueunits)
                         count = count + 1;
                        separatedT{count} =  channeltable(find([channeltable.unit] == uniqueunits(kkk)),:);     
                   end

                end

            end

            function separatedT = fromtxt(path_txt_or_pl2)

                rawT = readtable(path_txt_or_pl2);% sorted plus unsorted


                fid = fopen(path_txt_or_pl2);
                frewind(fid);
                first_line = Utl.deblankl(fgetl(fid)); % deblankl for removal of blanks, as VariableNames can't include blanks
                titlecell = split(first_line,',');
                lentitle = length(titlecell);
                rawT.Properties.VariableNames(1:lentitle) = titlecell;


                if ~isempty(find(strcmp('channelname',rawT.Properties.VariableNames)))
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
                    end

                elseif ~isempty(find(strcmp('channel',rawT.Properties.VariableNames))) % 如果导出格式不对，没有导出channelname而是channel时
                    channel = unique(rawT.('channel'));
                    separatedT = {};
                    count = 0;
                    for n = 1:length(channel)
                        temp = rawT(rawT.channel == channel(n),:);
                        unit_name = unique(temp.('unit'));
                        for m = 1:length(unit_name)
                            count = count + 1;
                            separatedT{count}= temp((temp.unit == unit_name(m)),:);
                        end
                    end

                end

                separatedT = Spike.rm0(separatedT);  % to get the unsorted spikes,what I need is to Extract 0

            end
        end
        
        function singleChannelT = extract_specific_channel(path_txt_or_pl2,channel_name)
            

            [~,~,ext] = fileparts(path_txt_or_pl2);
            switch ext
                case '.txt'
                    singleChannelT = extractfromtxt(path_txt_or_pl2,channel_name);
                case '.pl2'
                    singleChannelT = extractfrompl2(path_txt_or_pl2,channel_name);
            end



            function singleChannelT = extractfromtxt(path_txt_or_pl2,channel_name)
                rawT = readtable(path_txt_or_pl2);% sorted plus unsorted


                % These codes are for uniforming the column name of the rawT
                fid = fopen(path_txt_or_pl2);
                frewind(fid);
                first_line = Utl.deblankl(fgetl(fid)); % deblankl for removal of blanks, as VariableNames can't include blanks
                titlecell = split(first_line,',');
                lentitle = length(titlecell);
                rawT.Properties.VariableNames(1:lentitle) = titlecell;

                % Extract specific channel's rawT information as singleChannelT
                singleChannelT = rawT(startsWith(rawT.channelname,channel_name),:);
            end
            

            function singleChannelT = extractfrompl2(path_txt_or_pl2,channel_name)
                pl2info = PL2GetFileIndex(path_txt_or_pl2);

                for k = 1:length(pl2info.SpikeChannels)
                    spchnames{k} = pl2info.SpikeChannels{k}.Name;
                end

                spchnum = find(~cellfun(@isempty, regexp(spchnames,channel_name)));

                numunits = min(setdiff(find(pl2info.SpikeChannels{spchnum}.UnitCounts == 0),[1]))-1; % can't be the first unit (unit 1), because unsorted migght be zero
                chaname = regexp(pl2info.SpikeChannels{spchnum}.Name,'SPKC\d+','match');
                chaname = chaname{1};
                unitsnames = 0:1:numunits-1; % unsorted (==0) not included

                collectunit = {};
                wvswvs = {};
                tsts = {};
                for kk = 1:length(unitsnames)
                    wvs= PL2Waves(path_txt_or_pl2, pl2info.SpikeChannels{spchnum}.Name,unitsnames(kk));
                    wvswvs{kk} = wvs.Waves;
                    tsts{kk} = wvs.Ts;
                    collectunit{kk} = repmat(unitsnames(kk), size(wvs.Waves,1),1);
                end

                unit = vertcat(collectunit{:});
                timestamp = vertcat(tsts{:});
                allwvs = vertcat(wvswvs{:});
                [~,score,~] = pca(allwvs);

                pc1 = score(:,1);
                pc2 = score(:,2);
                pc3 = score(:,3);
                peak = max(allwvs,[],2);
                valley = min(allwvs,[],2);

                channelname = cellstr(repmat(convertCharsToStrings(chaname ),length(unit),1));
                channeltable = table(channelname,unit,timestamp,pc1,pc2,pc3,peak,valley);

                wvstable = array2table(allwvs);
                wvstable.Properties.VariableNames = ...% rename variable names
                    cellfun(@(x)strrep(x,'allwvs','Var'),wvstable.Properties.VariableNames,'Uni',0);

                singleChannelT = horzcat(channeltable ,wvstable);


            end



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
        
        
        function zeroOnlyT = extract0(separatedT) % Extract unit 0, which are unsorted
            
            is0 = [];
            k = 0;
            for idx = 1: length(separatedT)
            
               if unique(separatedT{idx}.unit)== 0
                   k = k + 1;
                   is0(k) = idx;
               end
            end
            
            zeroOnlyT  = separatedT(is0);  % 
                   
        end
        
          
        function checkWaveform(path_txt)
            % Quickly check the waveform based on txt file, to decide
            % whether the neuron should be resorted
            separatedT = Spike.split(path_txt);
            for k = 1:length(separatedT)
                spikes = Spike(separatedT{k});
                spikes.mean_waveform;
            end
        end
            
    end
    
    
end
