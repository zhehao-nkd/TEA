% a class to extract trigger info

classdef Trigger < handle
    
    properties
        equipment
        info
        plxname
        raw
        ephys_fs
    end
    
    methods
        
        function t = Trigger(path_plx)
            
            [~,t.plxname,~] = fileparts(path_plx);
            path_plx = convertStringsToChars(path_plx);
            plexonfile = readPLXFileC(path_plx,'events'); % read .plx files as a structure
            
            PULSE_LOCATION = 0;
            %
            % Find the trigger
            
            % Find the trigger plexonfile.EventChannels(1).Timestamps
            for n = 1:length(plexonfile.EventChannels)
                
                if  strcmp(plexonfile.EventChannels(n).Name, 'DIG01')
                    PULSE_LOCATION = n;
                    t.equipment = 'ZEUS';
                    t.raw = plexonfile.EventChannels(PULSE_LOCATION);
                    
                    % to judge whether to use the binary-code decoder or
                    % dig signal extractior
%                     diff_raw = diff(t.raw.Timestamps);
%                     
%                    
%                     diff_num = histcounts(diff_raw,unique(diff_raw));
%                     [diff_num,~] = sort(diff_num);
%                     sort(unique(diff_raw))
                    newdiffs = [];
                    diffs = diff(t.raw.Timestamps);
                    
                    for mm = 1: length(diffs)
                        newdiffs(mm) = round( double(diffs(mm)),-2 );
                    end
                    
                    method_choose_thres = 5000;
                    
                    newdiffs = newdiffs(newdiffs < method_choose_thres );
                    
                    if length(unique(newdiffs))== 1 % actually this is very dangerous !!!!!!!!!!!!!!!
                        t.info = t.digextract; % digital
                        clear plexonfile
                        return           
                    elseif length(unique(newdiffs))== 2
                        
                         t.info = t.digbidecode; % digital
                        clear plexonfile
                        return                     
                    end        
                end
                
            end
            
            plexonfile = readPLXFileC(path_plx,'continuous');
            
            for m = 1:length(plexonfile.ContinuousChannels)
                
                if  strcmp(plexonfile.ContinuousChannels(m).Name, 'AI01')
                    PULSE_LOCATION = m;
                    t.equipment = 'PLEXON';
                    t.raw = plexonfile.ContinuousChannels(PULSE_LOCATION);
                    t.info = t.extract; % not digital
                    clear plexonfile
                    return
                end
            end
            
            
            
        end   
        
        function time = getTime(t, pulse_number)
            if exist('pulse_number','var')
                time = t.info(pulse_number).time;
            else
                time = {t.info.time}.';
            end
        end
  
        function properties = getsingle(t, which)
            properties.equipment = t.equipment;
            properties.info = t.info(which);
            properties.plxname = t.plxname;
        end
        
        function info = digextract(t)
            
            
            fs = 30000;
            t.ephys_fs = fs;
            pulse_initials = t.raw.Timestamps;
            
            
            
            % figure; line([pulse_initials,pulse_initials],[0,5]); %for debug
            
            m = 1;
            groups ={};
            coder.varsize(groups);
            
            % As different stimuli have different number of square pulses, this works for grouping pulse clusters
            
           % WINDOW = fs; % The length of searching window, 30000samples means 1 second for Zeus
            % 07/11 2021 dynamic window
            intervals = pulse_initials(2:end) - pulse_initials(1:end-1);
            unique_intervals = unique(intervals,'sorted');
            WINDOW = unique_intervals(1) + 100; %% 100 is the redundancy for possible error  !!!! 不知道对不对
            
            for n = 1: length(pulse_initials)
                
                if n == length(pulse_initials)
                    groups{m}= [groups{m},pulse_initials(n)];

                else
                    
                    if pulse_initials(n+1)-pulse_initials(n)<WINDOW
                        if isempty(groups)
                            groups{m} = pulse_initials(n);
                        else
                            groups{m}=[ groups{m},pulse_initials(n)];
                        end
                        
                    else
                        if isempty(groups)
                            groups{m}= pulse_initials(n);
                            m = m+1;
                            groups{m} =[];
                        else
                            groups{m}= [groups{m},pulse_initials(n)];
                            m = m+1;
                            groups{m} =[];
                        end
                    end
                    
                    
                    
                end
                
            end
            
            
            
            % For collecting .time and .number
            detection.samplepoint = [];
            detection.name= [];
            coder.varsize(detection.name,detection.samplepoint);
            
            for o = 1:length(groups)
                
                detection.samplepoint(o) = groups{o}(1);
                detection.name(o) = length(groups{o});
            end
            
            detection.time = detection.samplepoint/fs;
            
            
            name = unique(detection.name);
            
            %sort triggers with same number of pulses
            
            for ii = 1:length(name) % for write name into the struct 'info'
                info(ii).name = name(ii);
            end
            
            for nn = 1:length(name)
                info(nn).repeats = length(find(detection.name==name(nn)));
                info(nn).time = detection.time(find(detection.name==name(nn)));
            end
            
            
        end
        
        function info = extract(t)
            
            fs = t.raw.ADFrequency; % avoid hard-coding;
            t.ephys_fs = fs;
            %dbstop if error
            pulse_channel = t.raw.Values;
            %detecting triggers and recording their time
            
            % debug pulse_channel = trigger.Values;
            
            MIN_THRESHOLD = 5000; % May vary in different situations
            MAX_THRESHOLD = 8500; % May vary in different situations
            
            
            pulse_channel(pulse_channel < MIN_THRESHOLD | pulse_channel > MAX_THRESHOLD) = 0; %Convert singals of pulse channel to 0 and 1
            %figure; plot(pulse_channel);
            
            pulse_channel(pulse_channel ~= 0) = 1; % based on the threshold
            %figure; plot(pulse_channel);
            
            
            
            k = find(pulse_channel);   % Find all the 1
            
            if k(1) == 1 % incase that k(1) = 1, when there are noises interrrput the pulse channel
                tempinitial = find(pulse_channel == 0);
                newinitial = tempinitial(1);
                pulse_channel = pulse_channel(newinitial:end);
                k = find(pulse_channel);
            end
            
            %    difference =  k(2:end)- k(1:end-1) ;
            %    idxDif = find(difference ~= 1);
            
            
            pulse_initials = k(pulse_channel(k - 1)==0);      % Find the initial of all the pulses
            pulse_ends = k(pulse_channel(k + 1)==0);          % Find the ends of all the pulses
            
            
            
            % Remove noises
            for n = 1:length(pulse_initials)
                
                if pulse_ends(n) - pulse_initials(n) < 10e-3 / 2 * fs % 20000 is the fs of AI channel
                    pulse_initials(n) = NaN;
                    pulse_ends(n) = NaN;
                end
                
            end
            
            
            pulse_initials (isnan(pulse_initials)) = [];
            pulse_ends (isnan(pulse_ends)) = [];
            
            
            % line([pulse_initials,pulse_initials],[0,5]); %for debug
            
            m = 1;
            groups ={};
            coder.varsize(groups);
            
            % As different stimuli have different number of square pulses, this works for grouping pulse clusters
            
            WINDOW = fs; % The length of searching window, 20000samples means 1 second
            
            
            for n = 1: length(pulse_initials)
                
                if n == length(pulse_initials)
                    groups{m}= [groups{m},pulse_initials(n)];
                    
                    
                    
                else
                    
                    if pulse_initials(n+1)-pulse_initials(n)<WINDOW
                        if isempty(groups)
                            groups{m} = pulse_initials(n);
                        else
                            groups{m}=[ groups{m},pulse_initials(n)];
                        end
                        
                    else
                        if isempty(groups)
                            groups{m}= pulse_initials(n);
                            m = m+1;
                            groups{m} =[];
                        else
                            groups{m}= [groups{m},pulse_initials(n)];
                            m = m+1;
                            groups{m} =[];
                        end
                    end
                    
                    
                    
                end
                
            end
            
            
            
            % For collecting .time and .number
            detection.samplepoint = [];
            detection.name= [];
            coder.varsize(detection.name,detection.samplepoint);
            
            for o = 1:length(groups)
                
                detection.samplepoint(o) = groups{o}(1);
                detection.name(o) = length(groups{o});
            end
            
            detection.time = detection.samplepoint/fs;
            
            
            name = unique(detection.name);
            
            %sort triggers with same number of pulses
            
            for n = 1:length(name) % for write name into the struct 'info'
                info(n).name = name(n);
            end
            
            for n = 1:length(name)
                info(n).repeats = length(find(detection.name==name(n)));
                info(n).time = detection.time(find(detection.name==name(n)));
            end
            
            
        end
        
        function info = digbidecode(t) % decode binary trigger information, which is compatible with the encoder
            
            
            fs = 30000;
            
           
            group_split_thres = 5000; % 0.1667 seconds, to split dig timestmaps to different groups for different stimuli
            raws = t.raw.Timestamps;
            diff_between_pulse = diff(raws);
            
            split_ids = find( diff_between_pulse >  group_split_thres);
            
            split_ids = [0; split_ids;length(raws)];
            
            groups = {};
            for g = 1: length(split_ids)-1
                the_first =  split_ids(g)+ 1;
                the_last = split_ids(g+ 1);
                groups{g} = raws(the_first:the_last);
            end
            
            temp = cellfun(@diff,groups,'UniformOutput',0);
            diff_collect_for_find_bi =arrayfun(@(x) round(x,-1), double(vertcat(temp{:})) );
            zero_value = min(  diff_collect_for_find_bi);
            one_value = max( diff_collect_for_find_bi);
            
            for u = 1: length(groups)
                inner_group_diff = round(diff(double(groups{u})),-1);
                inner_group_diff(inner_group_diff== zero_value) = 0;
                inner_group_diff(inner_group_diff== one_value) = 1;
                binary_code =  flip(inner_group_diff.'); % binary code
                if isempty(binary_code)
                    disp('Your encoder is wrong!!')
                end
                decimal_id(u) = bi2de(flip(binary_code));  % flip is necessar, as bi2de read in opposite direction
            end
            
            
            tris = unique(decimal_id);
            info = struct;
            
            for e = 1: length(tris)
                info(e).name = tris(e)+ 1;
                info(e).repeats = length(find(decimal_id == tris(e)));
                tri_id = find(decimal_id == tris(e));
                collects = {groups{tri_id}};
                for yy = 1: length(collects)
                    % important note here, here collects return int32, but what
                    % is necessary is double
                     info(e).time(yy) = double(collects{yy}(1))/fs;   % fs equals to 30000
                end
            end
            
        end
        
        
        
    end
    
    
    methods(Static)
        
        function binary(outdir) % add trigger channel to stimuli - Binary code
            
            dbstop if error
            
            pre_length = 0; post_length = 3; % 224
            GAP_DURATION = 1*8e-3; PULSE_DURATION = 5*8e-3; AMPLITUDE = 1;
            from_which = 1;
                
            dirpath = uigetdir();
            files = extract.filename(dirpath, '*.wav');
            files = flip(files,1);
            %outdir = 'AAAAAA'
            mkdir (outdir);
            
            for n = 1:length(files)
                
                [y,fs] = audioread(files{n}); % y means original y
                binary_code = de2bi(n); 
                zero_pulse = [zeros(GAP_DURATION*fs,1);AMPLITUDE*ones(GAP_DURATION*fs,1)];
                one_pulse = [zeros(GAP_DURATION*fs,1);zeros(GAP_DURATION*fs,1);AMPLITUDE*ones(GAP_DURATION*fs,1)];
                
                % write data channel
                yD = [zeros(pre_length*fs,1);y;zeros(post_length*fs,1)];
                
                % write trigger channel
                pulses = AMPLITUDE*ones(GAP_DURATION*fs,1); % Here the initial state of variable-pulses is a pulse
                coder.varsize(pulses); % yT means y Trigger channel
                
                for m = 1:length(binary_code)
                    if binary_code(m) == 0
                        pulses = [pulses;zero_pulse];
                    elseif binary_code(m) == 1
                        pulses = [pulses;one_pulse];
                    end
                end
                
                yT = [zeros(pre_length*fs,1);pulses;zeros(length(yD) - length(pulses) - pre_length*fs,1)]; % padding zeros
                yOut = [yT,yD];
                
                [~,name,~] = fileparts(files{n});
                audiowrite(sprintf('%s/%s-%uPulses.wav',outdir,name,from_which),yOut,fs);
                from_which = from_which + 1;
                
            end
            
        end
        
        function classic(outdir)  % add trigger channel to stimuli - Classical code
            
            %  This is a script for processing normalizede stimuli series.
            % It works by adding silent paddings together with trigger pulses to
            % another channel
            % By Zhehao Cheng, 0820,2020

            dbstop if error
            pre_length = 1; post_length = 3; % 224
            GAP_DURATION = 8e-3; PULSE_DURATION = 8e-3; AMPLITUDE = 1;
            setinitial = 1;
            initial = setinitial-1;
            
            dirpath = uigetdir();
            files = extract.filename(dirpath, '*.wav');
            files = flip(files,1);
            mkdir (outdir);
            
            for n = 1:length(files)
                [y,fs] = audioread(files{n}); % y means original y
                single_pulse = [AMPLITUDE*ones(GAP_DURATION*fs,1);zeros(GAP_DURATION*fs,1)];
                
                % write data channel
                yD = [zeros(pre_length*fs,1);y;zeros(post_length*fs,1)];
                
                % write trigger channel
                pulses = []; coder.varsize(pulses); % yT means y Trigger channel
                for m = 1:n+initial
                    pulses = [single_pulse;pulses];
                end
                
                yT = [zeros(pre_length*fs,1);pulses;zeros(length(yD) - length(pulses) - pre_length*fs,1)]; % padding zeros
                yOut = [yT,yD];
                
                [~,name,~] = fileparts(files{n});
                audiowrite(sprintf('%s/%s-%uPulses.wav',outdir,name,n+initial),yOut,fs);
                
            end

        end
        
        
    end
    
    
    
end




