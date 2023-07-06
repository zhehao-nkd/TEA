% a class to Extract trigger info

classdef Trigger < handle
    % To Extract and process trigger signals from .pl2 files
    % Trigger signal用来标定stimuli出现的时间
    properties
        equipment
        info
        pl2name
        raw
        trigger_fs
        path_pl2 % input path
        recording_time %记录时间
        
    end
    
    methods
        
        function t = Trigger(path_pl2)
            
            t.path_pl2 = path_pl2;
            [~,t.pl2name,~] = fileparts(path_pl2);
            path_pl2 = convertStringsToChars(path_pl2);
            try
                pl2info = PL2GetFileIndex(path_pl2);
            catch
                pl2info =  PL2GetFileIndex(strrep(path_pl2,'F:','D:'));
            end
            %t.recording_time = datetime(plexonfile.Date,'ConvertFrom','datenum'); 
            PULSE_LOCATION = 0;
            % find the strat and stop field
            
            if ~isempty(regexp(pl2info.CreatorSoftwareName,'Zeus')) || ~isempty(find(~arrayfun(@isempty, regexp(pl2info.FilePath,'Z\d{2}')))) % zeus
                if  isempty(vertcat(pl2info.EventChannels{:}))||isempty(find(~cellfun(@isempty, regexp( {vertcat(pl2info.EventChannels{:}).Name}.','DIG01'))))
                    code_method = 2;
                else
                    code_method = 1;
                end
            elseif ~isempty(regexp(pl2info.CreatorSoftwareName,'OmniPlex')) || ~isempty(find(~arrayfun(@isempty, regexp(pl2info.FilePath,'P\d{2}')))) % plexon
                code_method = 3;   
            end
            
            switch code_method
                case 1
                    t.coreTriggerZeusDig(path_pl2);
                case 2
                    t.coreTriggerZeusAna(path_pl2);
                case 3
                    t.coreTriggerPlexAna(path_pl2);
            end
            
        end
        
        function time = Frozen_getTime(t, pulse_number)
            if exist('pulse_number','var')
                time = t.info(pulse_number).time;
            else
                time = {t.info.time}.';
            end
        end
        
        function properties = getsingle(t, which)
            % 找到某一个stimuli对应的trigger标记的信息
            properties.equipment = t.equipment;
            properties.info = t.info(which);
            properties.pl2name = t.pl2name;
        end
        
    end
    
    methods(Access = private)
        
        function info = digextract(t)
            

            fs = 30000;
            t.trigger_fs = fs;
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
            WINDOW = unique_intervals(1) + 100/fs; %% 100 is the redundancy for possible error  !!!! 不知道对不对
            
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
                        
                    else  % add one more group
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
            
            detection.time = detection.samplepoint;
            
            
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
        
        function suminfo = Extract(t)

            [adfreq, ~, ts, frag_nums, ad] = plx_ad(t.path_pl2,'AI01'); % figure; plot(ad);
            if length(ad) <= 1 % 当AI01 channel不存在
                [adfreq, ~, ts, frag_nums, ad] = plx_ad(t.path_pl2,'AIN02');
            end

            if frag_nums(length(frag_nums))/adfreq < 45 % 有些时候frag_nums会出现问题,最后一段有一个极短的section, 比如小于三十秒
                frag_nums = frag_nums(1:end-1);
            end
            fs = adfreq; % avoid hard-coding;
            t.trigger_fs = fs;

            info_collect = {};
            counter = 0;

            edited_frag_nums = [0;frag_nums(1:end-1)]; % 暂时想不到好的命名方法了

            for k = 1:length(frag_nums) % fn 是每一个有data的section
                initial_timepoint = ts(k)*adfreq;
                values = ad(sum(edited_frag_nums(1:k)) + 1:sum(frag_nums(1:k)));
                % figure; plot(values);

                if max(values) > 0
                    pulse_channel = values;
                elseif max(values) < 0
                    pulse_channel = - values;
                end
                % figure; plot(pulse_channel); % for debug


                MIN_THRESHOLD = 25000; % May vary in different situations
                MAX_THRESHOLD = 40000; % May vary in different situations

%                 if max(pulse_channel) < MIN_THRESHOLD % Sarah's case
%                     MAX_THRESHOLD = max(pulse_channel)*1.3;
%                     MIN_THRESHOLD = max(pulse_channel)*0.7;
%                 end

                pulse_channel(pulse_channel < MIN_THRESHOLD | pulse_channel > MAX_THRESHOLD) = 0; %Convert singals of pulse channel to 0 and 1
                %figure; plot(pulse_channel);
                pulse_channel(pulse_channel ~= 0) = 1; % based on the threshold
                %figure; plot(pulse_channel);



                all_the_ones = find(pulse_channel);   % Find all the 1
                if all_the_ones(1) == 1 % incase that k(1) = 1, when there are noises interrrput the pulse channel
                    tempinitial = find(pulse_channel == 0);
                    newinitial = tempinitial(1);
                    pulse_channel = pulse_channel(newinitial:end);
                    all_the_ones = find(pulse_channel);
                end
                pulse_initials = all_the_ones(pulse_channel(all_the_ones - 1)==0);      % Find the initial of all the pulses
                pulse_ends = all_the_ones(pulse_channel(all_the_ones + 1)==0);          % Find the ends of all the pulses
                % Remove noises  这个不一定对， 可能 10e-3 要改
                for n = 1:length(pulse_initials)
                    if pulse_ends(n) - pulse_initials(n) < 10e-3 / 2 * fs % 20000 is the fs of AI channel
                        pulse_initials(n) = NaN;
                        pulse_ends(n) = NaN;
                    end
                end

                pulse_initials (isnan(pulse_initials)) = [];
                pulse_ends (isnan(pulse_ends)) = [];
                % 根据 diff的种类选择 decode 的 类型
                ini_diffs =  round(diff(pulse_initials),-1);
                group_split_thres = 3000;
                to_judge_decode_method = ini_diffs(ini_diffs < group_split_thres);


                if length(unique(to_judge_decode_method)) == 1 % use the direct decode method
                    info = struct; % 初始化info
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
                    detection.time = (detection.samplepoint + initial_timepoint)/fs;
                    name = unique(detection.name);

                    %sort triggers with same number of pulses
                    for n = 1:length(name) % for write name into the struct 'info'
                        info(n).name = name(n);
                        info(n).repeats = length(find(detection.name==name(n))); % 这三行是否适合合在一起有待商榷
                        info(n).time = detection.time(find(detection.name==name(n)));
                    end

                    counter = counter + 1;
                    info_collect{counter} = info;


                else  % use the bi 2 bi decode method 二进制

                    info = struct; % 初始化info
                    raws = pulse_initials;
                    diff_between_pulse  = round(diff(raws),-1);
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
                    decimal_id = []; % added 2022.12.11 initialize,, because there is a bug in Line 315: collects = {groups{tri_id}};
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
                        info(e).name = tris(e);
                        %info(e).name = tris(e)+ 1;
                        info(e).repeats = length(find(decimal_id == tris(e)));
                        tri_id = find(decimal_id == tris(e));
                        collects = {groups{tri_id}};
                        for yy = 1: length(collects)
                            % important note here, here collects return int32, but what
                            % is necessary is double
                            info(e).time(yy) = double(collects{yy}(1) +  initial_timepoint)/fs;   % fs equals to 30000 % 为什么我这里会写fs = 30000？ 搞不好这就能解释为什么会出现bug
                        end
                    end


                    counter = counter + 1;
                    info_collect{counter} = info;

                end


            end


            suminfo = horzcat(info_collect{:}); % 最后把所有的info拼接在一起
            % line([pulse_initials,pulse_initials],[0,5]); %for debug

        end
        
        function info = digbidecode(t) % decode binary trigger information, which is compatible with the encoder
            
            
            fs = 30000;
            t.trigger_fs = fs;
            group_split_thres = 0.1667; % 5000 samples,0.1667 seconds, to split dig timestmaps to different groups for different stimuli
            raws = t.raw.Timestamps;
            diff_between_pulse = diff(raws);
            
            split_ids = find( diff_between_pulse >  group_split_thres);
            
%             diff_between_pulse ( diff_between_pulse >  group_split_thres)
            split_ids = [0; split_ids;length(raws)];
            
            groups = {};
            for g = 1: length(split_ids)-1 % 把不同组的trigger分隔开
                the_first =  split_ids(g)+ 1;
                the_last = split_ids(g+ 1);
                groups{g} = raws(the_first:the_last);
            end
            
            temp = cellfun(@diff,groups,'UniformOutput',0);
            diff_collect_for_find_bi =arrayfun(@(x) round(x,3), double(vertcat(temp{:})) );
            zero_value = min(  diff_collect_for_find_bi);
            one_value = max( diff_collect_for_find_bi);
            
            for u = 1: length(groups)
                inner_group_diff = round(diff(double(groups{u})),3);
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
                info(e).name = tris(e);
                %info(e).name = tris(e)+ 1;
                info(e).repeats = length(find(decimal_id == tris(e)));
                tri_id = find(decimal_id == tris(e));
                collects = {groups{tri_id}};
                for yy = 1: length(collects)
                    % important note here, here collects return int32, but what
                    % is necessary is double
                    info(e).time(yy) = double(collects{yy}(1))%/fs;   % fs equals to 30000
                end
            end
            
        end
        
        function coreTriggerZeusDig(t,path_pl2)
            
            t.equipment = 'ZEUS';
           % t.raw = plexonfile.EventChannels(PULSE_LOCATION);
            t.raw.Timestamps = PL2EventTs(path_pl2,'DIG01').Ts;   
            raw_to_analysis = t.raw.Timestamps;
            
            % to judge whether to use the binary-code decoder or
            % dig signal extractior
            %                     diff_raw = diff(t.raw.Timestamps);
            %                     diff_num = histcounts(diff_raw,unique(diff_raw));
            %                     [diff_num,~] = sort(diff_num);
            %                     sort(unique(diff_raw))
            newdiffs = [];
            diffs = diff(raw_to_analysis);
            
            for mm = 1: length(diffs)
                newdiffs(mm) = round( double(diffs(mm)),3);
            end
            
            method_choose_thres = 0.3;
            newdiffs = newdiffs(newdiffs < method_choose_thres );
            
            if length(unique(newdiffs))== 1 % actually this is very dangerous !!!!!!!!!!!!!!!
                t.info = t.digextract; % digital
                clear plexonfile
            elseif length(unique(newdiffs))== 2
                t.info = t.digbidecode; % digital
                clear plexonfile % if system is zeus, the function will return here
            else

                warning('Trigger signal failed to be detected!!!! 异常');
                pause;
            end
            
        end
        
        function coreTriggerZeusAna(t,path_pl2)
            
            %             plexonfile = readpl2FileC(path_pl2,'continuous');
            %             PULSE_LOCATION = find(~cellfun(@isempty, regexp({plexonfile.ContinuousChannels.Name}.', 'AIN02')));
            %             t.equipment = 'ZEUS';
            %             t.raw = plexonfile.ContinuousChannels(PULSE_LOCATION);
            %             t.raw.Values = -t.raw.Values;
            %             t.info = t.Extract; % not digital

%             if isempty(t.raw.Values)
%                 t.raw = PL2Ad(path_pl2,'AIN02');
%             end

            
            t.equipment = 'ZEUS';

%             t.raw = PL2Ad(path_pl2,'AI01');
%             if isempty(t.raw.Values)
%                 t.raw = PL2Ad(path_pl2,'AIN02');
%             end
%             t.raw.Values = -t.raw.Values;
            t.info = t.Extract; % not digital
        end
        
        function coreTriggerPlexAna(t,path_pl2)
            
            t.equipment = 'PLEXON';
%             t.raw = PL2Ad(path_pl2,'AI01');

           
%             [gap_adfreq, gap_n, gap_ts, gap_fn] = plx_ad_gap_info(path_pl2,'AI01')
            %t.raw = plexonfile.ContinuousChannels(PULSE_LOCATION);
            t.info = t.Extract; % not digital
            
        end
        
    end
    
    
    methods(Static)
        
        
        
        function binary(outdir,from_which,pre,post) % add trigger channel to stimuli - Binary code
            
            dbstop if error
            
            if exist('pre','var')
                pre_length = pre;
            else
                pre_length = 0;
            end
            
            if exist('post','var')
                post_length = post;
            else
                post_length = 3; % default
            end
            
            
            GAP_DURATION = 1*8e-3; PULSE_DURATION = 5*8e-3; AMPLITUDE = 1;
            
            if exist('from_which','var')
                
            else
                from_which = 1;
            end
            
            dirpath = uigetdir();
            files = Extract.filename(dirpath, '*.wav');
            files = flip(files,1);
            %outdir = 'AAAAAA'
            mkdir (outdir);
            
            for n = 1:length(files)
                
                [y,fs] = audioread(files{n}); % y means original y
                binary_code = de2bi(n + from_which - 1);
                zero_pulse = [zeros(GAP_DURATION*fs,1);AMPLITUDE*ones(GAP_DURATION*fs,1)];  % 0-1
                one_pulse = [zeros(GAP_DURATION*fs,1);zeros(GAP_DURATION*fs,1);AMPLITUDE*ones(GAP_DURATION*fs,1)]; % 0-0-1
                
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
                audiowrite(sprintf('%s/%s-%02uPulses.wav',outdir,name,n + from_which - 1),yOut,fs);
                %from_which = from_which + 1;
                
            end
            
        end
        
        
        function adaptTime(outdir,from_which) % add trigger channel to stimuli - Binary code
            %如果是con song 那么就是和con song的长度一致（或者稍微再长一些）
            %如果是 deg或者repla，那么和con song的长度一致
            %如果是replacements，设为1
            %如果是individual syllables,设为1
            
            dbstop if error
            

            GAP_DURATION = 1*8e-3; PULSE_DURATION = 5*8e-3; AMPLITUDE = 1;
            
            if exist('from_which','var')
                
            else
                from_which = 1;
            end
            
            dirpath = uigetdir();
            files = Extract.filename(dirpath, '*.wav');
            files = flip(files,1);
            %outdir = 'AAAAAA'
            mkdir (outdir);
            
            for n = 1:length(files)

            %    [~,filename,~] = fileparts(files{n});

                
                [y,fs] = audioread(files{n}); % y means original y
                
                dur_y = length(y)/fs;
                
                if dur_y < 1
                    pre_length = 0.8; % minimum set to 0.8
                    post_length = 0.8;
                else
                    pre_length = dur_y;
                    post_length = dur_y;
                end
                
                binary_code = de2bi(n + from_which - 1);
                zero_pulse = [zeros(GAP_DURATION*fs,1);AMPLITUDE*ones(GAP_DURATION*fs,1)];  % 0-1
                one_pulse = [zeros(GAP_DURATION*fs,1);zeros(GAP_DURATION*fs,1);AMPLITUDE*ones(GAP_DURATION*fs,1)]; % 0-0-1
                
                % write data channel
                yD = [zeros(round(pre_length*fs),1);y;zeros(round(post_length*fs),1)];
                
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
                
                yT = [zeros(round(pre_length*fs),1);pulses;zeros(length(yD) - length(pulses) - round(pre_length*fs),1)]; % padding zeros
                yOut = [yT,yD];
                
                [~,name,~] = fileparts(files{n});
                audiowrite(sprintf('%s/%s-%02uPulses.wav',outdir,name,n + from_which - 1),yOut,fs);
                %from_which = from_which + 1;
                
            end
            
        end
        
        
        function binary_but_rand_gap(outdir,from_which,range)
            dbstop if error
            if ~exist('range','var')
                range = [0.5,0.75]; % unit is second
            end
            %pre_length = 0;
            
            GAP_DURATION = 1*8e-3; PULSE_DURATION = 5*8e-3; AMPLITUDE = 1;
            
            if exist('from_which','var')
                
            else
                from_which = 1;
            end
            
            dirpath = uigetdir();
            files = Extract.filename(dirpath, '*.wav');
            files = flip(files,1);
            %outdir = 'AAAAAA'
            mkdir (outdir);
            
            for n = 1:length(files)
                
                post_length = (range(2)-range(1))*rand + range(1);
                pre_length = (range(2)-range(1))*rand + range(1);
                
                [y,fs] = audioread(files{n}); % y means original y
                binary_code = de2bi(n + from_which - 1);
                zero_pulse = [zeros(GAP_DURATION*fs,1);AMPLITUDE*ones(GAP_DURATION*fs,1)];  % 0-1
                one_pulse = [zeros(GAP_DURATION*fs,1);zeros(GAP_DURATION*fs,1);AMPLITUDE*ones(GAP_DURATION*fs,1)]; % 0-0-1
                
                % write data channel
                yD = [zeros(round((pre_length*fs)),1);y;zeros(round(post_length*fs),1)];
                
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
                
                audiowrite(sprintf('%s/%s-%02uPulses.wav',outdir,name,n + from_which - 1),yOut,fs);
                %from_which = from_which + 1;
                
            end
            
        end
        
        function classic(outdir,pre,post)  % add trigger channel to stimuli - Classical code
            
            %  This is a script for processing normalizede stimuli series.
            % It works by adding silent paddings together with trigger pulses to
            % another channel
            % By Zhehao Cheng, 0820,2020
            
            dbstop if error
            
            if exist('pre','var')
                pre_length = pre;
            else
                pre_length = 0;
            end
            
            if exist('post','var')
                post_length = post;
            else
                post_length = 4;
            end
            
            GAP_DURATION = 8e-3; PULSE_DURATION = 8e-3; AMPLITUDE = 1;
            setinitial = 1;
            initial = setinitial-1;
            
            dirpath = uigetdir();
            files = Extract.filename(dirpath, '*.wav');
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
        
        function tobject = get_it(path_pl2)
            
            tobject = Trigger(path_pl2);
            
        end
        
        function t = shiftTrigger(path_pl2,shift_value) % shift_value in seconds
            % 此函数移动trigger signal以纠偏
            t = Trigger(path_pl2)
            % Edit t.info
            for k = 1: length(t.info)
                t.info(k).time = t.info(k).time + shift_value;
            end
            
            % Edit t.raw
            t.raw.Timestamps = t.raw.Timestamps + shift_value*t.trigger_fs;
            
        end
        
        function drawContinuousChannel(pl2_path)
            
            plexonfile = readpl2FileC(convertStringsToChars(path_pl2),'all');
            
            SPKC16 = plexonfile.ContinuousChannels(32).Values;
            
            triggerevent = plexonfile.EventChannels(1).Timestamps;
            
            % split the plexon file for drawing?
        end
        
    end
    
    
    
    
    
end




