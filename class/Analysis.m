classdef Analysis < handle
    % The class to do all the analysis for a single neuron
    %   Detailed explanation goes here
    
    properties % 必须
        neurons
        list_sum
        group
        
        source    % read xlxs into struct
        list
        
        fragnames % a cell of unique frag names
        degnames
        replanames
        
        normlist % list of norm-songs
        all_eleinf
        conspe_eleinf
        
        targets % unique(targetnames)
        
        neuinfo
        uniqueid
        birdid
        %formated_imagename
        
        % ids of neuros corresponding to different stimuli
        song_id_redundant
        song_id
        song_only_id % the neuron which has only song stimuli
        frag_id
        deg_id
        repla_id
        plx_data_fs
        other_id
        
        birdname
        zpid
        channelname
        unitname
        formated_imagename
        
        insonglist
        
        fig_size1 % size of figure ( and also position)
    end
    
    properties % 可选
        figdata
    end
    
    methods % 核心方法
        
        function a = Analysis(neurons)
            % 构造方法
            if exist('neurons','var')
                
                if isa(neurons,'Neuron'); a.neurons{1} = neurons;else; a.neurons = neurons; end % 输入可以是单个Neuron或者多个neuron组成的cell
                
                a.getNeuronInfo;
                a.setStimuliCorrespondingNeuronId;
                
                % if norm stimuli exist in multiple Neuron files, then
                % delect the first one which has only norm stimuli, as it is
                % far away from others in time
                a.updatelist;
                %a.normlist = a.neurons{a.song_id}.todisplay;%a.sort;
                
                temp = regexp(a.neurons{1}.plxname,'[RBOYRG]\d{3}','match');
                % set birdid uniqueid and formated_imagename
                if ~isempty(temp)
                    a.birdid = temp{1};
                end
                a.zpid = regexp(a.neurons{1}.plxname,'[ZP]\d{2}','match');
                a.channelname = a.neurons{1}.channelname;
                a.unitname = a.neurons{1}.unitname;
                try
                    a.formated_imagename = sprintf('%s_%s_%s_%u',a.birdid,a.zpid{1},a.channelname,a.unitname);
                catch ME
                end
                a.fig_size1 = [2091 -14 755 620];
                %a.judgeFragResp;
                %a.judgeConResp;
                %a.writeFigdata;
            end
            
        end
        
        function a = updatelist(a)
            % regenerate A.list
            to_remove_id = intersect(a.song_only_id,setdiff(a.song_id_redundant,a.song_id));%???
            to_calculate = setdiff(1: length(a.neurons),to_remove_id);
            whether_update_figure_or_not = 1;
            for k = 1: length(to_calculate)
                
                templist = a.neurons{to_calculate(k)}.todisplay(whether_update_figure_or_not);
                [templist.whichNeuron] = deal(k);
                lists{k} = templist;
            end
            a.list = horzcat(lists{:});
        end
        
        function a = getNeuronInfo(a)
            %  To generate info of each member of an Analysis object
            
            a.neuinfo = struct;
            
            for k = 1: length(a.neurons)
                a.neuinfo(k).uniqueid = a.neurons{k}.uniqueid;
                a.neuinfo(k).neuronname = a.neurons{k}.neuronname;
                a.neuinfo(k).keywords = [];
                
                if length(find(~cellfun(@isempty,regexp(cellstr({a.neurons{k}.slist.name}.'),'norm|Norm|song|Song')))) > 16
                    a.neuinfo(k).keywords = [a.neuinfo(k).keywords,"song"];
                end
                
                if ~isempty(find(~cellfun(@isempty,regexp(cellstr({a.neurons{k}.slist.name}.'),'syl|Syl|Ele|ele|frag|Frag'))))
                    a.neuinfo(k).keywords = [a.neuinfo(k).keywords,"frag"];
                end
                
                if ~isempty(find(~cellfun(@isempty,regexp(cellstr({a.neurons{k}.slist.name}.'),'deg|Deg'))))
                    a.neuinfo(k).keywords = [a.neuinfo(k).keywords,"deg"];
                end
                
                
                if ~isempty(find(~cellfun(@isempty,regexp(cellstr({a.neurons{k}.slist.name}.'),'repla|Repla|catego|Catego'))))
                    a.neuinfo(k).keywords = [a.neuinfo(k).keywords,"repla"];
                end
                
                if isempty(a.neuinfo(k).keywords)
                    a.neuinfo(k).keywords = [a.neuinfo(k).keywords,"other"];
                end
                
            end
        end
        
        function a = set_eleinf(a,eleinf)
            %从外部传入eleinf这个变量
            if isa(eleinf,'struct')
                a.all_eleinf = eleinf;
            elseif isa(eleinf,'string')|| isa(eleinf,'char')
                loaded = load(eleinf);
                a.all_eleinf = loaded.all_eleinf;
            end
            
            conspe_ids = find( ~cellfun(@isempty, regexp([a.all_eleinf.songname].','CON|SPE')) );
            a.conspe_eleinf = a.all_eleinf(conspe_ids);
        end
        
        function a = splitStimuliResponsePairsToDifferentTypes(a)
            % 从list提取出norm,frag,deg,repla等几个子集sublist
            % frag
            fragidx = find(~cellfun(@isempty, regexp(cellstr({a.list(:).stimuliname}.'),'Frag|syl|Syl|Ele|ele|sim|Sim')));
            fraglist = a.list(fragidx);
            for k = 1: length(fraglist)
                
                temp = strsplit(fraglist(k).stimuliname,'-');
                
                fragnames{k} = sprintf('%s-%s-%s',temp{1},temp{2},temp{3});
            end
            
            if exist('fragnames','var')
                a.fragnames = unique(fragnames);
            end
            
            
            % deg
            degidx = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'deg|Deg')));
            deglist = a.list(degidx);
            for k = 1: length(deglist)
                
                temp = strsplit(deglist(k).stimuliname,'-');
                
                degnames{k} = sprintf('%s-%s-%s',temp{1},temp{2},temp{3});
            end
            if exist('degnames','var')
                a.degnames = unique(degnames);
            end
            
            
            % repla
            replaidx = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'Repla|repla|catego|Catego')));
            replalist = a.list(replaidx);
            for k = 1: length(replalist)
                
                temp = strsplit(replalist(k).stimuliname,'-');
                
                replanames{k} = sprintf('%s-%s-%s',temp{1},temp{2},temp{3});
            end
            if exist('replanames','var')
                a.replanames = unique(replanames);
            end
            
            % norm
            normidx = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'norm')));
            a.normlist = a.list(normidx);
            
            % target
            targetidx = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'Repla')));
            targetlist = a.list(targetidx);
            for k = 1: length(targetlist)
                
                temp = strsplit(targetlist(k).stimuliname,'-');
                
                targetnames{k} = sprintf('%s-%s-%s',temp{6},temp{7},temp{8});
            end
            
            if exist('targetnames','var')
                a.targets = unique(targetnames);
            end
            
        end
        
        function a =setStimuliCorrespondingNeuronId(a)
            % 当analysis的neuron对应不同刺激类型，如frag，norm时，阐明各neuron对应的刺激类型
            a.song_id_redundant = [];
            a.frag_id = [];
            a.deg_id = [];
            a.repla_id = [];
            a.song_only_id = [];
            for k = 1: length(a.neuinfo)
                if ismember("song",a.neuinfo(k).keywords)
                    a.song_id_redundant = [a.song_id_redundant,k];
                end
                
                if ismember("frag",a.neuinfo(k).keywords)
                    a.frag_id = [a.frag_id,k];
                end
                
                if ismember("deg",a.neuinfo(k).keywords)
                    a.deg_id = [a.deg_id,k];
                end
                
                if ismember("repla",a.neuinfo(k).keywords)
                    a.repla_id = [a.repla_id,k];
                end
                
                if strcmp("song",a.neuinfo(k).keywords)
                    a.song_only_id = [a.song_only_id,k];
                end
                
                if strcmp("other",a.neuinfo(k).keywords)
                    a.song_only_id = [a.song_only_id,k];
                end
                
                
            end
            
            if length(a.song_id_redundant) == 1
                a.song_id = a.song_id_redundant;
            elseif length(a.song_id_redundant)~=0
                
                a.song_id = min(a.song_id_redundant);
                %[~,a.song_id] = max(cellfun(@length,{a.neuinfo(a.song_id_redundant).keywords}));
                % find the id which have max length of keywords, if the
                % result is multipl ids, then select the first one
                % But maybe the last one will be much proper!
            else
                a.song_id = 1; % Very bad meaningless code
            end
            
        end
        
        function a = judgeConResp_FR(a)
            % 判断对Cons是否反应，通过 Firing rate
            % firstly update e objectys 以后可以删掉这个部分
            for k = 1: length(a.neurons)
                for kk = 1: length( a.neurons{k}.e)
                    a.neurons{k}.e{kk}.setExtAndAllocate;
                end
            end
            a.updatelist;
            
            ids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'norm|deg|repla'))); % find all norms
            % ’syl'可以兼容旧的stimuli命名规则
            
            for n = 1: length(ids)
                thisi = ids(n);
                
                presdf = cal.sdf(a.list(thisi).prejudgerespsptimes,zeros(length(a.list(thisi).judgerespy),1),a.list(thisi).fs,0.001,0.02);
                sdf = cal.sdf(a.list(thisi).judgerespsptimes,a.list(thisi).judgerespy,a.list(thisi).fs,0.001,0.02); % 0.001,0.004
                [maxpresdf,~] = max(presdf);
                [maxsdf,maxidx] = max(sdf);
                
                pre_frs = cal.eachTrialFiringRate(a.list(thisi).prejudgerespsptimes,length(a.list(thisi).judgerespy)/a.list(thisi).fs);
                sti_frs = cal.eachTrialFiringRate(a.list(thisi).judgerespsptimes,length(a.list(thisi).judgerespy)/a.list(thisi).fs);
                [h,p] = ttest(sti_frs,pre_frs,'Tail','Right','Alpha',0.05);
                a.list(thisi).pvalue = p;
                a.list(thisi).label = 0; % 初始化
                if h == 1
                    a.list(thisi).label = 1;
                    %                     if (num_of_not_empty_trials/length(a.list(thisi).pltsptimes)<0.5)||(num_of_not_empty_trials<5) % if not 60% of the trails are not empty
                    %                         a.list(thisi).label = 0;
                    %                     end
                elseif maxsdf > 17 && maxsdf > maxpresdf % a rescue
                    a.list(thisi).label = 1;
                    
                end
            end
            
        end     
        
    end
    
    methods % 内部计算方法
        
        function a = calHarmRatio(a)
            % calculate harmonic noise ratio
            
            window_size = min([a.list.leny].');
            for k = 1:length(a.list)
                a.list(k).features.harmratio = harmonicRatio(a.list(k).y,a.list(k).fs,'Window',hamming(window_size,"periodic"),...
                    'OverlapLength',round(window_size*2/3) );
                %此处为照顾很短的frag改动了window size， 但或许更好的方法是放弃很短frag的数据
                a.list(k).meanfeatures.harmratio = mean(a.list(k).features.harmratio);
            end
        end
        
        function a = writeFigdata(a)
            % 生成这个Analysis的所有three plot的图片
            figmat = {};
            for k = 1: length(a.neurons)
                figmat{k} = a.neurons{k}.writeFigdata
            end
            a.figdata = horzcat(figmat{:});
        end
        
        function a = judgeFragResp(a)
            % 判断对frag 是否反应，通过自己定义的复杂的机制
            ids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'Frag|syl'))); % find all frags
            
            % ’syl'可以兼容旧的stimuli命名规则
            DUR = 0.2 ;% 200ms
            
            
            for n = 1: length(ids)
                thisi = ids(n);
                post_sptimes = extract.sptimes(a.list(thisi).rawsptimes, a.list( thisi).zpt, a.list( thisi).zpt + DUR);
                pre_sptimes = extract.sptimes(a.list(thisi).rawsptimes, a.list(thisi).zpt-DUR, a.list(thisi).zpt );
                post_mfr = length(vertcat(post_sptimes{:}))/DUR;
                pre_mfr = length(vertcat(pre_sptimes{:}))/DUR; % mfr: mean firing rate
                a.list(thisi).rs = post_mfr - pre_mfr; % response strength
                
                tempsum = cal.psth_frag(a.list(thisi).plty,a.list(thisi).fs,a.list(thisi).pltsptimes);
                halfsum = sum(tempsum(end/2:end));
                fullsum = sum(tempsum);
                maxvalue = max(cal.psth_frag(a.list(thisi).plty,a.list(thisi).fs,a.list(thisi).pltsptimes));
                a.list(thisi).maxvalue = maxvalue;
                a.list(thisi).halfsum = halfsum;
                a.list(thisi).fullsum = fullsum;
                %if maxvalue > 6 % here the threshold is very important % originally set as 8
                if a.list(thisi).rs > 51
                    a.list(thisi).label = 1;
                else
                    a.list(thisi).label = 0;
                end
                if isempty(find([a.list.label].' == 1)) || length(find([a.list.label].' == 1)) == 1
                    if a.list(n).rs > 24
                        a.list(n).label = 1;
                    else
                        a.list(n).label = 0;
                    end
                    
                end
            end
            
        end
        
        function a = judgeFragResp_FR(a)
            % 判断对frag 是否反应，通过 Firing rate
            for k = 1: length(a.neurons)
                for kk = 1: length( a.neurons{k}.e)
                    a.neurons{k}.e{kk}.setExtAndAllocate;
                end
            end
            a.updatelist;
            
            ids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'Frag|frag|syl|ele'))); % find all norms
            % ’syl'可以兼容旧的stimuli命名规则
            
            for n = 1: length(ids)
                thisi = ids(n);
                
%                 presdf = cal.sdf(a.list(thisi).prejudgerespsptimes,zeros(length(a.list(thisi).judgerespy),1),a.list(thisi).fs,0.001,0.02);
%                 sdf = cal.sdf(a.list(thisi).judgerespsptimes,a.list(thisi).judgerespy,a.list(thisi).fs,0.001,0.02); % 0.001,0.004
%                 
                presdf = cal.sdf(a.list(thisi).prejudgerespsptimes,zeros(length(a.list(thisi).judgerespy),1),a.list(thisi).fs,0.001,0.02);
                sdf = cal.sdf(a.list(thisi).judgerespsptimes,a.list(thisi).judgerespy,a.list(thisi).fs,0.001,0.02); % 0.001,0.004
                [maxpresdf,~] = max(presdf);
                [maxsdf,maxidx] = max(sdf);
                
                
                pre_frs = cal.eachTrialFiringRate(a.list(thisi).prejudgerespsptimes,length(a.list(thisi).judgerespy)/a.list(thisi).fs);
                sti_frs = cal.eachTrialFiringRate(a.list(thisi).judgerespsptimes,length(a.list(thisi).judgerespy)/a.list(thisi).fs);
                [h,p] = ttest(sti_frs,pre_frs,'Tail','Right','Alpha',0.05)
                a.list(thisi).pvalue = p;
                a.list(thisi).label = 0; % 初始化
                if h == 1
                    a.list(thisi).label = 1;
                    %                     if (num_of_not_empty_trials/length(a.list(thisi).pltsptimes)<0.5)||(num_of_not_empty_trials<5) % if not 60% of the trails are not empty
                    %                         a.list(thisi).label = 0;
                    %
                    
                end
            end
        end
        
        function a = judgeConResp(a)
            % 判断对Cons是否反应
            % firstly update e objectys 以后可以删掉这个部分
            for k = 1: length(a.neurons)
                for kk = 1: length( a.neurons{k}.e)
                    a.neurons{k}.e{kk}.setExtAndAllocate;
                end
            end
            a.updatelist;
            
            ids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'norm|deg'))); % find all norms
            % ’syl'可以兼容旧的stimuli命名规则
            
            for n = 1: length(ids)
                thisi = ids(n);
                
                presdf = cal.sdf(a.list(thisi).prejudgerespsptimes,zeros(length(a.list(thisi).judgerespy),1),a.list(thisi).fs,0.001,0.02);
                sdf = cal.sdf(a.list(thisi).judgerespsptimes,a.list(thisi).judgerespy,a.list(thisi).fs,0.001,0.02); % 0.001,0.004
                %figure; plot(sdf);
                % figure; draw.three(a.list(thisi).judgerespy,a.list(thisi).fs,a.list(thisi).judgerespsptimes);
                % figure; draw.three(a.list(thisi).plty,a.list(thisi).fs,a.list(thisi).pltsptimes);
                [maxpresdf,~] = max(presdf);
                [maxsdf,maxidx] = max(sdf);
                percentage_max = maxidx/length(sdf);
                time_max = length(a.list(thisi).judgerespy)/a.list(thisi).fs*percentage_max;
                % check whether the surroding are has spikes in most of the
                % trials
                extracted_sptimes = extract.sptimes(a.list(thisi).judgerespsptimes,time_max - 0.15, time_max + 0.15); % 前后 100ms
                num_of_not_empty_trials = length(find(~cellfun(@isempty, extracted_sptimes)));
                
                % mean_maxsdf = maxsdf/length(a.list(thisi).judgerespsptimes);
                
                a.list(thisi).label = 0; % 初始化
                if (maxsdf) > 17 && maxsdf > maxpresdf %如果是 time-locked response
                    a.list(thisi).label = 1;
                    
                    if (num_of_not_empty_trials/length(a.list(thisi).pltsptimes)<0.5)||(num_of_not_empty_trials<5) % if not 60% of the trails are not empty
                        a.list(thisi).label = 0;
                    end
                    
                elseif mean(sdf)> 9*mean(presdf) && mean(sdf)>0.6 % 如果不是 time-locked response
                    a.list(thisi).label = 1;   %  set to 2 ,biao ming shi fei time-locked response
                end
                
                
                
                
                
                %                 figure;
                %
                %                 draw.three(a.list(thisi).plty,a.list(thisi).fs,a.list(thisi).pltsptimes);
                %                 title(sprintf('Label is %u',a.list(thisi).label) );
                %                 close(gcf)
                
            end
            
        end
        
        function  a = judgeReplaResp(a)
            
            % evaluate the responsiveness of song replacements
            replaids = find( ~cellfun(@isempty, regexp([a.list.stimuliname].','Repla') ));
            
            % Find out that the repla stimuliname correspond to how many natrual songs
            bnames = unique(cellfun(@(x)convert.bid(x,2),[a.list(replaids).stimuliname].','Uni',0));
            
            findCorrespConID = @(x) intersect(find( ~cellfun(@isempty, regexp([a.list.stimuliname].','norm') )),...
                find( ~cellfun(@isempty, regexp([a.list.stimuliname].',x) )) ); % function handle to get corresponding norm ids
            
            findCorrespReplaID = @(x) intersect(find( ~cellfun(@isempty, regexp([a.list.stimuliname].','Repla') )),...
                find( ~cellfun(@isempty, regexp([a.list.stimuliname].',x) )) ); % function handle to get corresponding norm ids
            
            for k = 1:length(bnames)
                nid = findCorrespConID(bnames{k}); % 有可能得到多个nid
                % 如果有多个nid
                if length(nid) > 1
                    nid = nid(1);
                    disp('Dangerous Code: @Analysis.judgeReplaResp');
                end
                rids = findCorrespReplaID(bnames{k});
                
                for kk = 1:length(rids)
                    [Ini_y,Ini_replay] = Analysis.findConergentPointBetwenNormAndRepla(...
                        a.list(nid).y,...
                        a.list(rids(kk)).y );
                    a.list(rids(kk)).targety = a.list(rids(kk)).y(Ini_replay:Ini_replay + 0.1*32000); % 截取100ms
                    a.list(rids(kk)).targetsptimes = extract.sptimes_resetSP(a.list(rids(kk)).sptimes,Ini_replay,Ini_replay + 0.1*32000);
                    % corresponding pre (targety) data
                    a.list(rids(kk)).pretargety = zeros(length(a.list(rids(kk)).targety),1);%a.list(nid).prey(end - 0.1*32000:end); % 截取100ms
                    a.list(rids(kk)).pretargetsptimes = extract.sptimes_resetSP(...
                        a.list(rids(kk)).presptimes,length(a.list(rids(kk)).pretargety)/32000 -0.1*32000,length(a.list(rids(kk)).pretargety)/32000 );%length(a.list(rids(kk)).prey)/32000 -0.1*32000
                    %calculate and judege whether the neuron respond to the target area or not
                    [a.list(rids(kk)).replaresp,a.list(rids(kk)).replaPvalue] = Analysis.UseTtestToJudegeRespOrNot(...
                        a.list(rids(kk)).targety,a.list(rids(kk)).targetsptimes,...
                        a.list(rids(kk)).pretargety ,a.list(rids(kk)).pretargetsptimes ,32000);
                end
                a.list(nid).target = a.list(nid).y(Ini_y:Ini_y + 0.1*32000); % 截取100ms
                a.list(nid).targetsptimes = extract.sptimes_resetSP(a.list(nid).sptimes,Ini_y,Ini_y + 0.1*32000);
               
                
                
            end
            
            
        end
        
        function a = judgeDegResp(a)
            % evaluate the responsiveness of song degressive deletion
            
            % based on the name of replaced song, find out the corresponding % norm song
        end
    end
    
    methods % 外部计算方法
        
        % 计算部分
        function e_objects = getAllEphysObject(a)
            %提取Analysis所有的Ephys Object
            collects = {};
            for k = 1: length(a.neurons)
                collects{k} = a.neurons{k}.e(:);
            end
            e_objects = vertcat(collects{:});
        end
        
        function formated_imagename = neuronname(a)
            % 似乎没用
            formated_imagename = sprintf('%s_%u',a.birdid,a.uniqueid);
        end
        
        function reordered_ids = reorderListByResp(a)
            % reorder List By Response
            fragids = find(~cellfun(@isempty, regexp({a.list.stimuliname}.','Frag|syl|Syl'))); % find all frags,兼容 syl|Syl
            if ~isempty( fragids)
                fraglist = a.list(fragids);
                for n = 1: length(fraglist)
                    tempsum = cal.psth_frag(fraglist(n).rawy,fraglist(n).fs,fraglist(n).rawsptimes);
                    halfsum = sum(tempsum(end/2:end));
                    fullsum = sum(tempsum);
                    maxvalue = max(cal.psth_frag(fraglist(n).rawy,fraglist(n).fs,fraglist(n).rawsptimes));
                    fraglist(n).maxvalue = maxvalue;
                    fraglist(n).halfsum = halfsum;
                    fraglist(n).fullsum = fullsum;
                    if maxvalue > 6 % here the threshold is very important % originally set as 8
                        fraglist(n).label = 1;
                    else
                        fraglist(n).label = 0;
                    end
                end
                [~,fragindex] = sortrows( struct2table(fraglist) ,'maxvalue','descend');
                reordered_fragids = fragids(fragindex);
            else
                reordered_fragids = [];
            end
            
            replaids = find(~cellfun(@isempty, regexp({a.list.stimuliname}.','Repla|catego'))); % find all frags,兼容catego
            if ~isempty( replaids)
                replalist = a.list(replaids);
                for n = 1: length(replalist)
                    tempsum = cal.psth_frag(replalist(n).rawy,replalist(n).fs,replalist(n).rawsptimes);
                    halfsum = sum(tempsum(end/2:end));
                    fullsum = sum(tempsum);
                    maxvalue = max(cal.psth_frag(replalist(n).rawy,replalist(n).fs,replalist(n).rawsptimes));
                    replalist(n).maxvalue = maxvalue;
                    replalist(n).halfsum = halfsum;
                    replalist(n).fullsum = fullsum;
                    if maxvalue > 6 % here the threshold is very important % originally set as 8
                        replalist(n).label = 1;
                    else
                        replalist(n).label = 0;
                    end
                end
                [~,replaindex] = sortrows( struct2table(replalist) ,'maxvalue','descend');
                reordered_replaids = replaids(replaindex);
            else
                reordered_replaids = [];
            end
            
            normids = find(~cellfun(@isempty, regexp({a.list.stimuliname}.','norm'))); % find all frags
            if ~isempty( normids )
                normlist = a.list( normids);
                for n = 1: length( normlist)
                    tempsum = cal.psth_frag( normlist(n).rawy, normlist(n).fs, normlist(n).rawsptimes);
                    halfsum = sum(tempsum(end/2:end));
                    fullsum = sum(tempsum);
                    maxvalue = max(cal.psth_frag( normlist(n).rawy, normlist(n).fs, normlist(n).rawsptimes));
                    normlist(n).maxvalue = maxvalue;
                    normlist(n).halfsum = halfsum;
                    normlist(n).fullsum = fullsum;
                    if maxvalue > 6 % here the threshold is very important % originally set as 8
                        normlist(n).label = 1;
                    else
                        normlist(n).label = 0;
                    end
                end
                
                [~, normindex] = sortrows( struct2table( normlist) ,'maxvalue','descend');
                reordered_normids =  normids( normindex);
            else
                reordered_normids = [];
            end
            
            
            
            otherids = find(cellfun(@isempty, regexp({a.list.stimuliname}.','Frag|Repla|norm|catego|Syl|syl'))); % find all frags
            if ~isempty(otherids)
                otherlist = a.list(otherids);
                for n = 1: length(otherlist)
                    tempsum = cal.psth_frag( otherlist(n).rawy,otherlist(n).fs,otherlist(n).rawsptimes);
                    halfsum = sum(tempsum(end/2:end));
                    fullsum = sum(tempsum);
                    maxvalue = max(cal.psth_frag(otherlist(n).rawy,otherlist(n).fs,otherlist(n).rawsptimes));
                    otherlist(n).maxvalue = maxvalue;
                    otherlist(n).halfsum = halfsum;
                    otherlist(n).fullsum = fullsum;
                    if maxvalue > 6 % here the threshold is very important % originally set as 8
                        otherlist(n).label = 1;
                    else
                        otherlist(n).label = 0;
                    end
                end
                
                
                [~, otherindex] = sortrows( struct2table(otherlist) ,'maxvalue','descend');
                reordered_otherids =  otherids( otherindex);
                
            else
                reordered_otherids = [];
            end
            
            reordered_ids = [reordered_fragids;reordered_replaids;reordered_otherids;reordered_normids];
            
        end
        
        function meanWL = calMeanWaveLength(a)
            % 计算 meanWL，需要考虑到仪器的fs
            wavecollect = {};
            for s = 1: length(a.neurons)
                wavecollect{s} = a.neurons{s}.waveform;
            end
            waveforms = vertcat(wavecollect{:});
            %waveforms =  n.waveform;
            [~,troughstime] = min(waveforms,[],2);
            wavlen_units = [];
            
            for k = 1: size(waveforms,1) % length is dangerous!!!!!
                this_wf = waveforms(k,:);
                [~,wavlen_units(k)] =  max(this_wf (troughstime(k):end));
            end
            
            if ~isempty(regexp(a.zpid,'Z')) %zeus
                a.plx_data_fs = 30000; %hard code !!!!!! Dangerous
            elseif ~isempty(regexp(a.zpid,'P')) % plexon
                a.plx_data_fs = 40000;
            end
            
            meanWL =  mean(wavlen_units*(1/a.plx_data_fs)*1000); % ms
            
            
        end
        
        function [localSFR,h,p] = getSponFR(a,range)
            % calculate spontaneous firing rate
            sponFrInfo = struct;
            all_es = a.getAllEphysObject;
            for m = 1: length(all_es)
                
                
                % for prey
                sponFrInfo(m).triggerNum = all_es{m}.sound.trigger;
                sponFrInfo(m).presptimes = all_es{m}.presptimes
                sponFrInfo(m).preylen = length(all_es{m}.y)/all_es{m}.fs;
                sponFrInfo(m).repnum = size(all_es{m}.presptimes,2);
                temp = all_es{m}.presptimes.';
                sponFrInfo(m).localSpFr = length(find(vertcat(vertcat(temp{:}))))/(sponFrInfo(m).preylen*sponFrInfo(m).repnum);
                % for plty
                sponFrInfo(m).pltsptimes = all_es{m}.pltsptimes
                sponFrInfo(m).pltlen = length(all_es{m}.plty)/all_es{m}.fs;
                
            end
            
            
            localSFR = {};
            if exist('range','var')
                
                for k = 1: length(range)
                    
                    if k < length(range)
                        ids_in_range = intersect(find(range(k) <=[sponFrInfo.triggerNum].'), find( [sponFrInfo.triggerNum].' <range(k + 1)))
                    elseif k == length(range)
                        ids_in_range = find(range(k) <=[sponFrInfo.triggerNum].');
                    end
                    
                    
                    selected_sponFrInfo = sponFrInfo(ids_in_range);
                    
                    localSFR{k} = [selected_sponFrInfo.localSpFr].'
                    
                end
            end
            %             % for pre_y
            %             sponFrInfo(k).concat_pre_sptimes = concat_presptimes;
            %             sponFrInfo(k).concat_pre_len = sum_prelen;
            %             sponFrInfo(k).mean_pre_fr = length(concat_presptimes)/sum_prelen;
            %
            %             % for plt_y
            %             sponFrInfo(k).concat_plt_sptimes = concat_pltsptimes;
            %             sponFrInfo(k).concat_plt_len = sum_pltlen;
            %             sponFrInfo(k).mean_plt_fr = length(concat_pltsptimes)/sum_pltlen;
            [h,p]= ttest2(localSFR{1},localSFR{2});
            
        end
        
        
        function sponFrInfo = getSponFR_Sarah(a,range)
            % calculate spontaneous firing rate
            sponFrInfo = struct;
            all_es = a.getAllEphysObject;
            for m = 1: length(all_es)
                
                
                % for prey
                sponFrInfo(m).triggerNum = all_es{m}.sound.trigger;
                sponFrInfo(m).presptimes = all_es{m}.presptimes
                sponFrInfo(m).preylen = length(all_es{m}.y)/all_es{m}.fs;
                sponFrInfo(m).repnum = size(all_es{m}.presptimes,2);
                temp = all_es{m}.presptimes.';
                sponFrInfo(m).localSpFr = length(find(vertcat(vertcat(temp{:}))))/(sponFrInfo(m).preylen*sponFrInfo(m).repnum);
                % for plty
                sponFrInfo(m).pltsptimes = all_es{m}.pltsptimes
                sponFrInfo(m).pltlen = length(all_es{m}.plty)/all_es{m}.fs;
                
            end
            
            
            %             localSFR = {};
            %             if exist('range','var')
            %
            %                 for k = 1: length(range)
            %
            %                     if k < length(range)
            %                         ids_in_range = intersect(find(range(k) <=[sponFrInfo.triggerNum].'), find( [sponFrInfo.triggerNum].' <range(k + 1)))
            %                     elseif k == length(range)
            %                         ids_in_range = find(range(k) <=[sponFrInfo.triggerNum].');
            %                     end
            %
            %
            %                     selected_sponFrInfo = sponFrInfo(ids_in_range);
            %
            %                     localSFR{k} = [selected_sponFrInfo.localSpFr].'
            %
            %                 end
            %             end
            %             % for pre_y
            %             sponFrInfo(k).concat_pre_sptimes = concat_presptimes;
            %             sponFrInfo(k).concat_pre_len = sum_prelen;
            %             sponFrInfo(k).mean_pre_fr = length(concat_presptimes)/sum_prelen;
            %
            %             % for plt_y
            %             sponFrInfo(k).concat_plt_sptimes = concat_pltsptimes;
            %             sponFrInfo(k).concat_plt_len = sum_pltlen;
            %             sponFrInfo(k).mean_plt_fr = length(concat_pltsptimes)/sum_pltlen;
            % [h,p]= ttest2(localSFR{1},localSFR{2});
            
        end
        
        
        function fraglist = to2ndAcousticSpace(a)
            % 2nd acoustic space？
            ids = find(~cellfun(@isempty, regexp({a.list.stimuliname}.','repla'))); % find all frags
            
            fraglist = a.list(ids);
            
            for n = 1: length(fraglist)
                tempsum = cal.psth_frag(fraglist(n).rawy,fraglist(n).fs,fraglist(n).rawsptimes);
                range = 1 % the very fisrt 1 second
                beginmax = max(tempsum(1: ceil(length(tempsum)*range/(length(fraglist(n).rawy)/fraglist(n).fs)) ));% the maximum value of the begining 0.5 second
                halfsum = sum(tempsum(end/2:end));
                fullsum = sum(tempsum);
                maxvalue = max(cal.psth_frag(fraglist(n).rawy,fraglist(n).fs,fraglist(n).rawsptimes));
                fraglist(n).maxvalue = maxvalue;
                fraglist(n).halfsum = halfsum;
                fraglist(n).fullsum = fullsum;
                fraglist(n).beginmax = beginmax;
                if maxvalue > 8 % here the threshold is very important
                    fraglist(n).label = 1;
                else
                    fraglist(n).label = 0;
                end
            end
            
        end
        
        function split_fraginf = Deprecated_calResponseToWithinSongFragsFromEleinf(a)
            % deprecated
            latency = 50*0.001; % 50 ms
            
            split_fraginf = struct;
            
            counts = 0;
            if isempty(a.normlist)
                split_fraginf = struct([]);
                return
            end
            
            for k = 1: length(a.normlist)
                %songname = regexp(a.normlist(k).stimuliname,'(CON|SPE|norm)-([BRGOY]\d{3}|Fcall|Mcall|WNS|HET|Het)','match');
                songname = regexp(a.normlist(k).stimuliname,'([BRGOY]\d{3}|Fcall|Mcall|WNS|HET|Het)','match');
                
                % use autoseg to segment the elements
                %                [rawy,fiy,I,syledge,eleedge] = main(fiy,fs,birdid,CONFIG)
                %
                eleids = find(~cellfun(@isempty, regexp([a.conspe_eleinf.songname].',songname)));
                local_eleinf = a.conspe_eleinf(eleids);
                
                %                for w = 1: length(local_eleinf)
                %                    local_eleinf(w).yini = local_eleinf(w).initial - local_eleinf(1).initial;
                %                    local_eleinf(w).yter = local_eleinf(w).terminal - local_eleinf(1).initial;
                %                end
                
                % 为了方便命名和统一形式
                for w = 1: length(local_eleinf)
                    local_eleinf(w).yini = local_eleinf(w).initial - 0;
                    local_eleinf(w).yter = local_eleinf(w).terminal - 0;
                end
                
                
                for m = 1:length(local_eleinf)
                    counts = counts + 1;
                    split_fraginf(counts).y = local_eleinf(m).y;
                    split_fraginf(counts).fs = local_eleinf(m).fs;
                    split_fraginf(counts).padded_y = [local_eleinf(m).y;zeros(latency*local_eleinf(m).fs,1)]; % pad zeors with latency length
                    
                    split_fraginf(counts).initial = local_eleinf(m).yini/local_eleinf(m).fs;
                    split_fraginf(counts).terminal = local_eleinf(m).yter/local_eleinf(m).fs;
                    
                    split_fraginf(counts).padded_sptimes = extract.sptimes(a.normlist(k).rawsptimes, split_fraginf(counts).initial...
                        ,split_fraginf(counts).terminal + latency);
                    
                    front_percentage = local_eleinf(m).yini/length(a.normlist(k).rawy);
                    back_percentage = local_eleinf(m).yter/length(a.normlist(k).rawy);
                    
                    featurenames = fieldnames(a.normlist(k).rawfeatures); % iterate for each feature
                    featurenames = setdiff(featurenames,{'file_index','file_name'});
                    for omega = 1: length(featurenames)
                        lenfeature = length(a.normlist(k).rawfeatures.(featurenames{omega}));
                        if  round(front_percentage*lenfeature)== 0
                            
                            split_fraginf(counts).(featurenames{omega}) = a.normlist(k).rawfeatures.(featurenames{omega})(...
                                1:round(back_percentage*lenfeature));
                        else
                            split_fraginf(counts).(featurenames{omega}) = a.normlist(k).rawfeatures.(featurenames{omega})(...
                                round(front_percentage*lenfeature):round(back_percentage*lenfeature));
                        end
                        split_fraginf(counts).(sprintf('mean_%s',featurenames{omega})) = mean(split_fraginf(counts).(featurenames{omega}));
                        
                    end
                    
                end
                
            end
            
            for n = 1: length(split_fraginf)
                tempsum = cal.psth_frag(split_fraginf(n).padded_y,split_fraginf(n).fs,split_fraginf(n).padded_sptimes);
                if isempty(tempsum)
                    halfsum = 0;
                else
                    halfsum = sum(tempsum(end/2:end));
                end
                fullsum = sum(tempsum);
                maxvalue = max(cal.psth_frag(split_fraginf(n).padded_y,split_fraginf(n).fs,split_fraginf(n).padded_sptimes));
                split_fraginf(n).maxvalue = maxvalue;
                split_fraginf(n).halfsum = halfsum;
                split_fraginf(n).fullsum = fullsum;
                if maxvalue > 6 % here the threshold is very important % originally set as 8
                    split_fraginf(n).label = 1;
                else
                    split_fraginf(n).label = 0;
                end
            end
            
            
        end
        
        function insonglist = getInsongFragRespList(a)
            % get In-songFrag Response List
            cons_neuron = a.neurons{a.song_id};
            insonglist = cons_neuron.todisplayInsong;
            
        end
        
        function ainf = calPropertiesForClustering(a)
            % calculate lifetime sparseness,correlation index, number of
            % responsive songs, spontaneous firing rate, spike width
            % 目的是通过计算这些性质的值对neurons进行划分
            
            % number of responsive songs
            
            %             for kk = 1: length( a.neurons{a.song_id}.e)
            %                 a.neurons{a.song_id}.e{kk}.setExtAndAllocate;
            %             end
            %             a.updatelist;
            
            
            conkeywords = {'B346','B512','B521','B554','B606','G429','G506','G518','G548','G573',...
                'G578','O331','O507','O509','O540','Y515','Y606','Y616'};
            conids = [];
            for kk = 1: length(conkeywords)
                ids = find(~cellfun(@isempty,regexp(cellstr({a.list.stimuliname}.'),['norm\S+(?!(TUT|BOS))',conkeywords{kk},'(?!(TUT|BOS))'] )));
                if length(ids) == 1
                    conids(kk) = ids;
                elseif length(ids) >1
                    conids(kk) = ids(1);
                elseif length(ids) == 0
                    conids(kk) = nan;
                end
            end
            conids = rmmissing(conids);
            normids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'norm')));
            
            normlist = a.list(intersect(conids,normids));
            
            %thres = 0.001; % 1ms
            thres = 0.001;
            sum_prelen = 0; % summed prey length
            concat_presptimes = []; % concatenated prey sptimes
            
            sum_judgeresplen = 0; %summed prey( stimuli y, not plty or rawy) length
            concat_judgerespsptimes = []; %  % concatenated y sptimes
            
            sumNs = [];
            sdf_collect = {};
            for m = 1: length(normlist)
                if isempty(normlist(m).judgerespsptimes)
                    continue
                end
                
                jrsptimes = normlist(m).judgerespsptimes;
                all_spikes = vertcat(jrsptimes{:});
                all_Ns = length(find(abs(cal.allPairDiff(all_spikes))<thres));
                same_trail_Ns = [];
                for k = 1: length(jrsptimes)
                    same_trail_Ns(k) = length(find(abs(cal.allPairDiff(jrsptimes{k}))<thres));
                end
                Ns = all_Ns - sum(same_trail_Ns);
                sumNs(m) = Ns;
                M = length(jrsptimes); % number of presentation
                D = length(normlist(m).judgerespy)/normlist(m).fs; % stimulus duration
                r = length(all_spikes)/length(normlist(m).judgerespy);  % average firing rate
                omega = thres; % coincidence window
                normalization_factor = M*(M-1)*(r.^2)*omega*D;
                normlist(m).CI = Ns/normalization_factor;
                ainf.eachCI(m) = Ns/normalization_factor;
                
                sdf = cal.sdf(jrsptimes,normlist(m).judgerespy,normlist(m).fs,0.001,0.004);
                sdf_collect{m} = sdf;
                minsdf =min(sdf);
                maxsdf = max(sdf);
                hundredthres = linspace(minsdf,maxsdf,50);  % or 102
                fraction_above = [];
                for k = 1: length(hundredthres)
                    fraction_above(k) = length(find(sdf>hundredthres(k)))/length(sdf);
                end
                Avalue = trapz(hundredthres,fraction_above);
                normlist(m).sparseness = 1 - 2*Avalue;
                
                
                % for prey
                ainf.presptimes{m} = normlist(m).prejudgerespsptimes;
                ainf.preylen{m} = length(normlist(m).y)/normlist(m).fs;
                ainf.repnum{m} = size(normlist(m).prejudgerespsptimes,2);
                temp = normlist(m).prejudgerespsptimes.';
                concat_presptimes = [concat_presptimes;vertcat(vertcat(temp{:}))+ sum_prelen];
                sum_prelen = sum_prelen +  ainf.preylen{m};
                
                % for plty
                ainf.judgerespsptimes{m} = normlist(m).judgerespsptimes
                ainf.judgeresplen{m} = length(normlist(m).judgerespy)/normlist(m).fs;
                temp = normlist(m).judgerespsptimes.';
                %concat_judgerespsptimes = [concat_judgerespsptimes;vertcat(vertcat(temp{:}))+ sum_judgeresplen];
                concat_judgerespsptimes = [concat_judgerespsptimes; cellfun(@(x) x+sum_judgeresplen,temp,'Uni',0) ];
                sum_judgeresplen = sum_judgeresplen +  ainf.judgeresplen{m};
                
                
                % judge whether significant repsonse or not
                
                presdf = cal.sdf(normlist(m).prejudgerespsptimes,zeros(length(normlist(m).judgerespy),1),normlist(m).fs,0.001,0.02);
                sdf = cal.sdf(normlist(m).judgerespsptimes,normlist(m).judgerespy,normlist(m).fs,0.001,0.02); % 0.001,0.004
                [maxpresdf,~] = max(presdf);
                [maxsdf,~] = max(sdf);
                [minsdf,~] = min(sdf);
                
                % calculate sparseness for each song
                hundredthres = linspace(minsdf,maxsdf,100);  % or 102
                fraction_above = [];
                for k = 1: length(hundredthres)
                    fraction_above(k) = length(find(sdf>hundredthres(k)))/length(sdf);
                end
                
                
                
                Avalue = trapz(hundredthres,fraction_above)/100;
                ainf.eachsparseness(m) = 1 - 2*Avalue;
                
                
                
                
                pre_frs = cal.eachTrialFiringRate(normlist(m).prejudgerespsptimes,length(normlist(m).judgerespy)/normlist(m).fs);
                sti_frs = cal.eachTrialFiringRate(normlist(m).judgerespsptimes,length(normlist(m).judgerespy)/normlist(m).fs);
                [h,p] = ttest(sti_frs,pre_frs,'Tail','Right','Alpha',0.05)
                normlist(m).pvalue = p;
                normlist(m).label = 0; % 初始化
                if h == 1
                    normlist(m).label = 1;
                elseif maxsdf > 17 && maxsdf > maxpresdf % a rescue
                    normlist(m).label = 1;
                    
                end
                
                
            end
            
            
            allsdf = horzcat(sdf_collect{:});
            ainf.allsdf = allsdf;
            is1ids = find([normlist.label].' == 1);
            ainf.numrespsong = length(is1ids);
            label1normlist = normlist(is1ids);
            
            ainf.neuronname = a.formated_imagename;
            ainf.meanCI = mean([label1normlist.CI].');
            ainf.maxCI = max([label1normlist.CI].');
            ainf.minCI = min([label1normlist.CI].');
            ainf.meansparseness = mean([label1normlist.sparseness].');
            ainf.maxsparseness = max([normlist.sparseness].');
            ainf.forpca = horzcat(ainf.eachsparseness,ainf.eachCI);
            %
            
            % for norm songs, degressive songs, and detailed
            % frags/replas, calculate the concatenated firing rate one
            % by one
            % for pre_y
            ainf.concat_pre_sptimes = concat_presptimes;
            ainf.concat_pre_len = sum_prelen;
            ainf.mean_pre_fr = length(concat_presptimes)/sum_prelen;
            
            % for plt_y
            ainf.concat_judgeresp_sptimes = vertcat(concat_judgerespsptimes{:});
            ainf.concat_judgeresp_len = sum_judgeresplen;
            ainf.mean_judgeresp_fr = length(ainf.concat_judgeresp_sptimes)/sum_judgeresplen;
            ainf.meanWL = a.calMeanWaveLength;
            
            sumM = max([length(normlist(1).sptimes),length(normlist(2).sptimes),length(normlist(3).sptimes)]); % number of presentation
            sumD = sum_judgeresplen; % stimulus duration
            sum_allspikes = ainf.concat_judgeresp_sptimes;
            sumr = length(sum_allspikes)/sum_judgeresplen;  % average firing rate
            spikenum = length(sum_allspikes);
            ainf.spikenum = spikenum;
            ainf.afr =sumr;% average firing rate
            sumomega = thres; % coincidence window
            sumnormalization_factor = sumM*(sumM-1)*(sumr.^2)*sumomega*sumD;
            ainf.sumCI = sum(sumNs)/sumnormalization_factor;
            
            
            % calculate sum sdf
            %             sumsdf = cal.sdf(concat_judgerespsptimes,zeros(sum_judgeresplen*normlist(m).fs,1),normlist(m).fs,0.001,0.004);
            %             [sumsdf,~] = histcounts(ainf.concat_judgeresp_sptimes ,round(sumD/0.001));
            %
            sumsdf = allsdf;
            minsumsdf =min(sumsdf);
            maxsumsdf = max(sumsdf);
            sumhundredthres = linspace(minsumsdf,maxsumsdf,100);  % or 102
            sumfraction_above = [];
            for k = 1: length(sumhundredthres)
                sumfraction_above(k) = length(find(sumsdf>sumhundredthres(k)))/length(sumsdf);
            end
            
            devisdf = std(sumsdf);
            avgsdf = mean(sumsdf);
            
            summer = [];
            
            for k = 1: length(sumsdf)
                summer(k) = ((sumsdf(k) - avgsdf)/devisdf)^4;
            end
            
            
            ainf.kurtosis = sum(summer)/length(sumsdf) -3;
            
            sumAvalue = trapz(sumhundredthres,sumfraction_above)/100;
            ainf.sumsparseness = 1 - 2*sumAvalue;
            
            
        end
        
        function fr_info = multiRepeatsFiringRate(a)
            % 计算各种定义下，平均所有stimuli后的firing rate
            fr_info = struct;
            fr_info.neuronname = a.formated_imagename;
            
            sum_prelen = 0; % summed prey length
            concat_presptimes = []; % concatenated prey sptimes
            
            sum_pltlen = 0; %summed prey( stimuli y, not plty or rawy) length
            sum_ylen = 0;
            concat_pltsptimes = []; %  % concatenated y sptimes
            concat_ysptimes = [];
            all_es = a.getAllEphysObject;
            
            for m = 1: length(all_es)
                % for prey
                fr_info.presptimes{m} = all_es{m}.presptimes;
                fr_info.preylen{m} = length(all_es{m}.y)/all_es{m}.fs;
                fr_info.repnum{m} = size(all_es{m}.presptimes,2);
                temp = all_es{m}.presptimes.';
                concat_presptimes = [concat_presptimes;vertcat(vertcat(temp{:}))+ sum_prelen];
                sum_prelen = sum_prelen +  fr_info.preylen{m};
                
                % for plty
                fr_info.pltsptimes{m} = all_es{m}.pltsptimes
                fr_info.pltlen{m} = length(all_es{m}.plty)/all_es{m}.fs;
                temp = all_es{m}.pltsptimes.';
                concat_pltsptimes = [concat_pltsptimes;vertcat(vertcat(temp{:}))+ sum_pltlen];
                sum_pltlen = sum_pltlen +  fr_info.pltlen{m};
                
                % for y only
                
                fr_info.ysptimes{m} = all_es{m}.sptimes
                fr_info.ylen{m} = length(all_es{m}.y)/all_es{m}.fs;
                temp = all_es{m}.sptimes.';
                concat_ysptimes = [concat_ysptimes;vertcat(vertcat(temp{:}))+ sum_ylen];
                sum_ylen = sum_ylen +  fr_info.ylen{m};
                
            end
            
            % for norm songs, degressive songs, adn detailed
            % frags/replas, calculate the concatenated firing rate one
            % by one
            % for pre_y
            fr_info.concat_pre_sptimes = concat_presptimes;
            fr_info.concat_pre_len = sum_prelen;
            fr_info.mean_pre_fr = length(concat_presptimes)/sum_prelen;
            
            % for plt_y
            fr_info.concat_plt_sptimes = concat_pltsptimes;
            fr_info.concat_plt_len = sum_pltlen;
            fr_info.mean_plt_fr = length(concat_pltsptimes)/sum_pltlen;
            
            % for y only
            fr_info.concat_y_sptimes = concat_ysptimes;
            fr_info.concat_y_len = sum_ylen;
            fr_info.mean_y_fr = length(concat_ysptimes)/sum_ylen;
            
        end
        
        function [dists1,dists0,featurename] = calCumulativeFeatureDiff(a,featurename)
            %fraglist = a.judgeFragResponse;
            a.judgeFragResp_FR;
            fragids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'Frag|frag|syl|ele')));
            fraglist = a.list(fragids);
            labeled_frags = fraglist([fraglist.label].' == 1);
            fe1 = []; for k = 1:length(labeled_frags); eval(sprintf('fe1(k) = labeled_frags(k).meanfeatures.%s',featurename)); end
            if isempty(fe1)
                dists1 = [];
                dists0 = [];
                featurename = featurename;
                return
            end
            pairs1 = (nchoosek(fe1,2));
            dists1 = abs(pairs1(:,1)-pairs1(:,2));
            
            
            unlabeled_frags = fraglist([fraglist.label].' == 0);
            fe0 = [];  error = struct;
            for k = 1:length(unlabeled_frags)
                try
                    eval(sprintf('fe0(k) = unlabeled_frags(k).meanfeatures.%s',featurename));
                catch ME
                    error(k).ME = ME;
                    fe0(k) = nan;
                    disp('Feature Info Missing!!')
                end
            end
            fe0 = rmmissing(fe0);
            pairs0 = (nchoosek(fe0,2));
            dists0 = abs(pairs0(:,1)-pairs0(:,2));
            % fig = figure('Position',a.fig_size1,'Color','w');  hold on;
            
            %             if isempty(dists1)
            %                 frame = getframe(gcf);
            %                 img = frame.cdata;
            %                 close(gcf);
            %                 return
            %             end
            %             cdfplot(dists1);
            %             cdfplot(dists0);
            %             legend('Response-eliciting elements','Not eliciting elements','Location','best')
            %             xlabel(sprintf('Difference between Mean %s (Hz)',featurename),'interpreter', 'none');
            %             ylabel('Cumulative %');
            %             hold off
            %             %             [f0,x0]= ecdf(dists0)
            %             %             [f1,x1]= ecdf(dists1)
            %             [h,p] = kstest2(dists1,dists0);
            %             %set(fig,'defaultTextInterpreter','none')
            %             title(sprintf('%s P-value : %.8f',a.formated_imagename,p),'interpreter', 'none');
            %             %saveas(gcf,sprintf('CDF-%s.png',a.formated_imagename));
            %             frame = getframe(gcf);
            %             img = frame.cdata;
            %             close(gcf);
            
        end
        
        function deg_songnames = which_songs_are_targeted_for_further_modification(a)
            
            degids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'deg|Deg')));
            
            deglist = a.list(degids);
            
            deg_songnames = unique(cellstr(regexp([deglist.stimuliname].','[OGBYR]\d{3}','match')));
            
            
        end
    end
    
    methods % 作图方法
        
        % 只作图不保存
        
        function drawAllFigures(a)
            
            img1 = A.Three;
            img2 = A.saveDrawAlignedConsDegs;
            img3 = A.saveDrawAlignedConsFrag;
            img4 = A.saveDrawAlignedConsReplas;
            img5 = A.saveDrawSortedRespToFrags;
            img234 = padFigures({img2,img3,img4});
            allimg = padFigures({img1,img234,img5});
            imwrite(allimg,sprintf('%s.png',a.formated_imagename));
            
            % 首先生成一张大图
            % 之后可以生成一个报告
            %生成的A文件需要有自我更新的能力
            %             A.drawAllWaveform;
            %             A.drawFragScatter
            % probably draw it into one big figure
            
        end
        function img = Three(a,subfileid)
            % draw three plot, using plt-data
            es = getAllEphysObject(a);
            for idx = 1: length(es)
                es{idx}.pltthree;
                frame = getframe(gcf);
                I{idx} = frame.cdata;
                close(gcf);
            end
            
            a.drawFirstWaveform;     % draw waveform
            frame = getframe(gcf);
            I{length(I)+ 1} = frame.cdata;
            close(gcf);
            
            % draw blank white
            lieshu = 12;
            hangshu = ceil(length(I)/lieshu);
            rest = lieshu*hangshu - length(I);
            white = uint8(255*ones(size(I{1})));
            
            if rest > 0
                for k = 1:rest
                    I = [I,white];
                end
            end
            
            reshapedI = reshape(I, lieshu,[])';
            clear I
            img = cell2mat(reshapedI);
            imwrite(img,sprintf('Three_%s.png',a.formated_imagename));
            
        end
        
        
        function drawFirstWaveform(a)
            
            % temporialriy a.neurons{1}
            waveforms = a.neurons{1}.waveform;
            figure('Color','w','Position',PM.size3plots);
            hold on
            plot(waveforms.',':','Color',[.5,.5,.5]);
            plot(max(waveforms),'--','Color','blue');
            plot(min(waveforms),'--','Color','blue');
            plot(mean(waveforms),'Color','red');
            
        end
        
        function drawAllWaveform(a)
            % this function is used to draw waveform
            % what this works for concatenating all neurons together???
            for k = 1: length(a.neurons)
                waveforms{k} = a.neurons{k}.waveform;
            end
            concat_waveforms = vertcat(waveforms{:});
            figure('Color','w');
            hold on
            plot(concat_waveforms.',':','Color',[.5,.5,.5]);
            plot(max(concat_waveforms),'--','Color','blue');
            plot(min(concat_waveforms),'--','Color','blue');
            plot(mean(concat_waveforms),'Color','red');
            
        end
        
        function drawSeparatedWaveform(a)
            
            fig = figure('Color','w');
            cmap = colormap(flip(hsv(5)));
            for k = 1:length(a.neurons)
                local_waveform = a.neurons{k}.waveform;
                hold on
                plot(local_waveform.',':','Color',[cmap(k,:),0.3]);
                plot(max(local_waveform),'--','Color',cmap(k,:)*0.7);
                plot(min(local_waveform),'--','Color',cmap(k,:)*0.7);
                plot(mean(local_waveform),'-*','Color',cmap(k,:)*0.5);
                
            end
            colorbar('southoutside')
            
            title(sprintf('%u plexon files',length(a.neurons)));
            saveas(gcf,sprintf('WaveformsSeparated-%s.png', a.formated_imagename));
        end
        
        function drawFragScatter(a,not_tested_handle)
            % not_tested_handle = 1 means draw
            
            if ~exist('not_tested_handle','var') % to judge whether to draw not tested elements or not
                not_tested_handle = 0;
            end
            
            global_eleinf = a.all_eleinf;
            fraglist = a.judgeFragResponse;
            
            
            for k = 1: length(global_eleinf)
                uniqueid = sprintf('-%s-%u-',global_eleinf(k).songname, global_eleinf(k).fragid);
                
                if ~isempty (find(~cellfun(@isempty, regexp({fraglist.stimuliname}.',uniqueid)), 1)) % if this element was tested
                    id_in_fraglist = find(~cellfun(@isempty, regexp({fraglist.stimuliname}.',uniqueid)));
                    if fraglist(id_in_fraglist).label == 0 % if the tested ele does not trigger reponse
                        
                        global_eleinf(k).scatter = 0;
                    elseif fraglist(id_in_fraglist).label == 1 % if the tested ele trigger response
                        global_eleinf(k).scatter = 1;
                    end
                else
                    global_eleinf(k).scatter = -1;
                end
            end
            
            % section for drawing
            figure
            hold on
            for k = 1: length(global_eleinf)
                if global_eleinf(k).scatter == -1
                    if not_tested_handle == 1
                        scatter(global_eleinf(k).coor_1,global_eleinf(k).coor_2,[],'k','filled'); % black for not tested
                    elseif not_tested_handle == 0 % if not draw test handle, do nothing
                    end
                elseif global_eleinf(k).scatter == 0
                    scatter(global_eleinf(k).coor_1,global_eleinf(k).coor_2,[],'g','filled');   % green for tested but not response-eliciting
                    % text(double(global_eleinf(k).coor_1),double(global_eleinf(k).coor_2),sprintf('%s-%u',global_eleinf(k).songname, global_eleinf(k).fragid) )
                elseif global_eleinf(k).scatter == 1
                    scatter(global_eleinf(k).coor_1,global_eleinf(k).coor_2,[],'r','filled');
                    % text(double(global_eleinf(k).coor_1),double(global_eleinf(k).coor_2),sprintf('%s-%u',global_eleinf(k).songname, global_eleinf(k).fragid) )
                end
                
                
            end
            
            
            % label the targets
            global_names = [global_eleinf.songname].';
            global_fragids = [global_eleinf.fragid].';
            
            global_merged = {};
            for w = 1: length(global_names)
                global_merged{w} = sprintf('%s-%u',global_names(w),global_fragids(w));
            end
            global_merged = global_merged.';
            
            
            for u = 1: length(a.targets)
                [~,beta] = ismember( a.targets{u}, global_merged );
                scatter( global_eleinf(beta).coor_1,global_eleinf(beta).coor_2,[],'k','h');
            end
            hold off
            
            
        end
        
        function drawPCABasedOnAllFeatures(a)
            
            fraglist = a.judgeFragResponse;
            if isempty(fraglist)
                return;
            end
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).maxvalue;
            end
            
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).amplitude = fraglist(u).meanfeatures.amplitude;
                toshow(u).pitch = fraglist(u).meanfeatures.pitch;
                toshow(u).AM = fraglist(u).meanfeatures.AM;
                toshow(u).mean_frequency = fraglist(u).meanfeatures.mean_frequency;
                toshow(u).FM = fraglist(u).meanfeatures.FM;
                toshow(u).peak_frequency = fraglist(u).meanfeatures.peak_frequency;
                toshow(u).entropy = fraglist(u).meanfeatures.entropy;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            
            toshow  = table2struct(sortrows(struct2table(toshow), 'resp','ascend'));
            
            
            num_toshow = [toshow.amplitude; toshow.pitch; toshow.AM; toshow.mean_frequency; toshow.FM; ...
                toshow.peak_frequency; toshow.entropy].';
            
            [coe,score,~] = pca(num_toshow);
            
            cmap = 1- rescale( [toshow.resp].');
            % 1- to plot the data
            figure('Color','white');
            hold on
            scatter(score(:,1),score(:,2),[],repmat(cmap,1,3),'filled');
            %colorbar;
            xlabel('PCA Dim-1');
            ylabel('PCA Dim-2');
            %title('PCA data');
            
        end
        
        function saveDrawSSIMSimlarityMatrix(a)
            % This function firstly order the frags by response strength
            % then measure the pairwise similarity between elements
            fraglist = a.judgeFragResponse;
            sorted_fraglist = table2struct(sortrows( struct2table(fraglist) ,'halfsum','descend'));  %此处用 maxvalue或许不太对
            
            for k = 1: length(sorted_fraglist)
                fiy = bandpass(sorted_fraglist(k).y,[900 6000],sorted_fraglist(k).fs); % preiously 5000
                
                envy = rescale(smooth(abs(fiy),150)); % amplitude envelope of y
                %powery = downsample(fiy.^2/length(fiy),fs/1000);
                %downy = downsample(abs(fiy),fs/1000);
                I = cal.spec(fiy,sorted_fraglist(k).fs); % I is the image of the whole song
                sorted_fraglist(k).normalized_img = imresize(I,[257,50]);
            end
            
            for m = 1:length(sorted_fraglist)
                parfor p = 1:length(sorted_fraglist)
                    sim(m,p) = ssim(sorted_fraglist(m).normalized_img,sorted_fraglist(p).normalized_img);
                end
            end
            
            figure;
            imagesc(sim);
            
            saveas(gcf,sprintf('SSIMSimlarityMatrix-%s.png',a.formated_imagename));
            close(gcf);
        end
        
        
        function saveDrawMeanFeatureVsResp(a)
            % draw the distribution of mean features
            
            %             fragids = find(~cellfun(@isempty, regexp({a.list.stimuliname}.','frag|Frag|syl|Syl')));
            %
            %             if isempty(fragids)
            %                 return
            %             end
            %             fraglist = a.list(fragids);
            
            fraglist = a.judgeFragResponse;
            if isempty(fraglist)
                return;
            end
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).rs; % Response Strength
            end
            fraglist = table2struct(sortrows(struct2table(fraglist), 'responseMeasure'));
            
            figure;
            subplot(5,2,1)
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).value = fraglist(u).meanfeatures.amplitude;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            toshow  = table2struct(sortrows(struct2table(toshow), 'value'));
            plot([toshow.value].',[toshow.resp].');
            xlabel('Amplitude')
            ylabel('Maximum FR');
            
            subplot(5,2,2)
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).value = fraglist(u).meanfeatures.pitch;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            toshow  = table2struct(sortrows(struct2table(toshow), 'value'));
            plot([toshow.value].',[toshow.resp].');
            xlabel('Pitch')
            ylabel('Maximum FR');
            
            subplot(5,2,3)
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).value = fraglist(u).meanfeatures.mean_frequency;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            toshow  = table2struct(sortrows(struct2table(toshow), 'value'));
            plot([toshow.value].',[toshow.resp].');
            xlabel('Mean Frequency')
            ylabel('Maximum FR');
            
            subplot(5,2,4)
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).value = fraglist(u).meanfeatures.FM;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            toshow  = table2struct(sortrows(struct2table(toshow), 'value'));
            plot([toshow.value].',[toshow.resp].');
            xlabel('FM')
            ylabel('Maximum FR');
            
            subplot(5,2,5)
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).value = fraglist(u).meanfeatures.AM;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            toshow  = table2struct(sortrows(struct2table(toshow), 'value'));
            plot([toshow.value].',[toshow.resp].');
            xlabel('AM')
            ylabel('Maximum FR');
            
            subplot(5,2,6)
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).value = fraglist(u).meanfeatures.goodness;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            toshow  = table2struct(sortrows(struct2table(toshow), 'value'));
            plot([toshow.value].',[toshow.resp].');
            xlabel('Goodness')
            ylabel('Maximum FR');
            
            subplot(5,2,7)
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).value = fraglist(u).meanfeatures.entropy;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            toshow  = table2struct(sortrows(struct2table(toshow), 'value'));
            plot([toshow.value].',[toshow.resp].');
            xlabel('Entropy')
            ylabel('Maximum FR');
            
            subplot(5,2,8)
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).value = fraglist(u).meanfeatures.peak_frequency;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            toshow  = table2struct(sortrows(struct2table(toshow), 'value'));
            plot([toshow.value].',[toshow.resp].');
            xlabel('Peak Frequency')
            ylabel('Maximum FR');
            
            subplot(5,2,9)
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).value = fraglist(u).meanfeatures.continuity_t;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            toshow  = table2struct(sortrows(struct2table(toshow), 'value'));
            plot([toshow.value].',[toshow.resp].');
            xlabel('Time continuity')
            ylabel('Maximum FR');
            
            subplot(5,2,10)
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).value = fraglist(u).meanfeatures.continuity_f;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            toshow  = table2struct(sortrows(struct2table(toshow), 'value'));
            plot([toshow.value].',[toshow.resp].');
            xlabel('frequency continuity')
            ylabel('Maximum FR');
            
            
            saveas(gcf, sprintf('SeparatedFeatureVsResponse_%s.png',a.formated_imagename));
            close(gcf);
            
            
        end
        
        function saveDrawMeanFeaturesVsRespAsLineChart(a)
            % draw the distribution of mean features
            a.calHarmRatio;
            fraglist = a.judgeFragResponse;
            if isempty(fraglist)
                return;
            end
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).rs;
                % fraglist(k).responseMeasure = fraglist(k).maxvalue;
            end
            
            toshow = struct; error = struct;
            for u = 1: length(fraglist)
                try
                    toshow(u).amplitude = fraglist(u).meanfeatures.amplitude;
                    toshow(u).pitch = fraglist(u).meanfeatures.pitch;
                    toshow(u).AM = fraglist(u).meanfeatures.AM;
                    toshow(u).mean_frequency = fraglist(u).meanfeatures.mean_frequency;
                    toshow(u).FM = fraglist(u).meanfeatures.FM;
                    toshow(u).peak_frequency = fraglist(u).meanfeatures.peak_frequency;
                    toshow(u).entropy = fraglist(u).meanfeatures.entropy;
                    toshow(u).harmratio = fraglist(u).meanfeatures.harmratio;
                    toshow(u).resp = fraglist(u).responseMeasure;
                catch ME
                    %toshow(u) = [];
                    error(u).ME = ME;
                end
            end
            
            toshow = table2struct(rmmissing(struct2table(toshow)));
            toshow = toshow(find(~cellfun(@isempty,{toshow.amplitude}.')));
            toshow  = table2struct(sortrows(struct2table(toshow), 'resp','ascend'));
            
            
            num_toshow = [toshow.amplitude; toshow.pitch; toshow.AM; toshow.mean_frequency; toshow.FM; ...
                toshow.peak_frequency; toshow.entropy; toshow.harmratio].';
            
            znum_toshow = zscore(num_toshow,0,1);
            %
            %            figure
            %            plot(znum_toshow.');
            
            figure('Position',[1997 233 1388 658],'Color','w');
            hold on
            c = 1- rescale([toshow.resp].',0.1,1);
            for r = 1: size(znum_toshow,1)
                plot(znum_toshow(r,:),'Color',repmat(c(r),3,1));
                %drawnow
                % pause(0.5)
            end
            colormap(flip(repmat(unique(c),1,3)))
            colorbar
            xlim([0,9]);
            xticks([0 1 2 3 4 5 6 7 8,9])
            xticklabels({'','amplitude','pitch','AM','mean_frequency','FM','peak-frequency','entropy','Harmonic ratio',''});
            title(sprintf('Totally %u song elements',length(fraglist)));
            ylabel('Zscored Feature(averaged)');
            
            saveas(gcf,sprintf('LineChartMeanFeaturesVsResp-%s.png',a.formated_imagename));
            close(gcf);
        end
        
        function drawSelectedStimuliResp(a,stimulinames,range)
            % range是从plty的起始为始，以秒计
            dbstop if error; tic
            if isa(stimulinames,'cell')
                selectedids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),strjoin(stimulinames,'|') )));
            else
                selectedids = stimulinames;
            end
            
            
            songids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'norm|spe') ));
            songids = intersect(songids,selectedids);
            songlist = a.list(songids);
            %             [~,postunique] = unique(cellstr(cellfun(@convert.bid,{songlist.stimuliname}.','Uni',0)))
            %             songlist = songlist(postunique);
            RONGYU = 0.5;
            range = range + RONGYU;
            for k = 1: length(songlist)
                songlist(k).pady = [zeros(RONGYU*songlist(k).fs,1);songlist(k).plty];
                songlist(k).padsptimes = cellfun( @(x) x + RONGYU, songlist(k).pltsptimes,'uni',0);
            end
            
            % About Frag
            fragids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'frag|Frag|syl|Syl') ));
            fragids = intersect(fragids,selectedids);
            if ~isempty(fragids)
                fraglist = a.list(fragids);
                for m = 1: length(fraglist)
                    birdid = convert.bid(fraglist(m).stimuliname);
                    ids_norm = find(~cellfun(@isempty, regexp(cellstr({songlist.stimuliname}.'),birdid) ) );
                    if ~isempty(ids_norm)& length(ids_norm) == 1
                        fraglist(m).sylIni = Analysis.findIni(songlist(ids_norm).pady,fraglist(m).y);
                        fraglist(m).pady = [zeros([fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext,1]);fraglist(m).plty;zeros(length(songlist(ids_norm).plty)...
                            - (fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext) - length(fraglist(m).plty),1)];
                        fraglist(m).padsptimes = cellfun( @(x) x + (fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext)/fraglist(m).fs, fraglist(m).pltsptimes,'uni',0);
                    end
                end
            end
            
            % About Deg
            degids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'Deg|deg') ));
            degids = intersect(degids,selectedids);
            if ~isempty(degids)
                deglist = a.list(degids);
                for m = 1: length(deglist)
                    birdid = convert.bid(deglist(m).stimuliname);
                    ids_norm = find(~cellfun(@isempty, regexp(cellstr({songlist.stimuliname}.'),birdid) ) );
                    if ~isempty(ids_norm)& length(ids_norm) == 1
                        [deglist(m).sylIni,trump_diffvalue] = Analysis.findIni(songlist(ids_norm).pady,deglist(m).y);
                        fprintf('%s has the DIFFVALUE as: %u and INI as: %f \n ',deglist(m).stimuliname,trump_diffvalue,deglist(m).sylIni/32000);
                        deglist(m).pady = [zeros([deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext,1]);deglist(m).plty;zeros(length(songlist(ids_norm).plty)...
                            - (deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext) - length(deglist(m).plty),1)];
                        deglist(m).padsptimes = cellfun( @(x) x + (deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext)/deglist(m).fs, deglist(m).pltsptimes,'uni',0);
                    end
                end
            end
            
            % About Repla
            replaids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'Repla|repla|catego|Catego') ));
            replaids = intersect(replaids,selectedids);
            if ~isempty(replaids)
                replalist = a.list(replaids);
                ids_norm_collecting_box = [];
                for m = 1: length(replalist)
                    
                    afterBefore = regexp(replalist(m).stimuliname,'(?<=before-)\S*','match');
                    afterBefore = afterBefore{1};
                    birdid = convert.bid(afterBefore);
                    ids_norm = find(~cellfun(@isempty, regexp(cellstr({songlist.stimuliname}.'),birdid) ) );
                    
                    if ~isempty(ids_norm)%& length(ids_norm) == 1
                        ids_norm_1st = ids_norm(1);
                        ids_norm_collecting_box = [ids_norm_collecting_box, ids_norm];
                        afterpad_length = length(songlist(ids_norm_1st).plty) + RONGYU*songlist(ids_norm_1st).fs;  % +0.5s
                        
                        replalist(m).pady = [zeros(afterpad_length- length(replalist(m).plty),1);replalist(m).plty];
                        replalist(m).padsptimes = cellfun( @(x) x + length(zeros(afterpad_length- length(replalist(m).plty),1))/fraglist(m).fs,replalist(m).pltsptimes,'uni',0);
                    end
                end
            end
            
            % merge the Cons,Degs,Frags,Replas lists together
            for w = 1: length(songlist)
                
                if ~isempty(degids)
                    birdid = convert.bid(songlist(w).stimuliname);
                    ids_indeg = find(~cellfun(@isempty, regexp(cellstr({deglist.stimuliname}.'),birdid) ) );
                    selected_deglist = deglist(ids_indeg);
                    [~,temp_index] = sortrows([selected_deglist.sylIni].');
                    selected_deglist = selected_deglist(temp_index);
                else
                    selected_deglist = [];
                end
                
                if ~isempty(fragids)
                    birdid = convert.bid(songlist(w).stimuliname);
                    ids_infrag = find(~cellfun(@isempty, regexp(cellstr({fraglist.stimuliname}.'),birdid) ) );
                    selected_fraglist = fraglist(ids_infrag);
                    [~,temp_index] = sortrows([selected_fraglist.sylIni].');
                    selected_fraglist = selected_fraglist(temp_index);
                else
                    selected_fraglist = [];
                end
                
                if ~isempty(replaids)
                    birdid = convert.bid(songlist(w).stimuliname);
                    ids_inrepla = find(~cellfun(@isempty, regexp(cellstr({replalist.stimuliname}.'),birdid) ) );
                    selected_replalist = replalist(ids_inrepla);
                else
                    selected_replalist = [];
                end
                
                % draw figures
                if ~isempty(selected_deglist)
                    selected_deglist = rmfield(selected_deglist,'sylIni');
                end
                
                if ~isempty(selected_fraglist)
                    selected_fraglist = rmfield(selected_fraglist,'sylIni');
                end
                alllist = horzcat(songlist(w),selected_deglist,selected_fraglist,selected_replalist);
                len = length(alllist);
                figure('Position',[1935 -207 520 1458/9*len],'Color','none');
                
                ax = tight_subplot(2*len, 1, 0.002, 0.02, 0);
                for k = 1: len
                    
                    axes(ax(2*(k-1)+ 1)); % draw.two(,,);
                    %ax(2*(k-1)+ 1).Position(4) =  ax(2*(k-1)+ 1).Position(4);
                    if exist('range','var')
                        truncated_y = alllist(k).pady(range(1)*alllist(k).fs:range(2)*alllist(k).fs);
                        draw.spec(truncated_y,alllist(k).fs);
                    else
                        draw.spec(alllist(k).pady,alllist(k).fs);
                    end
                    xlabel('')
                    ylabel('')
                    
                    set(gca,'TickLength',[0 .01])
                    set(gca,'Yticklabel',[])
                    set(gca,'Xticklabel',[])
                    
                    axes(ax(2*k));
                    
                    if exist('range','var')
                        truncated_sptimes = extract.sptimes_resetSP(alllist(k).padsptimes, range(1), range(2));
                        truncated_y = alllist(k).pady(range(1)*alllist(k).fs:range(2)*alllist(k).fs);
                        %draw.raster(truncated_sptimes,truncated_y,alllist(k).fs);
                        draw.rasterBeta(truncated_sptimes,truncated_y,alllist(k).fs,2.8,'k');
                    else
                        %draw.raster(alllist(k).padsptimes,alllist(k).pady,alllist(k).fs);
                        draw.rasterBeta(alllist(k).padsptimes,alllist(k).pady,alllist(k).fs,2.8,'k');
                    end
                    
                    ylabel('')
                    set(gca,'TickLength',[0 .01])
                    
                    
                    
                    %                         set(gca,'Yticklabel',[])
                    %                         set(gca,'Xticklabel',[])
                    %                     else
                    %                         xlabel(alllist(k).stimuliname);
                    %                     end
                    
                end
                
                for k = 1:len
                    ax(2*k-1).Position(4) =  ax(2*k-1).Position(4) -0.008;
                    ax(2*k-1).Position(2) =  ax(2*k-1).Position(2)+0.006;
                    box(ax(2*k-1),'off');
                    set(ax(2*k-1),'XColor','none','YColor','none')
                    %set(ax(2*k-1), 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
                end
                
                for k = 1:len
                    ax(2*k).Position(2) =  ax(2*k).Position(2)+0.000;
                    ax(2*k).Position(4) =  ax(2*k).Position(4)+0.002;
                    box(ax(2*k),'off');
                    set(ax(2*k),'XColor','none','YColor','none')
                    %set(ax(2*k), 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
                end
                disp('pause here')
                
                
                
                
            end
        end
        
        function drawStimuliOnly(a,name_or_id,range)
            % range是从plty的起始为始，以秒计
            dbstop if error; tic
            if isnumeric(name_or_id)
                selectedid = name_or_id;
            else
                selectedid = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),name_or_id )));
            end
            
            k = selectedid;
            if exist('range','var')
                truncated_y = a.list(k).plty(range(1)*a.list(k).fs:range(2)*a.list(k).fs);
                figure('Position',[-36 437 length(truncated_y)/18 528]);
                
                
            else
                figure('Position',[-36 437 length(a.list(k).plty)/18 528]);
                draw.spec(a.list(k).plty,a.list(k).fs);
            end
            xlabel('')
            ylabel('')
            
            set(gca,'TickLength',[0 .01])
            set(gca,'Yticklabel',[])
            set(gca,'Xticklabel',[])
            set(gca,'XColor','none','YColor','none')
            box off
            disp('pause')
            
            %                         set(gca,'Yticklabel',[])
            %                         set(gca,'Xticklabel',[])
            
            
        end
        
        % 作图并导出
        
        function img = exportConRevMir(a)
            dbstop if error
            reverseids = find(~cellfun(@isempty,regexp(cellstr({a.list.stimuliname}.'),'reverse')));
            normids = find(~cellfun(@isempty,regexp(cellstr({a.list.stimuliname}.'),'norm')));
            mirrorids = find(~cellfun(@isempty,regexp(cellstr({a.list.stimuliname}.'),'mirror')));
            
            % find out the songs that are tested with the reversed version
            bids = cellfun(@(x) convert.bid(x),cellstr({a.list(reverseids).stimuliname}.'),'Uni',0);
            
            IMG = {};
            
            
            if isempty(bids)
                figure('Position',[552 -116 2410 1221],'Color','w');;
                %draw.three(a.list(sameFile_mirror_ids).plty,a.list(sameFile_mirror_ids).fs,a.list(sameFile_mirror_ids).pltsptimes);
                frame = getframe(gcf);
                IMG{1,1} = frame.cdata;
                IMG{2,1} = frame.cdata;
                IMG{3,1} = frame.cdata;
                close(gcf);
            end
            
            
            
            for k = 1:length(bids)
                
                samebirdids = find(~cellfun(@isempty,regexp(cellstr({a.list.stimuliname}.'),bids{k})));
                sb_reverse_ids = intersect(reverseids,samebirdids); % sab means same bird
                if length(sb_reverse_ids) > 1
                    sb_reverse_ids =   sb_reverse_ids(1); % 如果多项，暂时的对策是取第一项
                end
                sb_mirror_ids = intersect(mirrorids,samebirdids);
                sb_norm_ids = intersect(normids,samebirdids);
                
                samefileids  = find(~cellfun(@isempty,regexp(cellstr({a.list.Fid}.'),a.list(sb_reverse_ids).Fid)));
                sameFile_norm_ids = intersect(samefileids,sb_norm_ids);
                sameFile_reverse_ids = intersect(samefileids,sb_reverse_ids);
                sameFile_mirror_ids = intersect(samefileids,sb_mirror_ids);
                
                fig1 = figure('Position',[552 -116 2410 1221],'Color','w');
                draw.three(a.list(sameFile_norm_ids).y,a.list(sameFile_norm_ids).fs,a.list(sameFile_norm_ids).sptimes); % plty
                frame = getframe(fig1);
                IMG{1,k} = frame.cdata;
                close(fig1)
                
                fig2 = figure('Position',[552 -116 2410 1221],'Color','w');
                draw.three(a.list(sameFile_reverse_ids).y,a.list(sameFile_reverse_ids).fs,a.list(sameFile_reverse_ids).sptimes);
                frame = getframe(fig2);
                IMG{2,k} = frame.cdata;
                close(fig2)
                
                fig3 = figure('Position',[552 -116 2410 1221],'Color','w');
                draw.three(a.list(sameFile_mirror_ids).y,a.list(sameFile_mirror_ids).fs,a.list(sameFile_mirror_ids).sptimes);
                frame = getframe(fig3);
                IMG{3,k} = frame.cdata;
                close(fig3);
                
            end
            
            img = cell2mat(IMG);
            
            img = convert.colorEdge(img,'r');
            % I don't know why here the output could be unit8 or double
            % randomly
            % Now the temporal solution is
            img = uint8(img);
            
        end
        function img = exportDrawPitchHarmoLineChart(a)
            % draw the distribution of mean features
            a.calHarmRatio;
            fraglist = a.judgeFragResponse;
            if isempty(fraglist)
                return;
            end
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).maxvalue;
            end
            
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).pitch = fraglist(u).meanfeatures.pitch;
                toshow(u).harmratio = fraglist(u).meanfeatures.harmratio;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            
            toshow  = table2struct(sortrows(struct2table(toshow), 'resp','ascend'));
            
            
            num_toshow = [ toshow.pitch; toshow.harmratio].';
            
            znum_toshow = zscore(num_toshow,0,1);
            %
            %            figure
            %            plot(znum_toshow.');
            
            figure('Position',a.fig_size1,'Color','w');
            hold on
            c = 1- rescale([toshow.resp].',0.1,1);
            for r = 1: size(znum_toshow,1)
                plot(znum_toshow(r,:),'Color',repmat(c(r),3,1));
                %drawnow
                % pause(0.5)
            end
            colormap(flip(repmat(unique(c),1,3)))
            colorbar
            xlim([0,3]);
            xticks([0 1 2 3])
            xticklabels({'','Pitch','Harmonic ratio',''});
            title(sprintf('%s---%u song elements',a.formated_imagename,length(fraglist)),'interpreter','none');
            ylabel('Zscored Feature(averaged)');
            
            saveas(gcf,sprintf('PitchHarmoLineChart-%s.png',a.formated_imagename));
            
            frame = getframe(gcf);
            img = frame.cdata;
            close(gcf);
        end
        
        function img = exportDrawInsong_drawPitchHarmoLineChart(a)
            % draw the distribution of mean features
            fraglist = a.getInsongFragRespList;
            if isempty(fraglist)
                return;
            end
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).label;
            end
            
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).pitch = fraglist(u).meanfeatures.pitch;
                toshow(u).harmratio = fraglist(u).meanfeatures.harmratio;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            
            toshow  = table2struct(sortrows(struct2table(toshow), 'resp','ascend'));
            
            [~,index] = sortrows([toshow.resp].'); toshow = toshow(index); clear index
            num_toshow = [ toshow.pitch; toshow.harmratio].';
            
            znum_toshow = zscore(num_toshow,0,1);
            %
            %            figure
            %            plot(znum_toshow.');
            
            figure('Position',a.fig_size1,'Color','w');
            hold on
            c = 1- rescale([toshow.resp].',0.1,1);
            for r = 1: size(znum_toshow,1)
                plot(znum_toshow(r,:),'Color',repmat(c(r),3,1));
                %drawnow
                % pause(0.5)
            end
            colormap(flip(repmat(unique(c),1,3)))
            colorbar
            xlim([0,3]);
            xticks([0 1 2 3])
            xticklabels({'','Pitch','Harmonic ratio',''});
            title(sprintf('%s---%u song elements',a.formated_imagename,length(fraglist)),'interpreter','none');
            ylabel('Zscored Feature(averaged)');
            
            saveas(gcf,sprintf('Insong-PitchHarmoLineChart-%s.png',a.formated_imagename));
            
            frame = getframe(gcf);
            img = frame.cdata;
            close(gcf);
        end
        
        function img = exportDrawdraw2DPitchVsHarmo(a)
            
            fraglist = a.judgeFragResponse;
            if isempty(fraglist)||isempty(fraglist(1).meanfeatures)
                return
            end
            
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).rs;
            end
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).pitch = fraglist(u).meanfeatures.pitch;
                toshow(u).harmratio = fraglist(u).meanfeatures.harmratio;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            [~,index] = sortrows([toshow.resp].'); toshow = toshow(index); clear index
            
            name = {'pitch','harmratio'};
            ncomb = nchoosek(name,2);
            
            I = {};
            for idx = 1:size(ncomb,1) % 此处可用parfor
                
                %subplot(hangshu,lieshu,idx);
                figure('Position',a.fig_size1,'Color','w');
                hold on
                cmap = 1 - rescale([toshow.resp].');
                for  md = 1: length(toshow)
                    scatter(toshow(md).(ncomb{idx,1})', toshow(md).(ncomb{idx,2}),[],repmat(cmap(md,:),1,3),'filled');
                end
                xlabel(replace(ncomb{idx,1},'_','-'));
                ylabel(replace(ncomb{idx,2},'_','-'));
                title(a.formated_imagename,'interpreter', 'none')
                frame = getframe(gcf);
                I{idx} = frame.cdata;
                close(gcf);
                
            end
            
            reshapedI = reshape(I, 1,[])';
            clear I
            img = cell2mat(reshapedI);
            imwrite(img,sprintf('2DPitchHarmo_%s.png',a.formated_imagename));
            
        end
        
        function img = exportDrawInsong_draw2DPitchVsHarmo(a)
            
            fraglist = a.judgeFragResponse;
            if isempty(fraglist)||isempty(fraglist(1).meanfeatures)
                return
            end
            
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).label; % 二分法
            end
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).pitch = fraglist(u).meanfeatures.pitch;
                toshow(u).harmratio = fraglist(u).meanfeatures.harmratio;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            
            [~,index] = sortrows([toshow.resp].'); toshow = toshow(index); clear index
            name = {'pitch','harmratio'};
            ncomb = nchoosek(name,2);
            
            I = {};
            for idx = 1:size(ncomb,1) % 此处可用parfor
                
                %subplot(hangshu,lieshu,idx);
                figure('Position',a.fig_size1,'Color','w');
                hold on
                cmap = 0.9 - rescale([toshow.resp].')*0.9;
                for  md = 1: length(toshow)
                    scatter(toshow(md).(ncomb{idx,1})', toshow(md).(ncomb{idx,2}),[],repmat(cmap(md,:),1,3),'filled');
                end
                xlabel(replace(ncomb{idx,1},'_','-'));
                ylabel(replace(ncomb{idx,2},'_','-'));
                title(a.formated_imagename,'interpreter', 'none')
                frame = getframe(gcf);
                I{idx} = frame.cdata;
                close(gcf);
                
            end
            
            reshapedI = reshape(I, 1,[])';
            clear I
            img = cell2mat(reshapedI);
            imwrite(img,sprintf('Insong_2DPitchHarmo_%s.png',a.formated_imagename));
            
        end
        
        function img = exportDrawCumulativePitchDiff(a)
            fraglist = a.judgeFragResponse;
            labeled_frags = fraglist([fraglist.label].' == 1);
            pitch1 = []; for k = 1:length(labeled_frags); pitch1(k) = labeled_frags(k).meanfeatures.pitch; end
            pairs1 = (nchoosek(pitch1,2));
            dists1 = abs(pairs1(:,1)-pairs1(:,2));
            
            
            unlabeled_frags = fraglist([fraglist.label].' == 0);
            pitch0 = []; for k = 1:length(unlabeled_frags); pitch0(k) = unlabeled_frags(k).meanfeatures.pitch; end
            pairs0 = (nchoosek(pitch0,2));
            dists0 = abs(pairs0(:,1)-pairs0(:,2));
            fig = figure('Position',a.fig_size1,'Color','w');  hold on;
            
            if isempty(dists1)
                frame = getframe(gcf);
                img = frame.cdata;
                close(gcf);
                return
            end
            cdfplot(dists1);
            cdfplot(dists0);
            legend('Response-eliciting elements','Not eliciting elements','Location','best')
            xlabel('Difference between Mean Pitch (Hz)');
            ylabel('Cumulative %');
            hold off
            %             [f0,x0]= ecdf(dists0)
            %             [f1,x1]= ecdf(dists1)
            [h,p] = kstest2(dists1,dists0);
            %set(fig,'defaultTextInterpreter','none')
            title(sprintf('%s P-value : %.8f',a.formated_imagename,p),'interpreter', 'none');
            saveas(gcf,sprintf('CDF-%s.png',a.formated_imagename));
            frame = getframe(gcf);
            img = frame.cdata;
            close(gcf);
            
        end
        
        function img = exportdrawCumulativeFeatureDiff(a,featurename)
            fraglist = a.judgeFragResponse;
            labeled_frags = fraglist([fraglist.label].' == 1);
            fe1 = []; for k = 1:length(labeled_frags); eval(sprintf('fe1(k) = labeled_frags(k).meanfeatures.%s',featurename)); end
            pairs1 = (nchoosek(fe1,2));
            dists1 = abs(pairs1(:,1)-pairs1(:,2));
            
            
            unlabeled_frags = fraglist([fraglist.label].' == 0);
            fe0 = []; for k = 1:length(unlabeled_frags); eval(sprintf('fe0(k) = unlabeled_frags(k).meanfeatures.%s',featurename)); end
            pairs0 = (nchoosek(fe0,2));
            dists0 = abs(pairs0(:,1)-pairs0(:,2));
            fig = figure('Position',a.fig_size1,'Color','w');  hold on;
            
            if isempty(dists1)
                frame = getframe(gcf);
                img = frame.cdata;
                close(gcf);
                return
            end
            cdfplot(dists1);
            cdfplot(dists0);
            legend('Response-eliciting elements','Not eliciting elements','Location','best')
            xlabel(sprintf('Difference between Mean %s (Hz)',featurename),'interpreter', 'none');
            ylabel('Cumulative %');
            hold off
            %             [f0,x0]= ecdf(dists0)
            %             [f1,x1]= ecdf(dists1)
            [h,p] = kstest2(dists1,dists0);
            %set(fig,'defaultTextInterpreter','none')
            title(sprintf('%s P-value : %.8f',a.formated_imagename,p),'interpreter', 'none');
            %saveas(gcf,sprintf('CDF-%s.png',a.formated_imagename));
            frame = getframe(gcf);
            img = frame.cdata;
            close(gcf);
            
        end
        
        function img = exportdrawInsongCumulativePitchDiff(a)
            fraglist = a.getInsongFragRespList;
            labeled_frags = fraglist([fraglist.label].' == 1);
            pitch1 = []; for k = 1:length(labeled_frags); pitch1(k) = labeled_frags(k).meanfeatures.pitch; end
            if isempty(pitch1)||length(pitch1) == 1
                figure('Position',[1999 84 872 891],'Color','w');  hold on;
                title('malfunction')
                frame = getframe(gcf);
                img = frame.cdata;
                close(gcf);
                return
            end
            pairs1 = (nchoosek(pitch1,2));
            dists1 = abs(pairs1(:,1)-pairs1(:,2));
            
            
            unlabeled_frags = fraglist([fraglist.label].' == 0);
            pitch0 = []; for k = 1:length(unlabeled_frags); pitch0(k) = unlabeled_frags(k).meanfeatures.pitch; end
            pairs0 = (nchoosek(pitch0,2));
            dists0 = abs(pairs0(:,1)-pairs0(:,2));
            fig = figure('Position',a.fig_size1,'Color','w');  hold on;
            
            if isempty(dists1)
                frame = getframe(gcf);
                img = frame.cdata;
                close(gcf);
                return
            end
            cdfplot(dists1);
            cdfplot(dists0);
            legend('Response-eliciting elements','Not eliciting elements','Location','best')
            xlabel('Pairwise Mean Pitch Difference (Hz)');
            ylabel('Cumulative %');
            hold off
            %             [f0,x0]= ecdf(dists0)
            %             [f1,x1]= ecdf(dists1)
            [h,p] = kstest2(dists1,dists0);
            %set(fig,'defaultTextInterpreter','none')
            title(sprintf('%s P-value : %.8f',a.formated_imagename,p),'interpreter', 'none');
            saveas(gcf,sprintf('Insong-CDF-%s.png',a.formated_imagename));
            frame = getframe(gcf);
            img = frame.cdata;
            close(gcf);
            
        end
        
        function img = saveDrawSortedRespToFrags(a) %非常难以找到这个function
            %非常难以找到这个function
            fraglist =  judgeFragResponse(a);
            if isempty(fraglist)
                return
            end
            %sorted_fraglist = table2struct(sortrows( struct2table(fraglist) ,'maxvalue','descend'));
            sorted_fraglist = table2struct(sortrows( struct2table(fraglist) ,'rs','descend'));
            
            
            I = {}; % collection of frag-response-three images
            for k = 1: length(sorted_fraglist)
                
                if sorted_fraglist(k).label == 0
                    h =  figure('Position',[681 403 523 696],'Color','w');
                elseif sorted_fraglist(k).label == 1
                    h =  figure('Position',[681 403 523 696],'Color','y');
                end
                
                %h.WindowState = 'maximized';
                draw.two(sorted_fraglist(k).plty,sorted_fraglist(k).fs,sorted_fraglist(k).pltsptimes);
                xlabel(sprintf('%s-RS: %f',sorted_fraglist(k).stimuliname,sorted_fraglist(k).rs));
                temp = getframe(gcf);
                I{k} = temp.cdata;
                
                
                close(h)
            end
            
            %             figure;
            %             n.draw_waveform;     % draw waveform
            %             frame = getframe(gcf);
            %             I{length(I)+ 1} = frame.cdata;
            %             close(gcf);
            
            % draw blank white
            lieshu = 10;
            hangshu = ceil(length(I)/lieshu);
            rest = lieshu*hangshu - length(I);
            white = uint8(255*ones(size(I{1})));
            
            if rest > 0
                for k = 1:rest
                    I = [I,white];
                    %                     ax = gcf;
                    %                     ax.Position(3) = 560;
                    %                     ax.Position(4) = 420;
                end
            end
            
            reshapedI = reshape(I, lieshu,[])';
            clear I
            img = cell2mat(reshapedI);
            imwrite(img,sprintf('晋RespToFrags_%s.png',a.formated_imagename));
            
        end
        
        % 作图且保存
        function IMG = saveDrawlineChart_2D_Cumulative(a) % 秦
            
            img{1} = a.drawPitchHarmoLineChart;
            img{2} = a.draw2DPitchVsHarmo;
            img{3} = a.drawCumulativeFeatureDiff('pitch');
            img{4} = a.drawCumulativeFeatureDiff('mean_frequency');
            img{5} = a.drawCumulativeFeatureDiff('peak_frequency');
            img{6} = a.drawCumulativeFeatureDiff('FM');
            img{7} = a.drawCumulativeFeatureDiff('goodness');
            img{8} = a.drawCumulativeFeatureDiff('entropy');
            img{9} = a.drawCumulativeFeatureDiff('harmratio');
            img{10} = a.drawCumulativeFeatureDiff('amplitude');
            img{11} = a.drawCumulativeFeatureDiff('AM');
            
            % draw blank white
            hangshu = 3;
            lieshu = ceil(length(img)/hangshu);
            rest = lieshu*hangshu - length(img);
            white = uint8(255*ones(size(img{1})));
            
            if rest > 0
                for k = 1:rest;img = [img,white];end
            end
            reshapedI = reshape(img, lieshu,[])';
            clear img
            IMG = cell2mat(reshapedI);
            imwrite(IMG,sprintf('秦lineChart_2D_Cumulative_%s.png',a.formated_imagename));
            %imwrite(img,sprintf('新lineChart_2D_Cumulative_%s.png',a.formated_imagename));
        end
        
        function saveDrawInsong_lineChart_2D_Cumulative(a)
            img1 = a.Insong_drawPitchHarmoLineChart;
            img2 = a.Insong_draw2DPitchVsHarmo;
            img3 = a.Insong_drawCumulativePitchDiff;
            img = horzcat(img1,img2,img3);
            imwrite(img,sprintf('汉Insong_lineChart_2D_Cumulative_%s.png',a.formated_imagename));
        end
        
        function BinaryThresholdDrawMeanFeaturesVsRespAsLineChart(a) % draw the distribution of mean features
            
            fraglist = a.judgeFragResponse;
            if isempty(fraglist)
                return;
            end
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).maxvalue;
            end
            
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).amplitude = fraglist(u).meanfeatures.amplitude;
                toshow(u).pitch = fraglist(u).meanfeatures.pitch;
                toshow(u).AM = fraglist(u).meanfeatures.AM;
                toshow(u).mean_frequency = fraglist(u).meanfeatures.mean_frequency;
                toshow(u).FM = fraglist(u).meanfeatures.FM;
                toshow(u).peak_frequency = fraglist(u).meanfeatures.peak_frequency;
                toshow(u).entropy = fraglist(u).meanfeatures.entropy;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            
            toshow  = table2struct(sortrows(struct2table(toshow), 'resp','ascend'));
            
            
            num_toshow = [toshow.amplitude; toshow.pitch; toshow.AM; toshow.mean_frequency; toshow.FM; ...
                toshow.peak_frequency; toshow.entropy].';
            
            znum_toshow = zscore(num_toshow,0,1);
            %
            %            figure
            %            plot(znum_toshow.');
            
            figure('Position',[1997 233 1388 658],'Color','w');
            hold on
            % c = 1- rescale([toshow.resp].',0.1,1);
            
            
            % very artifical bad code
            switch a.formated_imagename
                case 'O686_34'
                    thres = 2;
                case 'R677_55'
                    thres = 6;
                case 'Y661_5'
                    thres = 6;
                case 'Y661_8'
                    thres = 2;
                case 'Y675_22'
                    thres = 3;
                case 'R677_59'
                    thres = 3;
                otherwise
                    thres = 5;
            end
            
            kemal = 0;
            for r = 1: size(znum_toshow,1)
                if toshow(r).resp > thres
                    kemal = kemal + 1;
                    defaultLW = get(gca,'LineWidth');
                    plot(znum_toshow(r,:),'Color',[1 0 0],'LineWidth',defaultLW*1.8);
                else
                    plot(znum_toshow(r,:),'Color',[0.6 0.6 0.6]);
                end
            end
            
            xlim([0,8]);
            xticks([0 1 2 3 4 5 6 7 8])
            xticklabels({'','amplitude','pitch','AM','mean_frequency','FM','peak-frequency','entropy',''});
            title(sprintf('%u of  %u song elements',kemal,length(fraglist)));
            ylabel('Zscored Feature(averaged)');
            
            saveas(gcf,sprintf('New_BinaThres_LineChartMeanFeaturesVsResp-%s.png',a.formated_imagename));
            close(gcf);
        end
        
        function saveDrawPairwiseFragmentsMeanFeaturesDistribution(a)
            
            fraglist = a.judgeFragResponse;
            if isempty(fraglist)
                return
            end
            if isempty(fraglist(1).meanfeatures) %%% Too bad !!
                return
            end
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).maxvalue;
            end
            
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).amplitude = fraglist(u).meanfeatures.amplitude;
                toshow(u).pitch = fraglist(u).meanfeatures.pitch;
                toshow(u).AM = fraglist(u).meanfeatures.AM;
                toshow(u).mean_frequency = fraglist(u).meanfeatures.mean_frequency;
                toshow(u).FM = fraglist(u).meanfeatures.FM;
                toshow(u).peak_frequency = fraglist(u).meanfeatures.peak_frequency;
                toshow(u).entropy = fraglist(u).meanfeatures.entropy;
                toshow(u).harmratio = fraglist(u).meanfeatures.harmratio;
                toshow(u).resp = fraglist(u).responseMeasure;
                
            end
            
            name = {'amplitude','pitch','AM','mean_frequency','FM','peak_frequency','entropy','harmratio'};
            ncomb = nchoosek(name,2);
            
            I = {};
            for idx = 1:size(ncomb,1) % 此处可用parfor
                
                %subplot(hangshu,lieshu,idx);
                figure('Color','w');
                hold on
                cmap = 1 - rescale([toshow.resp].');
                for  md = 1: length(toshow)
                    scatter(toshow(md).(ncomb{idx,1})', toshow(md).(ncomb{idx,2}),[],repmat(cmap(md,:),1,3),'filled');
                end
                xlabel(replace(ncomb{idx,1},'_','-'));
                ylabel(replace(ncomb{idx,2},'_','-'));
                
                frame = getframe(gcf);
                I{idx} = frame.cdata;
                close(gcf);
                
            end
            
            % PCA plot
            a.drawPCABasedOnAllFeatures;
            frame = getframe(gcf);
            I{length(I)+ 1} = frame.cdata;
            close(gcf);
            
            
            lieshu = 7;
            hangshu = ceil(length(I)/lieshu);
            
            rest = lieshu*hangshu - length(I);
            white = uint8(255*ones(size(I{1})));
            
            if rest > 0
                for k = 1:rest
                    I = [I,white];
                end
            end
            
            reshapedI = reshape(I, lieshu,[])';
            clear I
            IMG = cell2mat(reshapedI);
            imwrite(IMG,sprintf('PairwiseFragmentsMeanFeaturesDistribution_%s.png',a.formated_imagename));
            
        end
        
        function saveDrawDTWSimilarityMatrix(a)
            
            fraglist = a.judgeFragResponse;
            if isempty(fraglist)
                return;
            end
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).maxvalue;
            end
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).amplitude = fraglist(u).features.amplitude;
                toshow(u).pitch = fraglist(u).features.pitch;
                toshow(u).AM = fraglist(u).features.AM;
                toshow(u).mean_frequency = fraglist(u).features.mean_frequency;
                toshow(u).FM = fraglist(u).features.FM;
                toshow(u).peak_frequency = fraglist(u).features.peak_frequency;
                toshow(u).entropy = fraglist(u).features.entropy;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            toshow  = table2struct(sortrows(struct2table(toshow), 'resp','ascend'));
            
            
            names = {'amplitude','pitch','AM','FM','mean_frequency','peak_frequency','entropy'};
            
            I = {}; % inatialize the Image collection
            f = waitbar(0,'Start');
            
            for na = 1: length(names)
                waitbar(na/length(names),f,replace(sprintf('Generating %s similarity matrix ...',names{na}),'_','-'));
                figure('Color','w');
                sim = [];
                for s = 1: length(toshow)
                    parfor b = 1: length(toshow)
                        sim(s,b) = 1 - dtw(toshow(s).(names{na}),toshow(b).(names{na}));
                    end
                end
                imagesc(sim);
                xlabel(sprintf('%s-Similarity',names{na}));
                frame = getframe(gcf);
                I{na} = frame.cdata;
                close(gcf);
            end
            close(f);
            
            
            lieshu = 4;
            hangshu = ceil(length(I)/lieshu);
            rest = lieshu*hangshu - length(I);
            white = uint8(255*ones(size(I{1})));
            
            if rest > 0
                for k = 1:rest
                    I = [I,white];
                end
            end
            
            reshapedI = reshape(I, lieshu,[])';
            clear I
            IMG = cell2mat(reshapedI);
            imwrite(IMG,sprintf('DTWSimilarityMatrix_%s.png',a.formated_imagename));
            
        end
        
        function saveDrawDTWSimilarityMatrixBasedOnZscoredData(a)
            
            fraglist = a.judgeFragResponse;
            if isempty(fraglist)
                return;
            end
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).maxvalue;
            end
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).amplitude = zscore(fraglist(u).features.amplitude);
                toshow(u).pitch = zscore(fraglist(u).features.pitch);
                toshow(u).AM = zscore(fraglist(u).features.AM);
                toshow(u).mean_frequency = zscore(fraglist(u).features.mean_frequency);
                toshow(u).FM = zscore(fraglist(u).features.FM);
                toshow(u).peak_frequency = zscore(fraglist(u).features.peak_frequency);
                toshow(u).entropy = zscore(fraglist(u).features.entropy);
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            toshow  = table2struct(sortrows(struct2table(toshow), 'resp','ascend'));
            
            
            names = {'amplitude','pitch','AM','FM','mean_frequency','peak_frequency','entropy'};
            
            I = {}; % inatialize the Image collection
            f = waitbar(0,'Start');
            
            for na = 1: length(names)
                waitbar(na/length(names),f,replace(sprintf('Generating %s similarity matrix ...',names{na}),'_','-'));
                figure('Color','w');
                sim = [];
                for s = 1: length(toshow)
                    parfor b = 1: length(toshow)
                        sim(s,b) = 1 - dtw(toshow(s).(names{na}),toshow(b).(names{na}));
                    end
                end
                imagesc(sim);
                xlabel(sprintf('Zscored-DTWSimilarityMatrix-%s',names{na}));
                frame = getframe(gcf);
                I{na} = frame.cdata;
                close(gcf);
            end
            close(f);
            
            
            lieshu = 4;
            hangshu = ceil(length(I)/lieshu);
            rest = lieshu*hangshu - length(I);
            white = uint8(255*ones(size(I{1})));
            
            if rest > 0
                for k = 1:rest
                    I = [I,white];
                end
            end
            
            reshapedI = reshape(I, lieshu,[])';
            clear I
            IMG = cell2mat(reshapedI);
            imwrite(IMG,sprintf('ZscoredFeatureSimilarity_%s.png',a.formated_imagename));
            
        end
        
        function saveDrawCoeffOfFeaturesLinearFit(a)
            
            fraglist = a.judgeFragResponse;
            if isempty(fraglist)
                return
            end
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).maxvalue;
            end
            toshow = struct;
            for u = 1: length(fraglist)
                % for polyfits
                commonx = 1: length(fraglist(u).features.amplitude);
                toshow(u).amplitude_coef = polyfit(commonx, fraglist(u).features.amplitude,1);
                toshow(u).pitch_coef = polyfit(commonx, fraglist(u).features.pitch,1);
                toshow(u).AM_coef = polyfit(commonx, fraglist(u).features.AM,1);
                toshow(u).mean_frequency_coef = polyfit(commonx, fraglist(u).features.mean_frequency,1);
                toshow(u).FM_coef = polyfit(commonx, fraglist(u).features.FM,1);
                toshow(u).peak_frequency_coef = polyfit(commonx, fraglist(u).features.peak_frequency,1);
                toshow(u).entropy_coef = polyfit(commonx, fraglist(u).features.entropy,1);
                toshow(u).resp = fraglist(u).responseMeasure;
                
            end
            toshow  = table2struct(sortrows(struct2table(toshow), 'resp','ascend'));
            
            
            names = {'amplitude_coef','pitch_coef','AM_coef','FM_coef','mean_frequency_coef','peak_frequency_coef','entropy_coef'};
            I = {}; % inatialize the Image collection
            f = waitbar(0,'Start');
            
            for na = 1: length(names)
                waitbar(na/length(names),f,replace(sprintf('Creating %s scatter plot ...',names{na}),'_','-'));
                figure('Color','w');
                hold on
                feature_coefs = {toshow.(names{na})}.';
                feature_coefs = vertcat(feature_coefs{:});
                resp = [toshow.resp].';
                
                scatter(feature_coefs(:,1),feature_coefs(:,2),[],resp,'filled');
                xline(0,'--r');
                
                xlabel('1-th order coefficient','Interpreter','none');
                ylabel('0-th order coefficient','Interpreter','none');
                
                title(replace(sprintf('%s-1 and 0 coefficient',names{na}),'-','_'));
                frame = getframe(gcf);
                I{na} = frame.cdata;
                close(gcf);
            end
            close(f);
            
            
            lieshu = 4;
            hangshu = ceil(length(I)/lieshu);
            rest = lieshu*hangshu - length(I);
            white = uint8(255*ones(size(I{1})));
            
            if rest > 0
                for k = 1:rest
                    I = [I,white];
                end
            end
            
            reshapedI = reshape(I, lieshu,[])';
            clear I
            IMG = cell2mat(reshapedI);
            imwrite(IMG,sprintf('N-order-coefficient_%s.png',a.formated_imagename));
            
            
            
        end
        
        function saveDrawMeanFeaturesInSongVsRespAsLineChart(a) % draw the distribution of mean features
            
            fraglist = a.calResponseToWithinSongFragsFromEleinf;
            if isempty(fraglist)
                return;
            end
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).maxvalue;
            end
            
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).amplitude = fraglist(u).mean_amplitude;
                toshow(u).pitch = fraglist(u).mean_pitch;
                toshow(u).AM = fraglist(u).mean_AM;
                toshow(u).mean_frequency = fraglist(u).mean_mean_frequency;
                toshow(u).FM = fraglist(u).mean_FM;
                toshow(u).peak_frequency = fraglist(u).mean_peak_frequency;
                toshow(u).entropy = fraglist(u).mean_entropy;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            
            toshow  = table2struct(sortrows(struct2table(toshow), 'resp','ascend'));
            
            
            num_toshow = [toshow.amplitude; toshow.pitch; toshow.AM; toshow.mean_frequency; toshow.FM; ...
                toshow.peak_frequency; toshow.entropy].';
            
            znum_toshow = zscore(num_toshow,0,1);
            %            figure
            %            plot(znum_toshow.');
            
            figure('Position',[1997 233 1388 658],'Color','w');
            hold on
            c = 1- rescale([toshow.resp].',0.1,1);
            for r = 1: length(znum_toshow)
                plot(znum_toshow(r,:),'Color',repmat(c(r),3,1));
            end
            
            colormap(flip(repmat(unique(c),1,3)))
            colorbar
            xlim([0,8]);
            xticks([0 1 2 3 4 5 6 7 8])
            xticklabels({'','amplitude','pitch','AM','mean_frequency','FM','peak-frequency','entropy',''});
            
            ylabel('Zscored Feature(averaged)');
            title(sprintf('Totally %u song elements',length(fraglist)));
            saveas(gcf,sprintf('V1-WithinSongsLineChartMeanFeaturesVsResp-%s.png',a.formated_imagename));
            close(gcf);
        end
        
        function Iall = saveDrawAlignedConsFrag(a,songnames)
            % songnames 指定后则只会显示对应于songnames的fragments
            dbstop if error
            tic
            
            fragids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'frag|Frag|syl')));
            
            fraglist = a.list(fragids);
            frag_Fid = unique({fraglist.Fid}.');
            % normlist = Analysis(a.neurons{a.song_id}).normlist;
            
            subfile_frag_ids = find(~cellfun(@isempty,regexp(cellstr({a.list.Fid}.'),strjoin(frag_Fid,'|'))));
            fuckids = subfile_frag_ids;
            if exist('songnames','var')
                songnameids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),strjoin(songnames,'|'))));
                fuckids = intersect(subfile_frag_ids,songnameids);
            end
            
            fucklist = a.list(fuckids);
            
            normlist = fucklist(find(~cellfun(@isempty, regexp(cellstr({fucklist.stimuliname}.'),'norm'))));
            
            % about Norms % Redudant code
            %             normlist = Analysis(a.neurons{a.song_id}).normlist;
            %             [~,postunique] = unique(cellfun(@convert.bid,cellstr({normlist.stimuliname}.'),'Uni',0))
            %             normlist = normlist(postunique);
            
            
            % About Frag
            fragids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'frag|Frag|syl|Syl') ));
            if ~isempty(fragids)
                
                fraglist = a.list(fragids);
                
                for m = 1: length(fraglist)
                    
                    birdid = convert.bid(fraglist(m).stimuliname);
                    ids_norm = find(~cellfun(@isempty, regexp(cellstr({normlist.stimuliname}.'),birdid) ) );
                    
                    if ~isempty(ids_norm)& length(ids_norm) == 1
                        fraglist(m).sylIni = Analysis.findIni(normlist(ids_norm).plty,fraglist(m).y);
                        fraglist(m).pady = [zeros([fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext,1]);fraglist(m).plty;zeros(length(normlist(ids_norm).plty)...
                            - (fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext) - length(fraglist(m).plty),1)];
                        fraglist(m).padsptimes = cellfun( @(x) x + (fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext)/fraglist(m).fs, fraglist(m).pltsptimes,'uni',0);
                    end
                end
            end
            
            %             % About Deg
            %             fragids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'Deg|deg') ));
            %             if ~isempty(fragids)
            %
            %                 deglist = a.list(fragids);
            %                 for m = 1: length(deglist)
            %                     birdid = convert.bid(deglist(m).stimuliname);
            %                     ids_norm = find(~cellfun(@isempty, regexp(cellstr({normlist.stimuliname}.'),birdid) ) );
            %                     if ~isempty(ids_norm)& length(ids_norm) == 1
            %                         [deglist(m).sylIni,trump_diffvalue] = Analysis.findIni(normlist(ids_norm).plty,deglist(m).y);
            %                         fprintf('%s has the DIFFVALUE as: %u and INI as: %f \n ',deglist(m).stimuliname,trump_diffvalue,deglist(m).sylIni/32000);
            %                         deglist(m).pady = [zeros([deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext,1]);deglist(m).plty;zeros(length(normlist(ids_norm).plty)...
            %                             - (deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext) - length(deglist(m).plty),1)];
            %                         deglist(m).padsptimes = cellfun( @(x) x + (deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext)/deglist(m).fs, deglist(m).pltsptimes,'uni',0);
            %                     end
            %                 end
            %             end
            
            % merge the new fraglist and the deglist with the normlist
            I_of_each_column = {};
            for w = 1: length(normlist)
                
                %                 if ~isempty(fragids)
                %                     birdid = convert.bid(normlist(w).stimuliname);
                %                     ids_indeg = find(~cellfun(@isempty, regexp(cellstr({deglist.stimuliname}.'),birdid) ) );
                %                     selected_deglist = deglist(ids_indeg);
                %                     [~,temp_index] = sortrows([selected_deglist.sylIni].');
                %                     selected_deglist = selected_deglist(temp_index);
                %                 end
                
                if ~isempty(fragids)
                    birdid = convert.bid(normlist(w).stimuliname);
                    ids_infrag = find(~cellfun(@isempty, regexp(cellstr({fraglist.stimuliname}.'),birdid) ) );
                    selected_fraglist = fraglist(ids_infrag);
                    [~,temp_index] = sortrows([selected_fraglist.sylIni].');
                    selected_fraglist = selected_fraglist(temp_index);
                end
                
                
                % draw the basic figure
                Icollect = {};
                figure('Color','w');
                
                draw.two(normlist(w).plty,normlist(w).fs,a.normlist(w).pltsptimes);
                xlabel(normlist(w).stimuliname);
                frame = getframe(gcf);
                Icollect{1} = frame.cdata;
                close(gcf)
                
                if ~isempty(fragids)
                    for hh = 1: length(selected_deglist)
                        
                        figure('Color','w');
                        draw.two(selected_deglist(hh).pady,selected_deglist(hh).fs,selected_deglist(hh).padsptimes);
                        xlabel(selected_deglist(hh).stimuliname);
                        frame = getframe(gcf);
                        Icollect{1 + hh} = frame.cdata;
                        close(gcf);
                        
                    end
                end
                
                frozen_Icollect_len = length(Icollect);
                
                if ~isempty(fragids)
                    for bb = 1: length(selected_fraglist)
                        
                        figure('Color','w');
                        draw.two(selected_fraglist(bb).pady,selected_fraglist(bb).fs,selected_fraglist(bb).padsptimes);
                        xlabel(selected_fraglist(bb).stimuliname);
                        frame = getframe(gcf);
                        Icollect{frozen_Icollect_len + bb} = frame.cdata;
                        close(gcf);
                        
                    end
                end
                
                I_of_each_column{w} = vertcat(Icollect{:});
            end
            
            % this part draw individual elements from OTE songs
            if exist('fraglist','var')
                oteids = find(~cellfun(@isempty, regexp(cellstr({fraglist.stimuliname}.'),'OTE|ote|Ote') ));
                if ~isempty(oteids)
                    otelist = fraglist(oteids);
                    ote_Icollect = {};
                    for k = 1: ceil(length(otelist)/4)
                        % draw a figure for every four ote fragments
                        figure('Color','w');
                        for kk = flip([0:3])
                            if 4*k - kk <length(otelist)
                                subplot(2,4,4-kk);
                                
                                draw.spec(otelist(4*k - kk).plty,otelist(4*k - kk).fs);
                                subplot(2,4,4-kk + 4 );
                                draw.raster(otelist(4*k - kk).pltsptimes,otelist(4*k - kk).plty,otelist(4*k - kk).fs);
                                ylabel('')
                            end
                        end
                        frame = getframe(gcf);
                        ote_Icollect{k} = frame.cdata;
                        close(gcf)
                        
                    end
                    ote_img = vertcat( ote_Icollect{:});
                    I_of_each_column{length(I_of_each_column)+ 1} = ote_img;
                    
                end
            end
            
            a.drawSeparatedWaveform;
            temp = getframe(gcf);
            w_img = temp.cdata;
            I_of_each_column{length(I_of_each_column)+ 1} = w_img;
            
            % padding each I based on the maximum size of local I
            size1 = [];
            
            for oo = 1: length(I_of_each_column)
                size1(oo) = size(I_of_each_column{oo},1);
            end
            
            [max_size1,max_oo] = max(size1);
            
            Ipad = {};
            for oo = 1: length(I_of_each_column)
                localI = I_of_each_column{oo};
                Ibase= uint8(256*ones(size(I_of_each_column{max_oo})));
                Ibase(1:size(localI,1),1:size(localI,2),1:size(localI,3)) = localI;
                Ipad{oo} = Ibase;
            end
            
            Iall = horzcat(Ipad{:});
            
            imwrite(Iall,sprintf('Aligned_Normfrag_%s.png',a.formated_imagename));
            toc
            
        end
        
        function Iall = saveDrawAlignedConsDegs(a,songnames)
            dbstop if error
            
            %只有一个模式： 只针对二次播放里包含的norm songs进行degs的对齐
            tic
            
            degids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'deg|Deg')));
            
            deglist = a.list(degids);
            deg_Fid = unique({deglist.Fid}.');
            % normlist = Analysis(a.neurons{a.song_id}).normlist;
            
            
            
            subfile_deg_ids = find(~cellfun(@isempty,regexp(cellstr({a.list.Fid}.'),strjoin(deg_Fid,'|'))));
            fuckids = subfile_deg_ids;
            if exist('songnames','var')
                songnameids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),strjoin(songnames,'|'))));
                fuckids = intersect(subfile_deg_ids,songnameids);
            end
            
            fucklist = a.list(fuckids);
            
            normlist = fucklist(find(~cellfun(@isempty, regexp(cellstr({fucklist.stimuliname}.'),'norm'))));
            %             [~,postunique] = unique(cellfun(@convert.bid,cellstr({fucklist.stimuliname}.'),'Uni',0));
            %             normlist = normlist(postunique);
            
            % About Deg
            degids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'Deg|deg') ));
            if ~isempty(degids)
                
                deglist = a.list(degids);
                for m = 1: length(deglist)
                    birdid = convert.bid(deglist(m).stimuliname);
                    ids_norm = find(~cellfun(@isempty, regexp(cellstr({normlist.stimuliname}.'),birdid) ) );
                    if ~isempty(ids_norm)& length(ids_norm) == 1
                        [deglist(m).sylIni,trump_diffvalue] = Analysis.findIni(normlist(ids_norm).plty,deglist(m).y);
                        fprintf('%s has the DIFFVALUE as: %u and INI as: %f \n ',deglist(m).stimuliname,trump_diffvalue,deglist(m).sylIni/32000);
                        deglist(m).pady = [zeros([deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext,1]);deglist(m).plty;zeros(length(normlist(ids_norm).plty)...
                            - (deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext) - length(deglist(m).plty),1)];
                        deglist(m).padsptimes = cellfun( @(x) x + (deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext)/deglist(m).fs, deglist(m).pltsptimes,'uni',0);
                    end
                end
            end
            
            % merge the new fraglist and the deglist with the normlist
            I_of_each_column = {};
            for w = 1: length(normlist)
                
                if ~isempty(degids)
                    birdid = convert.bid(normlist(w).stimuliname);
                    ids_indeg = find(~cellfun(@isempty, regexp(cellstr({deglist.stimuliname}.'),birdid) ) );
                    selected_deglist = deglist(ids_indeg);
                    [~,temp_index] = sortrows([selected_deglist.sylIni].');
                    selected_deglist = selected_deglist(temp_index);
                end
                
                
                % draw the basic figure
                Icollect = {};
                figure('Color','w','Position',PM.size_wide);
                
                draw.two(normlist(w).plty,normlist(w).fs,normlist(w).pltsptimes);
                xlabel(normlist(w).stimuliname);
                frame = getframe(gcf);
                Icollect{1} = frame.cdata;
                close(gcf)
                
                
                
                if ~isempty(degids)
                    for hh = 1: length(selected_deglist)
                        figure('Color','w','Position',PM.size_wide);
                        draw.two(selected_deglist(hh).pady,selected_deglist(hh).fs,selected_deglist(hh).padsptimes);
                        xlabel(selected_deglist(hh).stimuliname);
                        frame = getframe(gcf);
                        Icollect{1 + hh} = frame.cdata;
                        close(gcf);
                    end
                end
                
                frozen_Icollect_len = length(Icollect);
                
                
                I_of_each_column{w} = vertcat(Icollect{:});
            end
            
            %             a.drawFirstWaveform;
            %             temp = getframe(gcf);
            %             w_img = temp.cdata;
            %             I_of_each_column{length(I_of_each_column)+ 1} = w_img;
            
            % padding each I based on the maximum size of local I
            size1 = [];
            
            for oo = 1: length(I_of_each_column)
                size1(oo) = size(I_of_each_column{oo},1);
            end
            
            [max_size1,max_oo] = max(size1);
            
            Ipad = {};
            for oo = 1: length(I_of_each_column)
                localI = I_of_each_column{oo};
                Ibase= uint8(256*ones(size(I_of_each_column{max_oo})));
                Ibase(1:size(localI,1),1:size(localI,2),1:size(localI,3)) = localI;
                Ipad{oo} = Ibase;
            end
            
            Iall = horzcat(Ipad{:});
            
            % imwrite(Iall,sprintf('Aligned_ConsDegs_%s.png',a.neurons{1}.neuronname));
            imwrite(Iall,sprintf('Aligned_ConsDegs_%s.png',a.formated_imagename));
            toc
            
        end
        
        function Iall = saveDrawAlignedConsReplas(a,songnames) % align replas and draw
            dbstop if error
            tic
            RONGYU = 0.5;
            
            tic
            
            replaids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'repla|Repla')));
            replalist = a.list(replaids);
            repla_Fid = unique({replalist.Fid}.');
            subfile_repla_ids = find(~cellfun(@isempty,regexp(cellstr({a.list.Fid}.'),strjoin(repla_Fid,'|'))));
            fuckids = subfile_repla_ids;
            if exist('songnames','var')
                songnameids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),strjoin(songnames,'|'))));
                fuckids = intersect(subfile_repla_ids,songnameids);
            end
            fucklist = a.list(fuckids);
            normlist = fucklist(find(~cellfun(@isempty, regexp(cellstr({fucklist.stimuliname}.'),'(?!repla|Repla)norm|Norm'))));
            
            
            
            % This is the new version 04.06.2022
            %             normids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'(?!repla|Repla)norm|Norm') ));
            %             normlist = a.list(normids);
            
            for k = 1: length(normlist)
                normlist(k).pady = [zeros(RONGYU*normlist(k).fs,1);normlist(k).plty];
                normlist(k).padsptimes = cellfun( @(x) x + RONGYU, normlist(k).pltsptimes,'uni',0);
            end
            
            % About Repla
            replaids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'Repla|repla|catego|Catego') ));
            if isempty(replaids); disp('No replas'); return; end
            if ~isempty(replaids)
                fraglist = a.list(replaids);
                ids_norm_collecting_box = [];
                for m = 1: length(fraglist)
                    
                    afterBefore = regexp(fraglist(m).stimuliname,'(?<=before-)\S*','match');
                    afterBefore = afterBefore{1};
                    birdid = convert.bid(afterBefore);
                    ids_norm = find(~cellfun(@isempty, regexp(cellstr({normlist.stimuliname}.'),birdid) ) );
                    
                    if ~isempty(ids_norm)%& length(ids_norm) == 1
                        ids_norm_1st = ids_norm(1);
                        ids_norm_collecting_box = [ids_norm_collecting_box, ids_norm];
                        afterpad_length = length(normlist(ids_norm_1st).plty) + RONGYU*normlist(ids_norm_1st).fs;  % +0.5s
                        
                        fraglist(m).pady = [zeros(afterpad_length- length(fraglist(m).plty),1);fraglist(m).plty];
                        fraglist(m).padsptimes = cellfun( @(x) x + length(zeros(afterpad_length- length(fraglist(m).plty),1))/fraglist(m).fs, fraglist(m).pltsptimes,'uni',0);
                    end
                end
            else
                return
            end
            
            ids_norm_collecting_box = unique(ids_norm_collecting_box); % Very important!
            
            
            % merge the new fraglist and the deglist with the selected_normlist
            selected_normlist = normlist(ids_norm_collecting_box);
            
            for w = 1: length(selected_normlist)
                
                if ~isempty(replaids)
                    birdid = convert.bid(selected_normlist(w).stimuliname);
                    ids_infrag = find(~cellfun(@isempty, regexp(cellstr({fraglist.stimuliname}.'),birdid) ) );
                    selected_fraglist = fraglist(ids_infrag);
                end
                
                
                % draw the basic figure
                Icollect = {};
                figure('Color','w','Position',[406 675 1378 420]);
                
                draw.two(selected_normlist(w).pady,selected_normlist(w).fs,selected_normlist(w).padsptimes);
                xlabel(selected_normlist(w).stimuliname);
                frame = getframe(gcf);
                Icollect{1} = frame.cdata;
                close(gcf)
                
                
                frozen_Icollect_len = length(Icollect);
                
                if ~isempty(replaids)
                    for bb = 1: length(selected_fraglist)
                        
                        figure('Color','w','Position',[406 675 1378 420]);
                        draw.two(selected_fraglist(bb).pady,selected_fraglist(bb).fs,selected_fraglist(bb).padsptimes);
                        xlabel(selected_fraglist(bb).stimuliname);
                        frame = getframe(gcf);
                        Icollect{frozen_Icollect_len + bb} = frame.cdata;
                        close(gcf);
                        
                    end
                end
                
                I_of_each_column{w} = vertcat(Icollect{:});
            end
            
            a.drawFirstWaveform;
            temp = getframe(gcf);
            w_img = temp.cdata;
            I_of_each_column{length(I_of_each_column)+ 1} = w_img;
            
            % padding each I based on the maximum size of local I
            size1 = [];
            
            for oo = 1: length(I_of_each_column)
                size1(oo) = size(I_of_each_column{oo},1);
            end
            
            [max_size1,max_oo] = max(size1);
            
            Ipad = {};
            for oo = 1: length(I_of_each_column)
                localI = I_of_each_column{oo};
                Ibase= uint8(256*ones(size(I_of_each_column{max_oo})));
                Ibase(1:size(localI,1),1:size(localI,2),1:size(localI,3)) = localI;
                Ipad{oo} = Ibase;
            end
            
            Iall = horzcat(Ipad{:});
            
            imwrite(Iall,sprintf('Aligned_Repla_With_Norm_%s.png',a.formated_imagename));
            toc
            
        end
        
        function saveDrawAlignedConsDegsReplasFrags(a)  % Align all together and draw
            dbstop if error
            tic
            
            % about songs been presented % Redudant code
            songlist = Analysis(a.neurons{a.song_id}).normlist;
            [~,postunique] = unique(cellstr(cellfun(@convert.bid,{songlist.stimuliname}.','Uni',0)))
            songlist = songlist(postunique);
            
            
            % About Frag
            fragids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'frag|Frag|syl|Syl') ));
            if ~isempty(fragids)
                
                fraglist = a.list(fragids);
                %                 normlist = Analysis(a.neurons{a.song_id}).normlist;
                %                 [~,postunique] = unique(cellfun(@convert.bid,[normlist.stimuliname].','Uni',0))
                %                 normlist = normlist(postunique);
                
                for m = 1: length(fraglist)
                    
                    birdid = convert.bid(fraglist(m).stimuliname);
                    ids_norm = find(~cellfun(@isempty, regexp(cellstr({songlist.stimuliname}.'),birdid) ) );
                    
                    if ~isempty(ids_norm)& length(ids_norm) == 1
                        fraglist(m).sylIni = Analysis.findIni(songlist(ids_norm).plty,fraglist(m).y);
                        fraglist(m).pady = [zeros([fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext,1]);fraglist(m).plty;zeros(length(songlist(ids_norm).plty)...
                            - (fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext) - length(fraglist(m).plty),1)];
                        fraglist(m).padsptimes = cellfun( @(x) x + (fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext)/fraglist(m).fs, fraglist(m).pltsptimes,'uni',0);
                    end
                end
            end
            
            % About Deg
            degids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'Deg|deg') ));
            if ~isempty(degids)
                
                deglist = a.list(degids);
                for m = 1: length(deglist)
                    birdid = convert.bid(deglist(m).stimuliname);
                    ids_norm = find(~cellfun(@isempty, regexp(cellstr({songlist.stimuliname}.'),birdid) ) );
                    if ~isempty(ids_norm)& length(ids_norm) == 1
                        [deglist(m).sylIni,trump_diffvalue] = Analysis.findIni(songlist(ids_norm).plty,deglist(m).y);
                        fprintf('%s has the DIFFVALUE as: %u and INI as: %f \n ',deglist(m).stimuliname,trump_diffvalue,deglist(m).sylIni/32000);
                        deglist(m).pady = [zeros([deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext,1]);deglist(m).plty;zeros(length(songlist(ids_norm).plty)...
                            - (deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext) - length(deglist(m).plty),1)];
                        deglist(m).padsptimes = cellfun( @(x) x + (deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext)/deglist(m).fs, deglist(m).pltsptimes,'uni',0);
                    end
                end
            end
            
            
            
            % About Repla
            
            % merge the new fraglist and the deglist with the songlist
            I_of_each_column = {};
            for w = 1: length(songlist)
                
                if ~isempty(degids)
                    birdid = convert.bid(songlist(w).stimuliname);
                    ids_indeg = find(~cellfun(@isempty, regexp(cellstr({deglist.stimuliname}.'),birdid) ) );
                    selected_deglist = deglist(ids_indeg);
                    [~,temp_index] = sortrows([selected_deglist.sylIni].');
                    selected_deglist = selected_deglist(temp_index);
                end
                
                if ~isempty(fragids)
                    birdid = convert.bid(songlist(w).stimuliname);
                    ids_infrag = find(~cellfun(@isempty, regexp(cellstr({fraglist.stimuliname}.'),birdid) ) );
                    selected_fraglist = fraglist(ids_infrag);
                    [~,temp_index] = sortrows([selected_fraglist.sylIni].');
                    selected_fraglist = selected_fraglist(temp_index);
                end
                
                
                % draw the basic figure
                Icollect = {};
                figure('Color','w');
                
                draw.two(songlist(w).plty,songlist(w).fs,songlist(w).pltsptimes);
                xlabel(songlist(w).stimuliname);
                frame = getframe(gcf);
                Icollect{1} = frame.cdata;
                close(gcf)
                
                if ~isempty(degids)
                    for hh = 1: length(selected_deglist)
                        
                        figure('Color','w');
                        draw.two(selected_deglist(hh).pady,selected_deglist(hh).fs,selected_deglist(hh).padsptimes);
                        xlabel(selected_deglist(hh).stimuliname);
                        frame = getframe(gcf);
                        Icollect{1 + hh} = frame.cdata;
                        close(gcf);
                        
                    end
                end
                
                frozen_Icollect_len = length(Icollect);
                
                if ~isempty(fragids)
                    for bb = 1: length(selected_fraglist)
                        
                        figure('Color','w');
                        draw.two(selected_fraglist(bb).pady,selected_fraglist(bb).fs,selected_fraglist(bb).padsptimes);
                        xlabel(selected_fraglist(bb).stimuliname);
                        frame = getframe(gcf);
                        Icollect{frozen_Icollect_len + bb} = frame.cdata;
                        close(gcf);
                        
                    end
                end
                
                I_of_each_column{w} = vertcat(Icollect{:});
            end
            
            % this part draw individual elements from OTE songs
            if exist('fraglist','var')
                oteids = find(~cellfun(@isempty, regexp(cellstr({fraglist.stimuliname}.'),'OTE|ote|Ote') ));
                if ~isempty(oteids)
                    otelist = fraglist(oteids);
                    ote_Icollect = {};
                    for k = 1: ceil(length(otelist)/4)
                        % draw a figure for every four ote fragments
                        figure('Color','w');
                        for kk = flip([0:3])
                            if 4*k - kk <length(otelist)
                                subplot(2,4,4-kk);
                                
                                draw.spec(otelist(4*k - kk).plty,otelist(4*k - kk).fs);
                                subplot(2,4,4-kk + 4 );
                                draw.raster(otelist(4*k - kk).pltsptimes,otelist(4*k - kk).plty,otelist(4*k - kk).fs);
                                ylabel('')
                            end
                        end
                        frame = getframe(gcf);
                        ote_Icollect{k} = frame.cdata;
                        close(gcf)
                        
                    end
                    ote_img = vertcat( ote_Icollect{:});
                    I_of_each_column{length(I_of_each_column)+ 1} = ote_img;
                    
                end
            end
            
            a.drawFirstWaveform;
            temp = getframe(gcf);
            w_img = temp.cdata;
            I_of_each_column{length(I_of_each_column)+ 1} = w_img;
            
            % padding each I based on the maximum size of local I
            size1 = [];
            
            for oo = 1: length(I_of_each_column)
                size1(oo) = size(I_of_each_column{oo},1);
            end
            
            [max_size1,max_oo] = max(size1);
            
            Ipad = {};
            for oo = 1: length(I_of_each_column)
                localI = I_of_each_column{oo};
                Ibase= uint8(256*ones(size(I_of_each_column{max_oo})));
                Ibase(1:size(localI,1),1:size(localI,2),1:size(localI,3)) = localI;
                Ipad{oo} = Ibase;
            end
            
            Iall = horzcat(Ipad{:});
            
            imwrite(Iall,sprintf('Aligned_Normfrag_%s.png',a.formated_imagename));
            toc
            
            
        end
        
        function IMG = draw_pltthree_subfile(a,neuronid)
            % 画此神经元的subfile对应的three plots
            dbstop if error
            n = a.neurons{neuronid};
            
            for idx = 1: length(n.e)
                n.e{idx}.pltthree;
                frame = getframe(gcf);
                I{idx} = frame.cdata;
                close(gcf);
            end
            
            figure('Color','w');
            n.draw_waveform;     % draw waveform
            frame = getframe(gcf);
            I{length(I)+ 1} = frame.cdata;
            close(gcf);
            
            n.drawPC12;
            frame = getframe(gcf);
            I{length(I)+ 1} = frame.cdata;
            close(gcf);
            
            n.drawValleyVsTime;
            frame = getframe(gcf);
            I{length(I)+ 1} = frame.cdata;
            close(gcf);
            
            n.drawPC1VsTime;
            frame = getframe(gcf);
            I{length(I)+ 1} = frame.cdata;
            close(gcf);
            
            % draw blank white
            if length(I) > 30
                lieshu = 15;
            else
                lieshu = 4; % previously it was 6
            end
            hangshu = ceil(length(I)/lieshu);
            rest = lieshu*hangshu - length(I);
            white = uint8(255*ones(size(I{1})));
            
            if rest > 0
                for k = 1:rest
                    I = [I,white];
                end
            end
            
            reshapedI = reshape(I, lieshu,[])';
            %clear I
            IMG = cell2mat(reshapedI);
            imwrite(IMG,sprintf('RawThree_%s.png',a.formated_imagename));
            
            
        end
        
        function Whether_NeuResp_To_SinFrags_Consis_Or_affected_By_Previous(a)
            dbstop if error
            tic
            
            % 首先，定义songlist
            songlist = a.list(find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'norm'))));
            
            %songlist = songlist(postunique);
            degids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'Deg|deg') ));
            if isempty(degids); return; end
            deglist = a.list(degids);
            deg_bids = unique(cellfun(@convert.bid,cellstr({deglist.stimuliname}.'),'Uni',0));
            
            I_song = {};
            for w = 1: length(deg_bids)
                
                
                Icollect = {};
                degexist_norm_list  = songlist(find(~cellfun(@isempty, regexp(cellstr({songlist.stimuliname}.'),deg_bids{w}))));
                
                for k = 1; length(degexist_norm_list)
                    smallfig = figure('Color','w','Position',[406 675 1378 420]);
                    draw.two(degexist_norm_list(k).plty,degexist_norm_list(k).fs,degexist_norm_list(k).pltsptimes);
                    xlabel(degexist_norm_list(k).stimuliname);
                    frame = getframe(smallfig); Icollect{1} = frame.cdata;close(gcf)
                    
                end
                
                I_song{w} = vertcat(Icollect{:});
                
            end
            [~,postunique] = unique(cellstr(cellfun(@convert.bid,{songlist.stimuliname}.','Uni',0)));
            songlist = songlist(postunique);
            
            % 其二，找到deressive songs对应的birdid
            degids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'Deg|deg') ));
            if isempty(degids); return; end
            deglist = a.list(degids);
            for m = 1: length(deglist)
                birdid_collect{m} = convert.bid(deglist(m).stimuliname);
                ids_norm = find(~cellfun(@isempty, regexp(cellstr({songlist.stimuliname}.'),birdid_collect{m}) ) );
                if ~isempty(ids_norm)& length(ids_norm) == 1
                    [deglist(m).sylIni,trump_diffvalue] = Analysis.findIni(songlist(ids_norm).plty,deglist(m).y);
                    fprintf('%s has the DIFFVALUE as: %u and INI as: %f \n ',deglist(m).stimuliname,trump_diffvalue,deglist(m).sylIni/32000);
                    deglist(m).pady = [zeros([deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext,1]);deglist(m).plty;zeros(length(songlist(ids_norm).plty)...
                        - (deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext) - length(deglist(m).plty),1)];
                    deglist(m).padsptimes = cellfun( @(x) x + (deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext)/deglist(m).fs, deglist(m).pltsptimes,'uni',0);
                end
            end
            
            unique_bid = unique(cellstr(birdid_collect));
            normlist = songlist(find(~cellfun(@isempty, regexp(cellstr({songlist.stimuliname}.'),strjoin(unique_bid,'|')) )));
            
            I_Deg = {}; % 最初定义
            for w = 1: length(normlist)
                
                birdid = convert.bid(normlist(w).stimuliname);
                ids_indeg = find(~cellfun(@isempty, regexp(cellstr({deglist.stimuliname}.'),birdid) ) );
                selected_deglist = deglist(ids_indeg);
                [~,temp_index] = sortrows([selected_deglist.sylIni].');
                selected_deglist = selected_deglist(temp_index);
                
                % draw the norm figure
                Icollect = {};figure('Color','w','Position',[406 675 1378 420]);
                draw.two(normlist(w).plty,normlist(w).fs,normlist(w).pltsptimes);
                xlabel(normlist(w).stimuliname);
                frame = getframe(gcf); Icollect{1} = frame.cdata;close(gcf)
                
                % draw the deg figure
                for hh = 1: length(selected_deglist)
                    figure('Color','w','Position',[406 675 1378 420]);
                    draw.two(selected_deglist(hh).pady,selected_deglist(hh).fs,selected_deglist(hh).padsptimes);
                    xlabel(selected_deglist(hh).stimuliname);
                    frame = getframe(gcf); Icollect{1 + hh} = frame.cdata;close(gcf);
                end
                
                I_Deg{w} = vertcat(Icollect{:});
            end
            
            % 其三，找到对应birdid的frags
            fragids1 = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'frag|Frag|syl|Syl') ));
            unique_bid = unique(cellstr(birdid_collect));
            fragids2 = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),strjoin(unique_bid,'|')) ));
            fragids = intersect(fragids1,fragids2);
            if ~isempty(fragids)
                fraglist = a.list(fragids);
                for m = 1: length(fraglist)
                    
                    birdid = convert.bid(fraglist(m).stimuliname);
                    ids_norm = find(~cellfun(@isempty, regexp(cellstr({songlist.stimuliname}.'),birdid) ) );
                    
                    if ~isempty(ids_norm)& length(ids_norm) == 1
                        fraglist(m).sylIni = Analysis.findIni(songlist(ids_norm).plty,fraglist(m).y);
                        fraglist(m).pady = [zeros([fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext,1]);fraglist(m).plty;zeros(length(songlist(ids_norm).plty)...
                            - (fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext) - length(fraglist(m).plty),1)];
                        fraglist(m).padsptimes = cellfun( @(x) x + (fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext)/fraglist(m).fs, fraglist(m).pltsptimes,'uni',0);
                    end
                end
            end
            
            I_Frag = {};
            for w = 1: length(normlist)
                
                if ~isempty(fragids)
                    birdid = convert.bid(normlist(w).stimuliname);
                    ids_infrag = find(~cellfun(@isempty, regexp(cellstr({fraglist.stimuliname}.'),birdid) ) );
                    selected_fraglist = fraglist(ids_infrag);
                    [~,temp_index] = sortrows([selected_fraglist.sylIni].');
                    selected_fraglist = selected_fraglist(temp_index);
                end
                
                Icollect = {};figure('Color','w','Position',[406 675 1378 420]);
                draw.two(normlist(w).plty,songlist(w).fs,normlist(w).pltsptimes);  xlabel(normlist(w).stimuliname);
                frame = getframe(gcf);   Icollect{1} = frame.cdata; close(gcf)
                
                if ~isempty(fragids)
                    for bb = 1: length(selected_fraglist)
                        figure('Color','w','Position',[406 675 1378 420]);
                        draw.two(selected_fraglist(bb).pady,selected_fraglist(bb).fs,selected_fraglist(bb).padsptimes);
                        xlabel(selected_fraglist(bb).stimuliname);
                        frame = getframe(gcf);Icollect{1 + bb} = frame.cdata;close(gcf);
                    end
                end
                I_Frag{w} = vertcat(Icollect{:});
            end
            
            
            
            % 其四，找到对应birdid的replas
            RONGYU = 0.5;
            normids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'(?!repla|Repla)norm|Norm') ));
            normlist = a.list(normids);
            
            for k = 1: length(normlist)
                normlist(k).pady = [zeros(RONGYU*normlist(k).fs,1);normlist(k).plty];
                normlist(k).padsptimes = cellfun( @(x) x + RONGYU, normlist(k).pltsptimes,'uni',0);
            end
            
            replaids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'Repla|repla|catego|Catego') ));
            if ~isempty(replaids)
                
                replalist = a.list(replaids);
                ids_norm_collecting_box = [];
                for m = 1: length(replalist)
                    
                    afterBefore = regexp(replalist(m).stimuliname,'(?<=before-)\S*','match');
                    afterBefore = afterBefore{1};
                    birdid = convert.bid(afterBefore);
                    ids_norm = find(~cellfun(@isempty, regexp(cellstr({normlist.stimuliname}.'),birdid) ) );
                    
                    if ~isempty(ids_norm)%& length(ids_norm) == 1
                        ids_norm_1st = ids_norm(1);
                        ids_norm_collecting_box = [ids_norm_collecting_box, ids_norm];
                        afterpad_length = length(normlist(ids_norm_1st).plty) + RONGYU*normlist(ids_norm_1st).fs;  % +0.5s
                        
                        replalist(m).pady = [zeros(afterpad_length- length(replalist(m).plty),1);replalist(m).plty];
                        replalist(m).padsptimes = cellfun( @(x) x + length(zeros(afterpad_length- length(replalist(m).plty),1))/replalist(m).fs, replalist(m).pltsptimes,'uni',0);
                    end
                end
            end
            
            ids_norm_collecting_box = unique(ids_norm_collecting_box); % Very important!
            selected_normlist = normlist(ids_norm_collecting_box);
            
            I_Repla = {};
            for w = 1: length(selected_normlist)
                
                if ~isempty(replaids)
                    birdid = convert.bid(selected_normlist(w).stimuliname);
                    ids_inrepla = find(~cellfun(@isempty, regexp(cellstr({replalist.stimuliname}.'),birdid) ) );
                    selected_replalist = replalist(ids_inrepla);
                end
                
                
                % draw the basic figure
                Icollect = {}; figure('Color','w','Position',[406 675 1378 420]);
                draw.two(selected_normlist(w).pady,selected_normlist(w).fs,selected_normlist(w).padsptimes);
                xlabel(selected_normlist(w).stimuliname);
                frame = getframe(gcf);  Icollect{1} = frame.cdata;  close(gcf);
                
                if ~isempty(replaids)
                    for bb = 1: length(selected_replalist)
                        figure('Color','w','Position',[406 675 1378 420]);
                        draw.two(selected_replalist(bb).pady,selected_replalist(bb).fs,selected_replalist(bb).padsptimes);
                        xlabel(selected_replalist(bb).stimuliname);
                        frame = getframe(gcf); Icollect{1 + bb} = frame.cdata; close(gcf);
                    end
                end
                I_Repla{w} = vertcat(Icollect{:});
            end
            
            
            % 其末 draw and save
            a.drawFirstWaveform;
            temp = getframe(gcf);close(gcf);
            w_img = temp.cdata;
            I_WF{1} = w_img;
            
            I_of_each_column = horzcat(I_song,I_Deg,I_Frag,I_Repla,I_WF);
            % padding each I based on the maximum size of local I
            size1 = [];
            for oo = 1: length(I_of_each_column)
                size1(oo) = size(I_of_each_column{oo},1);
            end
            
            [max_size1,max_oo] = max(size1);
            
            Ipad = {};
            for oo = 1: length(I_of_each_column)
                localI = I_of_each_column{oo};
                Ibase= uint8(256*ones(size(I_of_each_column{max_oo})));
                Ibase(1:size(localI,1),1:size(localI,2),1:size(localI,3)) = localI;
                Ipad{oo} = Ibase;
            end
            
            Iall = horzcat(Ipad{:});
            
            imwrite(Iall,sprintf('Whether_NeuResp_To_SinFrags_Coms_Or_FragsResps_affected_By_Pres_%s.png',a.formated_imagename));
            toc
            
            
        end
        
        function How_Do_Neurons_Differentiate_Sibling_Songs(a)
            dbstop if error
            tic
            
            % 首先，定义songlist
            songlist = a.list(find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'norm')))); % find all conspecific songs
            
            all_degids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'Deg|deg') ));
            all_fragids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'frag|Frag|syl|Syl') ));
            if isempty(all_degids); return; end % If no Degs at all, no need to draw
            bname_degAttached = unique(cellfun(@convert.bid,cellstr({a.list(all_degids).stimuliname}.'),'Uni',0));
            I_of_each_birdname = {};
            
            for w = 1: length(bname_degAttached)
                
                tempCollect = {};
                sub_normlist_degAttached  = songlist(... % 有对应degressive song 存在的普通所有norm songs
                    find(~cellfun(@isempty, regexp(cellstr({songlist.stimuliname}.'),bname_degAttached{w}))));
                
                for k = 1: length(sub_normlist_degAttached)
                    tempfig = figure('Color','w','Position',[7 347 2960 714]);
                    draw.two(sub_normlist_degAttached(k).plty,sub_normlist_degAttached(k).fs,sub_normlist_degAttached(k).pltsptimes);
                    xlabel(sub_normlist_degAttached(k).stimuliname);
                    frame = getframe(tempfig); tempCollect{1} = frame.cdata;close(gcf)
                end
                Inorm = vertcat(tempCollect{:});
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                deglist = a.list(intersect(all_degids,...
                    find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),bname_degAttached{w}))))); % conatin only and all deg stimuli
                ids_norm = find(... % 假设不会存在一首song的Degs在不同plx文件里出现
                    strcmp(deglist(1).Fid,{sub_normlist_degAttached.Fid}.')); % To find the normsong which (1) same birdid (2) same Fid
                if isempty(ids_norm); ids_norm = 1; disp('Warning!!!! Norm absent'); end
                for m = 1: length(deglist)
                    %if no normsong with same Fid, then use the first normsong instead
                    if ~isempty(ids_norm)&& length(ids_norm) == 1 % Based on ids_norm , pad Zeros
                        [deglist(m).sylIni,diffvalue] = Analysis.findIni(sub_normlist_degAttached(ids_norm).plty,deglist(m).y); % reference is the plty
                        fprintf('%s has the DIFFVALUE as: %u and INI as: %f \n ',deglist(m).stimuliname,diffvalue,deglist(m).sylIni/32000);
                        deglist(m).pady = [zeros([deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext,1]);deglist(m).plty;zeros(length(sub_normlist_degAttached(ids_norm).plty)...
                            - (deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext) - length(deglist(m).plty),1)];
                        deglist(m).padsptimes = cellfun( @(x) x + (deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext)/deglist(m).fs, deglist(m).pltsptimes,'uni',0);
                    end
                end
                
                % draw the norm figure
                tempCollect = {};figure('Color','w','Position',[7 347 2960 714]);
                
                draw.two(sub_normlist_degAttached(ids_norm).plty,...
                    sub_normlist_degAttached(ids_norm).fs,sub_normlist_degAttached(ids_norm).pltsptimes);
                xlabel(sub_normlist_degAttached(ids_norm).stimuliname);
                frame = getframe(gcf); tempCollect{1} = frame.cdata;close(gcf);
                % draw the deg figure
                for hh = 1: length(deglist)
                    figure('Color','w','Position',[7 347 2960 714]);
                    draw.two(deglist(hh).pady,deglist(hh).fs,deglist(hh).padsptimes);
                    xlabel(deglist(hh).stimuliname);
                    frame = getframe(gcf); tempCollect{1 + hh} = frame.cdata;close(gcf);
                end
                I_Deg = vertcat(tempCollect{:});
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fraglist = a.list(intersect(all_fragids,...
                    find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),bname_degAttached{w}))))); % conatin only and all deg stimuli
                ids_norm = find(... % 假设不会存在一首song的Degs在不同plx文件里出现
                    strcmp(fraglist(1).Fid,{sub_normlist_degAttached.Fid}.')); % To find the normsong which (1) same birdid (2) same Fid
                if isempty(ids_norm); ids_norm = 1; disp('Warning!!!! Norm absent'); end
                for m = 1: length(fraglist)
                    %birdid_collect{m} = convert.bid(deglist(m).stimuliname);
                    %if no normsong with same Fid, then use the first normsong instead
                    if ~isempty(ids_norm)& length(ids_norm) == 1
                        fraglist(m).sylIni = Analysis.findIni(sub_normlist_degAttached(ids_norm).plty,fraglist(m).y);
                        fraglist(m).pady = [zeros([fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext,1]);fraglist(m).plty;zeros(length(sub_normlist_degAttached(ids_norm).plty)...
                            - (fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext) - length(fraglist(m).plty),1)];
                        fraglist(m).padsptimes = cellfun( @(x) x + (fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext)/fraglist(m).fs, fraglist(m).pltsptimes,'uni',0);
                    end
                end
                
                % draw the norm figure
                tempCollect = {};figure('Color','w','Position',[7 347 2960 714]);
                
                draw.two(sub_normlist_degAttached(ids_norm).plty,...
                    sub_normlist_degAttached(ids_norm).fs,sub_normlist_degAttached(ids_norm).pltsptimes);
                xlabel(sub_normlist_degAttached(ids_norm).stimuliname);
                frame = getframe(gcf); tempCollect{1} = frame.cdata;close(gcf);
                % draw the deg figure
                for hh = 1: length(fraglist)
                    figure('Color','w','Position',[7 347 2960 714]);
                    draw.two(fraglist(hh).pady,fraglist(hh).fs,fraglist(hh).padsptimes);
                    xlabel(fraglist(hh).stimuliname);
                    frame = getframe(gcf); tempCollect{1 + hh} = frame.cdata;close(gcf);
                end
                I_Frag = vertcat(tempCollect{:});
                
                I_of_each_birdname{w} = {Inorm,I_Deg,I_Frag};
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            end
            
            % 其末 draw and save
            a.drawFirstWaveform;
            temp = getframe(gcf);close(gcf);
            w_img = temp.cdata;
            I_WF{1} = w_img;
            
            I_of_each_column = horzcat(I_of_each_birdname{:},I_WF);
            % padding each I based on the maximum size of local I
            size1 = [];
            for oo = 1: length(I_of_each_column)
                size1(oo) = size(I_of_each_column{oo},1);
            end
            
            [max_size1,max_oo] = max(size1);
            
            Ipad = {};
            for oo = 1: length(I_of_each_column)
                localI = I_of_each_column{oo};
                Ibase= uint8(256*ones(size(I_of_each_column{max_oo})));
                Ibase(1:size(localI,1),1:size(localI,2),1:size(localI,3)) = localI;
                Ipad{oo} = Ibase;
            end
            
            Iall = horzcat(Ipad{:});
            
            imwrite(Iall,sprintf('Sibling_Songs_%s.png',a.formated_imagename));
            toc
            
        end
        
    end
    
    methods(Static) % 静态方法
        
        function check_detail_folder(dirpath,global_eleinf)
            % This function draw distribution of song elements in spectral
            % space ( scatter / small spectrogram-rasterPlot)
            
            % extract .wav file names in the target folder
            filenames = extract.filename(dirpath,'*.wav');
            nameonly = {};
            for k = 1: length(filenames)
                [~,nameonly{k},~] = fileparts(filenames{k});
            end
            nameonly = nameonly.';
            nameonly = [nameonly{:}].';
            
            % concatenate each songname and fragid for comparing
            global_merge = {};
            for w = 1: length(global_eleinf)
                global_merge{w} = sprintf('%s-%02u',global_eleinf(w).songname,global_eleinf(w).fragid);
            end
            
            
            % extract the target song-element from the name of the target
            % directory
            temp = split(dirpath,'\');
            temp = temp{end};
            temp = split(temp,'-');
            targetname = sprintf('%s-%s-%s',temp{2},temp{3},temp{4});
            [~,targetid] = ismember( targetname,global_merge);
            
            % extract element data of song element Fragments/Replacements
            fragids = find(~cellfun(@isempty, regexp(nameonly,'Frag')))
            fragnames = nameonly(fragids);
            fragname_remove_frag = split(fragnames,'Frag-');
            fragname_remove_frag = fragname_remove_frag(:,2);
            
            % label global_eleinf by frag data
            for k = 1: length(global_eleinf)
                global_eleinf(k).whether_frag_is_tested = 0;
            end
            
            for u = 1: length(fragname_remove_frag)
                [~,beta] = ismember( fragname_remove_frag{u},global_merge);
                global_eleinf(beta).whether_frag_is_tested = 1;
            end
            
            frag_1_eleinf = global_eleinf(find([global_eleinf.whether_frag_is_tested].' == 1));
            frag_0_eleinf = global_eleinf(find([global_eleinf.whether_frag_is_tested].' == 0));
            %--% draw
            figure;
            hold on
            scatter([frag_0_eleinf.coor_1].', [frag_0_eleinf.coor_2].','k','filled');
            scatter([frag_1_eleinf.coor_1].', [frag_1_eleinf.coor_2].','r','filled');
            scatter(global_eleinf(targetid).coor_1,global_eleinf(targetid ).coor_2,[],'g','h');
            title(sprintf('Frag-%s',targetname));
            
            
            % extract data of replaced song
            replaids = find(~cellfun(@isempty, regexp(nameonly,'Repla')));
            replanames = nameonly(replaids);
            for e = 1: length(replanames)
                temp = split(replanames{e},'-before-');
                replanames{e} = temp{1};
            end
            replaname_remove_repla = split(replanames,'Repla-');
            if isempty(replaname_remove_repla)
                return % if there is no repla, the just return
            end
            replaname_remove_repla = replaname_remove_repla(:,2);
            
            
            % label global_eleinf by Repla data
            for k = 1: length(global_eleinf)
                global_eleinf(k).whether_repla_is_tested = 0;
            end
            
            for u = 1: length(replaname_remove_repla)
                [~,beta] = ismember( replaname_remove_repla{u},global_merge);
                global_eleinf(beta).whether_repla_is_tested = 1;
            end
            
            repla_1_eleinf = global_eleinf(find([global_eleinf.whether_repla_is_tested].' == 1));
            repla_0_eleinf = global_eleinf(find([global_eleinf.whether_repla_is_tested].' == 0));
            %--% draw
            figure;
            hold on
            scatter([repla_0_eleinf.coor_1].', [repla_0_eleinf.coor_2].','k','filled');
            scatter([repla_1_eleinf.coor_1].', [repla_1_eleinf.coor_2].','r','filled');
            
            if global_eleinf(targetid).fragid ~=1
                scatter(global_eleinf(targetid).coor_1,global_eleinf(targetid-1).coor_2,[],'g','h');
            end
            title(sprintf('Repla-%s',targetname));
            
            
        end
        
        function [SynIni,diff_value] = findIni(y, yfrag)
            % find the correspoding initial timestamps by aliging two time series
            %y = B521; yfrag = B521_7;
            
            %              counts = 0;
            diff_info = struct;
            parfor k = 1: length(y) - length(yfrag)
                totest = y(k:k+ length(yfrag) -1);
                
                if sum(totest) ~= 0
                    
                    diff = sum(abs(totest - yfrag));
                    
                    diff_info(k).diff = diff;
                    diff_info(k).kvalue = k;
                    %                      end
                else
                    
                    diff_info(k).diff = Inf;
                    diff_info(k).kvalue = k;
                end
                
                
            end
            
            [~,min_ids] = min([diff_info.diff].');
            SynIni = diff_info(min_ids).kvalue;
            diff_value = diff_info(min_ids).diff;
            
        end
        
        function [ConvergentIndexY,ConvergentIndexReplaY] = findConergentPointBetwenNormAndRepla(y, yrepla) % find the correspoding initial timestamps by aliging two time series
            % 此趋同点指的是趋同数据点在yrepla中的次序 （从前往后数）
            % 之后或可考虑更复杂的序列对比算法，这样的算法应更普适一些
            % 但现如今应该做的是增加当y和yrepla后补了很多零的问题
            % input的y可以是y或者plty,yrepla也可以是pltreplay或者yrepla
            dbstop if error 
            
            y = y(1:find(y,1,'last')); % remove the padded zeros， 11111111000000000 to  11111111，% 有点危险
            yrepla = yrepla(1:find(yrepla,1,'last'));
           
            freezey = y; freezeyrepla = yrepla;
            [maxlength,maxid] = max([length(y),length(yrepla)]);
            
            if maxid == 1 % ru-guo-y-bi-jiao-chang
                yrepla = [zeros(length(y)-length(yrepla),1);yrepla];
            elseif maxid==2 % ru-guo-yrepla-bi-jiao-chang
                y =  [zeros(length(yrepla)-length(y),1);y];
            end
            
            temp_convergentpoint = find(flip(y-yrepla), 1 ); % 找到第一个非零元素
            convergentpoint = maxlength - temp_convergentpoint + 2;
            
            if maxid == 1 % y更长
                ConvergentIndexY  = convergentpoint;
                ConvergentIndexReplaY = convergentpoint-(length(freezey)-length(freezeyrepla));
            elseif maxid ==2 % replay 更长
                ConvergentIndexY  = convergentpoint-(length(freezeyrepla)-length(freezey));
                ConvergentIndexReplaY = convergentpoint;
            end
            
            %             diff_info = struct;
            %             parfor k = 1: length(y) - length(yrepla)
            %                 totest = y(k:k+ length(yrepla) -1);
            %
            %                 if sum(totest) ~= 0
            %
            %                     diff = sum(abs(totest - yrepla));
            %
            %                     diff_info(k).diff = diff;
            %                     diff_info(k).kvalue = k;
            %                     %                      end
            %                 else
            %                     diff_info(k).diff = Inf;
            %                     diff_info(k).kvalue = k;
            %                 end
            %             end
            %
            %             [~,min_ids] = min([diff_info.diff].');
            %             SynIni = diff_info(min_ids).kvalue;
            %             diff_value = diff_info(min_ids).diff;
            
        end
        
        
        function [answer,pvalue] = UseTtestToJudegeRespOrNot(y,sptimes,prey,presptimes,fs)
            dbstop if error
            presdf = cal.sdf(presptimes,prey,fs,0.001,0.02);
            sdf = cal.sdf(sptimes,y,fs,0.001,0.02); % 0.001,0.004
            [maxpresdf,~] = max(presdf);
            [maxsdf,maxidx] = max(sdf);
            
            pre_frs = cal.eachTrialFiringRate(presptimes,length(prey)/fs);
            sti_frs = cal.eachTrialFiringRate(sptimes,length(y)/fs);
            [h,p] = ttest(sti_frs,pre_frs,'Tail','Right','Alpha',0.05)
            if h == 1
                answer = 1;
            elseif h == 0||isnan(h)
                answer = 0;
            end
            pvalue = p;
            
        end
    end
    
    methods % 弃用方法
        
        function Deprecated_Batch2(a)
            %a.rawthree_NeuronClassVersion;
            %a.sort_frags_by_response_strength_and_then_draw;
            a.drawMeanFeaturesVsRespAsLineChart;
            %a.drawAlignedNormDegsTwoPlots;
            %a.AlignReplasWithNormsThenDraw;
            a.rawthree_NeuronClassVersion;
            a.sort_frags_by_response_strength_and_then_draw;
            a.drawMeanFeaturesVsRespAsLineChart;
            a.drawAlignedNormDegsTwoPlots;
            a.AlignReplasWithNormsThenDraw;
            a.drawPairwiseFragmentsMeanFeaturesDistribution;
        end
        
        function fraglist = Deprecated_judgeFragResponse(a) %%% To judge whether the neuron response to a frag or not
            
            ids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'Frag|syl'))); % find all frags
            
            % ’syl'可以兼容旧的stimuli命名规则
            DUR = 0.2 ;% 200ms
            
            fraglist = a.list(ids);
            
            for n = 1: length(fraglist)
                
                post_sptimes = extract.sptimes(fraglist(n).rawsptimes, fraglist(n).zpt, fraglist(n).zpt + DUR);
                pre_sptimes = extract.sptimes(fraglist(n).rawsptimes, fraglist(n).zpt-DUR, fraglist(n).zpt );
                post_mfr = length(vertcat(post_sptimes{:}))/DUR;
                pre_mfr = length(vertcat(pre_sptimes{:}))/DUR; % mfr: mean firing rate
                fraglist(n).rs = post_mfr - pre_mfr; % response strength
                
                
                tempsum = cal.psth_frag(fraglist(n).plty,fraglist(n).fs,fraglist(n).pltsptimes);
                halfsum = sum(tempsum(end/2:end));
                fullsum = sum(tempsum);
                maxvalue = max(cal.psth_frag(fraglist(n).plty,fraglist(n).fs,fraglist(n).pltsptimes));
                fraglist(n).maxvalue = maxvalue;
                fraglist(n).halfsum = halfsum;
                fraglist(n).fullsum = fullsum;
                %if maxvalue > 6 % here the threshold is very important % originally set as 8
                if fraglist(n).rs > 51
                    fraglist(n).label = 1;
                else
                    fraglist(n).label = 0;
                end
                if isempty(find([fraglist.label].' == 1)) || length(find([fraglist.label].' == 1)) == 1
                    if fraglist(n).rs > 24
                        fraglist(n).label = 1;
                    else
                        fraglist(n).label = 0;
                    end
                    
                end
            end
            
        end
        
        function Conlist = Deprecated_evaluateConResponse(a) % Eveluate Conspecific song response
            
            ids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'norm|deg|repla'))); % find all frags
            
            % ’syl'可以兼容旧的stimuli命名规则
            
            Conlist = a.list(ids);
            
            for n = 1: length(Conlist)
                sdf = cal.sdf(Conlist(n).pltsptimes,Conlist(n).plty,Conlist(n).fs,0.001,0.004);
                [maxsdf,maxidx] = max(sdf);
                
                percentage_max = maxidx/length(sdf);
                time_max = length(Conlist(n).plty)/Conlist(n).fs*percentage_max;
                % check whether the surroding are has spikes in most of the
                % trials
                extracted_sptimes = extract.sptimes(Conlist(n).pltsptimes,time_max - 0.15, time_max + 0.15); % 前后 100ms
                
                num_of_not_empty_trials = length(find(~cellfun(@isempty, extracted_sptimes)));
                
                mean_maxsdf = maxsdf/length(Conlist(n).pltsptimes);
                tempsum = cal.psth_frag(Conlist(n).plty,Conlist(n).fs,Conlist(n).pltsptimes);
                halfsum = sum(tempsum(end/2:end));
                fullsum = sum(tempsum);
                maxvalue = max(cal.sdf(Conlist(n).pltsptimes,Conlist(n).plty,Conlist(n).fs)); % I changed the value here from cal.psth_frag to cal.sdf
                Conlist(n).sdf = cal.sdf(Conlist(n).pltsptimes,Conlist(n).plty,Conlist(n).fs);
                Conlist(n).maxvalue = maxvalue;
                Conlist(n).mean_maxsdf = mean_maxsdf;
                Conlist(n).halfsum = halfsum;
                Conlist(n).fullsum = fullsum;
                Conlist(n).maxsdf = maxsdf;
                %if maxvalue > 20 % here the threshold is very important originally 24
                if (mean_maxsdf > 2.22)
                    Conlist(n).label = 1;
                else
                    Conlist(n).label = 0;
                end
                
                if fullsum > 40 % if not 60% of the trails are not empty
                    Conlist(n).label = 1;
                end
                
                
                if (num_of_not_empty_trials/length(Conlist(n).pltsptimes)<0.5)||(num_of_not_empty_trials<5) % if not 60% of the trails are not empty
                    Conlist(n).label = 0;
                end
                
                if fullsum > 48 % if not 60% of the trails are not empty
                    Conlist(n).label = 1;
                end
                
                % A rescue
                slim_extracted_sptimes = extract.sptimes(Conlist(n).pltsptimes,time_max - 0.8, time_max + 0.8);
                slim_trials = length(find(~cellfun(@isempty,slim_extracted_sptimes)));
                if slim_trials/length(Conlist(n).pltsptimes)> 0.5 % if not 60% of the trails are not empty
                    Conlist(n).label = 1;
                end
                
                
                if (maxsdf) > 19.5&& strcmp(Conlist(n).stimuliname,'norm-Y515A-21Pulses') && strcmp(Conlist(n).plxname,'Y661_Z17')
                    disp('Incredible bug in pltsptimes of function Analysis.evaluateConResponse !!');
                    Conlist(n).label = 1;
                end
                
            end
            
        end
        
        function Deprecated_Neurons_Respond_To_Single_Frags_Or_Combinations(a)
            dbstop if error
            tic
            
            % About Deg
            degids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'Deg|deg') ));
            if isempty(degids);return;end
            
            deglist = a.list(degids);
            birdids = {};
            for m = 1: length(deglist)
                birdids{m} = convert.bid(deglist(m).stimuliname);
                ids_norm = find(~cellfun(@isempty, regexp(cellstr({normlist.stimuliname}.'),birdid) ) );
                if ~isempty(ids_norm)&& length(ids_norm) == 1
                    [deglist(m).sylIni,trump_diffvalue] = Analysis.findIni(normlist(ids_norm).plty,deglist(m).y);
                    fprintf('%s has the DIFFVALUE as: %u and INI as: %f \n ',deglist(m).stimuliname,trump_diffvalue,deglist(m).sylIni/32000);
                    deglist(m).pady = [zeros([deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext,1]);deglist(m).plty;zeros(length(normlist(ids_norm).plty)...
                        - (deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext) - length(deglist(m).plty),1)];
                    deglist(m).padsptimes = cellfun( @(x) x + (deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext)/deglist(m).fs, deglist(m).pltsptimes,'uni',0);
                end
            end
            
            
            % about Norms % Redudant code
            normids = intersect(find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'norm') )),...
                find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),birdid) )));
            
            normlist = a.list(normids);
            
            norm_bids = {};
            for k = 1:length(normlist)
                norm_bids{k} = convert.bid(normlist(k).stimuliname);
            end
            unique_bids = unique(cellstr(norm_bids));
            
            % merge the new fraglist and the deglist with the normlist
            I_of_each_column = {};
            for w = 1: length(unique_bids)
                
                if ~isempty(degids)
                    bid_tosearch = convert.bid(unique_bids{w});
                    ids_indeg = find(~cellfun(@isempty, regexp(cellstr({deglist.stimuliname}.'),bid_tosearch) ) );
                    selected_deglist = deglist(ids_indeg);
                    [~,temp_index] = sortrows([selected_deglist.sylIni].');
                    selected_deglist = selected_deglist(temp_index);
                end
                
                
                % draw the basic figure
                Icollect = {};
                figure('Color','w','Position',[406 675 1378 420]);
                
                draw.two(normlist(w).plty,normlist(w).fs,normlist(w).pltsptimes);
                xlabel(normlist(w).stimuliname);
                frame = getframe(gcf);
                Icollect{1} = frame.cdata;
                close(gcf)
                
                
                
                if ~isempty(degids)
                    for hh = 1: length(selected_deglist)
                        figure('Color','w','Position',[406 675 1378 420]);
                        draw.two(selected_deglist(hh).pady,selected_deglist(hh).fs,selected_deglist(hh).padsptimes);
                        xlabel(selected_deglist(hh).stimuliname);
                        frame = getframe(gcf);
                        Icollect{1 + hh} = frame.cdata;
                        close(gcf);
                    end
                end
                
                frozen_Icollect_len = length(Icollect);
                I_of_each_column{w} = vertcat(Icollect{:});
            end
            
            a.drawFirstWaveform;
            temp = getframe(gcf);
            w_img = temp.cdata;
            I_of_each_column{length(I_of_each_column)+ 1} = w_img;
            
            % padding each I based on the maximum size of local I
            size1 = [];
            
            for oo = 1: length(I_of_each_column)
                size1(oo) = size(I_of_each_column{oo},1);
            end
            
            [max_size1,max_oo] = max(size1);
            
            Ipad = {};
            for oo = 1: length(I_of_each_column)
                localI = I_of_each_column{oo};
                Ibase= uint8(256*ones(size(I_of_each_column{max_oo})));
                Ibase(1:size(localI,1),1:size(localI,2),1:size(localI,3)) = localI;
                Ipad{oo} = Ibase;
            end
            
            Iall = horzcat(Ipad{:});
            
            % imwrite(Iall,sprintf('Aligned_ConsDegs_%s.png',a.neurons{1}.neuronname));
            imwrite(Iall,sprintf('Aligned_ConsDegs_%s.png',a.formated_imagename));
            toc
        end
        
        function Deprecated_drawfrag(a,keyword,rangeratio,ids)
            
            mergedeleinf = a.all_eleinf;
            d = Display(a.list);
            
            if exist('keyword','var')
                
                if exist('rangeratio','var')
                    if exist('ids','var')
                        d.showfrag(keyword,mergedeleinf,rangeratio,ids);
                    else
                        d.showfrag(keyword,mergedeleinf,rangeratio);
                    end
                    
                else
                    d.showfrag(keyword,mergedeleinf);
                end
            else
                for k = 1: length(a.fragnames)
                    d.showfrag(a.fragnames{k},mergedeleinf);
                end
            end
            
        end
        
        function Deprecated_ThreeWithPitch(a)
            % draw three plots
            e_objects = getAllEphysObject(a);
            
            for idx = 1: length(e_objects)
                e_objects{idx}.threeWithFreqRelatedFeatures; % newer version of threeplot drawing method
                frame = getframe(gcf);
                I{idx} = frame.cdata;
                close(gcf);
            end
            
            figure('Color','w','Position',PM.size3plots);
            a.neurons{a.song_id}.draw_waveform;     % draw waveform
            % a.draw_waveform; % this is the original draw function
            frame = getframe(gcf);
            I{length(I)+ 1} = frame.cdata;
            close(gcf);
            
            % draw blank white
            lieshu = 10;
            hangshu = ceil(length(I)/lieshu);
            rest = lieshu*hangshu - length(I);
            white = uint8(255*ones(size(I{1})));
            
            if rest > 0
                for k = 1:rest
                    I = [I,white];
                    %                     ax = gcf;
                    %                     ax.Position(3) = 560;
                    %                     ax.Position(4) = 420;
                end
            end
            
            reshapedI = reshape(I, lieshu,[])';
            clear I
            IMG = cell2mat(reshapedI);
            imwrite(IMG,sprintf('Three_%s.png',a.formated_imagename));
            
        end
        
        function Deprecated_saveDraw_replas_only(a)
            % function to draw replas, but not aligned with cons
            d = Display(a.list);
            repla = d.findrepla();
            
            if isempty(repla)
                return
            end
            
            normlist = d.findnorm();
            
            replanames = {};
            for k = 1:length(repla)
                temp = strsplit(repla(k).stimuliname,'before')
                temp = regexp(temp{2},'[BGYOR]\d{3}','match');
                replanames{k} = temp{1};
            end
            replanames = replanames.';
            unique_replanames = unique(replanames);
            
            norm_ids = [];
            for k = 1: length(unique_replanames)
                temp = find(~cellfun(@isempty, regexp(cellstr({normlist.stimuliname}.'),unique_replanames{k})));
                if length(temp)>1 % this code is dangerous
                    norm_ids(k) = temp(length(temp));
                end
                repla(length(repla) + 1) = normlist(norm_ids(k));
            end
            repla = Display.f0(repla);
            
            
            I = {};
            for idx = 1: length(repla)
                figure('Color','w','Position',[2031 536 732 281]);
                draw.three(repla(idx).f0y,repla(idx).fs,repla(idx).f0sptimes); % newer version of threeplot drawing method
                frame = getframe(gcf);
                I{idx} = frame.cdata;
                close(gcf);
            end
            
            
            lieshu = 3;
            hangshu = ceil(length(I)/lieshu);
            rest = lieshu*hangshu - length(I);
            white = uint8(255*ones(size(I{1})));
            
            if rest > 0
                for k = 1:rest
                    I = [I,white];
                    %                     ax = gcf;
                    %                     ax.Position(3) = 560;
                    %                     ax.Position(4) = 420;
                end
            end
            
            reshapedI = reshape(I, lieshu,[])';
            clear I
            IMG = cell2mat(reshapedI);
            
            imwrite(IMG,sprintf('ReplaThree-%s.png',a.formated_imagename));
            
            %deg = Display.descend(deg);
            
            %replalist = d.info(find(~cellfun(@isempty, regexp({d.info.stimuliname}.','Repla|catego'))));
            
        end
        
        function Deprecated_saveDraw_deg(a)  % 这个不太对,应该说非常不对
            
            dbstop if error
            
            d = Display(a.list);
            deg = d.finddeg();
            
            if isempty(deg)
                return
            end
            
            normlist = d.findnorm();
            
            deganames = {};
            for k = 1:length(deg)
                temp = regexp(deg(k).stimuliname,'[BGYOR]\d{3}|Fcall|Mcall|Het','match');
                deganames{k} = temp{1};
            end
            deganames = deganames.';
            unique_replanames = unique(deganames);
            
            norm_ids = [];
            for k = 1: length(unique_replanames)
                temp = find(~cellfun(@isempty, regexp(cellstr({normlist.stimuliname}.'),unique_replanames{k})));
                norm_ids(k) = temp(end); % if there are multiple values in temp, using the last one
                deg(length(deg) + 1) = normlist(norm_ids(k));
            end
            deg = Display.f0(deg);
            
            
            I = {};
            for idx = 1: length(deg)
                figure('Color','w','Position',[2031 536 732 281]);
                draw.three(deg(idx).f0y,deg(idx).fs,deg(idx).f0sptimes); % newer version of threeplot drawing method
                frame = getframe(gcf);
                I{idx} = frame.cdata;
                close(gcf);
            end
            
            
            lieshu = 3;
            hangshu = ceil(length(I)/lieshu);
            rest = lieshu*hangshu - length(I);
            white = uint8(255*ones(size(I{1})));
            
            if rest > 0
                for k = 1:rest
                    I = [I,white];
                    %                     ax = gcf;
                    %                     ax.Position(3) = 560;
                    %                     ax.Position(4) = 420;
                end
            end
            
            reshapedI = reshape(I,[],hangshu)';
            clear I
            IMG = cell2mat(reshapedI);
            
            imwrite(IMG,sprintf('DegThree-%s.png',a.formated_imagename));
            
            %deg = Display.descend(deg);
            
            %replalist = d.info(find(~cellfun(@isempty, regexp({d.info.stimuliname}.','Repla|catego'))));
            
        end
        
        function Deprecated_saveDrawAlignFragsWithSong(a)
            dbstop if error
            
            normlist = a.normlist;
            
            ids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'Frag'))); % find all frags
            fraglist = a.list(ids);
            
            
            % generate addlagy and addlagsptimes for plotting
            for k = 1:length(normlist)
                
                splited = split( normlist(k).stimuliname,'-');
                songname = splited{3}; % name of the song
                this_fragids = find( ~cellfun(@isempty, regexp(cellstr({fraglist.stimuliname}.'),songname)));
                lensong = length(normlist(k).plty);
                
                normlist(k).fragids = this_fragids;
                
                for f = 1: length(this_fragids)
                    timelag = finddelay(highpass(fraglist(this_fragids(f)).y,450,32000),normlist(k).plty); % Dangerous!!!!!!
                    disp('Dangerous code exist!!!')
                    fraglist(this_fragids(f)).addlagy = [zeros(timelag,1);fraglist(this_fragids(f)).y; zeros((lensong- timelag - length(fraglist(this_fragids(f)).y)),1)];
                    fraglist(this_fragids(f)).addlagsptimes = cellfun(@(x) x+timelag/32000-fraglist(this_fragids(f)).pltext, fraglist(this_fragids(f)).pltsptimes,'un',0); % every sptimes will be added a timelag time
                    fraglist(this_fragids(f)).timelag = timelag;
                end
                
            end
            
            
            % draw Two plots
            
            I = {};
            for k = 1: length(normlist)
                
                figure('Position',[282 759 1444 272],'Color','w');
                draw.two(normlist(k).plty,normlist(k).fs,normlist(k).pltsptimes);
                frame = getframe(gcf);
                jihe{1} = frame.cdata;
                close(gcf);
                
                for v = 1:length(normlist(k).fragids)
                    localids = normlist(k).fragids(v);
                    figure('Position',[282 759 1444 272],'Color','w');
                    draw.two(fraglist(localids ).addlagy,fraglist(localids ).fs,fraglist(localids ).addlagsptimes);
                    frame = getframe(gcf);
                    jihe{1 + v} = frame.cdata;
                    close(gcf);
                end
                
                
                I{k} = vertcat(jihe{:});
                
            end
            
            for w = 1: length(I)
                imwrite(I{w},sprintf('Three-%u.png',w));
            end
            
            %             n.draw_waveform;     % draw waveform
            %             frame = getframe(gcf);
            %             I{length(I)+ 1} = frame.cdata;
            %             close(gcf);
            %
            % draw blank white
            %             lieshu = 9;
            %             hangshu = ceil(length(I)/lieshu);
            %             rest = lieshu*hangshu - length(I);
            %             white = uint8(255*ones(size(I{1})));
            %
            %             if rest > 0
            %                 for k = 1:rest
            %                     I = [I,white];
            %                     %                     ax = gcf;
            %                     %                     ax.Position(3) = 560;
            %                     %                     ax.Position(4) = 420;
            %                 end
            %             end
            %
            %             reshapedI = reshape(I, lieshu,[])';
            %             clear I
            %             IMG = cell2mat(reshapedI);
            %             imwrite(IMG,sprintf('Three.png'));
            %
            
            
            
        end
        
    end
    
end

