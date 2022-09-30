classdef Neuron < handle & Analysis_basicDrawings %& Analysis_acousticFeatures  & Analysis_categorization
    % The class to do all the analysis for neu single neuron
    %   Detailed explanation goes here
    
    properties % 必须
        neurons % raw features
        %list_sum %???
        %group %???
        
        
        % Class
        
       % source    % read xlxs into struct
       norm
       
       deg % class to process degressive song
       
       frag
       
       repla
       
        list
        %replalist
        
%         fragnames % neu cell of unique frag names
%         degnames
%         replanames
%         
%         normlist % list of norm-songs
%         all_eleinf
%         conspe_eleinf
%         
%         targets % unique(targetnames)
%         
%         neuinfo
%         uniqueid
%         birdid
%         %formated_imagename
%         
%         % ids of neuros corresponding to different stimuli
%         song_id_redundant
%         song_id
%         song_only_id % the neuron which has only song stimuli
%         frag_id
%         deg_id
%         repla_id
%         plx_data_fs
%         other_id
        
        birdname
        zpid
        channelname
        unitname
        formated_imagename
        
        insonglist
        test
        fig_size1 % size of figure ( and also position)
    end
    
    properties(Access = private) % 可选
        figdata
    end
    
    methods % 核心方法
        
        function neu = Neuron(neurons)
            % 构造方法
            if exist('neurons','var')
                
                if isa(neurons,'Experiment'); neu.neurons{1} = neurons;else; neu.neurons = neurons; end % 输入可以是单个Neuron或者多个neuron组成的cell
                
                neu.getNeuronInfo;
                neu.setStimuliCorrespondingNeuronId;
                
                % if norm stimuli exist in multiple Experiment files, then
                % delect the first one which has only norm stimuli, as it is
                % far away from others in time
                neu.updatelist;
                %neu.normlist = neu.neurons{neu.song_id}.toList;%neu.sort;
                
                temp = regexp(neu.neurons{1}.plxname,'[RBOYRG]\d{3}','match');
                % set birdid uniqueid and formated_imagename
                if ~isempty(temp)
                    neu.birdid = temp{1};
                end
                neu.zpid = regexp(neu.neurons{1}.plxname,'[ZP]\d{2}','match');
                neu.channelname = neu.neurons{1}.channelname;
                neu.unitname = neu.neurons{1}.unitname;
                try
                    neu.formated_imagename = sprintf('%s_%s_%s_%u',neu.birdid,neu.zpid{1},neu.channelname,neu.unitname);
                catch ME
                end
                neu.fig_size1 = [2091 -14 755 620];
                %neu.judgeFragResp;
                %neu.judgeConResp;
                %neu.writeFigdata;
            end
            
            
            % 构造subclass about_repla
            neu.repla = Repla(list);
            neu.deg = Deg(list);
            neu.frag = Frag(list);
        end
        
        function neu = updatelist(neu)
            % regenerate A.list
            to_remove_id = intersect(neu.song_only_id,setdiff(neu.song_id_redundant,neu.song_id));%???
            to_calculate = setdiff(1: length(neu.neurons),to_remove_id);
            whether_update_figure_or_not = 1;
            for k = 1: length(to_calculate)
                
                templist = neu.neurons{to_calculate(k)}.toList(whether_update_figure_or_not);
                [templist.whichNeuron] = deal(k);
                lists{k} = templist;
            end
            neu.list = horzcat(lists{:});
        end
        
        function neu = getNeuronInfo(neu)
            %  To generate info of each member of an Neuron object
            
            neu.neuinfo = struct;
            
            for k = 1: length(neu.neurons)
                neu.neuinfo(k).uniqueid = neu.neurons{k}.uniqueid;
                neu.neuinfo(k).neuronname = neu.neurons{k}.neuronname;
                neu.neuinfo(k).keywords = [];
                
                if length(find(~cellfun(@isempty,regexp(cellstr({neu.neurons{k}.slist.name}.'),'norm|Norm|song|Song')))) > 16
                    neu.neuinfo(k).keywords = [neu.neuinfo(k).keywords,"song"];
                end
                
                if ~isempty(find(~cellfun(@isempty,regexp(cellstr({neu.neurons{k}.slist.name}.'),'syl|Syl|Ele|ele|frag|Frag'))))
                    neu.neuinfo(k).keywords = [neu.neuinfo(k).keywords,"frag"];
                end
                
                if ~isempty(find(~cellfun(@isempty,regexp(cellstr({neu.neurons{k}.slist.name}.'),'deg|Deg'))))
                    neu.neuinfo(k).keywords = [neu.neuinfo(k).keywords,"deg"];
                end
                
                if ~isempty(find(~cellfun(@isempty,regexp(cellstr({neu.neurons{k}.slist.name}.'),'repla|Repla|catego|Catego'))))
                    neu.neuinfo(k).keywords = [neu.neuinfo(k).keywords,"repla"];
                end
                
                if isempty(neu.neuinfo(k).keywords)
                    neu.neuinfo(k).keywords = [neu.neuinfo(k).keywords,"other"];
                end
                
            end
        end

    end
    
    methods(Access = private)% 内部计算方法
        
        function neu = calHarmRatio(neu)
            % calculate harmonic noise ratio
            
            window_size = min([neu.list.leny].');
            for k = 1:length(neu.list)
                neu.list(k).features.harmratio = harmonicRatio(neu.list(k).y,neu.list(k).fs,'Window',hamming(window_size,"periodic"),...
                    'OverlapLength',round(window_size*2/3) );
                %此处为照顾很短的frag改动了window size， 但或许更好的方法是放弃很短frag的数据
                neu.list(k).meanfeatures.harmratio = mean(neu.list(k).features.harmratio);
            end
        end
        
        function neu = writeFigdata(neu)
            % 生成这个Analysis的所有three plot的图片
            figmat = {};
            for k = 1: length(neu.neurons)
                figmat{k} = neu.neurons{k}.writeFigdata
            end
            neu.figdata = horzcat(figmat{:});
        end
        
        function neu = judgeFragResp(neu)
            % 判断对frag 是否反应，通过自己定义的复杂的机制
            ids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'Frag|syl'))); % find all frags
            
            % ’syl'可以兼容旧的stimuli命名规则
            DUR = 0.2 ;% 200ms
            
            
            for n = 1: length(ids)
                thisi = ids(n);
                post_sptimes = Extract.sptimes(neu.list(thisi).rawsptimes, neu.list( thisi).zpt, neu.list( thisi).zpt + DUR);
                pre_sptimes = Extract.sptimes(neu.list(thisi).rawsptimes, neu.list(thisi).zpt-DUR, neu.list(thisi).zpt );
                post_mfr = length(vertcat(post_sptimes{:}))/DUR;
                pre_mfr = length(vertcat(pre_sptimes{:}))/DUR; % mfr: mean firing rate
                neu.list(thisi).rs = post_mfr - pre_mfr; % response strength
                
                tempsum = Cal.psth_frag(neu.list(thisi).plty,neu.list(thisi).fs,neu.list(thisi).pltsptimes);
                halfsum = sum(tempsum(end/2:end));
                fullsum = sum(tempsum);
                maxvalue = max(Cal.psth_frag(neu.list(thisi).plty,neu.list(thisi).fs,neu.list(thisi).pltsptimes));
                neu.list(thisi).maxvalue = maxvalue;
                neu.list(thisi).halfsum = halfsum;
                neu.list(thisi).fullsum = fullsum;
                %if maxvalue > 6 % here the threshold is very important % originally set as 8
                if neu.list(thisi).rs > 51
                    neu.list(thisi).label = 1;
                else
                    neu.list(thisi).label = 0;
                end
                if isempty(find([neu.list.label].' == 1)) || length(find([neu.list.label].' == 1)) == 1
                    if neu.list(n).rs > 24
                        neu.list(n).label = 1;
                    else
                        neu.list(n).label = 0;
                    end
                    
                end
            end
            
        end
        
        function neu = judgeFragResp_FR(neu)
            % 判断对frag 是否反应，通过 Firing rate
            for k = 1: length(neu.neurons)
                for kk = 1: length( neu.neurons{k}.e)
                    neu.neurons{k}.e{kk}.setExtAndAllocate;
                end
            end
            neu.updatelist;
            
            ids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'Frag|frag|syl|ele'))); % find all norms
            % ’syl'可以兼容旧的stimuli命名规则
            
            for n = 1: length(ids)
                thisi = ids(n);
                
%                 presdf = Cal.sdf(neu.list(thisi).prejudgerespsptimes,zeros(length(neu.list(thisi).judgerespy),1),neu.list(thisi).fs,0.001,0.02);
%                 sdf = Cal.sdf(neu.list(thisi).judgerespsptimes,neu.list(thisi).judgerespy,neu.list(thisi).fs,0.001,0.02); % 0.001,0.004
%                 
                presdf = Cal.sdf(neu.list(thisi).prejudgerespsptimes,zeros(length(neu.list(thisi).judgerespy),1),neu.list(thisi).fs,0.001,0.02);
                sdf = Cal.sdf(neu.list(thisi).judgerespsptimes,neu.list(thisi).judgerespy,neu.list(thisi).fs,0.001,0.02); % 0.001,0.004
                [maxpresdf,~] = max(presdf);
                [maxsdf,maxidx] = max(sdf);
                
                
                pre_frs = Cal.eachTrialFiringRate(neu.list(thisi).prejudgerespsptimes,length(neu.list(thisi).judgerespy)/neu.list(thisi).fs);
                sti_frs = Cal.eachTrialFiringRate(neu.list(thisi).judgerespsptimes,length(neu.list(thisi).judgerespy)/neu.list(thisi).fs);
                [h,p] = ttest(sti_frs,pre_frs,'Tail','Right','Alpha',0.05)
                neu.list(thisi).pvalue = p;
                neu.list(thisi).label = 0; % 初始化
                if h == 1
                    neu.list(thisi).label = 1;
                    %                     if (num_of_not_empty_trials/length(neu.list(thisi).pltsptimes)<0.5)||(num_of_not_empty_trials<5) % if not 60% of the trails are not empty
                    %                         neu.list(thisi).label = 0;
                    %
                    
                end
            end
        end
        
        function neu = judgeConResp(neu)
            % 判断对Cons是否反应
            % firstly update e objectys 以后可以删掉这个部分
            for k = 1: length(neu.neurons)
                for kk = 1: length( neu.neurons{k}.e)
                    neu.neurons{k}.e{kk}.setExtAndAllocate;
                end
            end
            neu.updatelist;
            
            ids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'norm|deg'))); % find all norms
            % ’syl'可以兼容旧的stimuli命名规则
            
            for n = 1: length(ids)
                thisi = ids(n);
                
                presdf = Cal.sdf(neu.list(thisi).prejudgerespsptimes,zeros(length(neu.list(thisi).judgerespy),1),neu.list(thisi).fs,0.001,0.02);
                sdf = Cal.sdf(neu.list(thisi).judgerespsptimes,neu.list(thisi).judgerespy,neu.list(thisi).fs,0.001,0.02); % 0.001,0.004
                %figure; plot(sdf);
                % figure; Draw.three(neu.list(thisi).judgerespy,neu.list(thisi).fs,neu.list(thisi).judgerespsptimes);
                % figure; Draw.three(neu.list(thisi).plty,neu.list(thisi).fs,neu.list(thisi).pltsptimes);
                [maxpresdf,~] = max(presdf);
                [maxsdf,maxidx] = max(sdf);
                percentage_max = maxidx/length(sdf);
                time_max = length(neu.list(thisi).judgerespy)/neu.list(thisi).fs*percentage_max;
                % check whether the surroding are has spikes in most of the
                % trials
                extracted_sptimes = Extract.sptimes(neu.list(thisi).judgerespsptimes,time_max - 0.15, time_max + 0.15); % 前后 100ms
                num_of_not_empty_trials = length(find(~cellfun(@isempty, extracted_sptimes)));
                
                % mean_maxsdf = maxsdf/length(neu.list(thisi).judgerespsptimes);
                
                neu.list(thisi).label = 0; % 初始化
                if (maxsdf) > 17 && maxsdf > maxpresdf %如果是 time-locked response
                    neu.list(thisi).label = 1;
                    
                    if (num_of_not_empty_trials/length(neu.list(thisi).pltsptimes)<0.5)||(num_of_not_empty_trials<5) % if not 60% of the trails are not empty
                        neu.list(thisi).label = 0;
                    end
                    
                elseif mean(sdf)> 9*mean(presdf) && mean(sdf)>0.6 % 如果不是 time-locked response
                    neu.list(thisi).label = 1;   %  set to 2 ,biao ming shi fei time-locked response
                end
                
                
                
                
                
                %                 figure;
                %
                %                 Draw.three(neu.list(thisi).plty,neu.list(thisi).fs,neu.list(thisi).pltsptimes);
                %                 title(sprintf('Label is %u',neu.list(thisi).label) );
                %                 close(gcf)
                
            end
            
        end

        function neu = judgeDegResp(neu)
            % evaluate the responsiveness of song degressive deletion
            
            % based on the name of replaced song, find out the corresponding % norm song
        end
        
        function neu = getReplalistFromA(neu)
            % 从 Neuron 文件中提取neuron对replas的反应，已经具体的replas的名称
            neu.judgeReplaResp;
            replalist = neu.list(find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'Repla'))));
            
            Replst_fnames = [replalist.stimuliname].';
            locY = @(x,y) x{y}; % get element in location y
            for k = 1:length(replalist)
                splited = split(Replst_fnames{k},{'-before-','-gapis-'});
                replalist(k).bname1 = Convert.bid(splited{1});
                replalist(k).fid1 = str2num(locY(split(splited{1},'-'),4));
                replalist(k).bname2 = Convert.bid(splited{2});
                replalist(k).fid2 = str2num(locY(split(splited{2},'-'),3));
                replalist(k).concat1 = sprintf('%s-%02u',replalist(k).bname1,replalist(k).fid1);
                replalist(k).concat2 = sprintf('%s-%02u',replalist(k).bname2,replalist(k).fid2);
                replalist(k).fullrepla = sprintf('%s-%02u-%s-%02u',replalist(k).bname1,replalist(k).fid1,replalist(k).bname2,replalist(k).fid2);
            end
            unique({replalist.concat1}.') % how many unique pre syllable
            unique({replalist.concat2}.') % how many unique post syllable
            unique(vertcat({replalist.concat1}.',{replalist.concat2}.'))  % how many syllable to be classified in total
            unique({replalist.fullrepla}.') % how many unique transition in neurons intotal?
            
            neu.replalist = replalist;
        end
                       
        function neu = set_eleinf(neu,eleinf)
            %从外部传入eleinf这个变量
            if isa(eleinf,'struct')
                neu.all_eleinf = eleinf;
            elseif isa(eleinf,'string')|| isa(eleinf,'char')
                loaded = load(eleinf);
                neu.all_eleinf = loaded.all_eleinf;
            end
            
            conspe_ids = find( ~cellfun(@isempty, regexp([neu.all_eleinf.songname].','CON|SPE')) );
            neu.conspe_eleinf = neu.all_eleinf(conspe_ids);
        end
        
        function neu = splitStimuliResponsePairsToDifferentTypes(neu)
            % 从list提取出norm,frag,deg,repla等几个子集sublist
            % frag
            fragidx = find(~cellfun(@isempty, regexp(cellstr({neu.list(:).stimuliname}.'),'Frag|syl|Syl|Ele|ele|sim|Sim')));
            fraglist = neu.list(fragidx);
            for k = 1: length(fraglist)
                
                temp = strsplit(fraglist(k).stimuliname,'-');
                
                fragnames{k} = sprintf('%s-%s-%s',temp{1},temp{2},temp{3});
            end
            
            if exist('fragnames','var')
                neu.fragnames = unique(fragnames);
            end
            
            
            % deg
            degidx = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'deg|Deg')));
            deglist = neu.list(degidx);
            for k = 1: length(deglist)
                
                temp = strsplit(deglist(k).stimuliname,'-');
                
                degnames{k} = sprintf('%s-%s-%s',temp{1},temp{2},temp{3});
            end
            if exist('degnames','var')
                neu.degnames = unique(degnames);
            end
            
            
            % repla
            replaidx = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'Repla|repla|catego|Catego')));
            replalist = neu.list(replaidx);
            for k = 1: length(replalist)
                
                temp = strsplit(replalist(k).stimuliname,'-');
                
                replanames{k} = sprintf('%s-%s-%s',temp{1},temp{2},temp{3});
            end
            if exist('replanames','var')
                neu.replanames = unique(replanames);
            end
            
            % norm
            normidx = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'norm')));
            neu.normlist = neu.list(normidx);
            
            % target
            targetidx = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'Repla')));
            targetlist = neu.list(targetidx);
            for k = 1: length(targetlist)
                
                temp = strsplit(targetlist(k).stimuliname,'-');
                
                targetnames{k} = sprintf('%s-%s-%s',temp{6},temp{7},temp{8});
            end
            
            if exist('targetnames','var')
                neu.targets = unique(targetnames);
            end
            
        end
        
        function neu =setStimuliCorrespondingNeuronId(neu)
            % 当analysis的neuron对应不同刺激类型，如frag，norm时，阐明各neuron对应的刺激类型
            neu.song_id_redundant = [];
            neu.frag_id = [];
            neu.deg_id = [];
            neu.repla_id = [];
            neu.song_only_id = [];
            for k = 1: length(neu.neuinfo)
                if ismember("song",neu.neuinfo(k).keywords)
                    neu.song_id_redundant = [neu.song_id_redundant,k];
                end
                
                if ismember("frag",neu.neuinfo(k).keywords)
                    neu.frag_id = [neu.frag_id,k];
                end
                
                if ismember("deg",neu.neuinfo(k).keywords)
                    neu.deg_id = [neu.deg_id,k];
                end
                
                if ismember("repla",neu.neuinfo(k).keywords)
                    neu.repla_id = [neu.repla_id,k];
                end
                
                if strcmp("song",neu.neuinfo(k).keywords)
                    neu.song_only_id = [neu.song_only_id,k];
                end
                
                if strcmp("other",neu.neuinfo(k).keywords)
                    neu.song_only_id = [neu.song_only_id,k];
                end
                
                
            end
            
            if length(neu.song_id_redundant) == 1
                neu.song_id = neu.song_id_redundant;
            elseif length(neu.song_id_redundant)~=0
                
                neu.song_id = min(neu.song_id_redundant);
                %[~,neu.song_id] = max(cellfun(@length,{neu.neuinfo(neu.song_id_redundant).keywords}));
                % find the id which have max length of keywords, if the
                % result is multipl ids, then select the first one
                % But maybe the last one will be much proper!
            else
                neu.song_id = 1; % Very bad meaningless code
            end
            
        end
        
        function neu = judgeConResp_FR(neu)
            % 判断对Cons是否反应，通过 Firing rate
            % firstly update e objectys 以后可以删掉这个部分
            for k = 1: length(neu.neurons)
                for kk = 1: length( neu.neurons{k}.e)
                    neu.neurons{k}.e{kk}.setExtAndAllocate;
                end
            end
            neu.updatelist;
            
            ids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'norm|deg|repla'))); % find all norms
            % ’syl'可以兼容旧的stimuli命名规则
            
            for n = 1: length(ids)
                thisi = ids(n);
                
                presdf = Cal.sdf(neu.list(thisi).prejudgerespsptimes,zeros(length(neu.list(thisi).judgerespy),1),neu.list(thisi).fs,0.001,0.02);
                sdf = Cal.sdf(neu.list(thisi).judgerespsptimes,neu.list(thisi).judgerespy,neu.list(thisi).fs,0.001,0.02); % 0.001,0.004
                [maxpresdf,~] = max(presdf);
                [maxsdf,maxidx] = max(sdf);
                
                pre_frs = Cal.eachTrialFiringRate(neu.list(thisi).prejudgerespsptimes,length(neu.list(thisi).judgerespy)/neu.list(thisi).fs);
                sti_frs = Cal.eachTrialFiringRate(neu.list(thisi).judgerespsptimes,length(neu.list(thisi).judgerespy)/neu.list(thisi).fs);
                [h,p] = ttest(sti_frs,pre_frs,'Tail','Right','Alpha',0.05);
                neu.list(thisi).pvalue = p;
                neu.list(thisi).label = 0; % 初始化
                if h == 1
                    neu.list(thisi).label = 1;
                    %                     if (num_of_not_empty_trials/length(neu.list(thisi).pltsptimes)<0.5)||(num_of_not_empty_trials<5) % if not 60% of the trails are not empty
                    %                         neu.list(thisi).label = 0;
                    %                     end
                elseif maxsdf > 17 && maxsdf > maxpresdf % neu rescue
                    neu.list(thisi).label = 1;
                    
                end
            end
            
        end
        
        function reordered_ids = reorderListByResp(neu)
            % reorder List By Response
            fragids = find(~cellfun(@isempty, regexp({neu.list.stimuliname}.','Frag|syl|Syl'))); % find all frags,兼容 syl|Syl
            if ~isempty( fragids)
                fraglist = neu.list(fragids);
                for n = 1: length(fraglist)
                    tempsum = Cal.psth_frag(fraglist(n).rawy,fraglist(n).fs,fraglist(n).rawsptimes);
                    halfsum = sum(tempsum(end/2:end));
                    fullsum = sum(tempsum);
                    maxvalue = max(Cal.psth_frag(fraglist(n).rawy,fraglist(n).fs,fraglist(n).rawsptimes));
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
            
            replaids = find(~cellfun(@isempty, regexp({neu.list.stimuliname}.','Repla|catego'))); % find all frags,兼容catego
            if ~isempty( replaids)
                replalist = neu.list(replaids);
                for n = 1: length(replalist)
                    tempsum = Cal.psth_frag(replalist(n).rawy,replalist(n).fs,replalist(n).rawsptimes);
                    halfsum = sum(tempsum(end/2:end));
                    fullsum = sum(tempsum);
                    maxvalue = max(Cal.psth_frag(replalist(n).rawy,replalist(n).fs,replalist(n).rawsptimes));
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
            
            normids = find(~cellfun(@isempty, regexp({neu.list.stimuliname}.','norm'))); % find all frags
            if ~isempty( normids )
                normlist = neu.list( normids);
                for n = 1: length( normlist)
                    tempsum = Cal.psth_frag( normlist(n).rawy, normlist(n).fs, normlist(n).rawsptimes);
                    halfsum = sum(tempsum(end/2:end));
                    fullsum = sum(tempsum);
                    maxvalue = max(Cal.psth_frag( normlist(n).rawy, normlist(n).fs, normlist(n).rawsptimes));
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
            
            
            
            otherids = find(cellfun(@isempty, regexp({neu.list.stimuliname}.','Frag|Repla|norm|catego|Syl|syl'))); % find all frags
            if ~isempty(otherids)
                otherlist = neu.list(otherids);
                for n = 1: length(otherlist)
                    tempsum = Cal.psth_frag( otherlist(n).rawy,otherlist(n).fs,otherlist(n).rawsptimes);
                    halfsum = sum(tempsum(end/2:end));
                    fullsum = sum(tempsum);
                    maxvalue = max(Cal.psth_frag(otherlist(n).rawy,otherlist(n).fs,otherlist(n).rawsptimes));
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
        
        
    end
    
    methods % 外部计算方法
        
        % 计算部分
        function e_objects = getAllEphysObject(neu)
            %提取Analysis所有的Ephys Object
            collects = {};
            for k = 1: length(neu.neurons)
                collects{k} = neu.neurons{k}.e(:);
            end
            e_objects = vertcat(collects{:});
        end
        
        function formated_imagename = neuronname(neu)
            % 似乎没用
            formated_imagename = sprintf('%s_%u',neu.birdid,neu.uniqueid);
        end
        
    
        function meanWL = calMeanWaveLength(neu)
            % 计算 meanWL，需要考虑到仪器的fs
            wavecollect = {};
            for s = 1: length(neu.neurons)
                wavecollect{s} = neu.neurons{s}.waveform;
            end
            waveforms = vertcat(wavecollect{:});
            %waveforms =  n.waveform;
            [~,troughstime] = min(waveforms,[],2);
            wavlen_units = [];
            
            for k = 1: size(waveforms,1) % length is dangerous!!!!!
                this_wf = waveforms(k,:);
                [~,wavlen_units(k)] =  max(this_wf (troughstime(k):end));
            end
            
            if ~isempty(regexp(neu.zpid,'Z')) %zeus
                neu.plx_data_fs = 30000; %hard code !!!!!! Dangerous
            elseif ~isempty(regexp(neu.zpid,'P')) % plexon
                neu.plx_data_fs = 40000;
            end
            
            meanWL =  mean(wavlen_units*(1/neu.plx_data_fs)*1000); % ms
            
            
        end
        
        function [localSFR,h,p] = getSponFR(neu,range)
            % calculate spontaneous firing rate
            sponFrInfo = struct;
            all_es = neu.getAllEphysObject;
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
        
        
        function sponFrInfo = getSponFR_Sarah(neu,range)
            % calculate spontaneous firing rate
            sponFrInfo = struct;
            all_es = neu.getAllEphysObject;
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
        
        
        function fraglist = to2ndAcousticSpace(neu)
            % 2nd acoustic space？
            ids = find(~cellfun(@isempty, regexp({neu.list.stimuliname}.','repla'))); % find all frags
            
            fraglist = neu.list(ids);
            
            for n = 1: length(fraglist)
                tempsum = Cal.psth_frag(fraglist(n).rawy,fraglist(n).fs,fraglist(n).rawsptimes);
                range = 1 % the very fisrt 1 second
                beginmax = max(tempsum(1: ceil(length(tempsum)*range/(length(fraglist(n).rawy)/fraglist(n).fs)) ));% the maximum value of the begining 0.5 second
                halfsum = sum(tempsum(end/2:end));
                fullsum = sum(tempsum);
                maxvalue = max(Cal.psth_frag(fraglist(n).rawy,fraglist(n).fs,fraglist(n).rawsptimes));
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
        
     
        function insonglist = getInsongFragRespList(neu)
            % get In-songFrag Response List
            cons_neuron = neu.neurons{neu.song_id};
            insonglist = cons_neuron.toList_Insong;
            
        end
        
        function ainf = calPropertiesForClustering(neu)
            % calculate lifetime sparseness,correlation index, number of
            % responsive songs, spontaneous firing rate, spike width
            % 目的是通过计算这些性质的值对neurons进行划分
            
            % number of responsive songs
            
            %             for kk = 1: length( neu.neurons{neu.song_id}.e)
            %                 neu.neurons{neu.song_id}.e{kk}.setExtAndAllocate;
            %             end
            %             neu.updatelist;
            
            
            conkeywords = {'B346','B512','B521','B554','B606','G429','G506','G518','G548','G573',...
                'G578','O331','O507','O509','O540','Y515','Y606','Y616'};
            conids = [];
            for kk = 1: length(conkeywords)
                ids = find(~cellfun(@isempty,regexp(cellstr({neu.list.stimuliname}.'),['norm\S+(?!(TUT|BOS))',conkeywords{kk},'(?!(TUT|BOS))'] )));
                if length(ids) == 1
                    conids(kk) = ids;
                elseif length(ids) >1
                    conids(kk) = ids(1);
                elseif length(ids) == 0
                    conids(kk) = nan;
                end
            end
            conids = rmmissing(conids);
            normids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'norm')));
            
            normlist = neu.list(intersect(conids,normids));
            
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
                all_Ns = length(find(abs(Cal.allPairDiff(all_spikes))<thres));
                same_trail_Ns = [];
                for k = 1: length(jrsptimes)
                    same_trail_Ns(k) = length(find(abs(Cal.allPairDiff(jrsptimes{k}))<thres));
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
                
                sdf = Cal.sdf(jrsptimes,normlist(m).judgerespy,normlist(m).fs,0.001,0.004);
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
                
                presdf = Cal.sdf(normlist(m).prejudgerespsptimes,zeros(length(normlist(m).judgerespy),1),normlist(m).fs,0.001,0.02);
                sdf = Cal.sdf(normlist(m).judgerespsptimes,normlist(m).judgerespy,normlist(m).fs,0.001,0.02); % 0.001,0.004
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
                
                
                
                
                pre_frs = Cal.eachTrialFiringRate(normlist(m).prejudgerespsptimes,length(normlist(m).judgerespy)/normlist(m).fs);
                sti_frs = Cal.eachTrialFiringRate(normlist(m).judgerespsptimes,length(normlist(m).judgerespy)/normlist(m).fs);
                [h,p] = ttest(sti_frs,pre_frs,'Tail','Right','Alpha',0.05)
                normlist(m).pvalue = p;
                normlist(m).label = 0; % 初始化
                if h == 1
                    normlist(m).label = 1;
                elseif maxsdf > 17 && maxsdf > maxpresdf % neu rescue
                    normlist(m).label = 1;
                    
                end
                
                
            end
            
            
            allsdf = horzcat(sdf_collect{:});
            ainf.allsdf = allsdf;
            is1ids = find([normlist.label].' == 1);
            ainf.numrespsong = length(is1ids);
            label1normlist = normlist(is1ids);
            
            ainf.neuronname = neu.formated_imagename;
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
            ainf.meanWL = neu.calMeanWaveLength;
            
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
            %             sumsdf = Cal.sdf(concat_judgerespsptimes,zeros(sum_judgeresplen*normlist(m).fs,1),normlist(m).fs,0.001,0.004);
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
        
   
   end
    
    methods % 作图方法
        
     
        % 作图并导出
        
        function img = exportConRevMir(neu)
            dbstop if error
            reverseids = find(~cellfun(@isempty,regexp(cellstr({neu.list.stimuliname}.'),'reverse')));
            normids = find(~cellfun(@isempty,regexp(cellstr({neu.list.stimuliname}.'),'norm')));
            mirrorids = find(~cellfun(@isempty,regexp(cellstr({neu.list.stimuliname}.'),'mirror')));
            
            % find out the songs that are tested with the reversed version
            bids = cellfun(@(x) Convert.bid(x),cellstr({neu.list(reverseids).stimuliname}.'),'Uni',0);
            
            IMG = {};
            
            
            if isempty(bids)
                figure('Position',[552 -116 2410 1221],'Color','w');;
                %Draw.three(neu.list(sameFile_mirror_ids).plty,neu.list(sameFile_mirror_ids).fs,neu.list(sameFile_mirror_ids).pltsptimes);
                frame = getframe(gcf);
                IMG{1,1} = frame.cdata;
                IMG{2,1} = frame.cdata;
                IMG{3,1} = frame.cdata;
                close(gcf);
            end
            
            
            
            for k = 1:length(bids)
                
                samebirdids = find(~cellfun(@isempty,regexp(cellstr({neu.list.stimuliname}.'),bids{k})));
                sb_reverse_ids = intersect(reverseids,samebirdids); % sab means same bird
                if length(sb_reverse_ids) > 1
                    sb_reverse_ids =   sb_reverse_ids(1); % 如果多项，暂时的对策是取第一项
                end
                sb_mirror_ids = intersect(mirrorids,samebirdids);
                sb_norm_ids = intersect(normids,samebirdids);
                
                samefileids  = find(~cellfun(@isempty,regexp(cellstr({neu.list.Fid}.'),neu.list(sb_reverse_ids).Fid)));
                sameFile_norm_ids = intersect(samefileids,sb_norm_ids);
                sameFile_reverse_ids = intersect(samefileids,sb_reverse_ids);
                sameFile_mirror_ids = intersect(samefileids,sb_mirror_ids);
                
                fig1 = figure('Position',[552 -116 2410 1221],'Color','w');
                Draw.three(neu.list(sameFile_norm_ids).y,neu.list(sameFile_norm_ids).fs,neu.list(sameFile_norm_ids).sptimes); % plty
                frame = getframe(fig1);
                IMG{1,k} = frame.cdata;
                close(fig1)
                
                fig2 = figure('Position',[552 -116 2410 1221],'Color','w');
                Draw.three(neu.list(sameFile_reverse_ids).y,neu.list(sameFile_reverse_ids).fs,neu.list(sameFile_reverse_ids).sptimes);
                frame = getframe(fig2);
                IMG{2,k} = frame.cdata;
                close(fig2)
                
                fig3 = figure('Position',[552 -116 2410 1221],'Color','w');
                Draw.three(neu.list(sameFile_mirror_ids).y,neu.list(sameFile_mirror_ids).fs,neu.list(sameFile_mirror_ids).sptimes);
                frame = getframe(fig3);
                IMG{3,k} = frame.cdata;
                close(fig3);
                
            end
            
            img = cell2mat(IMG);
            
            img = Convert.colorEdge(img,'r');
            % I don't know why here the output could be unit8 or double
            % randomly
            % Now the temporal solution is
            img = uint8(img);
            
        end
       
        function img = saveDrawSortedRespToFrags(neu) %非常难以找到这个function
            %非常难以找到这个function
            fraglist =  judgeFragResponse(neu);
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
                Draw.two(sorted_fraglist(k).plty,sorted_fraglist(k).fs,sorted_fraglist(k).pltsptimes);
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
            imwrite(img,sprintf('晋RespToFrags_%s.png',neu.formated_imagename));
            
        end
        
        % 作图且保存
         function Iall = saveDrawAlignedConsFrag(neu,songnames)
            % songnames 指定后则只会显示对应于songnames的fragments
            dbstop if error
            tic
            
            fragids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'frag|Frag|syl')));
            
            fraglist = neu.list(fragids);
            frag_Fid = unique({fraglist.Fid}.');
            % normlist = Neuron(neu.neurons{neu.song_id}).normlist;
            
            subfile_frag_ids = find(~cellfun(@isempty,regexp(cellstr({neu.list.Fid}.'),strjoin(frag_Fid,'|'))));
            hard_to_name_ids = subfile_frag_ids;
            if exist('songnames','var')
                songnameids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),strjoin(songnames,'|'))));
                hard_to_name_ids = intersect(subfile_frag_ids,songnameids);
            end
            
            fucklist = neu.list(hard_to_name_ids);
            
            normlist = fucklist(find(~cellfun(@isempty, regexp(cellstr({fucklist.stimuliname}.'),'norm'))));
            
            % about Norms % Redudant code
            %             normlist = Neuron(neu.neurons{neu.song_id}).normlist;
            %             [~,postunique] = unique(cellfun(@Convert.bid,cellstr({normlist.stimuliname}.'),'Uni',0))
            %             normlist = normlist(postunique);
            
            
            % About Frag
            fragids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'frag|Frag|syl|Syl') ));
            if ~isempty(fragids)
                
                fraglist = neu.list(fragids);
                
                for m = 1: length(fraglist)
                    
                    birdid = Convert.bid(fraglist(m).stimuliname);
                    ids_norm = find(~cellfun(@isempty, regexp(cellstr({normlist.stimuliname}.'),birdid) ) );
                    
                    if ~isempty(ids_norm)& length(ids_norm) == 1
                        fraglist(m).sylIni = Neuron.findIni(normlist(ids_norm).plty,fraglist(m).y);
                        fraglist(m).pady = [zeros([fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext,1]);fraglist(m).plty;zeros(length(normlist(ids_norm).plty)...
                            - (fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext) - length(fraglist(m).plty),1)];
                        fraglist(m).padsptimes = cellfun( @(x) x + (fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext)/fraglist(m).fs, fraglist(m).pltsptimes,'uni',0);
                    end
                end
            end
            
            %             % About Deg
            %             fragids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'Deg|deg') ));
            %             if ~isempty(fragids)
            %
            %                 deglist = neu.list(fragids);
            %                 for m = 1: length(deglist)
            %                     birdid = Convert.bid(deglist(m).stimuliname);
            %                     ids_norm = find(~cellfun(@isempty, regexp(cellstr({normlist.stimuliname}.'),birdid) ) );
            %                     if ~isempty(ids_norm)& length(ids_norm) == 1
            %                         [deglist(m).sylIni,trump_diffvalue] = Neuron.findIni(normlist(ids_norm).plty,deglist(m).y);
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
                %                     birdid = Convert.bid(normlist(w).stimuliname);
                %                     ids_indeg = find(~cellfun(@isempty, regexp(cellstr({deglist.stimuliname}.'),birdid) ) );
                %                     selected_deglist = deglist(ids_indeg);
                %                     [~,temp_index] = sortrows([selected_deglist.sylIni].');
                %                     selected_deglist = selected_deglist(temp_index);
                %                 end
                
                if ~isempty(fragids)
                    birdid = Convert.bid(normlist(w).stimuliname);
                    ids_infrag = find(~cellfun(@isempty, regexp(cellstr({fraglist.stimuliname}.'),birdid) ) );
                    selected_fraglist = fraglist(ids_infrag);
                    [~,temp_index] = sortrows([selected_fraglist.sylIni].');
                    selected_fraglist = selected_fraglist(temp_index);
                end
                
                
                % draw the basic figure
                Icollect = {};
                figure('Color','w');
                
                Draw.two(normlist(w).plty,normlist(w).fs,neu.normlist(w).pltsptimes);
                xlabel(normlist(w).stimuliname);
                frame = getframe(gcf);
                Icollect{1} = frame.cdata;
                close(gcf)
                
                if ~isempty(fragids)
                    for hh = 1: length(selected_deglist)
                        
                        figure('Color','w');
                        Draw.two(selected_deglist(hh).pady,selected_deglist(hh).fs,selected_deglist(hh).padsptimes);
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
                        Draw.two(selected_fraglist(bb).pady,selected_fraglist(bb).fs,selected_fraglist(bb).padsptimes);
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
                        % draw neu figure for every four ote fragments
                        figure('Color','w');
                        for kk = flip([0:3])
                            if 4*k - kk <length(otelist)
                                subplot(2,4,4-kk);
                                
                                Draw.spec(otelist(4*k - kk).plty,otelist(4*k - kk).fs);
                                subplot(2,4,4-kk + 4 );
                                Draw.raster(otelist(4*k - kk).pltsptimes,otelist(4*k - kk).plty,otelist(4*k - kk).fs);
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
            
            neu.drawSeparatedWaveform;
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
            
            imwrite(Iall,sprintf('Aligned_Normfrag_%s.png',neu.formated_imagename));
            toc
            
        end
        
        function Iall = saveDrawAlignedConsDegs(neu,songnames)
            dbstop if error
            
            %只有一个模式： 只针对二次播放里包含的norm songs进行degs的对齐
            tic
            
            degids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'deg|Deg')));
            
            deglist = neu.list(degids);
            deg_Fid = unique({deglist.Fid}.');
            % normlist = Neuron(neu.neurons{neu.song_id}).normlist;
            
            
            
            subfile_deg_ids = find(~cellfun(@isempty,regexp(cellstr({neu.list.Fid}.'),strjoin(deg_Fid,'|'))));
            hard_to_name_ids = subfile_deg_ids;
            if exist('songnames','var')
                songnameids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),strjoin(songnames,'|'))));
                hard_to_name_ids = intersect(subfile_deg_ids,songnameids);
            end
            
            fucklist = neu.list(hard_to_name_ids);
            
            normlist = fucklist(find(~cellfun(@isempty, regexp(cellstr({fucklist.stimuliname}.'),'norm'))));
            %             [~,postunique] = unique(cellfun(@Convert.bid,cellstr({fucklist.stimuliname}.'),'Uni',0));
            %             normlist = normlist(postunique);
            
            % About Deg
            degids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'Deg|deg') ));
            if ~isempty(degids)
                
                deglist = neu.list(degids);
                for m = 1: length(deglist)
                    birdid = Convert.bid(deglist(m).stimuliname);
                    ids_norm = find(~cellfun(@isempty, regexp(cellstr({normlist.stimuliname}.'),birdid) ) );
                    if ~isempty(ids_norm)& length(ids_norm) == 1
                        [deglist(m).sylIni,trump_diffvalue] = Neuron.findIni(normlist(ids_norm).plty,deglist(m).y);
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
                    birdid = Convert.bid(normlist(w).stimuliname);
                    ids_indeg = find(~cellfun(@isempty, regexp(cellstr({deglist.stimuliname}.'),birdid) ) );
                    selected_deglist = deglist(ids_indeg);
                    [~,temp_index] = sortrows([selected_deglist.sylIni].');
                    selected_deglist = selected_deglist(temp_index);
                end
                
                
                % draw the basic figure
                Icollect = {};
                figure('Color','w','Position',PM.size_wide);
                
                Draw.two(normlist(w).plty,normlist(w).fs,normlist(w).pltsptimes);
                xlabel(normlist(w).stimuliname);
                frame = getframe(gcf);
                Icollect{1} = frame.cdata;
                close(gcf)
                
                
                
                if ~isempty(degids)
                    for hh = 1: length(selected_deglist)
                        figure('Color','w','Position',PM.size_wide);
                        Draw.two(selected_deglist(hh).pady,selected_deglist(hh).fs,selected_deglist(hh).padsptimes);
                        xlabel(selected_deglist(hh).stimuliname);
                        frame = getframe(gcf);
                        Icollect{1 + hh} = frame.cdata;
                        close(gcf);
                    end
                end
                
                frozen_Icollect_len = length(Icollect);
                
                
                I_of_each_column{w} = vertcat(Icollect{:});
            end
            
            %             neu.drawFirstWaveform;
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
            
            % imwrite(Iall,sprintf('Aligned_ConsDegs_%s.png',neu.neurons{1}.neuronname));
            imwrite(Iall,sprintf('Aligned_ConsDegs_%s.png',neu.formated_imagename));
            toc
            
        end
        
        function Iall = saveDrawAlignedConsReplas(neu,songnames) % align replas and draw
            dbstop if error
            tic
            RONGYU = 0.5;
            
            tic
            
            replaids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'repla|Repla')));
            replalist = neu.list(replaids);
            repla_Fid = unique({replalist.Fid}.');
            subfile_repla_ids = find(~cellfun(@isempty,regexp(cellstr({neu.list.Fid}.'),strjoin(repla_Fid,'|'))));
            hard_to_name_ids = subfile_repla_ids;
            if exist('songnames','var')
                songnameids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),strjoin(songnames,'|'))));
                hard_to_name_ids = intersect(subfile_repla_ids,songnameids);
            end
            fucklist = neu.list(hard_to_name_ids);
            normlist = fucklist(find(~cellfun(@isempty, regexp(cellstr({fucklist.stimuliname}.'),'(?!repla|Repla)norm|Norm'))));
            
            
            
            % This is the new version 04.06.2022
            %             normids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'(?!repla|Repla)norm|Norm') ));
            %             normlist = neu.list(normids);
            
            for k = 1: length(normlist)
                normlist(k).pady = [zeros(RONGYU*normlist(k).fs,1);normlist(k).plty];
                normlist(k).padsptimes = cellfun( @(x) x + RONGYU, normlist(k).pltsptimes,'uni',0);
            end
            
            % About Repla
            replaids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'Repla|repla|catego|Catego') ));
            if isempty(replaids); disp('No replas'); return; end
            if ~isempty(replaids)
                fraglist = neu.list(replaids);
                ids_norm_collecting_box = [];
                for m = 1: length(fraglist)
                    
                    afterBefore = regexp(fraglist(m).stimuliname,'(?<=before-)\S*','match');
                    afterBefore = afterBefore{1};
                    birdid = Convert.bid(afterBefore);
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
                    birdid = Convert.bid(selected_normlist(w).stimuliname);
                    ids_infrag = find(~cellfun(@isempty, regexp(cellstr({fraglist.stimuliname}.'),birdid) ) );
                    selected_fraglist = fraglist(ids_infrag);
                end
                
                
                % draw the basic figure
                Icollect = {};
                figure('Color','w','Position',[406 675 1378 420]);
                
                Draw.two(selected_normlist(w).pady,selected_normlist(w).fs,selected_normlist(w).padsptimes);
                xlabel(selected_normlist(w).stimuliname);
                frame = getframe(gcf);
                Icollect{1} = frame.cdata;
                close(gcf)
                
                
                frozen_Icollect_len = length(Icollect);
                
                if ~isempty(replaids)
                    for bb = 1: length(selected_fraglist)
                        
                        figure('Color','w','Position',[406 675 1378 420]);
                        Draw.two(selected_fraglist(bb).pady,selected_fraglist(bb).fs,selected_fraglist(bb).padsptimes);
                        xlabel(selected_fraglist(bb).stimuliname);
                        frame = getframe(gcf);
                        Icollect{frozen_Icollect_len + bb} = frame.cdata;
                        close(gcf);
                        
                    end
                end
                
                I_of_each_column{w} = vertcat(Icollect{:});
            end
            
            neu.drawFirstWaveform;
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
            
            imwrite(Iall,sprintf('Aligned_Repla_With_Norm_%s.png',neu.formated_imagename));
            toc
            
        end
        

        function Whether_NeuResp_To_SinFrags_Consis_Or_affected_By_Previous(neu)
            dbstop if error
            tic
            
            % 首先，定义songlist
            songlist = neu.list(find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'norm'))));
            
            %songlist = songlist(postunique);
            degids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'Deg|deg') ));
            if isempty(degids); return; end
            deglist = neu.list(degids);
            deg_bids = unique(cellfun(@Convert.bid,cellstr({deglist.stimuliname}.'),'Uni',0));
            
            I_song = {};
            for w = 1: length(deg_bids)
                
                
                Icollect = {};
                degexist_norm_list  = songlist(find(~cellfun(@isempty, regexp(cellstr({songlist.stimuliname}.'),deg_bids{w}))));
                
                for k = 1; length(degexist_norm_list)
                    smallfig = figure('Color','w','Position',[406 675 1378 420]);
                    Draw.two(degexist_norm_list(k).plty,degexist_norm_list(k).fs,degexist_norm_list(k).pltsptimes);
                    xlabel(degexist_norm_list(k).stimuliname);
                    frame = getframe(smallfig); Icollect{1} = frame.cdata;close(gcf)
                    
                end
                
                I_song{w} = vertcat(Icollect{:});
                
            end
            [~,postunique] = unique(cellstr(cellfun(@Convert.bid,{songlist.stimuliname}.','Uni',0)));
            songlist = songlist(postunique);
            
            % 其二，找到deressive songs对应的birdid
            degids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'Deg|deg') ));
            if isempty(degids); return; end
            deglist = neu.list(degids);
            for m = 1: length(deglist)
                birdid_collect{m} = Convert.bid(deglist(m).stimuliname);
                ids_norm = find(~cellfun(@isempty, regexp(cellstr({songlist.stimuliname}.'),birdid_collect{m}) ) );
                if ~isempty(ids_norm)& length(ids_norm) == 1
                    [deglist(m).sylIni,trump_diffvalue] = Neuron.findIni(songlist(ids_norm).plty,deglist(m).y);
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
                
                birdid = Convert.bid(normlist(w).stimuliname);
                ids_indeg = find(~cellfun(@isempty, regexp(cellstr({deglist.stimuliname}.'),birdid) ) );
                selected_deglist = deglist(ids_indeg);
                [~,temp_index] = sortrows([selected_deglist.sylIni].');
                selected_deglist = selected_deglist(temp_index);
                
                % draw the norm figure
                Icollect = {};figure('Color','w','Position',[406 675 1378 420]);
                Draw.two(normlist(w).plty,normlist(w).fs,normlist(w).pltsptimes);
                xlabel(normlist(w).stimuliname);
                frame = getframe(gcf); Icollect{1} = frame.cdata;close(gcf)
                
                % draw the deg figure
                for hh = 1: length(selected_deglist)
                    figure('Color','w','Position',[406 675 1378 420]);
                    Draw.two(selected_deglist(hh).pady,selected_deglist(hh).fs,selected_deglist(hh).padsptimes);
                    xlabel(selected_deglist(hh).stimuliname);
                    frame = getframe(gcf); Icollect{1 + hh} = frame.cdata;close(gcf);
                end
                
                I_Deg{w} = vertcat(Icollect{:});
            end
            
            % 其三，找到对应birdid的frags
            fragids1 = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'frag|Frag|syl|Syl') ));
            unique_bid = unique(cellstr(birdid_collect));
            fragids2 = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),strjoin(unique_bid,'|')) ));
            fragids = intersect(fragids1,fragids2);
            if ~isempty(fragids)
                fraglist = neu.list(fragids);
                for m = 1: length(fraglist)
                    
                    birdid = Convert.bid(fraglist(m).stimuliname);
                    ids_norm = find(~cellfun(@isempty, regexp(cellstr({songlist.stimuliname}.'),birdid) ) );
                    
                    if ~isempty(ids_norm)& length(ids_norm) == 1
                        fraglist(m).sylIni = Neuron.findIni(songlist(ids_norm).plty,fraglist(m).y);
                        fraglist(m).pady = [zeros([fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext,1]);fraglist(m).plty;zeros(length(songlist(ids_norm).plty)...
                            - (fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext) - length(fraglist(m).plty),1)];
                        fraglist(m).padsptimes = cellfun( @(x) x + (fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext)/fraglist(m).fs, fraglist(m).pltsptimes,'uni',0);
                    end
                end
            end
            
            I_Frag = {};
            for w = 1: length(normlist)
                
                if ~isempty(fragids)
                    birdid = Convert.bid(normlist(w).stimuliname);
                    ids_infrag = find(~cellfun(@isempty, regexp(cellstr({fraglist.stimuliname}.'),birdid) ) );
                    selected_fraglist = fraglist(ids_infrag);
                    [~,temp_index] = sortrows([selected_fraglist.sylIni].');
                    selected_fraglist = selected_fraglist(temp_index);
                end
                
                Icollect = {};figure('Color','w','Position',[406 675 1378 420]);
                Draw.two(normlist(w).plty,songlist(w).fs,normlist(w).pltsptimes);  xlabel(normlist(w).stimuliname);
                frame = getframe(gcf);   Icollect{1} = frame.cdata; close(gcf)
                
                if ~isempty(fragids)
                    for bb = 1: length(selected_fraglist)
                        figure('Color','w','Position',[406 675 1378 420]);
                        Draw.two(selected_fraglist(bb).pady,selected_fraglist(bb).fs,selected_fraglist(bb).padsptimes);
                        xlabel(selected_fraglist(bb).stimuliname);
                        frame = getframe(gcf);Icollect{1 + bb} = frame.cdata;close(gcf);
                    end
                end
                I_Frag{w} = vertcat(Icollect{:});
            end
            
            
            
            % 其四，找到对应birdid的replas
            RONGYU = 0.5;
            normids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'(?!repla|Repla)norm|Norm') ));
            normlist = neu.list(normids);
            
            for k = 1: length(normlist)
                normlist(k).pady = [zeros(RONGYU*normlist(k).fs,1);normlist(k).plty];
                normlist(k).padsptimes = cellfun( @(x) x + RONGYU, normlist(k).pltsptimes,'uni',0);
            end
            
            replaids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'Repla|repla|catego|Catego') ));
            if ~isempty(replaids)
                
                replalist = neu.list(replaids);
                ids_norm_collecting_box = [];
                for m = 1: length(replalist)
                    
                    afterBefore = regexp(replalist(m).stimuliname,'(?<=before-)\S*','match');
                    afterBefore = afterBefore{1};
                    birdid = Convert.bid(afterBefore);
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
                    birdid = Convert.bid(selected_normlist(w).stimuliname);
                    ids_inrepla = find(~cellfun(@isempty, regexp(cellstr({replalist.stimuliname}.'),birdid) ) );
                    selected_replalist = replalist(ids_inrepla);
                end
                
                
                % draw the basic figure
                Icollect = {}; figure('Color','w','Position',[406 675 1378 420]);
                Draw.two(selected_normlist(w).pady,selected_normlist(w).fs,selected_normlist(w).padsptimes);
                xlabel(selected_normlist(w).stimuliname);
                frame = getframe(gcf);  Icollect{1} = frame.cdata;  close(gcf);
                
                if ~isempty(replaids)
                    for bb = 1: length(selected_replalist)
                        figure('Color','w','Position',[406 675 1378 420]);
                        Draw.two(selected_replalist(bb).pady,selected_replalist(bb).fs,selected_replalist(bb).padsptimes);
                        xlabel(selected_replalist(bb).stimuliname);
                        frame = getframe(gcf); Icollect{1 + bb} = frame.cdata; close(gcf);
                    end
                end
                I_Repla{w} = vertcat(Icollect{:});
            end
            
            
            % 其末 draw and save
            neu.drawFirstWaveform;
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
            
            imwrite(Iall,sprintf('Whether_NeuResp_To_SinFrags_Coms_Or_FragsResps_affected_By_Pres_%s.png',neu.formated_imagename));
            toc
            
            
        end
        
        function How_Do_Neurons_Differentiate_Sibling_Songs(neu)
            dbstop if error
            tic
            
            % 首先，定义songlist
            songlist = neu.list(find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'norm')))); % find all conspecific songs
            
            all_degids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'Deg|deg') ));
            all_fragids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'frag|Frag|syl|Syl') ));
            if isempty(all_degids); return; end % If no Degs at all, no need to draw
            bname_degAttached = unique(cellfun(@Convert.bid,cellstr({neu.list(all_degids).stimuliname}.'),'Uni',0));
            I_of_each_birdname = {};
            
            for w = 1: length(bname_degAttached)
                
                tempCollect = {};
                sub_normlist_degAttached  = songlist(... % 有对应degressive song 存在的普通所有norm songs
                    find(~cellfun(@isempty, regexp(cellstr({songlist.stimuliname}.'),bname_degAttached{w}))));
                
                for k = 1: length(sub_normlist_degAttached)
                    tempfig = figure('Color','w','Position',[7 347 2960 714]);
                    Draw.two(sub_normlist_degAttached(k).plty,sub_normlist_degAttached(k).fs,sub_normlist_degAttached(k).pltsptimes);
                    xlabel(sub_normlist_degAttached(k).stimuliname);
                    frame = getframe(tempfig); tempCollect{1} = frame.cdata;close(gcf)
                end
                Inorm = vertcat(tempCollect{:});
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                deglist = neu.list(intersect(all_degids,...
                    find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),bname_degAttached{w}))))); % conatin only and all deg stimuli
                ids_norm = find(... % 假设不会存在一首song的Degs在不同plx文件里出现
                    strcmp(deglist(1).Fid,{sub_normlist_degAttached.Fid}.')); % To find the normsong which (1) same birdid (2) same Fid
                if isempty(ids_norm); ids_norm = 1; disp('Warning!!!! Norm absent'); end
                for m = 1: length(deglist)
                    %if no normsong with same Fid, then use the first normsong instead
                    if ~isempty(ids_norm)&& length(ids_norm) == 1 % Based on ids_norm , pad Zeros
                        [deglist(m).sylIni,diffvalue] = Neuron.findIni(sub_normlist_degAttached(ids_norm).plty,deglist(m).y); % reference is the plty
                        fprintf('%s has the DIFFVALUE as: %u and INI as: %f \n ',deglist(m).stimuliname,diffvalue,deglist(m).sylIni/32000);
                        deglist(m).pady = [zeros([deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext,1]);deglist(m).plty;zeros(length(sub_normlist_degAttached(ids_norm).plty)...
                            - (deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext) - length(deglist(m).plty),1)];
                        deglist(m).padsptimes = cellfun( @(x) x + (deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext)/deglist(m).fs, deglist(m).pltsptimes,'uni',0);
                    end
                end
                
                % draw the norm figure
                tempCollect = {};figure('Color','w','Position',[7 347 2960 714]);
                
                Draw.two(sub_normlist_degAttached(ids_norm).plty,...
                    sub_normlist_degAttached(ids_norm).fs,sub_normlist_degAttached(ids_norm).pltsptimes);
                xlabel(sub_normlist_degAttached(ids_norm).stimuliname);
                frame = getframe(gcf); tempCollect{1} = frame.cdata;close(gcf);
                % draw the deg figure
                for hh = 1: length(deglist)
                    figure('Color','w','Position',[7 347 2960 714]);
                    Draw.two(deglist(hh).pady,deglist(hh).fs,deglist(hh).padsptimes);
                    xlabel(deglist(hh).stimuliname);
                    frame = getframe(gcf); tempCollect{1 + hh} = frame.cdata;close(gcf);
                end
                I_Deg = vertcat(tempCollect{:});
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fraglist = neu.list(intersect(all_fragids,...
                    find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),bname_degAttached{w}))))); % conatin only and all deg stimuli
                ids_norm = find(... % 假设不会存在一首song的Degs在不同plx文件里出现
                    strcmp(fraglist(1).Fid,{sub_normlist_degAttached.Fid}.')); % To find the normsong which (1) same birdid (2) same Fid
                if isempty(ids_norm); ids_norm = 1; disp('Warning!!!! Norm absent'); end
                for m = 1: length(fraglist)
                    %birdid_collect{m} = Convert.bid(deglist(m).stimuliname);
                    %if no normsong with same Fid, then use the first normsong instead
                    if ~isempty(ids_norm)& length(ids_norm) == 1
                        fraglist(m).sylIni = Neuron.findIni(sub_normlist_degAttached(ids_norm).plty,fraglist(m).y);
                        fraglist(m).pady = [zeros([fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext,1]);fraglist(m).plty;zeros(length(sub_normlist_degAttached(ids_norm).plty)...
                            - (fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext) - length(fraglist(m).plty),1)];
                        fraglist(m).padsptimes = cellfun( @(x) x + (fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext)/fraglist(m).fs, fraglist(m).pltsptimes,'uni',0);
                    end
                end
                
                % draw the norm figure
                tempCollect = {};figure('Color','w','Position',[7 347 2960 714]);
                
                Draw.two(sub_normlist_degAttached(ids_norm).plty,...
                    sub_normlist_degAttached(ids_norm).fs,sub_normlist_degAttached(ids_norm).pltsptimes);
                xlabel(sub_normlist_degAttached(ids_norm).stimuliname);
                frame = getframe(gcf); tempCollect{1} = frame.cdata;close(gcf);
                % draw the deg figure
                for hh = 1: length(fraglist)
                    figure('Color','w','Position',[7 347 2960 714]);
                    Draw.two(fraglist(hh).pady,fraglist(hh).fs,fraglist(hh).padsptimes);
                    xlabel(fraglist(hh).stimuliname);
                    frame = getframe(gcf); tempCollect{1 + hh} = frame.cdata;close(gcf);
                end
                I_Frag = vertcat(tempCollect{:});
                
                I_of_each_birdname{w} = {Inorm,I_Deg,I_Frag};
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            end
            
            % 其末 draw and save
            neu.drawFirstWaveform;
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
            
            imwrite(Iall,sprintf('Sibling_Songs_%s.png',neu.formated_imagename));
            toc
            
        end
        
    end
    
    methods(Static) % 静态方法
        
        function check_detail_folder(dirpath,global_eleinf)
            % This function draw distribution of song elements in spectral
            % space ( scatter / small spectrogram-rasterPlot)
            
            % Extract .wav file names in the target folder
            filenames = Extract.filename(dirpath,'*.wav');
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
            
            
            % Extract the target song-element from the name of the target
            % directory
            temp = split(dirpath,'\');
            temp = temp{end};
            temp = split(temp,'-');
            targetname = sprintf('%s-%s-%s',temp{2},temp{3},temp{4});
            [~,targetid] = ismember( targetname,global_merge);
            
            % Extract element data of song element Fragments/Replacements
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
            
            
            % Extract data of replaced song
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
            
            
        end
        
        
        function [answer,pvalue] = UseTtestToJudegeRespOrNot(y,sptimes,prey,presptimes,fs)
            dbstop if error
            presdf = Cal.sdf(presptimes,prey,fs,0.001,0.02);
            sdf = Cal.sdf(sptimes,y,fs,0.001,0.02); % 0.001,0.004
%             [maxpresdf,~] = max(presdf);
%             [maxsdf,maxidx] = max(sdf);
            
            pre_frs = Cal.eachTrialFiringRate(presptimes,length(prey)/fs);
            sti_frs = Cal.eachTrialFiringRate(sptimes,length(y)/fs);
            [h,p] = ttest(sti_frs,pre_frs,'Tail','Right','Alpha',0.05)
            if h == 1
                answer = 1;
            elseif h == 0||isnan(h)
                answer = 0;
            end
            pvalue = p;
            
        end
        
        
         function [answer,pvalue] = UseTtestToJudegeRespOrNot_Temp(y,sptimes,prey,presptimes,fs)
            dbstop if error
            presdf = Cal.sdf(presptimes,prey,fs,0.001,0.02);
            sdf = Cal.sdf(sptimes,y,fs,0.001,0.02); % 0.001,0.004
%             [maxpresdf,~] = max(presdf);
%             [maxsdf,maxidx] = max(sdf);
            
            pre_frs = Cal.eachTrialFiringRate(presptimes,length(prey)/fs);
            sti_frs = Cal.eachTrialFiringRate(sptimes,length(y)/fs);
            [h,p] = ttest(sti_frs,pre_frs,'Tail','Right','Alpha',0.05)
            if h == 1
                answer = 1;
            elseif h == 0||isnan(h)
                answer = 0;
            end
            pvalue = p;
            
         end
        
         
    end
    
    methods(Hidden = true) % 弃用方法
        
        function Deprecated_Batch2(neu)
            %neu.rawthree_NeuronClassVersion;
            %neu.sort_frags_by_response_strength_and_then_draw;
            neu.drawMeanFeaturesVsRespAsLineChart;
            %neu.drawAlignedNormDegsTwoPlots;
            %neu.AlignReplasWithNormsThenDraw;
            neu.rawthree_NeuronClassVersion;
            neu.sort_frags_by_response_strength_and_then_draw;
            neu.drawMeanFeaturesVsRespAsLineChart;
            neu.drawAlignedNormDegsTwoPlots;
            neu.AlignReplasWithNormsThenDraw;
            neu.drawPairwiseFragmentsMeanFeaturesDistribution;
        end
        
        function fraglist = Deprecated_judgeFragResponse(neu) %%% To judge whether the neuron response to neu frag or not
            
            ids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'Frag|syl'))); % find all frags
            
            % ’syl'可以兼容旧的stimuli命名规则
            DUR = 0.2 ;% 200ms
            
            fraglist = neu.list(ids);
            
            for n = 1: length(fraglist)
                
                post_sptimes = Extract.sptimes(fraglist(n).rawsptimes, fraglist(n).zpt, fraglist(n).zpt + DUR);
                pre_sptimes = Extract.sptimes(fraglist(n).rawsptimes, fraglist(n).zpt-DUR, fraglist(n).zpt );
                post_mfr = length(vertcat(post_sptimes{:}))/DUR;
                pre_mfr = length(vertcat(pre_sptimes{:}))/DUR; % mfr: mean firing rate
                fraglist(n).rs = post_mfr - pre_mfr; % response strength
                
                
                tempsum = Cal.psth_frag(fraglist(n).plty,fraglist(n).fs,fraglist(n).pltsptimes);
                halfsum = sum(tempsum(end/2:end));
                fullsum = sum(tempsum);
                maxvalue = max(Cal.psth_frag(fraglist(n).plty,fraglist(n).fs,fraglist(n).pltsptimes));
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
        
        function Conlist = Deprecated_evaluateConResponse(neu) % Eveluate Conspecific song response
            
            ids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'norm|deg|repla'))); % find all frags
            
            % ’syl'可以兼容旧的stimuli命名规则
            
            Conlist = neu.list(ids);
            
            for n = 1: length(Conlist)
                sdf = Cal.sdf(Conlist(n).pltsptimes,Conlist(n).plty,Conlist(n).fs,0.001,0.004);
                [maxsdf,maxidx] = max(sdf);
                
                percentage_max = maxidx/length(sdf);
                time_max = length(Conlist(n).plty)/Conlist(n).fs*percentage_max;
                % check whether the surroding are has spikes in most of the
                % trials
                extracted_sptimes = Extract.sptimes(Conlist(n).pltsptimes,time_max - 0.15, time_max + 0.15); % 前后 100ms
                
                num_of_not_empty_trials = length(find(~cellfun(@isempty, extracted_sptimes)));
                
                mean_maxsdf = maxsdf/length(Conlist(n).pltsptimes);
                tempsum = Cal.psth_frag(Conlist(n).plty,Conlist(n).fs,Conlist(n).pltsptimes);
                halfsum = sum(tempsum(end/2:end));
                fullsum = sum(tempsum);
                maxvalue = max(Cal.sdf(Conlist(n).pltsptimes,Conlist(n).plty,Conlist(n).fs)); % I changed the value here from Cal.psth_frag to Cal.sdf
                Conlist(n).sdf = Cal.sdf(Conlist(n).pltsptimes,Conlist(n).plty,Conlist(n).fs);
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
                slim_extracted_sptimes = Extract.sptimes(Conlist(n).pltsptimes,time_max - 0.8, time_max + 0.8);
                slim_trials = length(find(~cellfun(@isempty,slim_extracted_sptimes)));
                if slim_trials/length(Conlist(n).pltsptimes)> 0.5 % if not 60% of the trails are not empty
                    Conlist(n).label = 1;
                end
                
                
                if (maxsdf) > 19.5&& strcmp(Conlist(n).stimuliname,'norm-Y515A-21Pulses') && strcmp(Conlist(n).plxname,'Y661_Z17')
                    disp('Incredible bug in pltsptimes of function Neuron.evaluateConResponse !!');
                    Conlist(n).label = 1;
                end
                
            end
            
        end
        
        function Deprecated_Neurons_Respond_To_Single_Frags_Or_Combinations(neu)
            dbstop if error
            tic
            
            % About Deg
            degids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'Deg|deg') ));
            if isempty(degids);return;end
            
            deglist = neu.list(degids);
            birdids = {};
            for m = 1: length(deglist)
                birdids{m} = Convert.bid(deglist(m).stimuliname);
                ids_norm = find(~cellfun(@isempty, regexp(cellstr({normlist.stimuliname}.'),birdid) ) );
                if ~isempty(ids_norm)&& length(ids_norm) == 1
                    [deglist(m).sylIni,trump_diffvalue] = Neuron.findIni(normlist(ids_norm).plty,deglist(m).y);
                    fprintf('%s has the DIFFVALUE as: %u and INI as: %f \n ',deglist(m).stimuliname,trump_diffvalue,deglist(m).sylIni/32000);
                    deglist(m).pady = [zeros([deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext,1]);deglist(m).plty;zeros(length(normlist(ids_norm).plty)...
                        - (deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext) - length(deglist(m).plty),1)];
                    deglist(m).padsptimes = cellfun( @(x) x + (deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext)/deglist(m).fs, deglist(m).pltsptimes,'uni',0);
                end
            end
            
            
            % about Norms % Redudant code
            normids = intersect(find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'norm') )),...
                find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),birdid) )));
            
            normlist = neu.list(normids);
            
            norm_bids = {};
            for k = 1:length(normlist)
                norm_bids{k} = Convert.bid(normlist(k).stimuliname);
            end
            unique_bids = unique(cellstr(norm_bids));
            
            % merge the new fraglist and the deglist with the normlist
            I_of_each_column = {};
            for w = 1: length(unique_bids)
                
                if ~isempty(degids)
                    bid_tosearch = Convert.bid(unique_bids{w});
                    ids_indeg = find(~cellfun(@isempty, regexp(cellstr({deglist.stimuliname}.'),bid_tosearch) ) );
                    selected_deglist = deglist(ids_indeg);
                    [~,temp_index] = sortrows([selected_deglist.sylIni].');
                    selected_deglist = selected_deglist(temp_index);
                end
                
                
                % draw the basic figure
                Icollect = {};
                figure('Color','w','Position',[406 675 1378 420]);
                
                Draw.two(normlist(w).plty,normlist(w).fs,normlist(w).pltsptimes);
                xlabel(normlist(w).stimuliname);
                frame = getframe(gcf);
                Icollect{1} = frame.cdata;
                close(gcf)
                
                
                
                if ~isempty(degids)
                    for hh = 1: length(selected_deglist)
                        figure('Color','w','Position',[406 675 1378 420]);
                        Draw.two(selected_deglist(hh).pady,selected_deglist(hh).fs,selected_deglist(hh).padsptimes);
                        xlabel(selected_deglist(hh).stimuliname);
                        frame = getframe(gcf);
                        Icollect{1 + hh} = frame.cdata;
                        close(gcf);
                    end
                end
                
                frozen_Icollect_len = length(Icollect);
                I_of_each_column{w} = vertcat(Icollect{:});
            end
            
            neu.drawFirstWaveform;
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
            
            % imwrite(Iall,sprintf('Aligned_ConsDegs_%s.png',neu.neurons{1}.neuronname));
            imwrite(Iall,sprintf('Aligned_ConsDegs_%s.png',neu.formated_imagename));
            toc
        end
        
        function Deprecated_drawfrag(neu,keyword,rangeratio,ids)
            
            mergedeleinf = neu.all_eleinf;
            d = Display(neu.list);
            
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
                for k = 1: length(neu.fragnames)
                    d.showfrag(neu.fragnames{k},mergedeleinf);
                end
            end
            
        end
        
        function Deprecated_ThreeWithPitch(neu)
            % draw three plots
            e_objects = getAllEphysObject(neu);
            
            for idx = 1: length(e_objects)
                e_objects{idx}.threeWithFreqRelatedFeatures; % newer version of threeplot drawing method
                frame = getframe(gcf);
                I{idx} = frame.cdata;
                close(gcf);
            end
            
            figure('Color','w','Position',PM.size3plots);
            neu.neurons{neu.song_id}.draw_waveform;     % draw waveform
            % neu.draw_waveform; % this is the original draw function
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
            imwrite(IMG,sprintf('Three_%s.png',neu.formated_imagename));
            
        end
        
        function Deprecated_saveDraw_replas_only(neu)
            % function to draw replas, but not aligned with cons
            d = Display(neu.list);
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
                Draw.three(repla(idx).f0y,repla(idx).fs,repla(idx).f0sptimes); % newer version of threeplot drawing method
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
            
            imwrite(IMG,sprintf('ReplaThree-%s.png',neu.formated_imagename));
            
            %deg = Display.descend(deg);
            
            %replalist = d.info(find(~cellfun(@isempty, regexp({d.info.stimuliname}.','Repla|catego'))));
            
        end
        
        function Deprecated_saveDraw_deg(neu)  % 这个不太对,应该说非常不对
            
            dbstop if error
            
            d = Display(neu.list);
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
                Draw.three(deg(idx).f0y,deg(idx).fs,deg(idx).f0sptimes); % newer version of threeplot drawing method
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
            
            imwrite(IMG,sprintf('DegThree-%s.png',neu.formated_imagename));
            
            %deg = Display.descend(deg);
            
            %replalist = d.info(find(~cellfun(@isempty, regexp({d.info.stimuliname}.','Repla|catego'))));
            
        end
        
        function Deprecated_saveDrawAlignFragsWithSong(neu)
            dbstop if error
            
            normlist = neu.normlist;
            
            ids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'Frag'))); % find all frags
            fraglist = neu.list(ids);
            
            
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
                    fraglist(this_fragids(f)).addlagsptimes = cellfun(@(x) x+timelag/32000-fraglist(this_fragids(f)).pltext, fraglist(this_fragids(f)).pltsptimes,'un',0); % every sptimes will be added neu timelag time
                    fraglist(this_fragids(f)).timelag = timelag;
                end
                
            end
            
            
            % draw Two plots
            
            I = {};
            for k = 1: length(normlist)
                
                figure('Position',[282 759 1444 272],'Color','w');
                Draw.two(normlist(k).plty,normlist(k).fs,normlist(k).pltsptimes);
                frame = getframe(gcf);
                jihe{1} = frame.cdata;
                close(gcf);
                
                for v = 1:length(normlist(k).fragids)
                    localids = normlist(k).fragids(v);
                    figure('Position',[282 759 1444 272],'Color','w');
                    Draw.two(fraglist(localids ).addlagy,fraglist(localids ).fs,fraglist(localids ).addlagsptimes);
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



