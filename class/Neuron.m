classdef Neuron < handle & Neuron_basicDrawings & Neuron_CDF
    % The class to do all the analysis for neu single neuron
    %   Detailed explanation goes here

    properties % 践行 少用继承 多用组合

        experiments % raw features
        song
        repla
        deg % class to process degressive song
        frag
        simsong % class to process similar song
        simatrix % A class store frag-frag distanced measured for different features
        waveform
    end

    properties
        tags % 标签系统，判断神经元是NS还是BS， 是否反应，是否进行了所用测试，是否有 data mismatch云云,被并入了info里
        list
        info % contain all the necessary information
        extra % extra % 收集了一些后期可能会删掉的properties
        all_eleinf
        targets % unique(targetnames)
        insonglist
        fig_size1 % size of figure ( and also position)
        adfreq % analogous to digital frequency, 也即是 samplaing rate of the machine to capture SPKC data
        sectionsMFR % MeanFiringRate of each individual section
        figdata
    end

    methods % 核心方法
 
        function neu = Neuron(experiments,eleinf) % 暂时性的加上这一项
            % 构造方法
            if exist('experiments','var')

                if isa(experiments,'Experiment')
                    neu.experiments{1} = experiments;
                else
                    neu.experiments = experiments;
                end % 输入可以是单个Neuron或者多个neuron组成的cell
                neu.experiments = neu.experiments(~cellfun(@isempty,neu.experiments)); % 去除为空的experiments

                % 构造下属类 Song，Deg，Repla，Frag
                neu.getInfo_pre;
                neu.updatelist; % 生成 stimuli-response list

                %为list加上feature space上的坐标信息，这个feature space是由 VRAE（by
                %Dongqi）生成的
                if exist('eleinf','var')
                    for k = 1:length(neu.list)
                        hitbid = find(~cellfun(@isempty, regexp(cellstr({eleinf.bid}.'),Convert.bid(neu.list(k).stimuliname) )  ));
                        hitfragid = find( [eleinf.fragid].'== Convert.fragid(neu.list(k).stimuliname) );
                        hitrow = intersect(hitbid,hitfragid);
                        neu.list(k).coor = [eleinf(hitrow).coor_1,eleinf(hitrow).coor_2];
                    end
                end

%                 for k = 1:length(eleinf)
%                     eleinf(k).bid = Convert.bid(eleinf(k).songname);
%                 end

                % judge neuron's response
                neu.judgeConResp;
                neu.judgeFragResp_FR;

                neu.song = Song(neu.list); %如果大于19的实验数为零或者大于1，姑且让它报错
                neu.song.formated_name = neu.info.formated_name;


                neu.deg = Deg(neu.list);
                neu.deg.formated_name = neu.info.formated_name;
                neu.frag = Frag(neu.list);
                neu.frag.formated_name = neu.info.formated_name;
                neu.targets = neu.knowTarget;  %知道谁是target syllable
                neu.frag.knowCloseAndFarFrags(neu.targets);
                neu.simsong = Simsong(neu.list);
                neu.simsong.formated_name = neu.info.formated_name;
                latency = neu.calLatency;  % 计算 latency
                neu.repla = Repla(neu.list,latency.latency_halfweight);
                neu.repla.formated_name = neu.info.formated_name;

               
                 % 计算 tags
                neu.autoTag;

               
                %计算 multipleRepeats Meanfeatures
                neu.getInfo_post;
                neu.info = catstruct(neu.info,neu.getInfo_more);
                neu.getExtraInfo;
                neu.setStimuliCorrespondingNeuronId;
                neu.getInfo_last;
                neu.info.tags = neu.tags;

                % 生成experiments
                if length(neu.experiments)>1
                    input = {};
                    for k = 1:length(neu.experiments)
                        input{k} = neu.experiments{k}.waves.waveform;
                    end
                    neu.waveform = SpikeWaveform(input,'multiple');

                else
                    input = struct;
                    input.waveforms = neu.experiments{1}.waves.waveform;
                    input.times = neu.experiments{1}.waves.times;
                    input.timedurations = neu.experiments{1}.waves.times;

                    neu.waveform = SpikeWaveform(input,'single or merge');


                end
                

            end

        end

        function ZanShiDe_genExperiments(neu)

            if length(neu.experiments)>1
                input = {};
                for k = 1:length(neu.experiments)
                    input{k} = neu.experiments{k}.waves.waveform;
                end
                neu.waveform = SpikeWaveform(input,'multiple');

            else
                input = struct;
                input.waveforms = neu.experiments{1}.waves.waveform;
                input.times = neu.experiments{1}.waves.times;
                input.timedurations = neu.experiments{1}.timeSections;

                neu.waveform = SpikeWaveform(input,'single or merge');


            end

        end

        function neu = updatelist(neu)
            % firstly generate or regenerate A.list % 首次或者再次生成 list
            %to_remove_id = intersect(neu.song_only_id,setdiff(neu.song_id_redundant,neu.song_id));%???
            to_remove_id = [];
            to_calculate = setdiff(1: length(neu.experiments),to_remove_id);
            whether_update_figure_or_not = 1;
            for k = 1: length(to_calculate)

                templist = neu.experiments{to_calculate(k)}.toList(whether_update_figure_or_not);
                [templist.whichNeuron] = deal(k);
                lists{k} = templist;
            end
            neu.list = horzcat(lists{:});
        end

        function neu = autoTag(neu)
            % 自动注册标签系统
            % edited from the previous function named Deprecated_getNeuronInfo

            if isempty(neu.tags)
                neu.tags = {};
            end

            if ~isempty(find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'norm-'))))
                neu.tags = [neu.tags,'song'];
            end


            if ~isempty(find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'frag|Frag'))))
                neu.tags = [neu.tags,'frag'];
            end

            if ~isempty(find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'deg|Deg'))))
                neu.tags = [neu.tags,'deg'];
            end

            if ~isempty(find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'repla|Repla|catego|Catego'))))
                neu.tags = [neu.tags,'repla'];
            end

            if ~isempty(find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'R707|R705|R704|R703'))))
                neu.tags = [neu.tags,'sibling'];
            end


            % Whether the neuron is merged, single file, or several subfiles
            % check whether the respon is consistent

            if length(neu.experiments) > 1
                neu.tags = [neu.tags,'multiple'];

                neu.sectionsMFR = []; % 初始化 % mean firing rate
                for k = 1:length(neu.experiments)                 % 这里不一定对
                    try
                        neu.sectionsMFR(k) = length(neu.experiments{k}.sameChannelSpikes.timestamp)...
                            /(sum(neu.experiments{k}.nSections)/neu.experiments{k}.adfreq); % might not be correct !!!!!!! 这个有可能拯救shifted data
                        % 2022.12.13
                    catch

                        neu.sectionsMFR(k) = 0; % 非常差的权宜之计

                    end
                end

                if max(abs(diff(neu.sectionsMFR)))/max(neu.sectionsMFR) > 0.5 % max_rate_diff/max_firing_rate
                    neu.tags = [neu.tags,'posssible_signal_loss']  ;
                end
                % figure; histogram(A.experiments{1}.spikes.timestamp,1000)
            elseif length(neu.experiments{1}.timeSections) > 1

                neu.sectionsMFR = []; % 初始化
                neu.tags = [neu.tags,'merged'];

                for k = 1:length(neu.experiments{1}.timeSections)
                    subdur = neu.experiments{1}.nSections(k)/neu.experiments{1}.adfreq;
                    subspikenum = length(find(neu.experiments{1}.inputs.spikes.timestamp>=neu.experiments{1}.timeSections(k) &...
                        neu.experiments{1}.inputs.spikes.timestamp<neu.experiments{1}.timeSections(k) + subdur ));

                    neu.sectionsMFR(k) = subspikenum/subdur; % total mean fing rate of the section;

                end

                if max(abs(diff(neu.sectionsMFR)))/max(neu.sectionsMFR) > 0.5 % max_rate_diff/max_firing_rate
                    neu.tags = [neu.tags,'posssible_signal_loss']  ;
                end


            elseif length(neu.experiments{1}.timeSections) == 1
                neu.tags = [neu.tags,'single'];
            end


            if ~isempty( intersect(find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'frag|Frag'))),...
                    find([neu.list.label].' == 1)  ))
                neu.tags = [neu.tags,'resp2frag']; % means the neuron did response to at least one of the frag stimuli
            end

            if ~isempty( intersect(find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'deg|Deg'))),...
                    find([neu.list.label].' == 1)  ))
                neu.tags = [neu.tags,'resp2deg']; % means the neuron did response to at least one of the deg stimuli
            end

            if ~isempty( intersect(find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'repla|Repla'))),...
                    find([neu.list.label].' == 1)  ))
                neu.tags = [neu.tags,'resp2repla']; % means the neuron did response to at least one of the repla stimuli
            end

            if ~isempty( intersect(find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'norm|Norm'))),...
                    find([neu.list.label].' == 1)  ))
                neu.tags = [neu.tags,'resp2song']; % means the neuron did response to at least one of the norm(song) stimuli
            end




        end

        function neu = getExtraInfo(neu)
            %  To generate info of each member of an Neuron object

            neu.extra.expinfo = struct;

            for k = 1: length(neu.experiments)
                % neu.extra.expinfo(k).uniqueid = neu.experiments{k}.uniqueid;
                neu.extra.expinfo(k).neuronname = neu.experiments{k}.info.neuronname;
                neu.extra.expinfo(k).keywords = [];

                if length(find(~cellfun(@isempty,regexp(cellstr({neu.experiments{k}.slist.name}.'),'norm|Norm|song|Song')))) > 16
                    neu.extra.expinfo(k).keywords = [neu.extra.expinfo(k).keywords,"song"];
                end

                if ~isempty(find(~cellfun(@isempty,regexp(cellstr({neu.experiments{k}.slist.name}.'),'syl|Syl|Ele|ele|frag|Frag'))))
                    neu.extra.expinfo(k).keywords = [neu.extra.expinfo(k).keywords,"frag"];
                end

                if ~isempty(find(~cellfun(@isempty,regexp(cellstr({neu.experiments{k}.slist.name}.'),'deg|Deg'))))
                    neu.extra.expinfo(k).keywords = [neu.extra.expinfo(k).keywords,"deg"];
                end

                if ~isempty(find(~cellfun(@isempty,regexp(cellstr({neu.experiments{k}.slist.name}.'),'repla|Repla|catego|Catego'))))
                    neu.extra.expinfo(k).keywords = [neu.extra.expinfo(k).keywords,"repla"];
                end

                if isempty(neu.extra.expinfo(k).keywords)
                    neu.extra.expinfo(k).keywords = [neu.extra.expinfo(k).keywords,"other"];
                end

            end

        end

        function neu = getInfo_pre(neu)
            % 获取神经元信息，存入到matfile里，以便快速访问


            neu.info = struct;
            temp = regexp(neu.experiments{1, 1}.inputs.folder_wav,'[OGBYR]\d{3}','match','once');
            %temp = regexp(neu.experiments{1}.info.pl2name,'[RBOYRG]\d{3}','match');
            % set birdid uniqueid and formated_name
            if ~isempty(temp)
                neu.info.birdid = temp{1};
            end
            neu.info.zpid = regexp(neu.experiments{1}.info.pl2name,'[ZP]\d{2}','match');
            neu.info.channelname = neu.experiments{1}.info.channelname;
            neu.info.unitname = neu.experiments{1}.info.unitname;
            neu.info.formated_name = sprintf('%s_%s_%s_%u',neu.info.birdid,neu.info.zpid{1},neu.info.channelname,neu.info.unitname);


        
            %             neu.info.fr_info = neu.multiRepeatsFiringRate;

            % calculate number of responsive
            % songs!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            % 可能要从Sultan里面找，然后改对应的function

            if ~isempty(regexp(neu.info.zpid,'Z')) %zeus
                neu.info.pl2_data_fs = 30000; %hard code !!!!!! Dangerous
            elseif ~isempty(regexp(neu.info.zpid,'P')) % plexon
                neu.info.pl2_data_fs = 40000;
            end






            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % calculate fr_info = multiRepeatsFiringRate(s);
            %
            %                 all_info(k).neuronname = A.neuronname;
            %
            %                 sum_prelen = 0; % summed prey length
            %                 concat_presptimes = []; % concatenated prey sptimes
            %
            %                 sum_pltlen = 0; %summed prey( stimuli y, not plty or rawy) length
            %                 concat_pltsptimes = []; %  % concatenated y sptimes
            %
            %                 all_es = A.getAllEphysObject;

            %                 for m = 1: length(all_es)
            %
            %
            %                     % for prey
            %                     all_info(k).presptimes{m} = all_es{m}.presptimes
            %                     all_info(k).preylen{m} = length(all_es{m}.y)/all_es{m}.fs;
            %                     all_info(k).repnum{m} = size(all_es{m}.presptimes,2);
            %                     temp = all_es{m}.presptimes.';
            %                     concat_presptimes = [concat_presptimes;vertcat(vertcat(temp{:}))+ sum_prelen];
            %                     sum_prelen = sum_prelen +  all_info(k).preylen{m};
            %
            %                     % for plty
            %                     all_info(k).pltsptimes{m} = all_es{m}.pltsptimes
            %                     all_info(k).pltlen{m} = length(all_es{m}.plty)/all_es{m}.fs;
            %                     temp = all_es{m}.pltsptimes.';
            %                     concat_pltsptimes = [concat_pltsptimes;vertcat(vertcat(temp{:}))+ sum_pltlen];
            %                     sum_pltlen = sum_pltlen +  all_info(k).pltlen{m};
            %
            %                 end
            %                 for pre_y
            %                     all_info(k).concat_pre_sptimes = concat_presptimes;
            %                     all_info(k).concat_pre_len = sum_prelen;
            %                     all_info(k).mean_pre_fr = length(concat_presptimes)/sum_prelen;
            %
            %                     for plt_y
            %                         all_info(k).concat_plt_sptimes = concat_pltsptimes;
            %                         all_info(k).concat_plt_len = sum_pltlen;
            %                         all_info(k).mean_plt_fr = length(concat_pltsptimes)/sum_pltlen;

            %all_info = rmfield( all_info,'neuronname');
            %wlfr_info = table2struct([struct2table(wl_info),struct2table(fr_info)]);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%% For classify
            %                 classify_info = table2struct(readtable(path_classifyfile));
            %
            %                 Nid = regexp(regexp(ANAfiles{k},'_\d+','match'),'\d+','match');
            %                 Nid = Nid{1};
            %                 Nid = str2num(Nid{1});
            %
            %                 id_in_cla = find(ismember([classify_info.NeuronID].',Nid));
            %                 all_info(k).neurontype = classify_info(id_in_cla).Type;
            %                 all_info(k).reliableDeg = classify_info(id_in_cla).Deg;
            %                 all_info(k).reliableFrag = classify_info(id_in_cla).Frag;
            %                 all_info(k).reliableRepla = classify_info(id_in_cla).Repla;





        end


        function neu = getInfo_post(neu)
            %             Ninfo.num_respcons
            neu.info.wl = neu.calMeanWaveLength;  % calculate mean wavelength
            neu.info.num_resp_to_18 = neu.song.info.num_resp_to_18;
            neu.info.name_resp_to_18 = neu.song.info.name_resp_to_18;
            neu.info.fr_info = neu.multiRepeatsFiringRate;

        end


        function info = getInfo_more(neu)
            % calculate lifetime sparseness,correlation index, number of
            % responsive songs, spontaneous firing rate, spike width
            % 目的是通过计算这些性质的值对experiments进行划分

            normlist = neu.song.list18;

            %thres = 0.001; % 1ms
            thres = 0.001;
            sum_prelen = 0; % summed prey length
            concat_presptimes = []; % concatenated prey sptimes

            sum_judgeresplen = 0; %summed prey( stimuli y, not plty or rawy) length
            concat_judgerespsptimes = []; %  % concatenated y sptimes

            sumNs = [];
            sdf_collect = {};
            for m = 1: length(normlist) % 这一大段其实似乎只是在重复judgeNormResponse，只不过用了新的方法
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
                info.eachCI(m) = Ns/normalization_factor;

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
                info.presptimes{m} = normlist(m).prejudgerespsptimes;
                info.preylen{m} = length(normlist(m).y)/normlist(m).fs;
                info.repnum{m} = size(normlist(m).prejudgerespsptimes,2);
                temp = normlist(m).prejudgerespsptimes.';
                concat_presptimes = [concat_presptimes;vertcat(vertcat(temp{:}))+ sum_prelen];
                sum_prelen = sum_prelen +  info.preylen{m};

                % for plty
                info.judgerespsptimes{m} = normlist(m).judgerespsptimes
                info.judgeresplen{m} = length(normlist(m).judgerespy)/normlist(m).fs;
                temp = normlist(m).judgerespsptimes.';
                %concat_judgerespsptimes = [concat_judgerespsptimes;vertcat(vertcat(temp{:}))+ sum_judgeresplen];
                concat_judgerespsptimes = [concat_judgerespsptimes; cellfun(@(x) x+sum_judgeresplen,temp,'Uni',0) ];
                sum_judgeresplen = sum_judgeresplen +  info.judgeresplen{m};


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
                info.eachsparseness(m) = 1 - 2*Avalue;



% 
%                 pre_frs = Cal.eachTrialFiringRate(normlist(m).prejudgerespsptimes,length(normlist(m).judgerespy)/normlist(m).fs);
%                 sti_frs = Cal.eachTrialFiringRate(normlist(m).judgerespsptimes,length(normlist(m).judgerespy)/normlist(m).fs);
%                 [h,p] = ttest(sti_frs,pre_frs,'Tail','Right','Alpha',0.05)
%                 normlist(m).pvalue = p;
%                 normlist(m).label = 0; % 初始化
%                 if h == 1
%                     normlist(m).label = 1;
%                 elseif maxsdf > 17 && maxsdf > maxpresdf % neu rescue
%                     normlist(m).label = 1;
% 
%                 end


            end


            allsdf = horzcat(sdf_collect{:});
            info.allsdf = allsdf;
            is1ids = find([normlist.label].' == 1);
            info.numrespsong = length(is1ids);
            label1normlist = normlist(is1ids);

            if ~isempty(label1normlist)
                info.meanCI = mean([label1normlist.CI].');
                info.maxCI = max([label1normlist.CI].');
                info.minCI = min([label1normlist.CI].');
                info.meansparseness = mean([label1normlist.sparseness].');
                info.maxsparseness = max([normlist.sparseness].');
                info.forpca = horzcat(info.eachsparseness,info.eachCI);

            else

                info.meanCI =nan;
                info.maxCI = nan;
                info.minCI = nan;
                info.meansparseness = nan;
                info.maxsparseness = nan;
                info.forpca = nan;
            end

            %

            % for norm songs, degressive songs, and detailed
            % frags/replas, calculate the concatenated firing rate one
            % by one
            % for pre_y
            info.concat_pre_sptimes = concat_presptimes;
            info.concat_pre_len = sum_prelen;
            info.mean_pre_fr = length(concat_presptimes)/sum_prelen;

            % for plt_y
            if ~isempty(concat_judgerespsptimes)
                info.concat_judgeresp_sptimes = vertcat(concat_judgerespsptimes{:});
            else
                info.concat_judgeresp_sptimes =[];
            end
            info.concat_judgeresp_len = sum_judgeresplen;
            info.mean_judgeresp_fr = length(info.concat_judgeresp_sptimes)/sum_judgeresplen;
            info.meanWL = neu.calMeanWaveLength;

            if ~isempty(normlist)
                sumM = max([length(normlist(1).sptimes),length(normlist(2).sptimes),length(normlist(3).sptimes)]); % number of presentation
            else
                sumM = [];
            end
            sumD = sum_judgeresplen; % stimulus duration
            sum_allspikes = info.concat_judgeresp_sptimes;
            sumr = length(sum_allspikes)/sum_judgeresplen;  % average firing rate
            spikenum = length(sum_allspikes);
            info.spikenum = spikenum;
            info.afr =sumr;% average firing rate
            sumomega = thres; % coincidence window
            sumnormalization_factor = sumM*(sumM-1)*(sumr.^2)*sumomega*sumD;
            if ~isempty(sumNs)
                info.sumCI = sum(sumNs)/sumnormalization_factor;
            else
                info.sumCI = [];
            end


            % calculate sum sdf
            if ~isempty(normlist)  % 那么多的isempty应该有问题，后面需要修改一下 % 2022.12.11
                sumsdf = Cal.sdf({info.concat_judgeresp_sptimes},zeros(uint64(sum_judgeresplen*normlist(m).fs),1),normlist(m).fs,0.001,0.004);
                [sumsdf,~] = histcounts(info.concat_judgeresp_sptimes,round(sumD/0.001));


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


                info.kurtosis = sum(summer)/length(sumsdf) -3;

                sumAvalue = trapz(sumhundredthres,sumfraction_above)/100;
                info.sumsparseness = 1 - 2*sumAvalue;
            else

                info.kurtosis = [];

                info.sumsparseness = [];

            end


        end

        function neu = getInfo_last(neu)
            neu.info.percentages = neu.frag.percentageResponse;
            

        end

     

        function fr_info = multiRepeatsFiringRate(neu)

            % and also calculate the firing for each sub-file

            % 计算各种定义下，平均所有stimuli后的firing rate
            fr_info = struct;
            fr_info.neuronname = neu.info.formated_name;

            sum_prelen = 0; % summed prey length
            sum_preJlen = 0;
            concat_presptimes = []; % concatenated prey sptimes

            sum_pltlen = 0; %summed prey( stimuli y, not plty or rawy) length
            sum_ylen = 0;
            sum_Jlen = 0;

            concat_pltsptimes = []; %  % concatenated y sptimes
            concat_ysptimes = [];
            concat_preJsptimes = [];
            concat_Jsptimes = [];
            all_es = neu.getAllEphysObject;

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


                % for prejudgerespy
                fr_info.preJsptimes{m} = all_es{m}.prejudgerespsptimes
                fr_info.preJlen{m} = length(all_es{m}.judgerespy)/all_es{m}.fs;
                temp = all_es{m}.prejudgerespsptimes.';
                concat_preJsptimes = [concat_preJsptimes;vertcat(vertcat(temp{:}))+ sum_preJlen];
                sum_preJlen = sum_preJlen +  fr_info.preJlen{m};


                % for judgerespy
                fr_info.Jsptimes{m} = all_es{m}.judgerespsptimes;
                fr_info.Jlen{m} = length(all_es{m}.judgerespy)/all_es{m}.fs;
                temp = all_es{m}.judgerespsptimes.';
                concat_Jsptimes = [concat_Jsptimes;vertcat(vertcat(temp{:}))+ sum_Jlen];
                sum_Jlen = sum_Jlen +  fr_info.Jlen{m};

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

    end

    methods% 内部计算方法

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

        function neu = judgeFragResp_FR(neu)
            dbstop if error % 判断对frag 是否反应，通过 Firing rate

   

            ids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'Frag|frag|syl|ele'))); % find all norms
            % ’syl'可以兼容旧的stimuli命名规则

            if isempty(ids)
                return
            end
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
                [h,p] = ttest(sti_frs,pre_frs,'Tail','Right','Alpha',0.05); % should be 0.05!!!


                neu.list(thisi).pvalue = p;
                repeats = length(neu.list(thisi).judgerespsptimes);
                neu.list(thisi).sti_frs = length(vertcat(neu.list(thisi).judgerespsptimes{:}))...
                    /((length(neu.list(thisi).judgerespy)/neu.list(thisi).fs) *repeats);
                neu.list(thisi).maxsdf = maxsdf;
                neu.list(thisi).label = 0; % 初始化
                if h == 1 &&...
                        length(find(~cellfun(@isempty,neu.list(thisi).judgerespsptimes)))/length(neu.list(thisi).judgerespsptimes) >=0.5...
                        && maxsdf>13
                    neu.list(thisi).label = 1;
                end
                %                 disp([h,p,neu.list(thisi).label])  下面的三行是为测试之用
                %                 fig1 = figure('Position',[2407 186 529 403]); Draw.two(neu.list(thisi).judgerespy,32000,neu.list(thisi).judgerespsptimes);
                %                 fig2 = figure('Position',[1898 205 529 403]); Draw.two(zeros(length(neu.list(thisi).judgerespy),1),32000,neu.list(thisi).prejudgerespsptimes);
                %                 close(fig1); close(fig2);
            end

%             neu.frag = Frag(neu.list);


        end

        function neu = judgeConResp(neu,mode)

            if exist("mode",'var')  % 原来的备选方法

                % 判断对Cons是否反应，通过 Firing rate
                % firstly update e objectys 以后可以删掉这个部分
                %             for k = 1: length(neu.experiments)
                %                 for kk = 1: length( neu.experiments{k}.e)
                %                     neu.experiments{k}.e{kk}.setExtAndAllocate;
                %                 end
                %             end
                %             neu.updatelist;

                ids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'norm|deg|repla'))); % find all norms
                % ’syl'可以兼容旧的stimuli命名规则

                for n = 1: length(ids)
                    thisi = ids(n);

                    presdf = Cal.sdf(neu.list(thisi).prejudgerespsptimes,zeros(length(neu.list(thisi).judgerespy),1),neu.list(thisi).fs,0.001,0.02);
                    sdf = Cal.sdf(neu.list(thisi).judgerespsptimes,neu.list(thisi).judgerespy,neu.list(thisi).fs,0.001,0.02); % 0.001,0.004
                    %                 sdf = Cal.sdf(neu.list(thisi).judgerespsptimes,neu.list(thisi).judgerespy,neu.list(thisi).fs,0.001,0.02);
                    [maxpresdf,~] = max(presdf);
                    [maxsdf,maxidx] = max(sdf);

                    pre_frs = Cal.eachTrialFiringRate(neu.list(thisi).prejudgerespsptimes,length(neu.list(thisi).judgerespy)/neu.list(thisi).fs);
                    sti_frs = Cal.eachTrialFiringRate(neu.list(thisi).judgerespsptimes,length(neu.list(thisi).judgerespy)/neu.list(thisi).fs);

                    [h,p] = ttest(sti_frs,pre_frs,'Tail','Right','Alpha',0.05);
                    neu.list(thisi).pvalue = p;
                    neu.list(thisi).label = 0; % 初始化
                    neu.list(thisi).sti_frs = mean(sti_frs);
                    neu.list(thisi).pre_frs = mean(pre_frs);
                    if h == 1
                        neu.list(thisi).label = 1;
                        %                     if (num_of_not_empty_trials/length(neu.list(thisi).pltsptimes)<0.5)||(num_of_not_empty_trials<5) % if not 60% of the trails are not empty
                        %                         neu.list(thisi).label = 0;
                        %                     end
                    elseif maxsdf > 17 && maxsdf > maxpresdf % neu rescue
                        neu.list(thisi).label = 1;

                    end
                end

                return


            end

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

                % ttest2 to test whether sti_frs are significantly higher
                % than pre_frs
                pre_frs = Cal.eachTrialFiringRate(neu.list(thisi).prejudgerespsptimes,length(neu.list(thisi).judgerespy)/neu.list(thisi).fs);
                sti_frs = Cal.eachTrialFiringRate(neu.list(thisi).judgerespsptimes,length(neu.list(thisi).judgerespy)/neu.list(thisi).fs);

                [h,p] = ttest(sti_frs,pre_frs,'Tail','Right','Alpha',0.05);
                neu.list(thisi).pvalue = p;
                neu.list(thisi).label = 0; % 初始化


                % mean_maxsdf = maxsdf/length(neu.list(thisi).judgerespsptimes);
                neu.list(thisi).maxsdf = maxsdf;
                neu.list(thisi).label = 0; % 初始化
                if h == 1&&(maxsdf) > 17 && maxsdf > maxpresdf %如果是 time-locked response
                    neu.list(thisi).label = 1;

                    if (num_of_not_empty_trials/length(neu.list(thisi).pltsptimes)<0.5)||(num_of_not_empty_trials<5) % if not 60% of the trails are not empty
                        neu.list(thisi).label = 0;
                    end

                elseif mean(sdf)> 9*mean(presdf) && mean(sdf)>0.6 % 如果不是 time-locked response
                    neu.list(thisi).label = 1;   %  set to 2 ,biao ming shi fei time-locked response
                end

                neu.song = Song(neu.list); %不仅主list,而且subclass Song 的 norm list也要更新！

                %                 figure;  % 测试用代码
                %                 Draw.three(neu.list(thisi).plty,neu.list(thisi).fs,neu.list(thisi).pltsptimes);
                %                 title(sprintf('Label is %u',neu.list(thisi).label) );
                %                 close(gcf)

            end

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
            %neu.conspe_eleinf = neu.all_eleinf(conspe_ids);
        end

        function targets = knowTarget(neu)
            %找到哪个是目标syllable，目标syllable是生成deg stimuli和repla stimuli的基础
            %目前的版本只是通过repla判断，之后有时间也应该写通过deg判断的方法
            % 2022.12.21 追加通过degressive song 判断
            %当然，如果连degressive song 也没有，那就根本谈不上有targets
            %而为了通过degressive song 判断，首先需要解决的是JudgeDegResp


            % target
            targetidx = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'Repla')));
            targetlist = neu.list(targetidx);
            for k = 1: length(targetlist)

                temp = strsplit(targetlist(k).stimuliname,'-');

                targetnames{k} = sprintf('%s-%s-%s',temp{6},temp{7},temp{8});
            end

            if exist('targetnames','var')
                targets = unique(targetnames);
            else
                targets = {};
            end

        end
        
        function neu =setStimuliCorrespondingNeuronId(neu)
            % 当analysis的neuron对应不同刺激类型，如frag，norm时，阐明各neuron对应的刺激类型
            neu.extra.song_id_redundant = [];
            neu.extra.frag_id = [];
            neu.extra.deg_id = [];
            neu.extra.repla_id = [];
            neu.extra.song_only_id = [];
            for k = 1: length(neu.extra.expinfo)
                if ismember("song",neu.extra.expinfo(k).keywords)
                    neu.extra.song_id_redundant = [neu.extra.song_id_redundant,k];
                end

                if ismember("frag",neu.extra.expinfo(k).keywords)
                    neu.extra.frag_id = [neu.extra.frag_id,k];
                end

                if ismember("deg",neu.extra.expinfo(k).keywords)
                    neu.extra.deg_id = [neu.extra.deg_id,k];
                end

                if ismember("repla",neu.extra.expinfo(k).keywords)
                    neu.extra.repla_id = [neu.extra.repla_id,k];
                end

                if strcmp("song",neu.extra.expinfo(k).keywords)
                    neu.extra.song_only_id = [neu.extra.song_only_id,k];
                end

                if strcmp("other",neu.extra.expinfo(k).keywords)
                    neu.extra.song_only_id = [neu.extra.song_only_id,k];
                end


            end

            if length(neu.extra.song_id_redundant) == 1
                neu.extra.song_id = neu.extra.song_id_redundant;
            elseif length(neu.extra.song_id_redundant)~=0

                neu.extra.song_id = min(neu.extra.song_id_redundant);
                %[~,neu.song_id] = max(cellfun(@length,{neu.extra.expinfo(neu.song_id_redundant).keywords}));
                % find the id which have max length of keywords, if the
                % result is multipl ids, then select the first one
                % But maybe the last one will be much proper!
            else
                neu.extra.song_id = 1; % Very bad meaningless code
            end

        end

        function reordered_ids = Frozen_reorderListByResp(neu)
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
            for k = 1: length(neu.experiments)
                collects{k} = neu.experiments{k}.e(:);
            end
            e_objects = vertcat(collects{:});
        end

        function meanWL = calMeanWaveLength(neu)
            % 计算 meanWL，需要考虑到仪器的fs
            wavecollect = {};
            for s = 1: length(neu.experiments)
                wavecollect{s} = neu.experiments{s}.waves.waveform;
            end
            waveforms = vertcat(wavecollect{:});
            %waveforms =  n.waveform;
            [~,troughstime] = min(waveforms,[],2);
            wavlen_units = [];

            for k = 1: size(waveforms,1) % length is dangerous!!!!!
                this_wf = waveforms(k,:);
                [~,wavlen_units(k)] =  max(this_wf (troughstime(k):end));
            end



            meanWL =  mean(wavlen_units*(1/neu.info.pl2_data_fs)*1000); % ms


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


        function latency = calLatency(neu)

            latency = struct;

            Onelist = neu.frag.allfraglist_nonorm(find([neu.frag.allfraglist_nonorm.label].' == 1));
            % 第一种方法
            observed = [];
            spon = [];
            for k = 1:length(Onelist)
                %https://direct.mit.edu/neco/article/22/7/1675/7562/First-Spike-Latency-in-the-Presence-of-Spontaneous

                thisone = Onelist(k);
                temp1 = min(vertcat(thisone.sptimes{:}));
                temp2 = min(vertcat(thisone.presptimes{:}));
                if isempty(temp1)
                    temp1 = nan;
                end
                if isempty(temp2)
                    temp2 = nan;
                end

                obsevred(k) = temp1;
                spon(k) = temp2;

                %[observed(k),spon(k)] = e_objects{k}.calObservedFirstSpikeLatency;
            end

            % if isnan(mean(rmmissing(spon))) % 当pre duration里一个spike也没有时会出现这种情况
            spontanous_first_spike_latency = 0;


            % 现在这个方法依然不对，因为有可能会出现在pre stimuli的最后一秒，突然出现了first spike的情况
            latency.latency_firstspike = mean(rmmissing(observed))- spontanous_first_spike_latency; %

            % 第二种方法
            summer = [];
            for k = 1:length(Onelist)

                thisE = Onelist(k);
               % figure; Draw.three(thisE.plty,thisE.fs,thisE.pltsptimes) %测试用
             
                % 背景SDF的均值（baseline）
                prezpt_sptimes = Extract.sptimes_resetSP(thisE.rawsptimes, 0, thisE.zpt);
                baseline = mean(Cal.sdf( prezpt_sptimes,thisE.zpt*thisE.fs,thisE.fs));

                endtime = length(thisE.rawy)/thisE.fs;
                postzpt_sptimes = Extract.sptimes_resetSP(thisE.rawsptimes, thisE.zpt,endtime);
                sdfs = Cal.sdf(postzpt_sptimes,(endtime-thisE.zpt)*thisE.fs,thisE.fs);

                maxsdfs = max(sdfs);
                timeindex = find(sdfs > baseline + (maxsdfs - baseline)/2,1); % 第一次超过half weight的时刻
                summer(k) = timeindex/length(sdfs)*(endtime -  thisE.zpt);

            end

            latency.latency_halfweight = mean(summer);

        end

    end

    methods % 作图方法


        % 作图并导出

        function saveDrawAlignedConsDegsReplasFrags(neu)  % Align all together and draw
 
            % 我需要写一个方法得到每个neuron播放的stimuli的原始目的，比如fragments的构成 2022.12.21
            dbstop if error
            tic

            % about songs been presented % Redudant code
            songlist = neu.song.normlist;%Neuron(neu.neurons{neu.song_id}).normlist;
            [~,postunique] = unique(cellstr(cellfun(@Convert.bid,{songlist.stimuliname}.','Uni',0)))
            songlist = songlist(postunique);


            % About Frag
            fragids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'frag|Frag|syl|Syl') ));
            if ~isempty(fragids)

                fraglist = neu.list(fragids);
                %                 normlist = Neuron(neu.neurons{neu.song_id}).normlist;
                %                 [~,postunique] = unique(cellfun(@Convert.bid,[normlist.stimuliname].','Uni',0))
                %                 normlist = normlist(postunique);

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

            % About Deg
            degids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'Deg|deg') ));
            if ~isempty(degids)

                deglist = neu.list(degids);
                for m = 1: length(deglist)
                    birdid = Convert.bid(deglist(m).stimuliname);
                    ids_norm = find(~cellfun(@isempty, regexp(cellstr({songlist.stimuliname}.'),birdid) ) );
                    if ~isempty(ids_norm)& length(ids_norm) == 1
                        [deglist(m).sylIni,trump_diffvalue] = Neuron.findIni(songlist(ids_norm).plty,deglist(m).y);
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
                    birdid = Convert.bid(songlist(w).stimuliname);
                    ids_indeg = find(~cellfun(@isempty, regexp(cellstr({deglist.stimuliname}.'),birdid) ) );
                    selected_deglist = deglist(ids_indeg);
                    [~,temp_index] = sortrows([selected_deglist.sylIni].');
                    selected_deglist = selected_deglist(temp_index);
                end

                if ~isempty(fragids)
                    birdid = Convert.bid(songlist(w).stimuliname);
                    ids_infrag = find(~cellfun(@isempty, regexp(cellstr({fraglist.stimuliname}.'),birdid) ) );
                    selected_fraglist = fraglist(ids_infrag);
                    [~,temp_index] = sortrows([selected_fraglist.sylIni].');
                    selected_fraglist = selected_fraglist(temp_index);
                end


                % draw the basic figure
                Icollect = {};
                figure('Color','w');

                Draw.two(songlist(w).plty,songlist(w).fs,songlist(w).pltsptimes);
                xlabel(songlist(w).stimuliname);
                frame = getframe(gcf);
                Icollect{1} = frame.cdata;
                close(gcf)

                if ~isempty(degids)
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

            imwrite(Iall,sprintf('Aligned_Normfrag_%s.png',neu.info.formated_name));
            toc


        end

        function img = exportConRevMir(neu)
            % 作出Cons Revs 和 Mirs对照的图
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
                if length(sb_mirror_ids) > 1
                    sb_mirror_ids = sb_mirror_ids(1); % 如果多项，暂时的对策是取第一项
                end
                sb_norm_ids = intersect(normids,samebirdids);
                if length(sb_norm_ids) > 1
                    sb_norm_ids =sb_norm_ids(1); % 如果多项，暂时的对策是取第一项
                end

                samefileids  = find(~cellfun(@isempty,regexp(cellstr({neu.list.Fid}.'),neu.list(sb_reverse_ids).Fid)));
                sameFile_norm_ids = intersect(samefileids,sb_norm_ids);
                sameFile_reverse_ids = intersect(samefileids,sb_reverse_ids);
                sameFile_mirror_ids = intersect(samefileids,sb_mirror_ids);

                fig1 = figure('Position',[277 330 731 607],'Color','w');
                Draw.three(neu.list(sameFile_norm_ids).y,neu.list(sameFile_norm_ids).fs,neu.list(sameFile_norm_ids).sptimes); % plty
                frame = getframe(fig1);
                IMG{1,k} = frame.cdata;
                close(fig1)

                fig2 = figure('Position',[277 330 731 607],'Color','w');
                Draw.three(neu.list(sameFile_reverse_ids).y,neu.list(sameFile_reverse_ids).fs,neu.list(sameFile_reverse_ids).sptimes);
                frame = getframe(fig2);
                IMG{2,k} = frame.cdata;
                close(fig2)

                fig3 = figure('Position',[277 330 731 607],'Color','w');
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

            imwrite(img,'卫—CONsRevsMirs.png')

        end

        % 作图且保存
        function Iall = saveDrawAlignedConsFrag(neu,songnames)
            % songnames 指定后则只会显示对应于songnames的fragments
            dbstop if error
            tic

            fragids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'frag|Frag|syl')));

            fraglist = neu.list(fragids);
            frag_Fid = unique({fraglist.Fid}.');
            % normlist = Neuron(neu.experiments{neu.song_id}).normlist;

            subfile_frag_ids = find(~cellfun(@isempty,regexp(cellstr({neu.list.Fid}.'),strjoin(frag_Fid,'|'))));
            hard_to_name_ids = subfile_frag_ids;
            if exist('songnames','var')
                songnameids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),strjoin(songnames,'|'))));
                hard_to_name_ids = intersect(subfile_frag_ids,songnameids);
                fucklist = neu.list(hard_to_name_ids);

                normlist = fucklist(find(~cellfun(@isempty, regexp(cellstr({fucklist.stimuliname}.'),'norm'))));
            end

            normlist = neu.list(find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'norm'))));

            % about Norms % Redudant code
            %             normlist = Neuron(neu.experiments{neu.song_id}).normlist;
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

                Draw.two(normlist(w).plty,normlist(w).fs,normlist(w).pltsptimes);
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

            %             neu.drawSeparatedWaveform;
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

            imwrite(Iall,sprintf('魏Aligned_Normfrag_%s.png',neu.info.formated_name));
            toc

        end

        function drawFragScatter(neu,eleinf_with_coors, not_tested_handle)
            %syllable response的散点图，每一个点代表一个syllable
            %这个图也会显示syllables是否使该神经元反应
            % not_tested_handle = 1 means draw

            neu.pandingTargetSyllable;

            if ~exist('not_tested_handle','var') % to judge whether to draw not tested elements or not
                not_tested_handle = 0;
            end

            %global_eleinf = a.all_eleinf;
            global_eleinf = eleinf_with_coors;
            neu.judgeFragResp_FR;
            fraglist = neu.frag.fraglist;


            for k = 1: length(global_eleinf)
                uniqueid = sprintf('-%s-%u-',global_eleinf(k).songname, global_eleinf(k).fragid);

                if ~isempty (find(~cellfun(@isempty, regexp(cellstr({fraglist.stimuliname}.'),uniqueid)), 1)) % if this element was tested
                    id_in_fraglist = find(~cellfun(@isempty, regexp(cellstr({fraglist.stimuliname}.'),uniqueid)));
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
            figure('Position',[530 146 1230 886],'Color','w');
            hold on

            minus1_ids = find([global_eleinf.scatter].' == -1);
            zero_ids = find([global_eleinf.scatter].' == 0);
            one_ids = find([global_eleinf.scatter].' == 1);

            coor1s = [global_eleinf.coor_1].';
            coor2s = [global_eleinf.coor_2].';

            if not_tested_handle == 1
                scatter(coor1s(minus1_ids),coor2s(minus1_ids),[],[.5, .5, .5],'filled'); % grey for not tested
            end
            scatter(coor1s(zero_ids),coor2s(zero_ids),[],'g','filled'); % green for tested, and not eliciting
            scatter(coor1s(one_ids),coor2s(one_ids),[],'r','filled'); % red for response-eliciting

            % label the targets
            global_names = [global_eleinf.songname].';
            global_fragids = [global_eleinf.fragid].';

            global_merged = {};
            for w = 1: length(global_names)
                global_merged{w} = sprintf('%s-%u',global_names(w),global_fragids(w));
            end
            global_merged = global_merged.';

            try
                for u = 1: length(neu.targets)
                    [~,beta] = ismember( neu.targets{u}, global_merged );
                    scatter( global_eleinf(beta).coor_1,global_eleinf(beta).coor_2,[],'k','h');
                end

            end
            xlabel('Dim-1');
            ylabel('Dim-2');
            title(neu.formated_name,'interpreter','none');
            hold off
            saveas(gcf,sprintf('Another_Eleinf_SyllableScatter_%s.png',neu.formated_name));
            close gcf

        end

        function Iall = saveDrawAlignedConsDegs(neu,songnames)
            dbstop if error

            %只有一个模式： 只针对二次播放里包含的norm songs进行degs的对齐
            tic

            degids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'deg|Deg')));

            deglist = neu.list(degids);
            deg_Fid = unique({deglist.Fid}.');
            % normlist = Neuron(neu.experiments{neu.song_id}).normlist;

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
            deglist = neu.list(degids);

            for m = 1: length(deglist)
                birdid = Convert.bid(deglist(m).stimuliname);
                ids_norm = find(~cellfun(@isempty, regexp(cellstr({normlist.stimuliname}.'),birdid))); % 这个deg stimuli对应的norm stimuli

                if isempty(ids_norm)
                    continue %如果找不到对应norm就跳过
                elseif length(ids_norm)>1
                    ids_norm = ids_norm(1);%如果有多个，就取第一个
                end
                [deglist(m).sylIni,trump_diffvalue] = Neuron.findIni(normlist(ids_norm).plty,deglist(m).y);%此处目的是找到deg与norm的分歧点，所以具体是哪个E的norm不重要
                fprintf('%s has the DIFFVALUE as: %u and INI as: %f \n ',deglist(m).stimuliname,trump_diffvalue,deglist(m).sylIni/32000);
                deglist(m).pady = [zeros([deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext,1]);deglist(m).plty;zeros(length(normlist(ids_norm).plty)...
                    - (deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext) - length(deglist(m).plty),1)];
                deglist(m).padsptimes = cellfun( @(x) x + (deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext)/deglist(m).fs, deglist(m).pltsptimes,'uni',0);

            end

            % merge the new fraglist and the deglist with the normlist
            I_of_each_column = {};
            counter = 0;
            for w = 1: length(normlist)

                if ~isempty(degids)
                    birdid = Convert.bid(normlist(w).stimuliname);
                    ids_indeg = find(~cellfun(@isempty, regexp(cellstr({deglist.stimuliname}.'),birdid) ) );
                    selected_deglist = deglist(ids_indeg);
                    if isempty(ids_indeg)
                        continue
                    end
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
                counter = counter + 1;
                I_of_each_column{counter} = vertcat(Icollect{:});
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
            imwrite(Iall,sprintf('燕Aligned_ConsDegs_%s.png',neu.info.formated_name));
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

            imwrite(Iall,sprintf('韩Aligned_Repla_With_Norm_%s.png',neu.info.formated_name));
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

            end


            % 其末 draw and save
            figure;
            neu.waveform.draw1st;
            temp = getframe(gcf);close(gcf);
            w_img = temp.cdata;
            I_WF{1} = w_img;

            try
                I_of_each_column = horzcat(I_song,I_Deg,I_Frag,I_Repla,I_WF);
            catch
                I_of_each_column = horzcat(I_song,I_Deg,I_Frag,I_WF);
            end
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

            imwrite(Iall,sprintf('楚Whether_NeuResp_To_SinFrags_Coms_Or_FragsResps_affected_By_Pres_%s.png',neu.info.formated_name));
            toc


        end


    
        
    end

    methods(Static) % 静态方法



        function check_detail_folder(dirpath,global_eleinf)
            % 这个方法可以用来看所有 presented frag stimuli的分布
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
            temp = split(temp,'_');
            targetname = sprintf('%s-%s-%s',temp{1},temp{2},temp{3});
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
            frag_0_eleinf = global_eleinf([global_eleinf.whether_frag_is_tested].' == 0);
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

            %             diff_info = struct;

            [SynIni,~,diff_value] = findsignal(y,yfrag);

            %             parfor k = 1: length(y) - length(yfrag)
            %                 totest = y(k:k+ length(yfrag) -1);
            %
            %                 if sum(totest) ~= 0
            %
            %                     diff = sum(abs(totest - yfrag));
            %
            %                     diff_info(k).diff = diff;
            %                     diff_info(k).kvalue = k;
            %                     %                      end
            %                 else
            %
            %                     diff_info(k).diff = Inf;
            %                     diff_info(k).kvalue = k;
            %                 end
            %
            %
            %             end
            %
            %             [~,min_ids] = min([diff_info.diff].');
            %             SynIni = diff_info(min_ids).kvalue;
            %             diff_value = diff_info(min_ids).diff;

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

            %figure; subplot(311); plot(y);subplot(312); plot(yrepla);subplot(313); plot(yrepla-y);

            temp_convergentpoint = find(flip(y-yrepla), 1 ); % 找到第一个非零元素
            % 如果把y和yrepla倒过来之后第一个非零元素是第一位，说明y和yrepla不同，这stimuli有问题
            %             [coeffsy,~,~,~] = mfcc(flip(y),32000);
            %             [coeffsyrepla,~,~,~] = mfcc(flip(yrepla),32000);
            %             figure; plot(sum(coeffsy,2) - sum(coeffsyrepla,2));
            %             figure; subplot(311); Draw.spec(y,32000);subplot(312); Draw.spec(yrepla,32000);subplot(313);Draw.spec(yrepla-y,32000);
            %            figure; plot( flip(round(y,3)) - flip(round(yrepla,3)) )

            %figure; plot(flip(y-yrepla))
            convergentpoint = maxlength - temp_convergentpoint + 2;

            if maxid == 1 % y更长
                ConvergentIndexY  = convergentpoint;
                ConvergentIndexReplaY = convergentpoint-(length(freezey)-length(freezeyrepla));
            elseif maxid ==2 % replay 更长
                ConvergentIndexY  = convergentpoint-(length(freezeyrepla)-length(freezey));
                ConvergentIndexReplaY = convergentpoint;
            end


        end


        function [answer,pvalue] = UseTtestToJudgeRespOrNot(y,sptimes,prey,presptimes,fs)
            dbstop if error
            presdf = Cal.sdf(presptimes,prey,fs,0.001,0.02);
            sdf = Cal.sdf(sptimes,y,fs,0.001,0.02); % 0.001,0.004
            %             [maxpresdf,~] = max(presdf);
            %             [maxsdf,maxidx] = max(sdf);

            pre_frs = Cal.eachTrialFiringRate(presptimes,length(prey)/fs);
            sti_frs = Cal.eachTrialFiringRate(sptimes,length(y)/fs);
            [h,p] = ttest(sti_frs,pre_frs,'Tail','Right','Alpha',0.05);
            if h == 1
                answer = 1;
            elseif h == 0||isnan(h)
                answer = 0;
            end
            pvalue = p;

        end


        function result = UseMaxSdfToJudgeRespOrNot(y,sptimes,prey,presptimes,fs)
            dbstop if error
            presdf = Cal.sdf(presptimes,prey,fs,0.001,0.02);
            sdf = Cal.sdf(sptimes,y,fs,0.001,0.02); % 0.001,0.004


            % for separated sdfs
            for k = 1:length(presptimes)

            end
            
            %             [maxpresdf,~] = max(presdf);
            [maxpresdf,~] = max(presdf);
            [maxsdf,maxidx] = max(sdf);

            pre_frs = Cal.eachTrialFiringRate(presptimes,length(prey)/fs);
            sti_frs = Cal.eachTrialFiringRate(sptimes,length(y)/fs);

            [h,p] = ttest(sti_frs,pre_frs,'Tail','Right','Alpha',0.05);
            result.pvalue = p;
            result.label = 0; % 初始化
            result.sti_frs = mean(sti_frs);
            result.pre_frs = mean(pre_frs);
            result.maxsdf = maxsdf;
            if h == 1
                result.label = 1;
            
                %                     if (num_of_not_empty_trials/length(neu.list(thisi).pltsptimes)<0.5)||(num_of_not_empty_trials<5) % if not 60% of the trails are not empty
                %                         neu.list(thisi).label = 0;
                %                     end
            elseif maxsdf > 18 && maxsdf > maxpresdf % neu rescue
                result.label = 1;

            end

        end


    end

    methods(Hidden = true) % 弃用方法

        function sponFrInfo = Deprecated_getSponFR_Sarah(neu,range)
            % calculate spontaneous firing rate
            sponFrInfo = struct;
            all_es = neu.getAllEphysObject;
            for m = 1: length(all_es)


                % for prey
                sponFrInfo(m).triggerNum = all_es{m}.sound.trigger;
                sponFrInfo(m).presptimes = all_es{m}.presptimes;
                sponFrInfo(m).preylen = length(all_es{m}.y)/all_es{m}.fs;
                sponFrInfo(m).repnum = size(all_es{m}.presptimes,2);
                temp = all_es{m}.presptimes.';
                sponFrInfo(m).localSpFr = length(find(vertcat(vertcat(temp{:}))))/(sponFrInfo(m).preylen*sponFrInfo(m).repnum);
                % for plty
                sponFrInfo(m).pltsptimes = all_es{m}.pltsptimes;
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

        function neu = Deprecated_writeFigdata(neu)
            % 生成这个Analysis的所有three plot的图片
            figmat = {};
            for k = 1: length(neu.experiments)
                figmat{k} = neu.experiments{k}.writeFigdata
            end
            neu.figdata = horzcat(figmat{:});
        end

        function neu = Deprecated_splitStimuliResponsePairsToDifferentTypes(neu)
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


                if (maxsdf) > 19.5&& strcmp(Conlist(n).stimuliname,'norm-Y515A-21Pulses') && strcmp(Conlist(n).pl2name,'Y661_Z17')
                    disp('Incredible bug in pltsptimes of function Neuron.evaluateConResponse !!');
                    Conlist(n).label = 1;
                end

            end

        end

        function Deprecated_experiments_Respond_To_Single_Frags_Or_Combinations(neu)
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

            % imwrite(Iall,sprintf('Aligned_ConsDegs_%s.png',neu.experiments{1}.neuronname));
            imwrite(Iall,sprintf('Aligned_ConsDegs_%s.png',neu.formated_name));
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
            neu.experiments{neu.song_id}.draw_waveform;     % draw waveform
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
            imwrite(IMG,sprintf('Three_%s.png',neu.formated_name));

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

            imwrite(IMG,sprintf('ReplaThree-%s.png',neu.formated_name));

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

            imwrite(IMG,sprintf('DegThree-%s.png',neu.formated_name));

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

        function neu = Deprecated_judgeFragResp(neu)
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

        function split_fraginf = Deprecated_calResponseToWithinSongFragsFromEleinf(neu)
            % deprecated
            latency = 50*0.001; % 50 ms

            split_fraginf = struct;

            counts = 0;
            if isempty(neu.normlist)
                split_fraginf = struct([]);
                return
            end

            for k = 1: length(neu.normlist)
                %songname = regexp(neu.normlist(k).stimuliname,'(CON|SPE|norm)-([BRGOY]\d{3}|Fcall|Mcall|WNS|HET|Het)','match');
                songname = regexp(neu.normlist(k).stimuliname,'([BRGOY]\d{3}|Fcall|Mcall|WNS|HET|Het)','match');

                % use autoseg to segment the elements
                %                [rawy,fiy,I,syledge,eleedge] = main(fiy,fs,birdid,CONFIG)
                %
                eleids = find(~cellfun(@isempty, regexp([neu.conspe_eleinf.songname].',songname)));
                local_eleinf = neu.conspe_eleinf(eleids);

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

                    split_fraginf(counts).padded_sptimes = Extract.sptimes(neu.normlist(k).rawsptimes, split_fraginf(counts).initial...
                        ,split_fraginf(counts).terminal + latency);

                    front_percentage = local_eleinf(m).yini/length(neu.normlist(k).rawy);
                    back_percentage = local_eleinf(m).yter/length(neu.normlist(k).rawy);

                    featurenames = fieldnames(neu.normlist(k).rawfeatures); % iterate for each feature
                    featurenames = setdiff(featurenames,{'file_index','file_name'});
                    for omega = 1: length(featurenames)
                        lenfeature = length(neu.normlist(k).rawfeatures.(featurenames{omega}));
                        if  round(front_percentage*lenfeature)== 0

                            split_fraginf(counts).(featurenames{omega}) = neu.normlist(k).rawfeatures.(featurenames{omega})(...
                                1:round(back_percentage*lenfeature));
                        else
                            split_fraginf(counts).(featurenames{omega}) = neu.normlist(k).rawfeatures.(featurenames{omega})(...
                                round(front_percentage*lenfeature):round(back_percentage*lenfeature));
                        end
                        split_fraginf(counts).(sprintf('mean_%s',featurenames{omega})) = mean(split_fraginf(counts).(featurenames{omega}));

                    end

                end

            end

            for n = 1: length(split_fraginf)
                tempsum = Cal.psth_frag(split_fraginf(n).padded_y,split_fraginf(n).fs,split_fraginf(n).padded_sptimes);
                if isempty(tempsum)
                    halfsum = 0;
                else
                    halfsum = sum(tempsum(end/2:end));
                end
                fullsum = sum(tempsum);
                maxvalue = max(Cal.psth_frag(split_fraginf(n).padded_y,split_fraginf(n).fs,split_fraginf(n).padded_sptimes));
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

    end

end



