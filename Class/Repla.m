classdef Repla < handle
    %@ 类说明 Member class of Neuron

    properties
        replalist
        list
        formated_name
        simatrix
    end

    methods

        function rp = Repla(list)
            % 从 Neuron 文件中提取neuron对replas的反应，已经具体的replas的名称
            rp.list = list;
            rp.genReplalist
            rp.judgeReplaResp; % judgeReplaResponse是在生成replalist之后的，所以应该比较方便修改
            
            for k = 1:length(rp.replalist)

                pseudo_list = struct('stimuliname',{rp.replalist{k}.stimuliname}.','y',{rp.replalist{k}.replaceparty}.',...
                    'fs',{rp.replalist{k}.fs}.','label',{rp.replalist{k}.label}.'); % 把replalist变成fraglist的格式，以便计算simatrix
                for kk = 1:length(pseudo_list)
                    tempsat = SAT_sound(pseudo_list(kk).y,pseudo_list(kk).fs);
                    pseudo_list(kk).features = tempsat.features;
                    pseudo_list(kk).meanfeatures = tempsat.meanfeatures;

                end
                rp.simatrix{k} = Simatrix(pseudo_list);
            end


        end

        function genReplalist(rp)
            % 从list里筛选出每组repla stimuli，加上一些提取的信息生成replalist；


            replaids = find(~cellfun(@isempty, regexp(cellstr({rp.list.stimuliname}.'),'repla|Repla') ));
            list_contain_replas = rp.list(replaids);
            unpack = @(x) x{1};
            second_syl = {};
            for k = 1:length(list_contain_replas)
                second_syl{k} = regexp(list_contain_replas(k).stimuliname,'(?<=before-)\S+(?=-gapis-)','match'); % combination中的第二位
                if iscell(second_syl{k})
                    second_syl{k} = unpack(second_syl{k});
                end
            end
            tempo = cellstr(second_syl.');
            second_syl = tempo;
            unique_2ndsyl = unique(second_syl);
            unique_2ndbird = cellfun(@(x) Convert.bid(x,1),unique_2ndsyl,'Uni',0);
            %repla_birdid = cellfun(@(x) Convert.bid(x,2), cellstr({list_contain_replas.stimuliname}.'),'Uni',0);
            %unique_repla_birdids = unique(repla_birdid);
            unpack = @(x) x{1};

            for k = 1:length(unique_2ndsyl) % 筛选出repla stimuli构成list
                correspids = strcmp(unique_2ndsyl{k},second_syl); % 找到有这个birdid的所有repla的ids,这ids是关于list_contain_replas的
                rp.replalist{k} = list_contain_replas(correspids);

                corresp_fid = unpack(unique({rp.replalist{k}.Fid}.')); % 找到repla stimuli对应的fid

                corresp_norm = mintersect( find(~cellfun(@isempty, regexp(cellstr({rp.list.stimuliname}.'),unique_2ndbird{k}))),...
                    find(~cellfun(@isempty, regexp(cellstr({rp.list.stimuliname}.'),'norm'))),...
                    find(~cellfun(@isempty, regexp(cellstr({rp.list.Fid}.'),corresp_fid))) );

                if isempty(corresp_norm ) %如果因为实验中的疏漏没有在repla环节播放norm song， % 那么取其他次播放的norm song
                    corresp_norm = min(intersect( find(~cellfun(@isempty, regexp(cellstr({rp.list.stimuliname}.'),unique_2ndbird{k}))),...
                        find(~cellfun(@isempty, regexp(cellstr({rp.list.stimuliname}.'),'norm'))))); % 如果对应了多个值，那么取最小的那个

                elseif length(corresp_norm)>1
                    % 如果因为一些bug或者其他什么原因导致依旧有多个norm song的话，取序数最大的
                    % （这个做法很武断，但是权宜之计，而如果merged file只有一个Fid的bug能够被成功修复的话，这种话情况不会发生）
                    corresp_norm = max(corresp_norm);
                end

                rp.replalist{k} = vertcat(rp.list(corresp_norm).',rp.replalist{k}.');

            end


            for a = 1:length(rp.replalist) % 补充replalist特有的一些信息
                thislist = rp.replalist{a};

                Replst_fnames = {thislist.stimuliname}.';
                locY = @(x,y) x{y}; % get element in location y
                k_of_norm = [];
                for k = 1:length(thislist)
                    splited = split(Replst_fnames{k},{'-before-','-gapis-'});

                    if length(splited) == 1 % 如果这个stimuli是norm的话,可以看成是原来的song和原来的song发生replacement
                       k_of_norm = k;
                        continue %暂时跳过，目的是利用后续行的fid填充norm song的fid 
                    else
                        thislist(k).bname1 = Convert.bid(splited{1});
                        thislist(k).fid1 = str2num(locY(split(splited{1},'-'),4));
                        thislist(k).bname2 = Convert.bid(splited{2});
                        thislist(k).fid2 = str2num(locY(split(splited{2},'-'),3));
                        thislist(k).concat1 = sprintf('%s-%02u',thislist(k).bname1,thislist(k).fid1);
                        thislist(k).concat2 = sprintf('%s-%02u',thislist(k).bname2,thislist(k).fid2);
                        thislist(k).fullrepla = sprintf('%s-%02u-%s-%02u',thislist(k).bname1,thislist(k).fid1,thislist(k).bname2,thislist(k).fid2);
                    end
                end

                % 最后填充norm song的信息
                splited = split(Replst_fnames{k_of_norm},{'-before-','-gapis-'});
                thislist(k_of_norm).bname1 = Convert.bid(splited{1});
                thislist(k_of_norm).fid1 =  thislist(k_of_norm+1).fid2-1;
                thislist(k_of_norm).bname2 = thislist(k_of_norm).bname1 ;
                thislist(k_of_norm).fid2 = thislist(k_of_norm+1).fid2;
                thislist(k_of_norm).concat1 = sprintf('%s-%02u',thislist(k_of_norm).bname1,thislist(k_of_norm).fid1);
                thislist(k_of_norm).concat2 = sprintf('%s-%02u',thislist(k_of_norm).bname2,thislist(k_of_norm).fid2);
                thislist(k_of_norm).fullrepla = sprintf('%s-%02u-%s-%02u',thislist(k_of_norm).bname1,thislist(k_of_norm).fid1,thislist(k_of_norm).bname2,thislist(k_of_norm).fid2);
                rp.replalist{a} = thislist;

            end




        end
        
        function  rp = judgeReplaResp(rp)
            % 为什么会有V1和V2，写的时候没有注意，导致有两个版本，有些许区别，但似乎V1更正确
            % 如果需要对比，可使用 online code difference comparision 的工具
            dbstop if error

            % evaluate the responsiveness of song replacements
            % Find out that the repla stimuliname correspond to how many natrual songs
            %             bnames = unique(cellfun(@(x)Convert.bid(x,2),[rp.list(replaids).stimuliname].','Uni',0));
            %
            %             findCorrespConID = @(x) intersect(find( ~cellfun(@isempty, regexp([rp.list.stimuliname].','norm') )),...
            %                 find( ~cellfun(@isempty, regexp([rp.list.stimuliname].',x) )) ); % function handle to get corresponding norm ids
            %
            %             findCorrespReplaID = @(x) intersect(find( ~cellfun(@isempty, regexp([rp.list.stimuliname].','Repla') )),...
            %                 find( ~cellfun(@isempty, regexp([rp.list.stimuliname].',x) )) ); % function handle to get corresponding norm ids


            for k = 1:length(rp.replalist) % 对于每组replacement来说！ 一组针对一个目标syllables

                sublist = rp.replalist{k};
                rids = find( ~cellfun(@isempty, regexp(cellstr({sublist.stimuliname}.'),'Repla|repla') ));
                nid = find( ~cellfun(@isempty, regexp(cellstr({sublist.stimuliname}.'),'norm') ));

                Ini_y_collect = [];
                allids = [rids;nid];
                for kk = 1:length(allids)

                    if allids(kk) == nid

                        convergentpoint = unique(Ini_y_collect); % 因为norm和本身没有分歧，所以他的分歧点要靠其他repla来确定
            
                        Ini_replay = unique(Ini_y_collect) + (find(sublist(nid).y(unique(Ini_y_collect):end),1)-1);  
                    else
                        %Ini_y是在y上的分歧点坐标，而Ini_replay是在replay上的分歧点坐标
                        [Ini_y,Ini_replay] = Neuron.findConergentPointBetwenNormAndRepla( sublist(nid).y,sublist(allids(kk)).y);
                        if Ini_replay==1 && Ini_y>1
                            %万一 Ini_replay为1，并且Ini_y大于1（说明不是norm song），说明原来song的syllable
                            %first恰好替换到了syllable second 前面，等于是deg songs，没有发生改变，
                            %此时，Ini_y和Ini_replay需要重新计算，一下是补丁算法：
                            Ini_y = unique(Ini_y_collect);
                            Ini_replay =length(sublist(allids(kk)).y) - (length(sublist(nid).y) - unique(Ini_y_collect));
                        end
                        
                        %    figure; subplot(211);Draw.spec(sublist(nid).y,32000);subplot(212); Draw.spec(sublist(allids(kk)).y,32000);
                        if Ini_y == 0
                            Ini_y = 1; % dangerous code
                        end
                        if isempty(Ini_y)|| Ini_y>length(sublist(nid).y)
                            continue  % another dangerous code
                        end

                        num_of_zeros = find(sublist(nid).y(Ini_y:end),1)-1; % 在分歧点后多长一段是0
                        convergentpoint = unique(Ini_y_collect);
                        Ini_y_collect(kk) = Ini_y; %+ num_of_zeros;  % y of the responsive parts
                        raw_Ini_replay = Ini_replay;
                        Ini_replay = Ini_replay + num_of_zeros;

                    end

                    try %不管当前是不是norm song都不影响的
                        sublist(allids(kk)).targety = sublist(allids(kk)).y(Ini_replay:Ini_replay + 0.2*sublist(allids(kk)).fs); % 截取200ms
                    catch % if the second index overceed 当剩下的distance小于200ms时
                        sublist(allids(kk)).targety = sublist(allids(kk)).y(Ini_replay:end); % 截取到 end
                        disp(' MEssage@Neuron.judgeReplaResp :Overceed');
                    end
               
                    if allids(kk) == nid % 如果当前是norm song的话
                        if length(sublist(nid).y(1:Ini_replay)) < 2 %如果太短，说明替换发生在第一位
                            sublist(allids(kk)).replaceparty = zeros(3000,1);%nan;
                        else
                            frags = Segment.segNormalizedSong(sublist(nid).y(1:Ini_replay),sublist(nid).fs);
                            last_initial = max([frags.initial].'); % 最末一个syllable的起始点
                            sublist(allids(kk)).replaceparty = sublist(allids(kk)).y(last_initial:convergentpoint); % 被替换的那一部分的y值
                            %
                        end

%                           figure;
%                         plot(sublist(allids(kk)).y(1:Ini_replay))

                    else
                        sublist(allids(kk)).replaceparty = sublist(allids(kk)).y(1:raw_Ini_replay); % 被替换的那一部分的y值

                    end
                  
                
                    sublist(allids(kk)).targetsptimes = Extract.sptimes_resetSP(sublist(allids(kk)).sptimes,Ini_replay/32000,(Ini_replay + 0.2*sublist(allids(kk)).fs)/sublist(allids(kk)).fs);
                    sublist(allids(kk)).replacepartsptimes = Extract.sptimes_resetSP(... % replacepart, sptimes
                        sublist(allids(kk)).sptimes,1,Ini_replay/sublist(allids(kk)).fs);
                    % corresponding pre (targety) data
                    sublist(allids(kk)).pretargety = zeros(length(sublist(allids(kk)).targety),1);%sublist(nid).prey(end - 0.1*32000:end); % 截取100ms
                    sublist(allids(kk)).pretargetsptimes = Extract.sptimes_resetSP(...
                        sublist(allids(kk)).presptimes,length(sublist(allids(kk)).pretargety)/sublist(allids(kk)).fs -0.2*sublist(allids(kk)).fs,length(sublist(allids(kk)).pretargety)/32000 );%length(sublist(rids(kk)).prey)/32000 -0.1*32000
                    %calculate and judege whether the neuron respond to the target area or not
                    temp = Neuron.UseMaxSdfToJudegeRespOrNot(sublist(allids(kk)).targety,sublist(allids(kk)).targetsptimes,...
                        sublist(allids(kk)).pretargety,sublist(allids(kk)).pretargetsptimes,sublist(allids(kk)).fs);
                    sublist(allids(kk)).label = temp.label;
                    sublist(allids(kk)).replaPvalue = temp.pvalue;
                    sublist(allids(kk)).maxsdf = temp.maxsdf;
                    %                     [sublist(rids(kk)).label,sublist(rids(kk)).replaPvalue] = Neuron.UseTtestToJudegeRespOrNot(...
                    %                         sublist(rids(kk)).targety,sublist(rids(kk)).targetsptimes,... % 查看对replaced song的反应是否明显大于对pre period的response
                    %                         sublist(rids(kk)).pretargety ,sublist(rids(kk)).pretargetsptimes ,sublist(rids(kk)).fs);
                    sublist(allids(kk)).targetfr = length(vertcat(sublist(allids(kk)).targetsptimes{:}))/200; % per milisecond



                end



                % do the same analysi for norm song
                rp.replalist{k} = sublist;

            end


        end

        function toCheck(rp)
            % 功能说明：把replalist的神经元反应结果作图，进而探明到底哪些replacement可以使神经元反应
            dbstop if error

            for index = 1:length(rp.replalist)

                current_replalist = rp.replalist{index}; % let us make replalist as an internal variable

                purerepla = current_replalist(2:end);
                purerepla = table2struct(sortrows(struct2table(purerepla),'targetfr','descend')); % rank by firing rate


                replafig = {};
                prefig = {};
                replacepartfig = {};
                for k = 1:length(purerepla)

                    figure('Position',[2189 255 560 1000])
                    Draw.two(purerepla(k).targety, 32000,purerepla(k).targetsptimes); % something wrong with the getReplallist
                    xlabel(sprintf(...
                        'mFR (target dur): %f, Result: % u ',purerepla(k).targetfr,purerepla(k).replaresp));%???? big problem here
                    replacepartfig {k} = getframe(gcf).cdata;
                    replafig{k} = getframe(gcf).cdata;

                    close(gcf)

                    figure('Position',[2189 255 560 1000])  %basline level prestimuli
                    Draw.two(purerepla(k).pretargety, 32000,purerepla(k).pretargetsptimes);
                    prefig{k} = getframe(gcf).cdata;
                    close(gcf)

                    figure('Position',[2189 255 400 1000])

                    Draw.two(purerepla(k).replaceparty, 32000,purerepla(k).replacepartsptimes);
                    replacepartfig {k} = getframe(gcf).cdata;
                    close(gcf)
                end

                sum_replafig = vertcat(replafig{:});
                sum_prefig = vertcat(prefig{:});
                sum_replacepartfig = vertcat(replacepartfig{:});
                sum_all = horzcat(sum_prefig,sum_replacepartfig,sum_replafig);

                %imwrite(sum_two,sprintf('%s_Sum_repla_fig.png',A.formated_name));

                imwrite(sum_all,sprintf('晋%s_Sum_preNreplaceNtarget.png',rp.formated_name));


            end


        end


        function weather_response_eliciting_syllables_are_similar

        end

        function replalist = does_neuron_respond_to_specific_transition(rp,transition_statics,syllable_theircatego)
            % 根据Catego Trans Probility，解答最终的问题，即是否不常见的combination不会使neuron反应，或者正相反
            replalist = rp.replalist;
            for k = 1:length(replalist)

                ids1 = find(~cellfun(@isempty, regexp({syllable_theircatego.unireplanames}.', replalist(k).concat1)));

                if ~isempty(ids1)
                    replalist(k).catego1 =  syllable_theircatego(ids1).catego;
                else
                    replalist(k).catego1 = 0;
                end

                ids2 = find(~cellfun(@isempty, regexp({syllable_theircatego.unireplanames}.', replalist(k).concat2)));
                if ~isempty(ids2)
                    replalist(k).catego2 =  syllable_theircatego(ids2).catego;
                else
                    replalist(k).catego2 = 0;
                end

                replalist(k).catego1catego2 = [replalist(k).catego1,replalist(k).catego2];

                transID = find(cellfun(@(x)isequal(replalist(k).catego1catego2,x),{transition_statics.firstsecond}.','Uni',1));
                if ~isempty(transID)
                    replalist(k).transrank = transition_statics(transID).rank;
                    replalist(k).count = transition_statics(transID).counts;
                    replalist(k).chance = transition_statics(transID).chance;
                else
                    replalist(k).transrank = 0;
                    replalist(k).count = 0;
                    replalist(k).chance = 0;
                end

            end

        end



        function a = getFinalFromReplalistAndTransstatics(a, transition_statics,classification_table)
            % 根据Catego Trans Probility，解答最终的问题，即是否不常见的combination不会使neuron反应，或者正相反
            %a.replalist = a.replalist;
            a.getReplalistFromA;

            for k = 1:length(a.replalist)

                ids1 = find(~cellfun(@isempty, regexp({classification_table.unifragnames}.', a.replalist(k).concat1)));

                if ~isempty(ids1)
                    a.replalist(k).catego1 =  classification_table(ids1).catego;
                else
                    a.replalist(k).catego1 = 0;
                end

                ids2 = find(~cellfun(@isempty, regexp({classification_table.unifragnames}.', a.replalist(k).concat2)));
                if ~isempty(ids2)
                    a.replalist(k).catego2 =  classification_table(ids2).catego;
                else
                    a.replalist(k).catego2 = 0;
                end

                a.replalist(k).catego1catego2 = [a.replalist(k).catego1,a.replalist(k).catego2];

                transID = find(cellfun(@(x)isequal(a.replalist(k).catego1catego2,x),{transition_statics.firstsecond}.','Uni',1));
                if ~isempty(transID)
                    a.replalist(k).transrank = transition_statics(transID).rank;
                    a.replalist(k).count = transition_statics(transID).counts;
                else
                    a.replalist(k).transrank = 0;
                    a.replalist(k).count = 0;
                end

            end

        end


    end

    methods % deprecated

        function  rp = Deprecated_judgeReplaRespV2(rp)
            dbstop if error

            % evaluate the responsiveness of song replacements
            replaids = find( ~cellfun(@isempty, regexp(cellstr({rp.list.stimuliname}.'),'Repla') ));

            if isempty(replaids)
                return
            end
            % Find out that the repla stimuliname correspond to how many natrual songs
            bnames = unique(cellfun(@(x)Convert.bid(x,2),[rp.list(replaids).stimuliname].','Uni',0));

            findCorrespConID = @(x) intersect(find( ~cellfun(@isempty, regexp([rp.list.stimuliname].','norm') )),...
                find( ~cellfun(@isempty, regexp([rp.list.stimuliname].',x) )) ); % function handle to get corresponding norm ids

            findCorrespReplaID = @(x) intersect(find( ~cellfun(@isempty, regexp([rp.list.stimuliname].','Repla') )),...
                find( ~cellfun(@isempty, regexp([rp.list.stimuliname].',x) )) ); % function handle to get corresponding norm ids

            for k = 1:length(bnames)
                nid = findCorrespConID(bnames{k}); % 有可能得到多个nid
                % 如果有多个nid
                if length(nid) > 1
                    nid = nid(1);
                    disp('Dangerous Code: @Neuron.judgeReplaResp');
                end
                rids = findCorrespReplaID(sprintf('(?<=before\\S+%s)',bnames{k}));

                for kk = 1:length(rids)
                    [Ini_y,Ini_replay] = Repla.findConergentPointBetwenNormAndRepla(...
                        rp.list(nid).y,...
                        rp.list(rids(kk)).y );

                    num_of_zeros = find(rp.list(nid).y(Ini_y:end),1)-1; % 在分歧点后多长一段是0
                    Ini_y = Ini_y + num_of_zeros;
                    Ini_replay = Ini_replay + num_of_zeros;
                    try
                        rp.list(rids(kk)).targety = rp.list(rids(kk)).y(Ini_replay:Ini_replay + 0.2*32000); % 截取200ms
                        rp.list(rids(kk)).replaceparty = rp.list(rids(kk)).y(1:Ini_replay ); % 被替换的那一部分的y值
                    catch % if the second index overceed
                        rp.list(rids(kk)).targety = rp.list(rids(kk)).y(Ini_replay:end); % 截取到 end
                        rp.list(rids(kk)).replaceparty = rp.list(rids(kk)).y(1:Ini_replay ); % 被替换的那一部分的y值
                        disp(' MEssage@Neuron.judgeReplaResp :Overceed');
                    end
                    rp.list(rids(kk)).targetsptimes = Extract.sptimes_resetSP(rp.list(rids(kk)).sptimes,Ini_replay/32000,(Ini_replay + 0.2*32000)/32000);
                    rp.list(rids(kk)).replacepartsptimes = Extract.sptimes_resetSP(... % replacepart, sptimes
                        rp.list(rids(kk)).sptimes,1,Ini_replay/32000);
                    % corresponding pre (targety) data
                    rp.list(rids(kk)).pretargety = zeros(length(rp.list(rids(kk)).targety),1);%rp.list(nid).prey(end - 0.1*32000:end); % 截取100ms
                    rp.list(rids(kk)).pretargetsptimes = Extract.sptimes_resetSP(...
                        rp.list(rids(kk)).presptimes,length(rp.list(rids(kk)).pretargety)/32000 -0.2*32000,length(rp.list(rids(kk)).pretargety)/32000 );%length(rp.list(rids(kk)).prey)/32000 -0.1*32000
                    %calculate and judege whether the neuron respond to the target area or not
                    [rp.list(rids(kk)).replaresp,rp.list(rids(kk)).replaPvalue] = Neuron.UseTtestToJudegeRespOrNot(...
                        rp.list(rids(kk)).targety,rp.list(rids(kk)).targetsptimes,...
                        rp.list(rids(kk)).pretargety ,rp.list(rids(kk)).pretargetsptimes ,32000);
                    rp.list(rids(kk)).targetfr = length(vertcat(rp.list(rids(kk)).targetsptimes{:}))/200; % per milisecond

                end


            end


        end


        function a = Deprecated_getReplalistFromAB(a)
            % 从 Neuron 文件中提取neuron对replas的反应，已经具体的replas的名称
            a.judgeReplaResp;
            replalist = a.list(find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'Repla'))));

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

            a.replalist = replalist;
            %             unique({replalist.concat1}.') % how many unique pre syllable
            %             unique({replalist.concat2}.') % how many unique post syllable
            %             unique(vertcat({replalist.concat1}.',{replalist.concat2}.'))  % how many syllable to be classified in total
            %             unique({replalist.fullrepla}.') % how many unique transition in neurons intotal?

        end


        function  a = Deprecated_judgeReplaRespV3(a)
            dbstop if error

            % evaluate the responsiveness of song replacements
            replaids = find( ~cellfun(@isempty, regexp([a.list.stimuliname].','Repla') ));

            % Find out that the repla stimuliname correspond to how many natrual songs
            bnames = unique(cellfun(@(x)Convert.bid(x,2),[a.list(replaids).stimuliname].','Uni',0));

            findCorrespConID = @(x) intersect(find( ~cellfun(@isempty, regexp([a.list.stimuliname].','norm') )),...
                find( ~cellfun(@isempty, regexp([a.list.stimuliname].',x) )) ); % function handle to get corresponding norm ids

            findCorrespReplaID = @(x) intersect(find( ~cellfun(@isempty, regexp([a.list.stimuliname].','Repla') )),...
                find( ~cellfun(@isempty, regexp([a.list.stimuliname].',x) )) ); % function handle to get corresponding norm ids

            for k = 1:length(bnames)
                nid = findCorrespConID(bnames{k}); % 有可能得到多个nid
                % 如果有多个nid
                if length(nid) > 1
                    nid = nid(1);
                    disp('Dangerous Code: @Neuron.judgeReplaResp');
                end
                rids = findCorrespReplaID(bnames{k});

                for kk = 1:length(rids)
                    [Ini_y,Ini_replay] = Neuron.findConergentPointBetwenNormAndRepla(...
                        a.list(nid).y,...
                        a.list(rids(kk)).y );
                    num_of_zeros = find(a.list(nid).y(Ini_y:end),1)-1; % 在分歧点后多长一段是0
                    Ini_y = Ini_y + num_of_zeros;
                    Ini_replay = Ini_replay + num_of_zeros;
                    try
                        a.list(rids(kk)).targety = a.list(rids(kk)).y(Ini_replay:Ini_replay + 0.2*32000); % 截取200ms
                        a.list(rids(kk)).replaceparty = a.list(rids(kk)).y(1:Ini_replay ); % 被替换的那一部分的y值
                    catch % if the second index overceed
                        a.list(rids(kk)).targety = a.list(rids(kk)).y(Ini_replay:end); % 截取200ms
                        disp(' MEssage@Neuron.judgeReplaResp :Overceed');
                    end
                    a.list(rids(kk)).targetsptimes = Extract.sptimes_resetSP(a.list(rids(kk)).sptimes,Ini_replay/32000,(Ini_replay + 0.2*32000)/32000);
                    a.list(rids(kk)).replacepartsptimes = Extract.sptimes_resetSP(... % replacepart, sptimes
                        a.list(rids(kk)).sptimes,1,Ini_replay/32000);
                    % corresponding pre (targety) data
                    a.list(rids(kk)).pretargety = zeros(length(a.list(rids(kk)).targety),1);%a.list(nid).prey(end - 0.1*32000:end); % 截取100ms
                    a.list(rids(kk)).pretargetsptimes = Extract.sptimes_resetSP(...
                        a.list(rids(kk)).presptimes,length(a.list(rids(kk)).pretargety)/32000 -0.2*32000,length(a.list(rids(kk)).pretargety)/32000 );%length(a.list(rids(kk)).prey)/32000 -0.1*32000
                    %calculate and judege whether the neuron respond to the target area or not
                    [a.list(rids(kk)).replaresp,a.list(rids(kk)).replaPvalue] = Neuron.UseTtestToJudegeRespOrNot(...
                        a.list(rids(kk)).targety,a.list(rids(kk)).targetsptimes,...
                        a.list(rids(kk)).pretargety ,a.list(rids(kk)).pretargetsptimes ,32000);
                    a.list(rids(kk)).targetfr = length(vertcat(a.list(rids(kk)).targetsptimes{:}))/200; % per milisecond
                end
                % a.list(nid).target = a.list(nid).y(Ini_y:Ini_y + 0.2*32000); % 截取200ms
                % a.list(nid).targetsptimes = Extract.sptimes_resetSP(a.list(nid).sptimes,Ini_y,Ini_y + 0.2*32000);
                %a.list(nid).targetfr = length(vertcat(a.list(nid).targetsptimes{:}))/0.2;


            end


        end


        function neu = Deprecated_getReplalistFromA(neu)
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
            unique({replalist.fullrepla}.') % how many unique transition in experiments intotal?

            neu.replalist = replalist;
        end


    end




    methods(Static)

        function finalreport = permutationTest_4responsive(A)

            multisims = A.repla.simatrix;
            localreport = {};

            for sim = 1:length(multisims)

                s = multisims{sim}; % s is the local simatrix, "单一的这个simatrix"

                if isempty(s.yesids) ||length(s.yesids) ==1||length(s.yesids)==length(s.fraglist)||length(s.yesids)==length(s.fraglist)-1 %如果根本没有syllables which trigger significant response
                    continue
                end

                distances = struct;

                %                 if length(s.yesids) == 0||length(s.yesids) == 1 %如果根本没有syllables which trigger significant response
                %                     finalreport = [];
                %                     return
                %                 end
                yespairs = nchoosek(s.yesids,2); % yes means response-eliciting
                fnames = fieldnames(s.matrix);

                for f = 1:length(fnames)

                    for k = 1:size(yespairs,1)
                        eval(['distances.yesdist.',fnames{f},'(k) = s.matrix.',fnames{f},'(yespairs(k,1),yespairs(k,2));']);
                    end

                    %总长度
                    eval(['distances.yessum.',fnames{f},' = sum(distances.yesdist.',fnames{f},');']);

                    for k = 1:length(s.compareids)
                        local_comparepairs = nchoosek(s.compareids{k},2);
                        local_distance = [];
                        for kk = 1:size(local_comparepairs,1)
                            eval(['local_distance(kk) = s.matrix.',fnames{f},'(local_comparepairs(kk,1),local_comparepairs(kk,2));']);
                        end
                        eval(['distances.comparedist.',fnames{f},'{k} = local_distance;']);
                        eval(['distances.comparesum.',fnames{f},'(k) = sum(local_distance);']);

                    end

                    % calculate all the distances between syllables
                    local_allpairs = nchoosek(1:length(s.fraglist),2);
                    local_alldistance = [];
                    for k = 1:length(local_allpairs)
                        eval(['local_alldistance(k) = s.matrix.',fnames{f},'(local_allpairs(k,1),local_allpairs(k,2));']);
                    end
                    eval(['distances.alldist.',fnames{f},' = local_alldistance;']);

                end




                I = {};

                summary = struct; % 总结p值大于0.05的值


                for k = 1:length(fnames)
                    %     subplot(2,length(fnames),k);
                    figure('Position',[2036 736 608 376],'Color','w');
                    summary(k).fname = fnames{k};

                    set(gca,'LineWidth',1.2);
                    set(gca,'TickLength',[0.001,0.001]);
                    thres = distances.yessum.(fnames{k});
                    allvalues = distances.comparesum.(fnames{k});
                    summary(k).pvalue = Draw.Permutation(allvalues,thres);
                    ax = gca;
                    ax.FontSize = 13;
                    xlabel('Sum of distances(normalized)','FontSize',15,'FontWeight','bold');
                    ylabel('Percentage','FontSize',15,'FontWeight','bold');
                    I{1,k} = getframe(gcf).cdata;
                    close(gcf)

                end

                for k = 1:length(fnames)
                    %     subplot(2,length(fnames),length(fnames) + k);
                    figure('Position',[2036 736 608 376],'Color','w');

                    set(gca,'LineWidth',1.2);
                    set(gca,'TickLength',[0.001,0.001]);
                    y1 = distances.yesdist.(fnames{k});
                    y2 = distances.alldist.(fnames{k});
                    Draw.swarmchart(y1,y2);
                    set(gca,'TickDir','out');
                    %ylim([-0.1,1])
                    xticks([1 2]);
                    set(gca,'xticklabels',{'Response-eliciting','All presented'});
                    ax = gca;
                    ax.FontSize = 13;
                    ax.XAxis.FontSize = 15;
                    set(get(gca, 'XAxis'), 'FontWeight', 'bold');
                    ylabel(sprintf('Distances in %s',fnames{k}),'fontsize',15,'FontWeight','bold');
                    set(ax,'defaultTextInterpreter','none');
                    I{2,k} = getframe(gcf).cdata;
                    close(gcf)

                end

                sumI = cell2mat(I);
                imwrite(sumI,sprintf('幽州ReplaElicit%s.png',A.info.formated_name));

                those_names = {summary.fname}.';
                localreport{sim}.neuronname = A.info.formated_name;
                localreport{sim}.simatrixid = sim;
                localreport{sim}.sigfeatures = setdiff(those_names([summary.pvalue].'>=0.95),'global');
                localreport{sim}.global = ismember('FM',those_names([summary.pvalue].'>=0.95));




            end

            finalreport = vertcat(localreport{:});

        end


        function finalreport = permutationTest_4not(A)

            multisims = A.repla.simatrix;

            localreport = {};

            for sim = 1:length(multisims)

                s = multisims{sim}; % s is the local simatrix, "单一的这个simatrix"

                noids = setdiff(1:length(s.fraglist),s.yesids);
                compareids = cell(500,1); % no ids 对应的compareids要重新取，因为 s.compareids是对应于yesids的
                for k = 1:500 % 姑且取500次
                    compareids{k} = randsample(1:length(s.fraglist),length(noids)); % 这个并不耗时
                end

                if isempty(noids)||length(noids) == 1  %如果根本没有syllables which trigger significant response
                    continue
                end

                distances = struct;
                nopairs = nchoosek(noids,2); % yes means response-eliciting
                fnames = fieldnames(s.matrix);

                for f = 1:length(fnames)

                    for k = 1:size(nopairs,1)
                        eval(['distances.nodist.',fnames{f},'(k) = s.matrix.',fnames{f},'(nopairs(k,1),nopairs(k,2));']);
                    end

                    %总长度
                    eval(['distances.nosum.',fnames{f},' = sum(distances.nodist.',fnames{f},');']);

                    for k = 1:length(compareids)
                        local_comparepairs = nchoosek(compareids{k},2);
                        local_distance = [];
                        for kk = 1:size(local_comparepairs,1)
                            eval(['local_distance(kk) = s.matrix.',fnames{f},'(local_comparepairs(kk,1),local_comparepairs(kk,2));']);
                        end
                        eval(['distances.comparedist.',fnames{f},'{k} = local_distance;']);
                        eval(['distances.comparesum.',fnames{f},'(k) = sum(local_distance);']);

                    end

                    % calculate all the distances between syllables
                    local_allpairs = nchoosek(1:length(s.fraglist),2);
                    local_alldistance = [];
                    for k = 1:length(local_allpairs)
                        eval(['local_alldistance(k) = s.matrix.',fnames{f},'(local_allpairs(k,1),local_allpairs(k,2));']);
                    end
                    eval(['distances.alldist.',fnames{f},' = local_alldistance;']);

                end




                I = {};

                summary = struct; % 总结p值大于0.05的值


                for k = 1:length(fnames)
                    %     subplot(2,length(fnames),k);
                    figure('Position',[2036 736 608 376],'Color','w');
                    summary(k).fname = fnames{k};

                    set(gca,'LineWidth',1.2);
                    set(gca,'TickLength',[0.001,0.001]);
                    thres = distances.nosum.(fnames{k});
                    allvalues = distances.comparesum.(fnames{k});
                    summary(k).pvalue = Draw.Permutation(allvalues,thres);
                    ax = gca;
                    ax.FontSize = 13;
                    xlabel('Sum of distances(normalized)','FontSize',15,'FontWeight','bold');
                    ylabel('Percentage','FontSize',15,'FontWeight','bold');
                    I{1,k} = getframe(gcf).cdata;
                    close(gcf)

                end

                for k = 1:length(fnames)
                    %     subplot(2,length(fnames),length(fnames) + k);
                    figure('Position',[2036 736 608 376],'Color','w');

                    set(gca,'LineWidth',1.2);
                    set(gca,'TickLength',[0.001,0.001]);
                    y1 = distances.nodist.(fnames{k});
                    y2 = distances.alldist.(fnames{k});
                    Draw.swarmchart(y1,y2);
                    set(gca,'TickDir','out');
                    %ylim([-0.1,1])
                    xticks([1 2]);
                    set(gca,'xticklabels',{'Response-prohibiting','All presented'});
                    ax = gca;
                    ax.FontSize = 13;
                    ax.XAxis.FontSize = 15;
                    set(get(gca, 'XAxis'), 'FontWeight', 'bold');
                    ylabel(sprintf('Distances in %s',fnames{k}),'fontsize',15,'FontWeight','bold');
                    set(ax,'defaultTextInterpreter','none');
                    I{2,k} = getframe(gcf).cdata;
                    close(gcf);

                end

                sumI = cell2mat(I);
                imwrite(sumI,sprintf('云州ReplaProhibit%s.png',A.info.formated_name));

                those_names = {summary.fname}.';
                localreport{sim}.neuronname = A.info.formated_name;
                localreport{sim}.simatrixid = sim;
                localreport{sim}.sigfeatures = setdiff(those_names([summary.pvalue].'>=0.95),'global');
                localreport{sim}.global = ismember('global',those_names([summary.pvalue].'>=0.95));


            end


            finalreport = vertcat(localreport{:});


        end

        function result = getNumTestedNumResponsive(A,bingyang)
            locallist = A.repla.replalist{1};%这里取{1}有些武断了
            normid = find(~cellfun(@isempty, regexp(cellstr({locallist.stimuliname}.'),'norm')));

            if sum(locallist(normid).replaceparty) == 0
                result = [];
                return

            end

            catego_norm = MetaStimuli.categorizeByBingyang(locallist(normid).replaceparty,32000, bingyang);

            purelocallist = locallist(find(cellfun(@isempty, regexp(cellstr({locallist.stimuliname}.'),'norm'))));% 去掉了norm项

            for k = 1:length(purelocallist)
                purelocallist(k).catego = MetaStimuli.categorizeByBingyang(purelocallist(k).replaceparty,32000,bingyang);
            end

          
            result.neuronname = A.info.formated_name;
            % 同组
            result.num_tested_simgroup = length(find([purelocallist.catego].'==catego_norm));
            result.num_responsive_simgroup = length(find([purelocallist.catego].'==catego_norm & [purelocallist.label].'==1));

            % 异组
            result.num_tested_difgroup = length(find([purelocallist.catego].'~=catego_norm));
            result.num_responsive_difgroup = length(find([purelocallist.catego].'~=catego_norm & [purelocallist.label].'==1));

            %所有组
            result.num_tested_allgroup = length(find([purelocallist.catego].'));
            result.num_responsive_allgroup = length(find([purelocallist.label].'==1));

            %比例
            result.ratio_sim = result.num_responsive_simgroup/result.num_tested_simgroup;
            result.ratio_dif = result.num_responsive_difgroup/result.num_tested_difgroup;
            result.ratio_all = result.num_responsive_allgroup/result.num_tested_allgroup;


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



        function namelist = get_all_elements_involved(rawdir_path)
            % 功能说明： 从dir里包含的所用音频文件里提取参与到replacement的音节的名称，称为namelist
            % from all ZP folder, know all the syllable transitions presented as repla stimuli
            filenames = cellfun(@Utl.fileparts,Extract.filesAllLevel(rawdir_path, '*.wav'),'Uni',0);
            repla_filenames = filenames(find(~cellfun(@isempty, regexp(filenames,'Repla|repla'))));
            %上面那行可以用来获取ZPdir里面的Individual replaments的filename，以期后续对其分类

            namelist = struct;
            locY = @(x,y) x{y}; % get element in location y
            for k = 1:length(repla_filenames)
                splited = split(repla_filenames{k},{'-before-','-gapis-'});
                namelist(k).bname1 = Convert.bid(splited{1});
                namelist(k).fid1 = str2num(locY(split(splited{1},'-'),4));
                namelist(k).bname2 = Convert.bid(splited{2});
                namelist(k).fid2 = str2num(locY(split(splited{2},'-'),3));
                namelist(k).concat1 = sprintf('%s-%02u',namelist(k).bname1,namelist(k).fid1);
                namelist(k).concat2 = sprintf('%s-%02u',namelist(k).bname2,namelist(k).fid2);
                namelist(k).fullrepla = sprintf('%s-%02u-%s-%02u',namelist(k).bname1,namelist(k).fid1,namelist(k).bname2,namelist(k).fid2);
            end

            [~,uniqueids,~] = unique({namelist.fullrepla}.'); % 去除重复的rows
            namelist = namelist(uniqueids);

        end


        function metadata = namelist2meta(namelist,refer_list)
            % @功能说明： convert namelist to standard metadata, based on refer_list
            unireplanames = unique(vertcat({namelist.concat1}.',{namelist.concat2}.')); % unique replanames in ZPdirs of repla stimuli
            %reference_eleinf = Utl.load('con_and_allspe_eleinf.mat');

            bnames = cellfun(@Convert.bid,{refer_list.fullname}.','UniformOutput',false);
            func1 = @(x) regexp(x,'(?<=-)\d+','match'); % reference_eleinf(1).fullname
            func2 = @(x) x{1};
            replaids = cellfun(@(x) func2(func1(x)),{refer_list.fullname}.','UniformOutput',false) ;
            bname_replaid = cellfun(@(x,y) horzcat(x,'-',y), bnames,replaids,'Uni',0); % reference_eleinf

            metadata = struct;
            for k = 1:length(unireplanames)
                ids = find(strcmp(unireplanames{k},bname_replaid));
                metadata(k).unireplanames = unireplanames{k};
                if ~isempty(ids)
                    metadata(k).y = refer_list(ids).y;
                    metadata(k).fs= refer_list(ids).fs;
                    metadata(k).fullname = refer_list(ids).fullname;
                else
                    metadata(k).y = [];
                    metadata(k).fs= [];
                    metadata(k).fullname = [];
                end
            end

            metadata = metadata(~cellfun(@isempty,{metadata.fullname}.'));

        end


        function [transition_statics,relation_map] = get_trans_statics(sequential_song_meta)
            % @功能说明： 生成transition statics
            % remove those with catego num as 0
            idsnot0 = find([sequential_song_meta.catego].' ~= 0);
            sequential_song_meta = sequential_song_meta(idsnot0);
            disp('Remove those with catego num as 0');
            splited = MetaStimuli.split(sequential_song_meta); % 把大的struct分解成对单个song的
            collect_transfering = {};
            for k = 1: length(splited) % 对每个song计算transfering的计数
                transfering = {};
                parfor kk = 1: length(splited{k})-1
                    transfering{kk} = [splited{k}(kk).catego,splited{k}(kk+1).catego];
                end
                collect_transfering{k} = transfering;
            end
            all_transfering = horzcat(collect_transfering{:}); % Concat这些song里面所有的transfering

            % 集合所有存在的transition事件
            transitioncollect = struct;
            for k = 1:length(all_transfering)
                transitioncollect(k).first = all_transfering{k}(1);
                transitioncollect(k).second = all_transfering{k}(2);
            end

            % 对每对可能的catego transition统计对应的频数
            all_categos = unique([sequential_song_meta .catego].');
            relation_map = zeros(length(all_categos),length(all_categos)); % initialization
            transition_statics = struct;

            count = 0;
            for k = 1:length(all_categos)
                for kk = 1:length(all_categos)
                    count = count + 1;
                    relation_map(k,kk) = length(intersect(...
                        find([transitioncollect.first] == all_categos(k)),...
                        find([transitioncollect.second] == all_categos(kk))));

                    transition_statics(count).firstsyllable = k;
                    transition_statics(count).secondsyllable = kk;
                    transition_statics(count).firstsecond = [k,kk];
                    transition_statics(count).counts = relation_map(k,kk);
                end
            end

            all_counts = sum(sum(relation_map));

            % sort the transition_statics and assign the rank
            transition_statics = table2struct(sortrows(struct2table(transition_statics),'counts','descend'));
            for k = 1:length(transition_statics)
                transition_statics(k).rank = k; % rank越小， transition越常见
                transition_statics(k).chance =  transition_statics(k).counts/all_counts;
            end

        end


        function drawreplaSpec(all_eleinf,mode)
            % A function to Draw spectrograms of already categorized replaments

            fs = 32000;
            img = {};

            switch mode
                case 'by categos' % by number of categories
                    categos = unique([all_eleinf.catego].'); % 计算categos的数量
                    concated_y = {};
                    for k = 1: length(categos)

                        ids = find([all_eleinf.catego].' == categos(k));
                        gap = zeros(0.3*fs,1);
                        categoys = cellfun(@(x) vertcat(x,gap),{all_eleinf(ids).y}.','Uni',0);
                        concated_y{k} = vertcat(categoys{:});

                    end

                    max_y_length = max(cellfun(@length,concated_y));


                    for k = 1: length(categos)


                        concated_y{k} = [concated_y{k};zeros(max_y_length - length(concated_y{k}),1)];

                        figure('Position',[8 310 6000 716]);
                        Draw.spec(concated_y{k},fs);
                        img{k} = getframe(gcf).cdata;
                        close(gcf);
                    end


                case 'by rows' % by number of rows

                    [~,index] = sortrows([all_eleinf.catego].'); all_eleinf = all_eleinf(index); clear index % sort all_eleinf
                    column_num = 10; % I should use fixed syllables for replacement
                    row_num = round(length(all_eleinf)/column_num);

                    for k = 1: row_num

                        ids = column_num*(k-1) + 1: min(column_num*k,length(all_eleinf));
                        gap = zeros(0.3*fs,1);
                        categoys = cellfun(@(x) vertcat(x,gap),{all_eleinf(ids).y}.','Uni',0);
                        concated_y = vertcat(categoys{:});

                        figure('Position',[14 829 1902 280]);
                        Draw.spec(concated_y,fs);
                        img{k} = getframe(gcf).cdata;
                        close(gcf);

                    end
            end

            whole_img = vertcat(img{:});
            imwrite(whole_img,'New_replas_separated_by_replas.png')
        end


        % from csv sort figures


        function csv2Sortfig(csv_path,figdir)
            %csvpath = "E:\My_Analysis_by_date\0925_CheckNeuralResponseToReplas\FigOTE.csv";

            otestruct = table2struct(readtable(csv_path))
            %figdir = "E:\My_Analysis_by_date\0925_CheckNeuralResponseToReplas\Figs_NewCon_OTE_For_TransStatastics\Figs_OTE";
            [~,figdirname,~]= fileparts(figdir)
            which_categos = unique([otestruct.catego].');
            outdir = sprintf('Sorted_%s_figs',figdirname);

            for k = 1:length(which_categos)
                subdirs{k} = sprintf('%s\\%u',outdir,which_categos(k));
                mkdir(subdirs{k});
            end

            if min(which_categos) == 0
                method = 0;
            else
                methods = 1;
            end


            for k = 1:length(otestruct)
                if method == 0
                    targetdir = subdirs{otestruct(k).catego + 1};

                else
                    targetdir = subdirs{otestruct(k).catego};
                end
                source = fullfile(figdir,otestruct(k).name);
                destiney = fullfile(targetdir,otestruct(k).name);
                copyfile(source,destiney);
            end

        end

        % from sorted fig regenrate csv
        function figdir2csv(figdir)
            %figdir = "E:\My_Analysis_by_date\0925_CheckNeuralResponseToReplas\Sorted_ote_repla_figs";

            subdirs = Extract.folder(figdir);

            figstruct = struct;
            count = 0;
            for k = 1:length(subdirs)

                [~,categoname,~] = fileparts(subdirs{k});
                categonum = int16(str2double(categoname));
                files = Extract.filename(subdirs{k},'*.png');
                for kk = 1:length(files)
                    count = count + 1;
                    [~,name,ext] = fileparts(files{kk});
                    figstruct(count).name = strcat(name, ext);
                    figstruct(count).catego = categonum;
                end

            end

            figtable = struct2table(figstruct);
            [~,figdirname,~] = fileparts(figdir);

            writetable(figtable,sprintf('%s.csv',figdirname));
        end


        function fraglist = to2ndAcousticSpace(rp)
            % compare the similarity between two-syllable sequence which
            % triggered the neuron to respond and those did not


            %             ids = find(~cellfun(@isempty, regexp({rp.list.stimuliname}.','repla'))); % find all frags
            %
            %             fraglist = neu.list(ids);

            replalist = rp.replalist;

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


    end
end

