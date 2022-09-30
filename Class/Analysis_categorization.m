classdef Analysis_categorization < handle
    %UNTITLED6 Summary of this class goes here
    properties
        replalist
    end
    
    methods (Access = protected)
        
        function a = getReplalistFromAB(a)
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
        
        
        function  a = judgeReplaResp(a)
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
        
    end
    
    
    methods
        
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
    
    
    
    
    methods(Static)
        
        
        
    end
end

