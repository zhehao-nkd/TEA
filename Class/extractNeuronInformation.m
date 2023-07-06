classdef extractNeuronInformation
    %这个Class的作用是从mat文件中提取出想要的信息，
    % input 应当是单一的matfile，以及想要的信息的名称，比如wavelength，或者是response218songs

    % 提取这些信息可能需要load不同的variable，比如wavelength可以从info中获得？ 真的需要这个class吗？？？

    properties
        Property1
    end

    methods(Static)

        function result = getAll(list)
            dbstop if error

            

            names19 = {'B346','B512','B521','B554','B606','G429','G506','G518','G548','G573',...
                'G578','O331','O507','O509','O540','Y515','Y606','Y616','G699'};


            wnslist = list(~cellfun(@isempty,regexp(cellstr({list.stimuliname}.'),'WNS')));
            wns_samefile_Fid = wnslist.Fid;
            samefile_ids = find(~cellfun(@isempty, regexp(cellstr({list.Fid}.'),wns_samefile_Fid)));
            all_norm_ids = find(~cellfun(@isempty, regexp(cellstr({list.stimuliname}.'),'norm')));
            temp_normlist = list(intersect(samefile_ids,all_norm_ids));
            temp_normlist = Song.judgeConResp(temp_normlist);




            normlist19 = temp_normlist(find(~cellfun(@isempty,regexp(cellstr({temp_normlist.stimuliname}.'),strjoin(names19,'|')))));
            %disp(k)
            result.numRepSong19 = sum([normlist19.label].');
            if sum([temp_normlist.label].') > 0
                result.RepSong19 = 1;  % 似乎 make sense
            else
                result.RepSong19 = 0;
            end

            names18 = {'B346','B512','B521','B554','B606','G429','G506','G518','G573',...
                'G578','O331','O507','O509','O540','Y515','Y606','Y616','G699'};

            normlist18 = temp_normlist(find(~cellfun(@isempty,regexp(cellstr({temp_normlist.stimuliname}.'),strjoin(names18,'|')))));
            %disp(k)
            result.numRepSong18 = sum([normlist18.label].');
            if sum([temp_normlist.label].') > 0
                result.RepSong18 = 1; % 
            else
                result.RepSong18 = 0;
            end

         

            % 80 frags就是那些既含有Type又含有Frag这两个关键词的
            criteria1 = find(~cellfun(@isempty,regexp(cellstr({list.stimuliname}.'),'Type')));
            criteria2 = find(~cellfun(@isempty,regexp(cellstr({list.stimuliname}.'),'Frag')));
            frag80 =  Frag.judgeFragResp(list(intersect(criteria1,criteria2)));
            if  length(frag80) >= 80
                result.frag80tested = 1;
            else
                result.frag80tested = 0;
            end
            result.numRespFrag80 = length(find([frag80.label].'==1));

         
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% About Degs
            
             % 这个Neuron 是否被测试了Degressive song？ 如果测了跟没测一样的话不计入测试
            degids = find(~cellfun(@isempty, regexp(cellstr({list.stimuliname}.'),'Deg|deg') ));

            if isempty(degids)

                result.degtested = 0;
            else
                list_contain_degs = list(degids);
                deg_birdids = cellfun(@Convert.bid, cellstr({list_contain_degs.stimuliname}.'),'Uni',0);
                unique_deg_birdids = unique(deg_birdids);
                unpack = @(x) x{1};
                summer = {};
                for k = 1:length(unique_deg_birdids)
                    correspids = strcmp(unique_deg_birdids{k},deg_birdids); % 找到有这个birdid的所有deg的ids,这ids是关于list_contain_degs的
                    deglist{k} = list_contain_degs(correspids);

                    corresp_fid = unpack(unique({deglist{k}.Fid}.')); % 找到deg stimuli对应的fid
                    corresp_normids = mintersect( find(~cellfun(@isempty, regexp(cellstr({list.stimuliname}.'),unique_deg_birdids{k}))),...
                        find(~cellfun(@isempty, regexp(cellstr({list.stimuliname}.'),'norm'))),...
                        find(~cellfun(@isempty, regexp(cellstr({list.Fid}.'),corresp_fid))) );
                    %                 if isempty(corresp_normids)
                    temp = intersect( find(~cellfun(@isempty, regexp(cellstr({list.stimuliname}.'),unique_deg_birdids{k}))),...
                        find(~cellfun(@isempty, regexp(cellstr({list.stimuliname}.'),'norm'))));
                    summer{k} = list(temp(1)); % 取第一位
                    deglist{k} = horzcat(list(corresp_normids),deglist{k});
                end
                corresp_normlist = vertcat(summer{:});
                if length(corresp_normlist) > 1
                    warning('Consistency Line 56 出现了错误！');
                end
                deg_answer = Song.judgeConResp(corresp_normlist).label;
                if deg_answer == 1
                    result.degtested = 1;
                else
                    result.degtested = 0.5;
                end

            end
       
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%About samseongfrags

             % 是否对任意一个single syllable of the song 反应？  这个Neuron 是否对单个elements （same song） 反应？如果是，有多少个？
            raw_samesongfrags = list(find(~cellfun(@isempty,regexp( cellstr({list.stimuliname}.'),'samesong'))));
            [counts, groupnames] = groupcounts({raw_samesongfrags.Fid}.');
            [maxvalue, maxid] = max(counts);
            if length(counts) > 1
                warning('Dangerous!!!!!!');
                %pause
            end
            selected_fid = groupnames(maxid);
            samesongfrags = raw_samesongfrags(find(~cellfun(@isempty,regexp( cellstr({raw_samesongfrags.Fid}.'),selected_fid))));
            judged_frags = Frag.judgeFragResp(samesongfrags);

            if any([judged_frags.label].'==1)
                result.resp2SameSongFrags = 1;
            else
                result.resp2SameSongFrags = 0;
            end

            result.numResp2SameSongFrags = length(find([judged_frags.label].'==1));



            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            targets = Consistency.knowTarget(list);%{1};
        
            

            if isempty(targets)
                result.replatested = 0;
            else
                

%                 songname = regexp(convertCharsToStrings(thetarget),'[OBGYR]\d+','match');
%                 fragname = regexp(convertCharsToStrings(thetarget),'[OBGYR]\d+-\d+','match');
%                 replaname2nd = fragname;
%                 parts = split(fragname,'-');
%                 replaname1st = strjoin({parts{1},num2str(str2num(parts{2})-1,'%02u')},'-');
% 
%                 norm_temp =   find(~cellfun(@isempty,regexp(cellstr({list.stimuliname}.'),'norm')));
%                 frag_temp =   find(~cellfun(@isempty,regexp(cellstr({list.stimuliname}.'),'frag|Frag')));
%                 repla_temp =   find(~cellfun(@isempty,regexp(cellstr({list.stimuliname}.'),'repla|Repla')));
% 
%                 songname_temp = find(~cellfun(@isempty,regexp(cellstr({list.stimuliname}.'),songname)));
%                 fragname_temp = find(~cellfun(@isempty,regexp(cellstr({list.stimuliname}.'),fragname)));
%                 replaname2nd_temp = find(~cellfun(@isempty,regexp(cellstr({list.stimuliname}.'),replaname2nd)));
%                 replaname1st_temp = find(~cellfun(@isempty,regexp(cellstr({list.stimuliname}.'),replaname1st)));
%                 result.songids = intersect(norm_temp,songname_temp);
%                 result.fragids = intersect(frag_temp,fragname_temp);
%                 result.replaids = mintersect(replaname1st_temp, replaname2nd_temp);
%                 allids = [result.songids;result.fragids;result.replaids];
%                 result.sublist = list(allids);


                %是否被测了 replas + 80 frags
              
                all_repla_ids = find(~cellfun(@isempty,regexp(cellstr({list.stimuliname}.'),'Repla')));
                if isempty(all_repla_ids)
                    result.replatested = 0;
                else

                    
                    repla_ids =   find(~cellfun(@isempty,regexp(cellstr({list.stimuliname}.'),'repla|Repla')));
                    repla_Fid = unique(unique({list(repla_ids).Fid}.'));
                    norm_ids =   find(~cellfun(@isempty,regexp(cellstr({list.stimuliname}.'),'norm')));
                    repla_samefile_ids =find(~cellfun(@isempty,regexp(cellstr({list.Fid}.'),repla_Fid)));
                    repla_corresp_norm = list(intersect(norm_ids,repla_samefile_ids));
                    if isempty(repla_corresp_norm)
                        result.replatested = 0.4;
                    else
                        repla_answer = Song.judgeConResp(repla_corresp_norm).label;
                        % figure; Draw.three(repla_corresp_norm.plty,repla_corresp_norm.fs,repla_corresp_norm.pltsptimes);

                        if repla_answer == 1
                            result.replatested = 1;
                        else
                            result.replatested = 0.5;
                        end

                    end
                end

                

            end

          


           

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            criteria1 = find(~cellfun(@isempty,regexp(cellstr({list.stimuliname}.'),'Type')));
            criteria2 = find(~cellfun(@isempty,regexp(cellstr({list.stimuliname}.'),'Frag')));
            frag80 =  Frag.judgeFragResp(list(intersect(criteria1,criteria2)));
            if  length(frag80) >= 80
                result.frag80tested = 1;
            else
                result.frag80tested = 0;
            end
            result.numRespFrag80 = length(find([frag80.label].'==1));



            %是否被测了 replas + 80 frags
            %consis_songs = result.sublist(find(~cellfun(@isempty,regexp(cellstr({result.sublist.stimuliname}.'),'norm'))));
            all_repla_ids = find(~cellfun(@isempty,regexp(cellstr({list.stimuliname}.'),'Repla')));
            if isempty(all_repla_ids)
                result.replatested = 0;
                result.comparable = 0;
                result.resp2target_replafile = 0;
                return

            end



            repla_Fid = list(unique(all_repla_ids)).Fid;
            repla_samefile_ids = find(~cellfun(@isempty,regexp(cellstr({list.Fid}.'),repla_Fid)))
            norm_ids = find(~cellfun(@isempty,regexp(cellstr({list.stimuliname}.'),'norm')));
            corresp_norm_ids = intersect(repla_samefile_ids,norm_ids);

            repla_corresp_norm = list(corresp_norm_ids);
            if isempty(repla_corresp_norm)
                result.replatested = 0.4;
            else
                repla_answer = Song.judgeConResp(repla_corresp_norm).label;
                % figure; Draw.three(repla_corresp_norm.plty,repla_corresp_norm.fs,repla_corresp_norm.pltsptimes);

                if repla_answer == 1
                    result.replatested = 1;
                else
                    result.replatested = 0.5;
                end

            end


           
            if strcmp(repla_Fid, frag80(1).Fid)
                result.comparable = 1;
            else
                result.comparable = 0;
            end


            if  result.comparable == 1

                temp = targets{1};
                thetarget = regexp(convertCharsToStrings(temp),'[OGBYR]\d+-\d+','match');
                repla_ids =   find(~cellfun(@isempty,regexp(cellstr({list.stimuliname}.'),'repla|Repla')));
                repla_Fid = unique(unique({list(repla_ids).Fid}.'));
                repla_samefile_ids = find(~cellfun(@isempty,regexp(cellstr({list.Fid}.'),repla_Fid)));
                fragids = find(~cellfun(@isempty,regexp(cellstr({list.stimuliname}.'),'frag|Frag')));
                targetids = find(~cellfun(@isempty,regexp(cellstr({list.stimuliname}.'),thetarget)));

                target_syllable_response = list(mintersect(repla_samefile_ids,targetids,fragids));
                target_syllable_response = Frag.judgeFragResp(target_syllable_response);
                if isempty(target_syllable_response)

                    result.resp2target_replafile = 0;
                elseif length(target_syllable_response) > 1
                    if any([target_syllable_response.label].' ==1)
                        result.resp2target_replafile = 1;

                    else
                        result.resp2target_replafile = 0.5;
                    end

                elseif target_syllable_response.label ==1
                    result.resp2target_replafile = 1;
                else
                    result.resp2target_replafile = 0.5;
                end

            else

                result.resp2target_replafile = 0;

            end

            

 


   

        end

    end


end
