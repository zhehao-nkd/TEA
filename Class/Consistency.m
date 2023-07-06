classdef Consistency
    %查看Neuron在经历了多个stimuli set之后，其反应是否持续到了最后

    properties
        songids
        fragids
        replaids
        sublist


        %是否对song反应
        RepSong18
        RepSong19
        numRepSong18
        numRepSong19

        frag80tested
        numRespFrag80
        degtested % 是否测了deg,如果根本没测,0,测了但同一个song不trigger response:0.5,测了且song trigger 反应, 1
        replatested % 是否测了repla,如果根本没测,0,测了但同一个song不trigger response:0.5,测了且song trigger 反应, 1
        comparable

        resp2SameSongFrags
        numResp2SameSongFrags
        formated_name
    end

    methods
        function cs = Consistency(list)
        

             %  这个Neuron是否对某个Conspecific Song起反应？
            separated ={};
            unique_fids = unique({list.Fid}.');
            for k = 1:length(unique_fids) %有多少file，对应多少experiment
                temp = list(find(~cellfun(@isempty, regexp(cellstr({list.Fid}.'),unique_fids{k}))));
                separated{k}= temp(find(~cellfun(@isempty, regexp(cellstr({temp.stimuliname}.'),'norm'))));
            end
            [~,maxloc] = max(cellfun(@length,separated));
            temp_normlist = separated{maxloc};
            temp_normlist = Song.judgeConResp(temp_normlist);

            names19 = {'B346','B512','B521','B554','B606','G429','G506','G518','G548','G573',...
                'G578','O331','O507','O509','O540','Y515','Y606','Y616','G699'};

            normlist19 = temp_normlist(find(~cellfun(@isempty,regexp(cellstr({temp_normlist.stimuliname}.'),strjoin(names19,'|')))));
            %disp(k)
            cs.numRepSong19 = sum([normlist19.label].');
            if sum([temp_normlist.label].') > 0
                cs.RepSong19 = 1;  % 似乎 make sense
            else
                cs.RepSong19 = 0;
            end

            names18 = {'B346','B512','B521','B554','B606','G429','G506','G518','G573',...
                'G578','O331','O507','O509','O540','Y515','Y606','Y616','G699'};

            normlist18 = temp_normlist(find(~cellfun(@isempty,regexp(cellstr({temp_normlist.stimuliname}.'),strjoin(names18,'|')))));
            %disp(k)
            cs.numRepSong18 = sum([normlist18.label].');
            if sum([temp_normlist.label].') > 0
                cs.RepSong18 = 1; % 
            else
                cs.RepSong18 = 0;
            end

         

            % 80 frags就是那些既含有Type又含有Frag这两个关键词的
            criteria1 = find(~cellfun(@isempty,regexp(cellstr({list.stimuliname}.'),'Type')));
            criteria2 = find(~cellfun(@isempty,regexp(cellstr({list.stimuliname}.'),'Frag')));
            frag80 =  Frag.judgeFragResp(list(intersect(criteria1,criteria2)));
            if  length(frag80) >= 80
                cs.frag80tested = 1;
            else
                cs.frag80tested = 0;
            end
            cs.numRespFrag80 = length(find([frag80.label].'==1));

            targets = Consistency.knowTarget(list);%{1};

            if isempty(targets)
                return
            end


             % 这个Neuron 是否被测试了Degressive song？ 如果测了跟没测一样的话不计入测试
            degids = find(~cellfun(@isempty, regexp(cellstr({list.stimuliname}.'),'Deg|deg') ));
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
                cs.degtested = 1;
            else
                cs.degtested = 0;
            end

            thetarget = targets{1};
         
            songname = regexp(convertCharsToStrings(thetarget),'[OBGYR]\d+','match');
            fragname = regexp(convertCharsToStrings(thetarget),'[OBGYR]\d+-\d+','match');
            replaname2nd = fragname;
            parts = split(fragname,'-');
            replaname1st = strjoin({parts{1},num2str(str2num(parts{2})-1,'%02u')},'-');

            norm_temp =   find(~cellfun(@isempty,regexp(cellstr({list.stimuliname}.'),'norm')));
            frag_temp =   find(~cellfun(@isempty,regexp(cellstr({list.stimuliname}.'),'frag|Frag')));
            repla_temp =   find(~cellfun(@isempty,regexp(cellstr({list.stimuliname}.'),'repla|Repla')));
            
            songname_temp = find(~cellfun(@isempty,regexp(cellstr({list.stimuliname}.'),songname)));

            fragname_temp = find(~cellfun(@isempty,regexp(cellstr({list.stimuliname}.'),fragname)));

            replaname2nd_temp = find(~cellfun(@isempty,regexp(cellstr({list.stimuliname}.'),replaname2nd)));

            replaname1st_temp = find(~cellfun(@isempty,regexp(cellstr({list.stimuliname}.'),replaname1st)));

            cs.songids = intersect(norm_temp,songname_temp);

            cs.fragids = intersect(frag_temp,fragname_temp);

            cs.replaids = mintersect(replaname1st_temp, replaname2nd_temp);

            allids = [cs.songids;cs.fragids;cs.replaids];
            cs.sublist = list(allids);





            %是否被测了 replas + 80 frags
            consis_songs = cs.sublist(find(~cellfun(@isempty,regexp(cellstr({cs.sublist.stimuliname}.'),'norm'))));
            all_repla_ids = find(~cellfun(@isempty,regexp(cellstr({list.stimuliname}.'),'Repla')));
            if isempty(all_repla_ids)
                cs.replatested = 0;

            else

                repla_Fid = cs.sublist(find(~cellfun(@isempty,regexp(cellstr({cs.sublist.stimuliname}.'),'Repla')))).Fid;
                repla_corresp_norm = consis_songs(find(~cellfun(@isempty,regexp(cellstr({consis_songs.Fid}.'),repla_Fid))));
                if isempty(repla_corresp_norm)
                    cs.replatested = 0.4;
                else
                    repla_answer = Song.judgeConResp(repla_corresp_norm).label;
                    % figure; Draw.three(repla_corresp_norm.plty,repla_corresp_norm.fs,repla_corresp_norm.pltsptimes);

                    if repla_answer == 1
                        cs.replatested = 1;
                    else
                        cs.replatested = 0.5;
                    end

                end
            end


            criteria1 = find(~cellfun(@isempty,regexp(cellstr({list.stimuliname}.'),'Type')));
            criteria2 = find(~cellfun(@isempty,regexp(cellstr({list.stimuliname}.'),'Frag')));
            frag80 =  Frag.judgeFragResp(list(intersect(criteria1,criteria2)));
           
            if strcmp(repla_Fid, frag80(1).Fid)
                cs.comparable = 1;
            else
                cs.comparable = 0;
            end


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
                cs.resp2SameSongFrags = 1;
            else
                cs.resp2SameSongFrags = 0;
            end

            cs.numResp2SameSongFrags = length(find([judged_frags.label].'==1));

       
         
        end


        function drawThree(cs,outdir)

            if ~exist('outdir','var')
                outdir = '.';
            end

            if isempty(cs.sublist)
                return
            end

            I = {};
            for k = 1:length(cs.sublist)
                figure('Position',PM.size1,'color','w');
                Draw.three(cs.sublist(k).plty,cs.sublist(k).fs,cs.sublist(k).pltsptimes);
                xlabel(sprintf('%s-%s',cs.sublist(k).stimuliname,cs.sublist(k).Fid));

                frame = getframe(gcf);
                I{k} = frame.cdata;
                close(gcf);
            end

            % draw blank white
            lieshu = 5;
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
            imwrite(img,sprintf('%s\\凉州Consistency_%s.png',outdir,cs.formated_name));


        end

    end

    methods(Static)

        function targets = knowTarget(list)
            %找到哪个是目标syllable，目标syllable是生成deg stimuli和repla stimuli的基础
            %目前的版本只是通过repla判断，之后有时间也应该写通过deg判断的方法
            % 2022.12.21 追加通过degressive song 判断
            %当然，如果连degressive song 也没有，那就根本谈不上有targets
            %而为了通过degressive song 判断，首先需要解决的是JudgeDegResp


            % target
            targetidx = find(~cellfun(@isempty, regexp(cellstr({list.stimuliname}.'),'Repla')));
            targetlist = list(targetidx);
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

    end
end