classdef Consistency
    %查看Neuron在经历了多个stimuli set之后，其反应是否持续到了最后

    properties
        songids
        fragids
        replaids
        sublist
        formated_name
    end

    methods
        function cs = Consistency(list)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here

            targets = Consistency.knowTarget(list);%{1};

            if isempty(targets)
                return
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

       
         
        end


        function drawThree(cs)

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
            imwrite(img,sprintf('凉州Consistency_%s.png',cs.formated_name));


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