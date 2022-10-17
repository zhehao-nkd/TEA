classdef Frag
    
    properties
        fraglist
        list
        formated_name
    end
    methods
        function fg = Frag(list)
            fg.list = list;
            fg.fraglist = list( find(~cellfun(@isempty, regexp(cellstr({fg.list.stimuliname}.'),'Frag|syl'))));
        end
    
    
    function fg = judgeFragResp(fg)
            % 判断对frag 是否反应，通过自己定义的复杂的机制
            ids = find(~cellfun(@isempty, regexp(cellstr({fg.list.stimuliname}.'),'Frag|syl'))); % find all frags
            
            % ’syl'可以兼容旧的stimuli命名规则
            DUR = 0.2 ;% 200ms
            
            
            for n = 1: length(ids)
                thisi = ids(n);
                post_sptimes = Extract.sptimes(fg.list(thisi).rawsptimes, fg.list( thisi).zpt, fg.list( thisi).zpt + DUR);
                pre_sptimes = Extract.sptimes(fg.list(thisi).rawsptimes, fg.list(thisi).zpt-DUR, fg.list(thisi).zpt );
                post_mfr = length(vertcat(post_sptimes{:}))/DUR;
                pre_mfr = length(vertcat(pre_sptimes{:}))/DUR; % mfr: mean firing rate
                fg.list(thisi).rs = post_mfr - pre_mfr; % response strength
                
                tempsum = Cal.psth_frag(fg.list(thisi).plty,fg.list(thisi).fs,fg.list(thisi).pltsptimes);
                halfsum = sum(tempsum(end/2:end));
                fullsum = sum(tempsum);
                maxvalue = max(Cal.psth_frag(fg.list(thisi).plty,fg.list(thisi).fs,fg.list(thisi).pltsptimes));
                fg.list(thisi).maxvalue = maxvalue;
                fg.list(thisi).halfsum = halfsum;
                fg.list(thisi).fullsum = fullsum;
                %if maxvalue > 6 % here the threshold is very important % originally set as 8
                if fg.list(thisi).rs > 51
                    fg.list(thisi).label = 1;
                else
                    fg.list(thisi).label = 0;
                end
                if isempty(find([fg.list.label].' == 1)) || length(find([fg.list.label].' == 1)) == 1
                    if fg.list(n).rs > 24
                        fg.list(n).label = 1;
                    else
                        fg.list(n).label = 0;
                    end
                    
                end
            end
            
        end
    
    
    end
    
end