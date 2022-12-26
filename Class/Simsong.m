classdef Simsong < handle % similar song
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        list
        simlist
        formated_name
    end

    methods
        function ss = Simsong(input_list)

            % insert names of similar (sibling/nonsibling) songs
            siblinggroups = {["R707","R703","R704","R705"],["Y515","G548"]};

            for k = 1:length(siblinggroups)

                jointed = strjoin( siblinggroups{k},'|');

                ss.simlist{k} = input_list(find(~cellfun(@isempty, regexp(cellstr({input_list.stimuliname}.'),jointed))));

            end


        end

        function outputArg = drawThree(ss)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            for k = 1:length(ss.simlist)
                locallist = ss.simlist{k}
                if length(locallist)~= 0
                    img = {};

                    for kk = 1:length(locallist)
                        figure;
                        Draw.twoForPoster(locallist(kk).y,locallist(kk).sptimes, locallist(kk).fs);

                        img{kk} = getframe(gcf).cdata;
                    end
                    sumimg = vertcat(img{:});
                    figure
                    imshow(sumimg)




                end

            end
           
        end
   
    
        function How_Do_experiments_Differentiate_Sibling_Songs(neu)
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
                ids_norm = find(... % 假设不会存在一首song的Degs在不同pl2文件里出现
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
                ids_norm = find(... % 假设不会存在一首song的Degs在不同pl2文件里出现
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

            imwrite(Iall,sprintf('Sibling_Songs_%s.png',neu.formated_name));
            toc

        end

    
    end
end