classdef autoseg
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dirpath
        names
        config
    end
    
    methods
        
        function a = autoseg(dirpath)
            
            a.dirpath = dirpath;
            
            % plot the mean and sd
            dbstop if error
            
            names = Extract.filename(dirpath,'*.wav')
            a.names = names;
            
        end
        
        function Deprecated_autoSegmenter(a)
            tic;
            % plot the mean and sd
            dbstop if error
            
            dirpath = a.dirpath;
            names = a.names;
            
            outdirname = 'SegDataIndi';
            outdir = sprintf('%s\\%s',dirpath,outdirname);
            mkdir(outdir);
            
            segdata = struct;
            for idx = 1: length(names)
                
                if ~isempty(regexp(names{idx},'WNS', 'once'))  % if the name contain the 'WNS'
                    continue
                end
                
                wavpath = names{idx};
                [y,fs] = audioread(wavpath);
                [~,birdid,~] = fileparts(wavpath);
                [rawy,~,I,syledge,eleedge] = autoseg.main(y,fs,birdid);
                segdata.rawy = rawy;  % this might be useful
                segdata.I = I;
                segdata.syledge = syledge;
                segdata.eleedge = eleedge;
                segdata.birdid = birdid;
                
                save(sprintf('%s\\%s-segdata.mat',outdir,segdata.birdid),'segdata');
                
            end
            
            toc;
            
        end  %Deprecated
        
        function standard(a,num_of_bouts ) % 此乃标准算法
            % if no num_of_bouts,then just no such restrictions
            dbstop if error    
            input_folder = a.dirpath;
            names = a.names;
            files = Extract.filename(input_folder,'*.wav');
            outdirname = 'SegData';% to create outdir
            outdir = sprintf('%s\\%s',input_folder,outdirname);
            mkdir(outdir);
            
            if exist('num_of_bouts','var') % Case-1 只有需要限制bout数量的时候才需要把所有的wav拼接在一起
                
                y = {};
                for h = 1: length(files)
                    [y{h},fs] = audioread(files{h}); % concatenate all y
                end
                sumy = vertcat(y{:});
                
                %通过低精度的分割找到大概有多少个bout，从而 restrict bout
                fiy = highpass(sumy, 500, 32000);
                img = Cal.spec(fiy,32000);
               
                %img = flip(log(1+abs(S)));
            
                proI = smooth(mean(img),960); % preocessed I
                used_for_seg = proI;
                seg_threshold = 0.06;
                used_for_seg(used_for_seg > seg_threshold) = 1;
                used_for_seg(used_for_seg <= seg_threshold) = 0;
                curve_diff = diff(used_for_seg);
                test_initials = find(curve_diff == 1) + 1;
                test_terminals = find(curve_diff == -1);
                
                temp = split(input_folder,'\');
                birdid = sprintf('OTE-%s',temp{length(temp)}); % OTE means songs never played as CONs
                
                if num_of_bouts > length(test_terminals)
                    end_terminal = test_terminals(end)
                else
                    end_terminal = test_terminals(num_of_bouts);  % end_terminal is the terminal of the 10th bouts
                end
                
                y_to_analyze = (end_terminal/1000 + 0.1)*fs; % 0.1 is the compensation
                
                a.standard_config('remove short elements'); % set configurtation
                [rawy,~,I,syledge,eleedge] = autoseg.main(sumy(1:y_to_analyze),fs,birdid,a.config); %%% for default,a.config can be removed
                
                % 组装计算的量到最终的数据结构中
                segdata.rawy = rawy;  % this might be useful
                segdata.fs = fs;
                segdata.I = I;
                segdata.syledge = syledge;
                segdata.eleedge = eleedge;
                segdata.birdid = birdid;
                %segdata.sumy = y_to_analyze;
                save(sprintf('%s\\%s-segdata.mat',outdir,segdata.birdid),'segdata');
            
                
                
            else  %如果不需要restrict bout的数量
                
                for idx = 1: length(names)
                    
                    if ~isempty(regexp(names{idx},'WNS', 'once'))  % if the name contain the 'WNS'
                        wavpath = names{idx};
                        [y,fs] = audioread(wavpath);
                        [~,birdid,~] = fileparts(wavpath);
                        [~,~,I,syledge,eleedge] = autoseg.main(y,fs,birdid);
                        segdata.rawy = y;  % this might be useful
                        segdata.fs = fs;
                        segdata.I = I;
                        segdata.syledge = syledge;
                        segdata.eleedge = eleedge;
                        segdata.birdid = birdid;
                        save(sprintf('%s\\%s-segdata.mat',outdir,segdata.birdid),'segdata');
                        continue
                    end
                    
                    wavpath = names{idx};
                    [y,fs] = audioread(wavpath);
                    [~,birdid,~] = fileparts(wavpath);
                    [rawy,~,I,syledge,eleedge] = autoseg.main(y,fs,birdid);
                    segdata.rawy = rawy;  % this might be useful
                    segdata.fs = fs;
                    segdata.I = I;
                    segdata.syledge = syledge;
                    segdata.eleedge = eleedge;
                    segdata.birdid = birdid;
                    save(sprintf('%s\\%s-segdata.mat',outdir,segdata.birdid),'segdata');
                end
               
            end
            
       
             
        end
        
        function eleinf = Deprecated_singleSegmenter(a)
            
            dirpath = a.dirpath;
            names = a.names;
            
            tic;
            % plot the mean and sd
            dbstop if error
            
            names = Extract.filename(dirpath,'*.wav')
            outdirname = 'SegData';
            outdir = sprintf('%s\\%s',dirpath,outdirname);
            mkdir(outdir);
            
            segdata = struct;
            for idx = 1: length(names)
                
                if ~isempty(regexp(names{idx},'WNS', 'once'))  % if the name contain the 'WNS'
                    continue
                end
                
                wavpath = names{idx};
                [y,fs] = audioread(wavpath);
                [~,birdid,~] = fileparts(wavpath);
                [rawy,~,I,syledge,eleedge] = autoseg.main(y,fs,birdid);
                segdata.rawy = rawy;  % this might be useful
                segdata.I = I;
                segdata.syledge = syledge;
                segdata.eleedge = eleedge;
                segdata.birdid = birdid;
                
                ele{idx}  = seg2ele(segdata);
                
            end
            
            eleinf = horzcat(ele{:});
            
            
            toc;
            
        end  %Deprecated
        
        function eleinf = Deprecated_singleSegmenterSyl(a)
            
            dirpath = a.dirpath;
            names = a.names;
            
            tic;
            % plot the mean and sd
            dbstop if error
            
            names = Extract.filename(dirpath,'*.wav')
            outdirname = 'SegData';
            outdir = sprintf('%s\\%s',dirpath,outdirname);
            mkdir(outdir);
            
            segdata = struct;
            for idx = 1: length(names)
                
                if ~isempty(regexp(names{idx},'WNS', 'once'))  % if the name contain the 'WNS'
                    continue
                end
                
                wavpath = names{idx};
                [y,fs] = audioread(wavpath);
                [~,birdid,~] = fileparts(wavpath);
                [rawy,~,I,syledge,eleedge] = autoseg.main(y,fs,birdid);
                segdata.rawy = rawy;  % this might be useful
                segdata.I = I;
                segdata.syledge = syledge;
                segdata.eleedge = eleedge;
                segdata.birdid = birdid;
                
                ele{idx}  = autoseg.seg2syl(segdata); % %%% Here is the difference between singleSegmenter and singleSegmenterSyl
                
            end
            
            eleinf = horzcat(ele{:});
            
            
            toc;
            
        end  %Deprecated
        
        function Deprecated_intoElements(a,num_of_bouts)
            
            
            dbstop if error
            input_folder = a.dirpath;
            names = a.names;
            files = Extract.filename(input_folder,'*.wav');
            % to create outdir
            outdirname = 'SegData';
            outdir = sprintf('%s\\%s',input_folder,outdirname);
            mkdir(outdir);
            
            if exist('num_of_bouts','var') % Case-1 只有需要限制bout数量的时候才需要把所有的wav拼接在一起
                
                y = {};
                for h = 1: length(files)
                    [y{h},fs] = audioread(files{h}); % concatenate all y
                end
                sumy = vertcat(y{:});
                
                % process sumy to get feature for segmentation
                fiy = highpass(sumy, 500, 32000);
                img = Cal.spec(fiy,32000);
                proI = smooth(mean(img),960); % preocessed I
                % figure; plot(proI);

                % Unprecise segmentation to restrict bout
                used_for_seg = proI;
                seg_threshold = 0.06;
                used_for_seg(used_for_seg > seg_threshold) = 1;
                used_for_seg(used_for_seg <= seg_threshold) = 0;
                curve_diff = diff(used_for_seg);
                test_initials = find(curve_diff == 1) + 1;
                test_terminals = find(curve_diff == -1);
                
                if exist( 'num_of_bouts','var')
                    if num_of_bouts > length(test_terminals)
                        end_terminal = test_terminals(end);
                    else
                        end_terminal = test_terminals(num_of_bouts);  % end_terminal is the terminal of the 10th bouts
                    end
                else
                    end_terminal = test_terminals(end);
                end
                
                y_to_analyze = (end_terminal/1000 + 0.1)*fs; % 0.1 is the compensation
                
                temp = split(input_folder,'\');
                birdid = sprintf('OTE-%s',temp{length(temp)}); % OTE means songs never played as CONs
                [rawy,~,I,syledge,eleedge] = autoseg.main(sumy(1:y_to_analyze),fs,birdid);
                
            else  %如果不需要restrict bout的数量
                
                
                for idx = 1: length(names)
                    
                    if ~isempty(regexp(names{idx},'WNS', 'once'))  % if the name contain the 'WNS'
                        continue
                    end
                    
                    wavpath = names{idx};
                    [y,fs] = audioread(wavpath);
                    [~,birdid,~] = fileparts(wavpath);
                    [rawy,~,I,syledge,eleedge] = autoseg.main(y,fs,birdid);
                    segdata.rawy = rawy;  % this might be useful
                    segdata.I = I;
                    segdata.syledge = syledge;
                    segdata.eleedge = eleedge;
                    segdata.birdid = birdid;
                    
                    save(sprintf('%s\\%s-segdata.mat',outdir,segdata.birdid),'segdata');
                    
                end
            end
                
                
            segdata.rawy = rawy;  % this might be useful
            segdata.fs = fs;
            segdata.I = I;
            segdata.syledge = syledge;
            segdata.eleedge = eleedge;
            segdata.birdid = birdid;
            %segdata.sumy = y_to_analyze;
            
            save(sprintf('%s\\%s-segdata.mat',outdir,segdata.birdid),'segdata');
            
            
            
        end
        
        function standard_config(a,name) % no output, this function is to update autoseg.config
            % input can be a string , that string represent a config group
            % finally send config into the main function
            switch name
                
                case 'normal' % normnal is the threshold I used previously to balance segmentation of syllables and elements
                    a.config.filterlow = 900;
                    a.config.filterhigh = 6000;
                    a.config.envy_smoothnum = 150;
                    a.config.bIthres = 0.25;
                    a.config.strict_bIthres = 0.45;
                    a.config.smooth_thres = 12;
                    a.config.syllable_thres = 0.10;
                    a.config.large_locs_MinPeakHeight = 0.25;
                    a.config.too_short_seglen = 20;
                    a.config.gap_too_short_thres = 5;
                    a.config.min_sim = 0.8;
                    a.config.e_max_MinPeakProminence = 0.08;
                    a.config.minus_e_max_MinPeakProminence = 0.05;
                    a.config.minus_e_max_MinPeakHeight = -0.7;
                    a.config.flip_emax_MinPeakProminence = 0.08;
                    a.config.minus_flip_emax_MinPeakProminence = 0.05;
                    a.config.minus_flip_emax_MinPeakHeight = -0.7;
                    a.config.diffy_threshold = 0.32;
                    a.config.tpoint_x_thres = 30;
                    a.config.max_pre_pks_end = 0.28;
                    a.config.tpoint_x_close_to_initial = 18;
                    a.config.flip_diffy_threshold = 0.33;
                    a.config.max_flip_pre_pks = 0.38;
                    a.config.too_close_to_end_flip_tpoint_x = 18;
                    a.config.ele_time_judge = 28;
                    
                case 'remove short elements'
%                     a.config.filterlow = 900;
%                     a.config.filterhigh = 6000;
%                     a.config.envy_smoothnum = 150;
%                     a.config.bIthres = 0.25;
%                     a.config.strict_bIthres = 0.45;
%                     a.config.smooth_thres = 12;
%                     a.config.syllable_thres = 0.10;
%                     a.config.large_locs_MinPeakHeight = 0.25;
%                     a.config.too_short_seglen = 20;
%                     a.config.gap_too_short_thres = 5;
%                     a.config.min_sim = 0.8;
%                     a.config.e_max_MinPeakProminence = 0.08;
%                     a.config.minus_e_max_MinPeakProminence = 0.05;
%                     a.config.minus_e_max_MinPeakHeight = -0.7;
%                     a.config.flip_emax_MinPeakProminence = 0.08;
%                     a.config.minus_flip_emax_MinPeakProminence = 0.05;
%                     a.config.minus_flip_emax_MinPeakHeight = -0.7;
%                     a.config.diffy_threshold = 0.32;
%                     a.config.tpoint_x_thres = 30;
%                     a.config.max_pre_pks_end = 0.28;
%                     a.config.tpoint_x_close_to_initial = 18;
%                     a.config.flip_diffy_threshold = 0.33;
%                     a.config.max_flip_pre_pks = 0.38;
%                     a.config.too_close_to_end_flip_tpoint_x = 18;
                    a.config.ele_time_judge = 35;
                    
                
            end
        end
        
        
    end
    
    methods(Static)
        
        % 想办法写一个可以更新threshold（configuration)的main function！！！！！！！！！！！！
        
        function [rawy,fiy,I,syledge,eleedge] = main(inputy,fs,birdid,CONFIG) % this the main function of autoseg
            
            rawy = inputy;
            fiy = bandpass(inputy,[900 6000],fs); % preiously 5000
            
            envy = rescale(smooth(abs(fiy),150)); % amplitude envelope of y
            %powery = downsample(fiy.^2/length(fiy),fs/1000);
            %downy = downsample(abs(fiy),fs/1000);
            I = Cal.spec(fiy,fs); % I is the image of the whole song   
            [S,F,T] = spectrogram(2*fiy,hamming(512),512-round(fs/1e3),512,fs);
            % if the input stimuli is white noise
            if regexp(birdid,'WNS|wns|Wns')
                eleedge = [];
                syledge = [0,size(I,2)];
                return
            end
            
            % @#%$%*() This is very dangerous step to rescale the I
            
            bI = imbinarize(I,0.25); % biI is the binarized I
            strict_bI = imbinarize(I,0.45);
            smooth_num = 12; % number of points for smoothing
            
            meanbI = rescale(smooth(mean(bI,1),smooth_num));
            stdbI =   rescale( smooth(std(bI,1),smooth_num ));
            sumbI =  rescale(smooth(sum(bI,1),smooth_num));
            
            strict_meanbI = rescale(smooth(mean(bI,1),smooth_num));
            strict_stdbI =   rescale( smooth(std(bI,1),smooth_num ));
            strict_sumbI =  rescale(smooth(sum(bI,1),smooth_num));
            
            
            Imean = rescale(mean(I,1));
            Istd =   rescale(std(I,1));
            Isum =  rescale(sum(I,1));
            Imax = rescale(max(I)); % maximum value for each column of I
            [~,maxloc] = max(I);
            maxloc = rescale(-maxloc); % frequency-domain location of the frequency with maximum intensity
            
            r1 = 1;
            r2 = 1;
            r3 = 1;
            Imixed =  rescale(Imean*r1 + Istd*r2 + Isum*r3);
            mixedbI = rescale(meanbI*r1 + stdbI*r2 + sumbI*r3);
            
            
            % segment song into syllables
            syllable_threshold = 0.10;
            
            autoseg.drawsylfig(I,meanbI,stdbI,sumbI,Imax,maxloc,syllable_threshold); % function
            close(gcf)
            % FIND INTERSECTIONS    <@> One way
            used_for_seg = Imax;
            
            % <@> Anothere way
            
            used_for_seg(used_for_seg > syllable_threshold) = 1;
            used_for_seg(used_for_seg <= syllable_threshold) = 0;
            
            curve_diff = diff(used_for_seg);
            
            test_initials = find(curve_diff == 1) + 1;
            test_terminals = find(curve_diff == -1);
            
            edges = sort([test_initials,test_terminals]);
            
            
            
            if used_for_seg(1) > syllable_threshold %if the terminal has  large intensity
                edges = [1;edges(:)];
            end
            
            if used_for_seg(length(used_for_seg)) > syllable_threshold %if the terminal has  large intensity
                edges = [edges(:);length(used_for_seg)];
            end
            
            
            [~,large_locs] = findpeaks(used_for_seg,'MinPeakHeight',0.25); % Within the segmentted duration, maxI must override 0.25
            
            for k = 1: length(edges)/2
                if isempty(find(( large_locs >= edges(2*k-1)  & large_locs <= edges(2*k)  ), 1))
                    edges(2*k-1) = nan;
                    edges(2*k) = nan;
                end
                
            end
            
            edges = rmmissing(edges);
            
            % remove segments which are too short
            seglen = edges(2:2:end)-edges(1:2:end); % length of segments
            too_short = find(seglen < 20);
            edges(2*too_short - 1) = nan;
            edges(2*too_short) = nan;
            edges = rmmissing(edges);
            
            
            % remove the intervals which are too short
            check_short_initial = edges(1:2:end);
            check_short_initial = check_short_initial(2:end);
            
            check_short_terminal= edges(2:2:end);
            check_short_terminal = check_short_terminal(1:end-1);
            
            interval_bet_seg = check_short_initial  - check_short_terminal;
            
            gap_too_short = find(interval_bet_seg < 5);  % Here 5 is the threshold for shortest element
            
            
            edges(2*gap_too_short + 1) = nan;
            edges(2*gap_too_short ) = nan;
            edges = rmmissing(edges);
            
            
            autoseg.draw_seg_result(edges,I,Imax,syllable_threshold) % % % % draw the result of primary segmentation
            close gcf
            
            initials = edges(1:2:end);
            ends = edges(2:2:end);
            
            % @-------------------@ This section is for finding the similar fragments
            % and only segment the first one of each cluster
            
            for u = 1: length(initials)
                frag_img{u} = imresize( I(:,initials(u):ends(u)),[257,50]);
            end
            
            img_sim = []; % image similarityu matrix
            for z = 1: length(frag_img)
                for s = 1: length(frag_img)
                    img_sim(z,s) = ssim(frag_img{z},frag_img{s});
                end
            end
            
            img_sim = (img_sim + img_sim.')/2; % make img_sim symmetric
            min_sim = 0.8;% minimum similarity to sort fragments together
            % figure; imagesc(img_sim); % plotting
            img_sim_tri = triu(img_sim,1); % triangular
            for c = 1:length(img_sim)
                group{c} = find(img_sim(c,:)>min_sim)
            end
            
            %<solution-2>
            num_of_catego = 0;
            not_accessed = linspace(1,length(group),length(group));
            categos = {};
            
            for k = 1 :length(group)
                if ismember(k,not_accessed)
                    
                    num_of_catego = num_of_catego + 1;
                    
                    categos{num_of_catego} =  group{k};
                    while sum(cellfun(@length,cellfun(@(x) intersect(categos{num_of_catego},x),group(not_accessed),'UniformOutput',0)))~= 0
                        ids= find(cellfun(@length,cellfun(@(x) intersect(categos{num_of_catego},x),group(not_accessed),'UniformOutput',0)));
                        categos{num_of_catego} = unique(horzcat(categos{num_of_catego}, horzcat( group{not_accessed(ids) })));
                        not_accessed(ids) = []; % delete accessed
                        disp('Loop!!!!!!!')
                    end
                    
                else
                    continue
                end
            end
            
            
            
            % segment song into element
            element_tp = {};
            
            container(1:length(categos)) = struct('data',[],'tag',0);  % contain segmentation information, "tag" is to check whether the catego has been accessed or not
            
            for o = 1: length(initials)
                
                
                eimg = I(:,initials(o):ends(o));
                eimg_std = Istd (initials(o):ends(o));
                eimg_mean = Imean (initials(o):ends(o));
                eimg_ysum = Isum (initials(o):ends(o));
                eimg_mixed = Imixed (initials(o):ends(o));
                ebi_ysum =  sumbI(initials(o):ends(o));
                e_max = Imax(initials(o):ends(o));
                e_maxloc = maxloc(initials(o):ends(o));
                
                
                ey = envy(initials(o)/1000*fs:ends(o)/1000*fs); % pre processed
                later_processed_y = rescale (abs(fiy(initials(o)/1000*fs:ends(o)/1000*fs)));
                
                
                [highy,highx] = findpeaks(e_max,'MinPeakProminence',0.08);
                [temp,lowx] = findpeaks(-e_max,'MinPeakProminence',0.05,'MinPeakHeight',-0.7);
                lowy = -temp;
                % % % % % % % % % %
                % % % % % % % % % %     %after flipping
                flip_emax = flip(e_max);
                [flip_highy,flip_highx] = findpeaks(flip_emax ,'MinPeakProminence',0.08);
                [flip_temp,flip_lowx] = findpeaks(-flip_emax ,'MinPeakProminence',0.05,'MinPeakHeight',-0.7);
                flip_lowy = -flip_temp;
                
                tsm_highx = length(e_max)+ 1 - flip_highx; % tsm means transformed: transform flipped coordinate to not flipped coordinate
                tsm_lowx = length(e_max)+ 1 - flip_lowx; % transformation
                
                
                
                
                
                
                %(***) judge whether the element from a new catego or not
                this_catego = find( cellfun(@sum, cellfun(@(x)ismember(o,x),categos,'uni',false)));
                
                if container(this_catego).tag == 1 % if  the same catego has been accessed
                    tpoint_x = container(this_catego).data;
                    element_tp{o} = tpoint_x + initials(o);
                    close(gcf)
                    continue
                end
                
                
                % if there is only one high peak, no segmentation is necessary
                
                if length(highx)== 1
                    close(gcf)
                    container(this_catego).tag = 1;
                    continue
                else
                    disp('Need segmentation!')
                    
                    % For normal peak/valley
                    xcollect = [1,highx,lowx,length(e_max)];
                    ycollect = [e_max(1),highy,lowy,e_max(end)];
                    [sortedx, sortIndex] = sort(xcollect);
                    sortedy = ycollect(sortIndex);
                    diffy = diff(sortedy);
                    diffy_threshold = 0.32;
                    transfer_point_idx = find(diffy >diffy_threshold);
                    tpoint_x = sortedx(transfer_point_idx);
                    tpoint_x(tpoint_x == 1) = [];   % remove the first idx: which indicating the first element
                    tpoint_y = sortedy(transfer_point_idx);
                    tpoint_y(tpoint_y == 1) = [];
                    % Find  the previous peak which is  cloest to the tpoint_x, if the
                    % differenc is too small, then remove this tpoint_x
                    for  s = 1:length(tpoint_x)
                        if tpoint_x < 30 % withe the case that the tpoint is too short
                            previous_max = e_max(1: tpoint_x(s));
                            [pre_pks,pre_locs] = findpeaks(previous_max);
                            if ~isempty(pre_locs)
                                if max(pre_pks(end))< 0.28 % if the difference between the previous peak and this valley is too small
                                    tpoint_x(s) = nan;
                                    tpoint_y(s) = nan;
                                end
                            end
                        end
                    end
                    
                    tpoint_x(tpoint_x < 18) = nan; % remove those too short???Maybe what I mean is too close to the initial
                    tpoint_x = rmmissing(tpoint_x);
                    
                    
                    
                    % For fipped peak/valley
                    flip_xcollect = [1,flip_highx,flip_lowx,length(flip_emax)];
                    flip_ycollect = [flip_emax(1),flip_highy,flip_lowy,flip_emax(end)];
                    [flip_sortedx, flip_sortIndex] = sort(flip_xcollect);
                    flip_sortedy = flip_ycollect(flip_sortIndex);
                    flip_diffy = diff(flip_sortedy);
                    flip_diffy_threshold = 0.33; % DIff in Y
                    flip_transfer_point_idx = find(flip_diffy >flip_diffy_threshold);
                    flip_tpoint_x = flip_sortedx(flip_transfer_point_idx);
                    flip_tpoint_x(flip_tpoint_x == 1) = [];   % remove the first idx: which indicating the first element
                    flip_tpoint_y = flip_sortedy(flip_transfer_point_idx);
                    flip_tpoint_y(flip_tpoint_y == 1) = [];
                    
                    
                    for  s = 1:length(flip_tpoint_x)
                        flip_previous_max = flip_emax(1: flip_tpoint_x(s));
                        [flip_pre_pks,flip_pre_locs] = findpeaks(flip_previous_max);
                        if ~isempty(flip_pre_locs)
                            if max(flip_pre_pks) < 0.38 % Which means the rest of the signal is too weak
                                %if flip_pre_pks(end)- flip_tpoint_y(s) < 0.2 % if the difference between the previous peak and this valley is too small
                                flip_tpoint_x(s) = nan;
                                flip_tpoint_y(s) = nan;
                            end
                        end
                    end
                    
                    
                    flip_tpoint_x( flip_tpoint_x < 18) =   nan; % remove those which are too close to the end
                    flip_tpoint_x = rmmissing(flip_tpoint_x);
                    
                    trans_flip_tpoint_x = length(flip_emax) + 1 - flip_tpoint_x;
                    
                    
                    %%% draw segmentation of element
                    autoseg.drawelefig(eimg,eimg_mean,eimg_std,eimg_ysum,eimg_mixed,e_max,e_maxloc,highx,highy...
                        ,lowx,lowy,tsm_highx,flip_highy,tsm_lowx,flip_lowy,tpoint_x,trans_flip_tpoint_x )
                    
                    
                    all_X = union(tpoint_x,trans_flip_tpoint_x);
                    
                    % remove some of the all_X if the distance between them are too short
                    ele_time_judge = [1,all_X(:).',length(e_max)];
                    
                    if exist( 'CONFIG','var')
                        short_ele = find(diff(ele_time_judge) < CONFIG.ele_time_judge);
                    else
                        short_ele = find(diff(ele_time_judge) < 28); % originallly set to 15
                    end
                    
                    
                    for yy = 1: length(short_ele)
                        if short_ele(yy) == 1
                            ele_time_judge(1) = nan; % remove the second for the first block
                        elseif short_ele(yy) == length(ele_time_judge) - 1
                            ele_time_judge(length(ele_time_judge) - 1) = nan; % remove the len-1 for the first block
                        else
                            
                            % judge which value is lower
                            left_point =   e_max(short_ele(yy));
                            right_point =   e_max(short_ele(yy)+ 1);
                            
                            if left_point < right_point % keep the lower one
                                ele_time_judge(short_ele(yy)+ 1) = nan;
                            else
                                ele_time_judge(short_ele(yy)) = nan;
                            end
                        end
                        
                    end
                    
                    ele_time_judge = rmmissing(ele_time_judge);
                    
                    all_X = ele_time_judge(2:end-1);
                    
                    
                    element_tp{o} = all_X + initials(o);
                    
                    %(***) update the state of catego
                    
                    if isempty(container(this_catego).data) % if sgementation from the same catego has been done
                        container(this_catego).data = all_X;
                        container(this_catego).tag = 1;
                    end
                    
                end
                
            end
            
            for tt = 1: length(element_tp)
                element_tp{tt} = element_tp{tt}(:); % force each element in element_tp to be a column vector
            end
            
            tp_sum = vertcat(element_tp{:}); % sum tps from all elements together
            
            % @_@_@_@_@_@_@
            if exist('birdid','var')
                autoseg.drawfinalfig(I,edges,tp_sum,birdid); % plot the final result
            end
            
            
            % assign syllable edges and also element edges
            syledge = T(edges);
            eleedge = T(tp_sum);
            
            
            disp('长风几万里，吹度玉门关')
            
        end
        
        function drawsylfig(I,meanbI,stdbI,sumbI,Imax,maxloc,syllable_threshold)
            
            sylfig = figure
            frame_h = get(handle(sylfig),'JavaFrame');
            set(frame_h,'Maximized',1);
            
            t = tiledlayout(sylfig,1,1);
            
            ax1 = axes(t);
            imagesc(ax1,I)
            colormap(gray)
            ax1.Color = 'none';ax1.Box = 'off';
            set(ax1, 'YTickLabel', [],'YTick',[]);
            
            % ax2 = axes(t);
            % plot(ax2,meanbI,'Color','g')
            % xlim([0 length(I)]);
            % ax2.Color = 'none';ax2.Box = 'off';
            %
            % ax3 = axes(t);
            % plot(ax3,stdbI,'Color','r')
            % xlim([0 length(I)]);
            % ax3.Color = 'none';ax3.Box = 'off';
            %
            % ax4 = axes(t);
            % plot(ax4,sumbI,'Color','y')
            % xlim([0 length(I)]);
            % ax4.Color = 'none';ax4.Box = 'off';
            
            ax5 = axes(t);
            plot(ax5,Imax,'Color','Cyan')
            xlim([0 length(Imax)]);
            ax5.Color = 'none';ax5.Box = 'off';
            
            % ax6 = axes(t);
            % plot(ax6,maxloc,'Color',[.6 .1 .9])
            % xlim([0 length(maxloc)]);
            % ax6.Color = 'none';ax6.Box = 'off';
            
            
            % % &&& To plot it the intersections
            ax5;
            %plot(linspace(1,length(mixedbI),length(mixedbI)),mixedbI);
            hold on
            plot(linspace(...
                1,length(Imax),length(Imax)),syllable_threshold*ones(1,length(Imax)));
            
            % delete some peaks
            
        end
        
        function drawelefig(eimg,eimg_mean,eimg_std,eimg_ysum,eimg_mixed,e_max,e_maxloc,highx,highy...
                ,lowx,lowy,tsm_highx,flip_highy,tsm_lowx,flip_lowy,tpoint_x,trans_flip_tpoint_x )
            
            figure('Position',[1308 123 598 934])
            %set(gcf, 'Position', get(0, 'Screensize'));
            t = tiledlayout(1,1);
            ax1 = axes(t);
            imagesc(ax1,eimg);
            colormap(gray)
            ax1.Color = 'none';
            ax1.Box = 'off';
            set(ax1, 'YTickLabel', [],'YTick',[]);
            
            ax2 = axes(t);
            plot(ax2,eimg_mean,'Color','g');
            xlim([0 size(eimg,2)]);
            ylim([0,1])
            ax2.Color = 'none';
            ax2.Box = 'off';
            
            ax3 = axes(t);
            plot(ax3,eimg_std,'Color','r');
            xlim([0 size(eimg,2)]);
            ylim([0,1])
            ax3.Color = 'none';
            ax3.Box = 'off';
            
            ax4 = axes(t);
            plot(ax4,eimg_ysum,'Color','y');
            xlim([0 size(eimg,2)]);
            ylim([0,1])
            ax4.Color = 'none';
            ax4.Box = 'off';
            
            ax5 = axes(t);
            plot(ax5,eimg_mixed,'Color','Cyan');
            xlim([0 size(eimg,2)]);
            ylim([0,1])
            ax5.Color = 'none';
            ax5.Box = 'off';
            
            ax6 = axes(t);
            plot(ax6,e_max,'Color','[.6 .6 .9]');
            xlim([0 length(e_max)]);
            ylim([0,1])
            ax6.Color = 'none';
            ax6.Box = 'off';
            
            ax7 = axes(t);
            plot(ax7,e_maxloc,'Color','[.9 .6 .1]');
            xlim([0 length(e_maxloc)]);
            ylim([0,1])
            ax7.Color = 'none';
            ax7.Box = 'off';
            
            
            axes(ax6)
            hold on
            % plot the nomral peak/valley
            scatter(ax6,highx,highy,60,'g','^');
            if ~isempty(lowx)
                hold on
                scatter(ax6,lowx,lowy,60,'g','o');
            end
            
            % plot the flipped peak/valley
            scatter(ax6,tsm_highx,flip_highy,60,'y','v');
            if ~isempty(tsm_lowx)
                hold on
                scatter(ax6,tsm_lowx,flip_lowy,60,'y','*');
            end
            
            % plot the segmentation line
            axes(ax6)
            for t = 1: length(tpoint_x)
                hold on
                line([tpoint_x(t),tpoint_x(t)],[0,1],'Color','w')
            end
            
            
            axes(ax6)
            for t = 1: length(trans_flip_tpoint_x)
                hold on
                line([trans_flip_tpoint_x(t),trans_flip_tpoint_x(t)],[0,1],'Color','Magenta','LineStyle',':')
            end
            hold off
            
            %close(gcf)
            
        end
        
        function drawfinalfig(I,edges,tp_sum,birdid)
            
            % plot the final result
            f2 = figure;
            frame_h = get(handle(f2),'JavaFrame');
            set(frame_h,'Maximized',1);
            yyaxis left
            imagesc(I)
            colormap(gray)
            
            for k = 1: length(edges)/2
                hold on
                q = area([edges(2*k-1),edges(2*k)],[257,257],'FaceColor','r','FaceAlpha',0.2);
            end
            
            for p = 1: length(tp_sum)
                yyaxis right
                hold on
                line([tp_sum(p),tp_sum(p)],[0,1],'Color','w');
            end
            
            title(sprintf('%s.png',birdid));
            
            figdir = 'Segmentation_of_CON_senatus';
            mkdir(figdir);
            
            if exist(sprintf('%s\\FinalResult-%s.png',figdir,birdid),'file')
                delete(sprintf('%s\\FinalResult-%s.png',figdir,birdid));
            end
            
            saveas(f2,sprintf('%s\\FinalResult-%s.png',figdir,birdid));
            
        end
        
        function draw_seg_result(edges,I,Imax,syllable_threshold)
            
            % draw the result of primary segmentation
            area_fig = figure
            frame_h = get(handle(area_fig),'JavaFrame');
            set(frame_h,'Maximized',1);
            
            t = tiledlayout(area_fig,1,1);
            ax1 = axes(t);
            
            imagesc(I)
            colormap(gray)
            
            for k = 1: length(edges)/2
                hold on
                q = area(ax1,[edges(2*k-1),edges(2*k)],[257,257],'FaceColor','r','FaceAlpha',0.2);
            end
            
            ax2 = axes(t);
            % plot(ax2,mixedbI,'Color','Cyan')
            % xlim([0 length(mixedbI)]);
            plot(ax2,Imax,'Color','Cyan')
            used_for_seg = Imax;
            xlim([0 length(used_for_seg)]);
            ax2.Color = 'none';ax2.Box = 'off';
            
            yline(syllable_threshold,'Color','Magenta')
            
            % % % % % % if exist(sprintf('SyllableSeg-%s.png',birdid),'file')
            % % % % % %     delete(sprintf('SyllableSeg-%s.png',birdid));
            % % % % % % end
            % % % % % % saveas(area_fig,sprintf('SyllableSeg-%s.png',birdid));
            
        end
        
        function eleinf = seg2ele(segdata)
            
            all_collect = {};
            song_collect = {};
            
            
            two_eleedge = repmat(segdata.eleedge,[2,1]);
            alledges = sort( vertcat(segdata.syledge(:),reshape(two_eleedge,[],1) ));
            
            % the following is a very bad temporily code
            if ~isfield(segdata,'fs') % if there is no such field named fs
                fs = 32000;
            else
                fs = segdata.fs;
            end
            
            initials = alledges(1:2:end);
            terminals = alledges(2:2:end);
            
            song_eleinf = struct;
            
            for w = 1: length(initials) % can add a par-for
                song_eleinf(w).initial = initials(w)*fs/1000;
                song_eleinf(w).terminal = terminals(w)*fs/1000;
                song_eleinf(w).songname = segdata.birdid;
                
                if isfield(segdata,'y')
                    song_eleinf(w).y = segdata.y(initials(w)*fs/1000:terminals(w)*fs/1000);  % very bad code, but meibanfa
                elseif isfield(segdata,'rawy')
                    song_eleinf(w).y = segdata.rawy(initials(w)*fs/1000:terminals(w)*fs/1000);
                end
                
                 fiy = highpass(song_eleinf(w).y, 500, 32000);
                 I = Cal.spec(fiy,fs); 
                 song_eleinf(w).fragI = imresize(I,[257,50]);
                
                song_eleinf(w).fs = fs;
            end
            
            parfor f = 1:length(song_eleinf)
                song_eleinf(f).fragid = f;
            end
            
            eleinf = song_eleinf;
            
        end
        
        function eleinf = seg2syl(segdata)
            
            all_collect = {};
            song_collect = {};
            
           
            alledges = sort( vertcat(segdata.syledge(:)));
            
            % the following is a very bad temporily code
            if ~isfield(segdata,'fs') % if there is no such field named fs
                fs = 32000;
            else
                fs = segdata.fs;
            end
            
            initials = alledges(1:2:end);
            terminals = alledges(2:2:end);
            
            song_eleinf = struct;
            
            for w = 1: length(initials) % can add a par-for
                song_eleinf(w).initial = initials(w)*fs/1000;
                song_eleinf(w).terminal = terminals(w)*fs/1000;
                song_eleinf(w).songname = segdata.birdid;
                
                if isfield(segdata,'y')
                    song_eleinf(w).y = segdata.y(initials(w)*fs/1000:terminals(w)*fs/1000);  % very bad code, but meibanfa
                elseif isfield(segdata,'rawy')
                    song_eleinf(w).y = segdata.rawy(initials(w)*fs/1000:terminals(w)*fs/1000);
                end
                
                 fiy = highpass(song_eleinf(w).y, 500, 32000);
                 I = Cal.spec(fiy,fs); 
                 song_eleinf(w).fragI = imresize(I,[257,50]);
                
                song_eleinf(w).fs = fs;
            end
            
            parfor f = 1:length(song_eleinf)
                song_eleinf(f).fragid = f;
            end
            
            eleinf = song_eleinf;
            
        end
        
    end
end

