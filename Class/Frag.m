classdef Frag < handle

    properties
        allfraglist % concatenated fragment list
        allfraglist_nonorm % contain all the fragment stimuli but no norm song stimuli
        close_fraglist %那些离target比较近的syllables
        far_fraglist  % 那些离target比较远的syllables
        fraglist
        list
        formated_name
        simatrix
        weird_fraglist %某些fraglist包含了一些经过了特殊modification的syllables
    end
    methods
        function fg = Frag(list)
            fg.list = list;

            include_weird_fragids = find(~cellfun(@isempty, regexp(cellstr({fg.list.stimuliname}.'),'frag|Frag') ));
            pure_fagids = find(~cellfun(@isempty, regexp(cellstr({fg.list.stimuliname}.'),'^Frag') ));
            fg.weird_fraglist = fg.list(setdiff(include_weird_fragids,pure_fagids));
            list_contain_frags = fg.list(pure_fagids);
            frag_birdids = cellfun(@Convert.bid, cellstr({list_contain_frags.stimuliname}.'),'Uni',0);
            unique_frag_birdids = unique(frag_birdids);
            unpack = @(x) x{1};
            for k = 1:length(unique_frag_birdids)
                correspids = strcmp(unique_frag_birdids{k},frag_birdids); % 找到有这个birdid的所有frag的ids,这ids是关于list_contain_frags的
                fg.fraglist{k} = list_contain_frags(correspids);

                corresp_fid = unpack(unique({fg.fraglist{k}.Fid}.')); % 找到frag stimuli对应的fid

                corresp_norm = mintersect( find(~cellfun(@isempty, regexp(cellstr({fg.list.stimuliname}.'),unique_frag_birdids{k}))),...
                    find(~cellfun(@isempty, regexp(cellstr({fg.list.stimuliname}.'),'norm'))),...
                    find(~cellfun(@isempty, regexp(cellstr({fg.list.Fid}.'),corresp_fid))) );  %这里找的是和fragment同时播放的norm song

                fg.fraglist{k} = horzcat(fg.list(corresp_norm),fg.fraglist{k});

            end
            if ~isempty(fg.fraglist)
                fg.allfraglist = horzcat(fg.fraglist{:});
                fg.allfraglist_nonorm = fg.allfraglist(find(~cellfun(@isempty, regexp(cellstr({fg.allfraglist.stimuliname}.'),'Frag|frag'))));


                %构造Simatrix

                %             fraglist = neu.list( find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'), 'Frag|frag') )));
                fg.simatrix = Simatrix(fg.allfraglist_nonorm);
                fg.simatrix.formated_name = fg.formated_name;
            end
            

        end


        function  knowCloseAndFarFrags(fg,targetname) % know which syllables is a par one , and which is a close one


            % 计算使神经元反应的syllables占所有播放的syllables的比例
            % the ratio of responded syllables compares with all tested
            % syllables
            if isempty(targetname)
                disp('Sadly that Targetname is empty')
                return
            end
            targetname = targetname{1}; % 这个targetname其实可以有很多，所以用{1}不太对啊
            
            local_fraglist = fg.allfraglist_nonorm;
            targetid = find(~cellfun(@isempty, regexp(cellstr({local_fraglist.stimuliname}.'), targetname)));


            key = 0;
            if key == 1 % 暂且不使用
                for k = 1:length(local_fraglist)
                    local_fraglist(k).sim2target = fg.simatrix.matrix.global(targetid,k);
                end

                threshold = 0.85; % similarity to judge very similar syllables and not very similar syllables
                % the list which contains the very similar syllables
                fg.close_fraglist = local_fraglist([local_fraglist.sim2target].'>threshold);
                % the list which contains the not very similar syllables
                fg.far_fraglist = local_fraglist([local_fraglist.sim2target].'<=threshold);
            end
            for k = 1:length(local_fraglist)
                try
                    local_fraglist(k).dis2target = norm(local_fraglist(targetid).coor-local_fraglist(k).coor);
                catch
                    local_fraglist(k).dis2target = nan;
                end
            end

            threshold = 0.85; % similarity to judge very similar syllables and not very similar syllables
            temp_tobe_deleted = sort([local_fraglist.dis2target].');
            threshold = max(temp_tobe_deleted(1:31)); % 权宜之计，这个后面必须修改
            % the list which contains the very similar syllables
            fg.close_fraglist = local_fraglist([local_fraglist.dis2target].'>threshold);
            % the list which contains the not very similar syllables
            fg.far_fraglist = local_fraglist([local_fraglist.dis2target].'<=threshold);

            


        end

        function percentage = percentageResponse(fg) % percentage of responsive syllables
            if isempty(fg.fraglist)
                percentage = [];
                return
            end

            percentage.numerator_all = length(find([fg.allfraglist_nonorm.label].' == 1));
            percentage.denominator_all = length(fg.allfraglist);
            percentage.ratio_all = percentage.numerator_all/percentage.denominator_all ;



            % 与target接近的syllable的反应/测试比例
            if isempty(fg.close_fraglist)
                percentage.numerator_close = 0;
            else
            percentage.numerator_close = length(find([fg.close_fraglist.label].' == 1));
            end
            percentage.denominator_close = length(fg.close_fraglist);
            percentage.ratio_close =  percentage.numerator_close /percentage.denominator_close ;

            % 与target不接近的syllable的反应/测试比例
            percentage.numerator_far = length(find([fg.allfraglist_nonorm.label].' == 1));
            percentage.denominator_far = length(fg.allfraglist);
            percentage.ratio_far = percentage.numerator_all/length(fg.allfraglist);




        end


        function fg = jufgeFragResp(fg) % 似乎暂时并未使用
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


        function saveDrawSSIMSimlarityMatrix(neu)
            % This function firstly order the frags by response strength
            % then measure the pairwise similarity between elements
            fraglist = neu.judgeFragResponse;
            sorted_fraglist = table2struct(sortrows( struct2table(fraglist) ,'halfsum','descend'));  %此处用 maxvalue或许不太对

            for k = 1: length(sorted_fraglist)
                fiy = bandpass(sorted_fraglist(k).y,[900 6000],sorted_fraglist(k).fs); % preiously 5000

                envy = rescale(smooth(abs(fiy),150)); % amplitude envelope of y
                %powery = downsample(fiy.^2/length(fiy),fs/1000);
                %downy = downsample(abs(fiy),fs/1000);
                I = Cal.spec(fiy,sorted_fraglist(k).fs); % I is the image of the whole song
                sorted_fraglist(k).normalized_img = imresize(I,[257,50]);
            end

            for m = 1:length(sorted_fraglist)
                parfor p = 1:length(sorted_fraglist)
                    sim(m,p) = ssim(sorted_fraglist(m).normalized_img,sorted_fraglist(p).normalized_img);
                end
            end

            figure;
            imagesc(sim);

            saveas(gcf,sprintf('SSIMSimlarityMatrix-%s.png',neu.formated_name));
            close(gcf);
        end

        function saveDrawMeanFeatureVsResp(neu)
            % draw the distribution of mean features

            %             fragids = find(~cellfun(@isempty, regexp({neu.list.stimuliname}.','frag|Frag|syl|Syl')));
            %
            %             if isempty(fragids)
            %                 return
            %             end
            %             fraglist = neu.list(fragids);

            fraglist = neu.judgeFragResponse;
            if isempty(fraglist)
                return;
            end
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).rs; % Response Strength
            end
            fraglist = table2struct(sortrows(struct2table(fraglist), 'responseMeasure'));

            figure;
            subplot(5,2,1)
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).value = fraglist(u).meanfeatures.amplitude;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            toshow  = table2struct(sortrows(struct2table(toshow), 'value'));
            plot([toshow.value].',[toshow.resp].');
            xlabel('Amplitude')
            ylabel('Maximum FR');

            subplot(5,2,2)
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).value = fraglist(u).meanfeatures.pitch;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            toshow  = table2struct(sortrows(struct2table(toshow), 'value'));
            plot([toshow.value].',[toshow.resp].');
            xlabel('Pitch')
            ylabel('Maximum FR');

            subplot(5,2,3)
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).value = fraglist(u).meanfeatures.mean_frequency;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            toshow  = table2struct(sortrows(struct2table(toshow), 'value'));
            plot([toshow.value].',[toshow.resp].');
            xlabel('Mean Frequency')
            ylabel('Maximum FR');

            subplot(5,2,4)
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).value = fraglist(u).meanfeatures.FM;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            toshow  = table2struct(sortrows(struct2table(toshow), 'value'));
            plot([toshow.value].',[toshow.resp].');
            xlabel('FM')
            ylabel('Maximum FR');

            subplot(5,2,5)
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).value = fraglist(u).meanfeatures.AM;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            toshow  = table2struct(sortrows(struct2table(toshow), 'value'));
            plot([toshow.value].',[toshow.resp].');
            xlabel('AM')
            ylabel('Maximum FR');

            subplot(5,2,6)
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).value = fraglist(u).meanfeatures.goodness;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            toshow  = table2struct(sortrows(struct2table(toshow), 'value'));
            plot([toshow.value].',[toshow.resp].');
            xlabel('Goodness')
            ylabel('Maximum FR');

            subplot(5,2,7)
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).value = fraglist(u).meanfeatures.entropy;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            toshow  = table2struct(sortrows(struct2table(toshow), 'value'));
            plot([toshow.value].',[toshow.resp].');
            xlabel('Entropy')
            ylabel('Maximum FR');

            subplot(5,2,8)
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).value = fraglist(u).meanfeatures.peak_frequency;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            toshow  = table2struct(sortrows(struct2table(toshow), 'value'));
            plot([toshow.value].',[toshow.resp].');
            xlabel('Peak Frequency')
            ylabel('Maximum FR');

            subplot(5,2,9)
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).value = fraglist(u).meanfeatures.continuity_t;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            toshow  = table2struct(sortrows(struct2table(toshow), 'value'));
            plot([toshow.value].',[toshow.resp].');
            xlabel('Time continuity')
            ylabel('Maximum FR');

            subplot(5,2,10)
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).value = fraglist(u).meanfeatures.continuity_f;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            toshow  = table2struct(sortrows(struct2table(toshow), 'value'));
            plot([toshow.value].',[toshow.resp].');
            xlabel('frequency continuity')
            ylabel('Maximum FR');


            saveas(gcf, sprintf('SeparatedFeatureVsResponse_%s.png',neu.formated_name));
            close(gcf);


        end

        function saveDrawMeanFeaturesVsRespAsLineChart(neu)
            % draw the distribution of mean features
            neu.calHarmRatio;
            fraglist = neu.judgeFragResponse;
            if isempty(fraglist)
                return;
            end
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).rs;
                % fraglist(k).responseMeasure = fraglist(k).maxvalue;
            end

            toshow = struct; error = struct;
            for u = 1: length(fraglist)
                try
                    toshow(u).amplitude = fraglist(u).meanfeatures.amplitude;
                    toshow(u).pitch = fraglist(u).meanfeatures.pitch;
                    toshow(u).AM = fraglist(u).meanfeatures.AM;
                    toshow(u).mean_frequency = fraglist(u).meanfeatures.mean_frequency;
                    toshow(u).FM = fraglist(u).meanfeatures.FM;
                    toshow(u).peak_frequency = fraglist(u).meanfeatures.peak_frequency;
                    toshow(u).entropy = fraglist(u).meanfeatures.entropy;
                    toshow(u).harmratio = fraglist(u).meanfeatures.harmratio;
                    toshow(u).resp = fraglist(u).responseMeasure;
                catch ME
                    %toshow(u) = [];
                    error(u).ME = ME;
                end
            end

            toshow = table2struct(rmmissing(struct2table(toshow)));
            toshow = toshow(find(~cellfun(@isempty,{toshow.amplitude}.')));
            toshow  = table2struct(sortrows(struct2table(toshow), 'resp','ascend'));


            num_toshow = [toshow.amplitude; toshow.pitch; toshow.AM; toshow.mean_frequency; toshow.FM; ...
                toshow.peak_frequency; toshow.entropy; toshow.harmratio].';

            znum_toshow = zscore(num_toshow,0,1);
            %
            %            figure
            %            plot(znum_toshow.');

            figure('Position',[1997 233 1388 658],'Color','w');
            hold on
            c = 1- rescale([toshow.resp].',0.1,1);
            for r = 1: size(znum_toshow,1)
                plot(znum_toshow(r,:),'Color',repmat(c(r),3,1));
                %drawnow
                % pause(0.5)
            end
            colormap(flip(repmat(unique(c),1,3)))
            colorbar
            xlim([0,9]);
            xticks([0 1 2 3 4 5 6 7 8,9])
            xticklabels({'','amplitude','pitch','AM','mean_frequency','FM','peak-frequency','entropy','Harmonic ratio',''});
            title(sprintf('Totally %u song elements',length(fraglist)));
            ylabel('Zscored Feature(averaged)');

            saveas(gcf,sprintf('LineChartMeanFeaturesVsResp-%s.png',neu.formated_name));
            close(gcf);
        end


        function IMG = saveDrawlineChart_2D_Cumulative(neu) % 秦

            img{1} = neu.drawPitchHarmoLineChart;
            img{2} = neu.draw2DPitchVsHarmo;
            img{3} = neu.drawCumulativeFeatureDiff('pitch');
            img{4} = neu.drawCumulativeFeatureDiff('mean_frequency');
            img{5} = neu.drawCumulativeFeatureDiff('peak_frequency');
            img{6} = neu.drawCumulativeFeatureDiff('FM');
            img{7} = neu.drawCumulativeFeatureDiff('goodness');
            img{8} = neu.drawCumulativeFeatureDiff('entropy');
            img{9} = neu.drawCumulativeFeatureDiff('harmratio');
            img{10} = neu.drawCumulativeFeatureDiff('amplitude');
            img{11} = neu.drawCumulativeFeatureDiff('AM');

            % draw blank white
            hangshu = 3;
            lieshu = ceil(length(img)/hangshu);
            rest = lieshu*hangshu - length(img);
            white = uint8(255*ones(size(img{1})));

            if rest > 0
                for k = 1:rest;img = [img,white];end
            end
            reshapedI = reshape(img, lieshu,[])';
            clear img
            IMG = cell2mat(reshapedI);
            imwrite(IMG,sprintf('秦lineChart_2D_Cumulative_%s.png',neu.formated_name));
            %imwrite(img,sprintf('新lineChart_2D_Cumulative_%s.png',neu.formated_name));
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


        function BinaryThresholdDrawMeanFeaturesVsRespAsLineChart(neu) % draw the distribution of mean features

            fraglist = neu.judgeFragResponse;
            if isempty(fraglist)
                return;
            end
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).maxvalue;
            end

            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).amplitude = fraglist(u).meanfeatures.amplitude;
                toshow(u).pitch = fraglist(u).meanfeatures.pitch;
                toshow(u).AM = fraglist(u).meanfeatures.AM;
                toshow(u).mean_frequency = fraglist(u).meanfeatures.mean_frequency;
                toshow(u).FM = fraglist(u).meanfeatures.FM;
                toshow(u).peak_frequency = fraglist(u).meanfeatures.peak_frequency;
                toshow(u).entropy = fraglist(u).meanfeatures.entropy;
                toshow(u).resp = fraglist(u).responseMeasure;
            end

            toshow  = table2struct(sortrows(struct2table(toshow), 'resp','ascend'));


            num_toshow = [toshow.amplitude; toshow.pitch; toshow.AM; toshow.mean_frequency; toshow.FM; ...
                toshow.peak_frequency; toshow.entropy].';

            znum_toshow = zscore(num_toshow,0,1);
            %
            %            figure
            %            plot(znum_toshow.');

            figure('Position',[1997 233 1388 658],'Color','w');
            hold on
            % c = 1- rescale([toshow.resp].',0.1,1);


            % very artifical bad code
            switch neu.formated_name
                case 'O686_34'
                    thres = 2;
                case 'R677_55'
                    thres = 6;
                case 'Y661_5'
                    thres = 6;
                case 'Y661_8'
                    thres = 2;
                case 'Y675_22'
                    thres = 3;
                case 'R677_59'
                    thres = 3;
                otherwise
                    thres = 5;
            end

            kemal = 0;
            for r = 1: size(znum_toshow,1)
                if toshow(r).resp > thres
                    kemal = kemal + 1;
                    defaultLW = get(gca,'LineWidth');
                    plot(znum_toshow(r,:),'Color',[1 0 0],'LineWidth',defaultLW*1.8);
                else
                    plot(znum_toshow(r,:),'Color',[0.6 0.6 0.6]);
                end
            end

            xlim([0,8]);
            xticks([0 1 2 3 4 5 6 7 8])
            xticklabels({'','amplitude','pitch','AM','mean_frequency','FM','peak-frequency','entropy',''});
            title(sprintf('%u of  %u song elements',kemal,length(fraglist)));
            ylabel('Zscored Feature(averaged)');

            saveas(gcf,sprintf('New_BinaThres_LineChartMeanFeaturesVsResp-%s.png',neu.formated_name));
            close(gcf);
        end

        function saveDrawPairwiseFragmentsMeanFeaturesDistribution(neu)

            fraglist = neu.judgeFragResponse;
            if isempty(fraglist)
                return
            end
            if isempty(fraglist(1).meanfeatures) %%% Too bad !!
                return
            end
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).maxvalue;
            end

            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).amplitude = fraglist(u).meanfeatures.amplitude;
                toshow(u).pitch = fraglist(u).meanfeatures.pitch;
                toshow(u).AM = fraglist(u).meanfeatures.AM;
                toshow(u).mean_frequency = fraglist(u).meanfeatures.mean_frequency;
                toshow(u).FM = fraglist(u).meanfeatures.FM;
                toshow(u).peak_frequency = fraglist(u).meanfeatures.peak_frequency;
                toshow(u).entropy = fraglist(u).meanfeatures.entropy;
                toshow(u).harmratio = fraglist(u).meanfeatures.harmratio;
                toshow(u).resp = fraglist(u).responseMeasure;

            end

            name = {'amplitude','pitch','AM','mean_frequency','FM','peak_frequency','entropy','harmratio'};
            ncomb = nchoosek(name,2);

            I = {};
            for idx = 1:size(ncomb,1) % 此处可用parfor

                %subplot(hangshu,lieshu,idx);
                figure('Color','w');
                hold on
                cmap = 1 - rescale([toshow.resp].');
                for  md = 1: length(toshow)
                    scatter(toshow(md).(ncomb{idx,1})', toshow(md).(ncomb{idx,2}),[],repmat(cmap(md,:),1,3),'filled');
                end
                xlabel(replace(ncomb{idx,1},'_','-'));
                ylabel(replace(ncomb{idx,2},'_','-'));

                frame = getframe(gcf);
                I{idx} = frame.cdata;
                close(gcf);

            end

            % PCA plot
            neu.drawPCABasedOnAllFeatures;
            frame = getframe(gcf);
            I{length(I)+ 1} = frame.cdata;
            close(gcf);


            lieshu = 7;
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
            IMG = cell2mat(reshapedI);
            imwrite(IMG,sprintf('PairwiseFragmentsMeanFeaturesDistribution_%s.png',neu.formated_name));

        end

        function saveDrawDTWSimilarityMatrix(neu)

            temp = neu.frag.fraglist;
            fraglist = horzcat(temp{:});
            if isempty(fraglist)
                return;
            end
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).maxsdf;
            end
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).amplitude = fraglist(u).features.amplitude;
                toshow(u).pitch = fraglist(u).features.pitch;
                toshow(u).AM = fraglist(u).features.AM;
                toshow(u).mean_frequency = fraglist(u).features.mean_frequency;
                toshow(u).FM = fraglist(u).features.FM;
                toshow(u).peak_frequency = fraglist(u).features.peak_frequency;
                toshow(u).entropy = fraglist(u).features.entropy;
                toshow(u).resp = fraglist(u).responseMeasure;
            end

            toshow = toshow(~cellfun(@isempty,{toshow.resp}.')); % only keep those not empty rows  !!!! 后面要调查一下为什么有些会是maxsdfempty？ 2022.12.11
            toshow  = table2struct(sortrows(struct2table(toshow), 'resp','ascend'));


            names = {'amplitude','pitch','AM','FM','mean_frequency','peak_frequency','entropy'};

            I = {}; % inatialize the Image collection
            f = waitbar(0,'Start');

            for na = 1: length(names)
                waitbar(na/length(names),f,replace(sprintf('Generating %s similarity matrix ...',names{na}),'_','-'));
                figure('Color','w');
                sim = [];
                for s = 1: length(toshow)
                    parfor b = 1: length(toshow)
                        sim(s,b) = 1 - dtw(toshow(s).(names{na}),toshow(b).(names{na}));
                    end
                end
                imagesc(sim);
                xlabel(sprintf('%s-Similarity',names{na}));
                frame = getframe(gcf);
                I{na} = frame.cdata;
                close(gcf);
            end
            close(f);


            lieshu = 4;
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
            IMG = cell2mat(reshapedI);
            imwrite(IMG,sprintf('DTWSimilarityMatrix_%s.png',neu.info.formated_name));

        end

        function saveDrawDTWSimilarityMatrixBasedOnZscoredData(neu)

            temp = neu.frag.fraglist;
            fraglist = horzcat(temp{:});
            if isempty(fraglist)
                return;
            end
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).maxsdf;
            end
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).amplitude = zscore(fraglist(u).features.amplitude);
                toshow(u).pitch = zscore(fraglist(u).features.pitch);
                toshow(u).AM = zscore(fraglist(u).features.AM);
                toshow(u).mean_frequency = zscore(fraglist(u).features.mean_frequency);
                toshow(u).FM = zscore(fraglist(u).features.FM);
                toshow(u).peak_frequency = zscore(fraglist(u).features.peak_frequency);
                toshow(u).entropy = zscore(fraglist(u).features.entropy);
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            toshow = toshow(~cellfun(@isempty,{toshow.resp}.')); % only keep those not empty rows  !!!! 后面要调查一下为什么有些会是maxsdfempty？ 2022.12.11

            toshow  = table2struct(sortrows(struct2table(toshow), 'resp','ascend'));


            names = {'amplitude','pitch','AM','FM','mean_frequency','peak_frequency','entropy'};

            I = {}; % inatialize the Image collection
            f = waitbar(0,'Start');

            for na = 1: length(names)
                waitbar(na/length(names),f,replace(sprintf('Generating %s similarity matrix ...',names{na}),'_','-'));
                figure('Color','w');
                sim = [];
                for s = 1: length(toshow)
                    parfor b = 1: length(toshow)
                        sim(s,b) = 1 - dtw(toshow(s).(names{na}),toshow(b).(names{na}));
                    end
                end
                imagesc(sim);
                xlabel(sprintf('Zscored-DTWSimilarityMatrix-%s',names{na}));
                frame = getframe(gcf);
                I{na} = frame.cdata;
                close(gcf);
            end
            close(f);


            lieshu = 4;
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
            IMG = cell2mat(reshapedI);
            imwrite(IMG,sprintf('ZscoredFeatureSimilarity_%s.png',neu.info.formated_name));

        end

        function saveDrawCoeffOfFeaturesLinearFit(neu)

            fraglist = neu.judgeFragResponse;
            if isempty(fraglist)
                return
            end
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).maxvalue;
            end
            toshow = struct;
            for u = 1: length(fraglist)
                % for polyfits
                commonx = 1: length(fraglist(u).features.amplitude);
                toshow(u).amplitude_coef = polyfit(commonx, fraglist(u).features.amplitude,1);
                toshow(u).pitch_coef = polyfit(commonx, fraglist(u).features.pitch,1);
                toshow(u).AM_coef = polyfit(commonx, fraglist(u).features.AM,1);
                toshow(u).mean_frequency_coef = polyfit(commonx, fraglist(u).features.mean_frequency,1);
                toshow(u).FM_coef = polyfit(commonx, fraglist(u).features.FM,1);
                toshow(u).peak_frequency_coef = polyfit(commonx, fraglist(u).features.peak_frequency,1);
                toshow(u).entropy_coef = polyfit(commonx, fraglist(u).features.entropy,1);
                toshow(u).resp = fraglist(u).responseMeasure;

            end
            toshow  = table2struct(sortrows(struct2table(toshow), 'resp','ascend'));


            names = {'amplitude_coef','pitch_coef','AM_coef','FM_coef','mean_frequency_coef','peak_frequency_coef','entropy_coef'};
            I = {}; % inatialize the Image collection
            f = waitbar(0,'Start');

            for na = 1: length(names)
                waitbar(na/length(names),f,replace(sprintf('Creating %s scatter plot ...',names{na}),'_','-'));
                figure('Color','w');
                hold on
                feature_coefs = {toshow.(names{na})}.';
                feature_coefs = vertcat(feature_coefs{:});
                resp = [toshow.resp].';

                scatter(feature_coefs(:,1),feature_coefs(:,2),[],resp,'filled');
                xline(0,'--r');

                xlabel('1-th order coefficient','Interpreter','none');
                ylabel('0-th order coefficient','Interpreter','none');

                title(replace(sprintf('%s-1 and 0 coefficient',names{na}),'-','_'));
                frame = getframe(gcf);
                I{na} = frame.cdata;
                close(gcf);
            end
            close(f);


            lieshu = 4;
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
            IMG = cell2mat(reshapedI);
            imwrite(IMG,sprintf('N-order-coefficient_%s.png',neu.formated_name));



        end


        function drawPCABasedOnAllFeatures(neu)

            fraglist = neu.judgeFragResponse;
            if isempty(fraglist)
                return;
            end
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).maxvalue;
            end

            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).amplitude = fraglist(u).meanfeatures.amplitude;
                toshow(u).pitch = fraglist(u).meanfeatures.pitch;
                toshow(u).AM = fraglist(u).meanfeatures.AM;
                toshow(u).mean_frequency = fraglist(u).meanfeatures.mean_frequency;
                toshow(u).FM = fraglist(u).meanfeatures.FM;
                toshow(u).peak_frequency = fraglist(u).meanfeatures.peak_frequency;
                toshow(u).entropy = fraglist(u).meanfeatures.entropy;
                toshow(u).resp = fraglist(u).responseMeasure;
            end

            toshow  = table2struct(sortrows(struct2table(toshow), 'resp','ascend'));


            num_toshow = [toshow.amplitude; toshow.pitch; toshow.AM; toshow.mean_frequency; toshow.FM; ...
                toshow.peak_frequency; toshow.entropy].';

            [coe,score,~] = pca(num_toshow);

            cmap = 1- rescale( [toshow.resp].');
            % 1- to plot the data
            figure('Color','white');
            hold on
            scatter(score(:,1),score(:,2),[],repmat(cmap,1,3),'filled');
            %colorbar;
            xlabel('PCA Dim-1');
            ylabel('PCA Dim-2');
            %title('PCA data');

        end


        function img = saveDrawSortedRespToFrags(neu) %非常难以找到这个function
            %非常难以找到这个function

            neu.judgeFragResp_FR;
            %             finalfig = getframe(gcf).cdata;
            %             close(gcf)
            ids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'Frag|frag|syl|ele')));
            fraglist = neu.list(ids);
            if isempty(fraglist)
                return
            end
            %sorted_fraglist = table2struct(sortrows( struct2table(fraglist) ,'maxvalue','descend'));
            sorted_fraglist = table2struct(sortrows( struct2table(fraglist) ,'maxsdf','descend'));


            I = {}; % collection of frag-response-three images
            for k = 1: length(sorted_fraglist)

                if sorted_fraglist(k).label == 0
                    h =  figure('Position',[681 403 523 696],'Color','w');
                elseif sorted_fraglist(k).label == 1
                    h =  figure('Position',[681 403 523 696],'Color','y');
                end

                %h.WindowState = 'maximized';
                Draw.two(sorted_fraglist(k).plty,sorted_fraglist(k).fs,sorted_fraglist(k).pltsptimes);
                xlabel(sprintf('%s-maxsdf: %f',sorted_fraglist(k).stimuliname,sorted_fraglist(k).maxsdf));
                temp = getframe(gcf);
                I{k} = temp.cdata;


                close(h)
            end

            %             I{k+1} = finalfig; % performance figure

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
                end
            end

            reshapedI = reshape(I, lieshu,[])';
            clear I;
            img = cell2mat(reshapedI);
            imwrite(img,sprintf('赵RespToFrags_%s.png',neu.formated_name));

        end



    end

    methods(Static)

        function result = getNumTestedNumResponsive(A,bingyang)
            temp = A.knowTarget;
            locallist  = A.frag.allfraglist_nonorm;
            if ~isempty(temp)
                targetname = temp{1};% 这个{1}也似乎比较武断
                % 找到targetid
                targetid = find(~cellfun(@isempty,regexp(cellstr({locallist.stimuliname}.'),targetname)));
            else
                [~,targetid] = max([locallist.maxsdf].');

            end

            catego_target = MetaStimuli.categorizeByBingyang(locallist( targetid).y,32000, bingyang);
            purelocallist = locallist; % 不去掉target本身

            for k = 1:length(purelocallist)
                purelocallist(k).catego = MetaStimuli.categorizeByBingyang(purelocallist(k).y,32000,bingyang);
            end

            % 同组
            result.num_tested_simgroup = length(find([purelocallist.catego].'==catego_target));
            result.num_responsive_simgroup = length(find([purelocallist.catego].'==catego_target & [purelocallist.label].'==1));

            % 异组
            result.num_tested_difgroup = length(find([purelocallist.catego].'~=catego_target));
            result.num_responsive_difgroup = length(find([purelocallist.catego].'~=catego_target & [purelocallist.label].'==1));

            %所有组
            result.num_tested_allgroup = length(find([purelocallist.catego].'));
            result.num_responsive_allgroup = length(find([purelocallist.label].'==1));

            %比例
            result.ratio_sim = result.num_responsive_simgroup/result.num_tested_simgroup;
            result.ratio_dif = result.num_responsive_difgroup/result.num_tested_difgroup;
            result.ratio_all = result.num_responsive_allgroup/result.num_tested_allgroup;


        end

        function sorted_presented_eleinf = findDist2Target(presented_eleinf,targetname,all_eleinf)
            % 根据与target之间的距离对eleinf进行排序
            
            % presented_eleinf = A.frag.fraglist;
            % targetname = A.knowTarget{1}
            % all_eleinf = all_eleinf

            %首先把all_eleinf里的坐标值全部赋予presented_eleinf
            for k = 1:length(presented_eleinf)
                temp = erase(presented_eleinf(k).stimuliname,'Frag-');
                purename = erase(temp,regexp(presented_eleinf(k).stimuliname,'-\d+Pulses','match'));
                songname = erase(purename,regexp(purename,'-\d+$','match'));
                fragid = str2num(regexp(purename,'(?<=-)\d+$','match'));

                hitsongname = find(strcmp(songname,[all_eleinf.songname].'));
                hitfragid = find(eq(fragid,[all_eleinf.fragid].'));
                hitfinal = intersect(hitsongname,hitfragid);

                presented_eleinf(k).coor = [all_eleinf(hitfinal).coor_1,all_eleinf(hitfinal).coor_2];

            end

            % 找到target并计算所有presented syllables与target之间的距离
            targetid = find(~cellfun(@isempty,regexp([presented_eleinf.stimuliname].',targetname)));
            if length(targetid)> 1
                disp('@Frag.findDist2Target---presented frag stimuli可能有多个')
                targetid = targetid(1); % targetid可能有多个，这时只取第一个
            end
            for k = 1:length(presented_eleinf)

                try
                    presented_eleinf(k).dist_with_target = norm(presented_eleinf(k).coor-presented_eleinf(targetid).coor);
                catch
                    presented_eleinf(k).dist_with_target = inf;  % 有的时候这个值不存在，那就设为inf
                end

               

            end


            [~,index] = sortrows([presented_eleinf.dist_with_target].');
            sorted_presented_eleinf = presented_eleinf(index); % close id 是从相似度由高到低排列的
            clear index


        end

    end

end