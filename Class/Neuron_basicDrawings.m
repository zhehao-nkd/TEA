classdef Neuron_basicDrawings <handle
    
    %Axillary functions of Neuron
    
    
    methods
        function drawAllFigures(a)
            
            img1 = A.Three;
            img2 = A.saveDrawAlignedConsDegs;
            img3 = A.saveDrawAlignedConsFrag;
            img4 = A.saveDrawAlignedConsReplas;
            img5 = A.saveDrawSortedRespToFrags;
            img234 = padFigures({img2,img3,img4});
            allimg = padFigures({img1,img234,img5});
            imwrite(allimg,sprintf('%s.png',a.formated_name));
            
            % 首先生成一张大图
            % 之后可以生成一个报告
            %生成的A文件需要有自我更新的能力
            %             A.drawAllWaveform;
            %             A.drawFragScatter
            % probably draw it into one big figure
            
        end
        function Three(a,keyword)
            
            % draw three plot, using plt-data
            es = getAllEphysObject(a);
            getSoundname = @(x) x.sound.name;
            names = {};
            for k = 1:length(es)
                names{k} = getSoundname(es{k});
            end

            if exist('keyword','var')
                validids = find(~cellfun(@isempty,regexp(cellstr(names),keyword)));
                valid_es = es(validids);
            else
                valid_es = es;
            end
            
            for idx = 1: length(valid_es)
                valid_es{idx}.pltthree;
                frame = getframe(gcf);
                I{idx} = frame.cdata;
                close(gcf);
            end

            a.drawFirstWaveform;     % draw waveform
            frame = getframe(gcf);
            I{length(I)+ 1} = frame.cdata;
            close(gcf);
            
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
            imwrite(img,sprintf('Three_%s.png',a.formated_name));
            
        end
        
        
        function drawFirstWaveform(a)
            
            % temporialriy a.neurons{1}
            waveforms = a.neurons{1}.waveform;
            figure('Color','w','Position',PM.size1);
            hold on
            plot(waveforms.',':','Color',[.5,.5,.5]);
            plot(max(waveforms),'--','Color','blue');
            plot(min(waveforms),'--','Color','blue');
            plot(mean(waveforms),'Color','red');
            
        end
        
        function drawAllWaveform(a)
            % this function is used to draw waveform
            % what this works for concatenating all neurons together???
            for k = 1: length(a.neurons)
                waveforms{k} = a.neurons{k}.waveform;
            end
            concat_waveforms = vertcat(waveforms{:});
            figure('Color','w');
            hold on
            plot(concat_waveforms.',':','Color',[.5,.5,.5]);
            plot(max(concat_waveforms),'--','Color','blue');
            plot(min(concat_waveforms),'--','Color','blue');
            plot(mean(concat_waveforms),'Color','red');
            
        end
        
        function drawSeparatedWaveform(a)
            
            fig = figure('Color','w');
            cmap = colormap(flip(hsv(5)));
            for k = 1:length(a.neurons)
                local_waveform = a.neurons{k}.waveform;
                hold on
                plot(local_waveform.',':','Color',[cmap(k,:),0.3]);
                plot(max(local_waveform),'--','Color',cmap(k,:)*0.7);
                plot(min(local_waveform),'--','Color',cmap(k,:)*0.7);
                plot(mean(local_waveform),'-*','Color',cmap(k,:)*0.5);
                
            end
            colorbar('southoutside')
            
            title(sprintf('%u plexon files',length(a.neurons)));
            saveas(gcf,sprintf('WaveformsSeparated-%s.png', a.formated_name));
        end
        
        function drawFragScatter(a,not_tested_handle)
            % not_tested_handle = 1 means draw
            
            if ~exist('not_tested_handle','var') % to judge whether to draw not tested elements or not
                not_tested_handle = 0;
            end
            
            global_eleinf = a.all_eleinf;
            fraglist = a.judgeFragResponse;
            
            
            for k = 1: length(global_eleinf)
                uniqueid = sprintf('-%s-%u-',global_eleinf(k).songname, global_eleinf(k).fragid);
                
                if ~isempty (find(~cellfun(@isempty, regexp({fraglist.stimuliname}.',uniqueid)), 1)) % if this element was tested
                    id_in_fraglist = find(~cellfun(@isempty, regexp({fraglist.stimuliname}.',uniqueid)));
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
            figure
            hold on
            for k = 1: length(global_eleinf)
                if global_eleinf(k).scatter == -1
                    if not_tested_handle == 1
                        scatter(global_eleinf(k).coor_1,global_eleinf(k).coor_2,[],'k','filled'); % black for not tested
                    elseif not_tested_handle == 0 % if not draw test handle, do nothing
                    end
                elseif global_eleinf(k).scatter == 0
                    scatter(global_eleinf(k).coor_1,global_eleinf(k).coor_2,[],'g','filled');   % green for tested but not response-eliciting
                    % text(double(global_eleinf(k).coor_1),double(global_eleinf(k).coor_2),sprintf('%s-%u',global_eleinf(k).songname, global_eleinf(k).fragid) )
                elseif global_eleinf(k).scatter == 1
                    scatter(global_eleinf(k).coor_1,global_eleinf(k).coor_2,[],'r','filled');
                    % text(double(global_eleinf(k).coor_1),double(global_eleinf(k).coor_2),sprintf('%s-%u',global_eleinf(k).songname, global_eleinf(k).fragid) )
                end
                
                
            end
            
            
            % label the targets
            global_names = [global_eleinf.songname].';
            global_fragids = [global_eleinf.fragid].';
            
            global_merged = {};
            for w = 1: length(global_names)
                global_merged{w} = sprintf('%s-%u',global_names(w),global_fragids(w));
            end
            global_merged = global_merged.';
            
            
            for u = 1: length(a.targets)
                [~,beta] = ismember( a.targets{u}, global_merged );
                scatter( global_eleinf(beta).coor_1,global_eleinf(beta).coor_2,[],'k','h');
            end
            hold off
            
            
        end
        
        function drawPCABasedOnAllFeatures(a)
            
            fraglist = a.judgeFragResponse;
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
        
        function saveDrawSSIMSimlarityMatrix(a)
            % This function firstly order the frags by response strength
            % then measure the pairwise similarity between elements
            fraglist = a.judgeFragResponse;
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
            
            saveas(gcf,sprintf('SSIMSimlarityMatrix-%s.png',a.formated_name));
            close(gcf);
        end
        
        function saveDrawMeanFeatureVsResp(a)
            % draw the distribution of mean features
            
            %             fragids = find(~cellfun(@isempty, regexp({a.list.stimuliname}.','frag|Frag|syl|Syl')));
            %
            %             if isempty(fragids)
            %                 return
            %             end
            %             fraglist = a.list(fragids);
            
            fraglist = a.judgeFragResponse;
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
            
            
            saveas(gcf, sprintf('SeparatedFeatureVsResponse_%s.png',a.formated_name));
            close(gcf);
            
            
        end
        
        function saveDrawMeanFeaturesVsRespAsLineChart(a)
            % draw the distribution of mean features
            a.calHarmRatio;
            fraglist = a.judgeFragResponse;
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
            
            saveas(gcf,sprintf('LineChartMeanFeaturesVsResp-%s.png',a.formated_name));
            close(gcf);
        end
        
        function drawSelectedStimuliResp(a,stimulinames,range)
            % range是从plty的起始为始，以秒计
            dbstop if error; tic
            if isa(stimulinames,'cell')
                selectedids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),strjoin(stimulinames,'|') )));
            else
                selectedids = stimulinames;
            end
            
            
            songids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'norm|spe') ));
            songids = intersect(songids,selectedids);
            songlist = a.list(songids);
            %             [~,postunique] = unique(cellstr(cellfun(@Convert.bid,{songlist.stimuliname}.','Uni',0)))
            %             songlist = songlist(postunique);
            RONGYU = 0.5;
            range = range + RONGYU;
            for k = 1: length(songlist)
                songlist(k).pady = [zeros(RONGYU*songlist(k).fs,1);songlist(k).plty];
                songlist(k).padsptimes = cellfun( @(x) x + RONGYU, songlist(k).pltsptimes,'uni',0);
            end
            
            % About Frag
            fragids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'frag|Frag|syl|Syl') ));
            fragids = intersect(fragids,selectedids);
            if ~isempty(fragids)
                fraglist = a.list(fragids);
                for m = 1: length(fraglist)
                    birdid = Convert.bid(fraglist(m).stimuliname);
                    ids_norm = find(~cellfun(@isempty, regexp(cellstr({songlist.stimuliname}.'),birdid) ) );
                    if ~isempty(ids_norm)& length(ids_norm) == 1
                        fraglist(m).sylIni = Neuron.findIni(songlist(ids_norm).pady,fraglist(m).y);
                        fraglist(m).pady = [zeros([fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext,1]);fraglist(m).plty;zeros(length(songlist(ids_norm).plty)...
                            - (fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext) - length(fraglist(m).plty),1)];
                        fraglist(m).padsptimes = cellfun( @(x) x + (fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext)/fraglist(m).fs, fraglist(m).pltsptimes,'uni',0);
                    end
                end
            end
            
            % About Deg
            degids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'Deg|deg') ));
            degids = intersect(degids,selectedids);
            if ~isempty(degids)
                deglist = a.list(degids);
                for m = 1: length(deglist)
                    birdid = Convert.bid(deglist(m).stimuliname);
                    ids_norm = find(~cellfun(@isempty, regexp(cellstr({songlist.stimuliname}.'),birdid) ) );
                    if ~isempty(ids_norm)& length(ids_norm) == 1
                        [deglist(m).sylIni,trump_diffvalue] = Neuron.findIni(songlist(ids_norm).pady,deglist(m).y);
                        fprintf('%s has the DIFFVALUE as: %u and INI as: %f \n ',deglist(m).stimuliname,trump_diffvalue,deglist(m).sylIni/32000);
                        deglist(m).pady = [zeros([deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext,1]);deglist(m).plty;zeros(length(songlist(ids_norm).plty)...
                            - (deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext) - length(deglist(m).plty),1)];
                        deglist(m).padsptimes = cellfun( @(x) x + (deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext)/deglist(m).fs, deglist(m).pltsptimes,'uni',0);
                    end
                end
            end
            
            % About Repla
            replaids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'Repla|repla|catego|Catego') ));
            replaids = intersect(replaids,selectedids);
            if ~isempty(replaids)
                replalist = a.list(replaids);
                ids_norm_collecting_box = [];
                for m = 1: length(replalist)
                    
                    afterBefore = regexp(replalist(m).stimuliname,'(?<=before-)\S*','match');
                    afterBefore = afterBefore{1};
                    birdid = Convert.bid(afterBefore);
                    ids_norm = find(~cellfun(@isempty, regexp(cellstr({songlist.stimuliname}.'),birdid) ) );
                    
                    if ~isempty(ids_norm)%& length(ids_norm) == 1
                        ids_norm_1st = ids_norm(1);
                        ids_norm_collecting_box = [ids_norm_collecting_box, ids_norm];
                        afterpad_length = length(songlist(ids_norm_1st).plty) + RONGYU*songlist(ids_norm_1st).fs;  % +0.5s
                        
                        replalist(m).pady = [zeros(afterpad_length- length(replalist(m).plty),1);replalist(m).plty];
                        replalist(m).padsptimes = cellfun( @(x) x + length(zeros(afterpad_length- length(replalist(m).plty),1))/fraglist(m).fs,replalist(m).pltsptimes,'uni',0);
                    end
                end
            end
            
            % merge the Cons,Degs,Frags,Replas lists together
            for w = 1: length(songlist)
                
                if ~isempty(degids)
                    birdid = Convert.bid(songlist(w).stimuliname);
                    ids_indeg = find(~cellfun(@isempty, regexp(cellstr({deglist.stimuliname}.'),birdid) ) );
                    selected_deglist = deglist(ids_indeg);
                    [~,temp_index] = sortrows([selected_deglist.sylIni].');
                    selected_deglist = selected_deglist(temp_index);
                else
                    selected_deglist = [];
                end
                
                if ~isempty(fragids)
                    birdid = Convert.bid(songlist(w).stimuliname);
                    ids_infrag = find(~cellfun(@isempty, regexp(cellstr({fraglist.stimuliname}.'),birdid) ) );
                    selected_fraglist = fraglist(ids_infrag);
                    [~,temp_index] = sortrows([selected_fraglist.sylIni].');
                    selected_fraglist = selected_fraglist(temp_index);
                else
                    selected_fraglist = [];
                end
                
                if ~isempty(replaids)
                    birdid = Convert.bid(songlist(w).stimuliname);
                    ids_inrepla = find(~cellfun(@isempty, regexp(cellstr({replalist.stimuliname}.'),birdid) ) );
                    selected_replalist = replalist(ids_inrepla);
                else
                    selected_replalist = [];
                end
                
                % draw figures
                if ~isempty(selected_deglist)
                    selected_deglist = rmfield(selected_deglist,'sylIni');
                end
                
                if ~isempty(selected_fraglist)
                    selected_fraglist = rmfield(selected_fraglist,'sylIni');
                end
                alllist = horzcat(songlist(w),selected_deglist,selected_fraglist,selected_replalist);
                len = length(alllist);
                figure('Position',[1935 -207 520 1458/9*len],'Color','none');
                
                ax = tight_subplot(2*len, 1, 0.002, 0.02, 0);
                for k = 1: len
                    
                    axes(ax(2*(k-1)+ 1)); % Draw.two(,,);
                    %ax(2*(k-1)+ 1).Position(4) =  ax(2*(k-1)+ 1).Position(4);
                    if exist('range','var')
                        truncated_y = alllist(k).pady(range(1)*alllist(k).fs:range(2)*alllist(k).fs);
                        Draw.spec(truncated_y,alllist(k).fs);
                    else
                        Draw.spec(alllist(k).pady,alllist(k).fs);
                    end
                    xlabel('')
                    ylabel('')
                    
                    set(gca,'TickLength',[0 .01])
                    set(gca,'Yticklabel',[])
                    set(gca,'Xticklabel',[])
                    
                    axes(ax(2*k));
                    
                    if exist('range','var')
                        truncated_sptimes = Extract.sptimes_resetSP(alllist(k).padsptimes, range(1), range(2));
                        truncated_y = alllist(k).pady(range(1)*alllist(k).fs:range(2)*alllist(k).fs);
                        %Draw.raster(truncated_sptimes,truncated_y,alllist(k).fs);
                        Draw.rasterBeta(truncated_sptimes,truncated_y,alllist(k).fs,2.8,'k');
                    else
                        %Draw.raster(alllist(k).padsptimes,alllist(k).pady,alllist(k).fs);
                        Draw.rasterBeta(alllist(k).padsptimes,alllist(k).pady,alllist(k).fs,2.8,'k');
                    end
                    
                    ylabel('')
                    set(gca,'TickLength',[0 .01])
                    
                    
                    
                    %                         set(gca,'Yticklabel',[])
                    %                         set(gca,'Xticklabel',[])
                    %                     else
                    %                         xlabel(alllist(k).stimuliname);
                    %                     end
                    
                end
                
                for k = 1:len
                    ax(2*k-1).Position(4) =  ax(2*k-1).Position(4) -0.008;
                    ax(2*k-1).Position(2) =  ax(2*k-1).Position(2)+0.006;
                    box(ax(2*k-1),'off');
                    set(ax(2*k-1),'XColor','none','YColor','none')
                    %set(ax(2*k-1), 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
                end
                
                for k = 1:len
                    ax(2*k).Position(2) =  ax(2*k).Position(2)+0.000;
                    ax(2*k).Position(4) =  ax(2*k).Position(4)+0.002;
                    box(ax(2*k),'off');
                    set(ax(2*k),'XColor','none','YColor','none')
                    %set(ax(2*k), 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
                end
                disp('pause here')
                
                
                
                
            end
        end
        
        function drawStimuliOnly(a,name_or_id,range)
            % range是从plty的起始为始，以秒计
            dbstop if error; tic
            if isnumeric(name_or_id)
                selectedid = name_or_id;
            else
                selectedid = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),name_or_id )));
            end
            
            k = selectedid;
            if exist('range','var')
                truncated_y = a.list(k).plty(range(1)*a.list(k).fs:range(2)*a.list(k).fs);
                figure('Position',[-36 437 length(truncated_y)/18 528]);
                
                
            else
                figure('Position',[-36 437 length(a.list(k).plty)/18 528]);
                Draw.spec(a.list(k).plty,a.list(k).fs);
            end
            xlabel('')
            ylabel('')
            
            set(gca,'TickLength',[0 .01])
            set(gca,'Yticklabel',[])
            set(gca,'Xticklabel',[])
            set(gca,'XColor','none','YColor','none')
            box off
            disp('pause')
            
            %                         set(gca,'Yticklabel',[])
            %                         set(gca,'Xticklabel',[])
            
            
        end
        
        function split_fraginf = Deprecated_calResponseToWithinSongFragsFromEleinf(a)
            % deprecated
            latency = 50*0.001; % 50 ms
            
            split_fraginf = struct;
            
            counts = 0;
            if isempty(a.normlist)
                split_fraginf = struct([]);
                return
            end
            
            for k = 1: length(a.normlist)
                %songname = regexp(a.normlist(k).stimuliname,'(CON|SPE|norm)-([BRGOY]\d{3}|Fcall|Mcall|WNS|HET|Het)','match');
                songname = regexp(a.normlist(k).stimuliname,'([BRGOY]\d{3}|Fcall|Mcall|WNS|HET|Het)','match');
                
                % use autoseg to segment the elements
                %                [rawy,fiy,I,syledge,eleedge] = main(fiy,fs,birdid,CONFIG)
                %
                eleids = find(~cellfun(@isempty, regexp([a.conspe_eleinf.songname].',songname)));
                local_eleinf = a.conspe_eleinf(eleids);
                
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
                    
                    split_fraginf(counts).padded_sptimes = Extract.sptimes(a.normlist(k).rawsptimes, split_fraginf(counts).initial...
                        ,split_fraginf(counts).terminal + latency);
                    
                    front_percentage = local_eleinf(m).yini/length(a.normlist(k).rawy);
                    back_percentage = local_eleinf(m).yter/length(a.normlist(k).rawy);
                    
                    featurenames = fieldnames(a.normlist(k).rawfeatures); % iterate for each feature
                    featurenames = setdiff(featurenames,{'file_index','file_name'});
                    for omega = 1: length(featurenames)
                        lenfeature = length(a.normlist(k).rawfeatures.(featurenames{omega}));
                        if  round(front_percentage*lenfeature)== 0
                            
                            split_fraginf(counts).(featurenames{omega}) = a.normlist(k).rawfeatures.(featurenames{omega})(...
                                1:round(back_percentage*lenfeature));
                        else
                            split_fraginf(counts).(featurenames{omega}) = a.normlist(k).rawfeatures.(featurenames{omega})(...
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
        
        function deg_songnames = which_songs_are_targeted_for_further_modification(a)
            
            degids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'deg|Deg')));
            
            deglist = a.list(degids);
            
            deg_songnames = unique(cellstr(regexp([deglist.stimuliname].','[OGBYR]\d{3}','match')));
            
            
        end
        
        
        function IMG = saveDrawlineChart_2D_Cumulative(a) % 秦
            
            img{1} = a.drawPitchHarmoLineChart;
            img{2} = a.draw2DPitchVsHarmo;
            img{3} = a.drawCumulativeFeatureDiff('pitch');
            img{4} = a.drawCumulativeFeatureDiff('mean_frequency');
            img{5} = a.drawCumulativeFeatureDiff('peak_frequency');
            img{6} = a.drawCumulativeFeatureDiff('FM');
            img{7} = a.drawCumulativeFeatureDiff('goodness');
            img{8} = a.drawCumulativeFeatureDiff('entropy');
            img{9} = a.drawCumulativeFeatureDiff('harmratio');
            img{10} = a.drawCumulativeFeatureDiff('amplitude');
            img{11} = a.drawCumulativeFeatureDiff('AM');
            
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
            imwrite(IMG,sprintf('秦lineChart_2D_Cumulative_%s.png',a.formated_name));
            %imwrite(img,sprintf('新lineChart_2D_Cumulative_%s.png',a.formated_name));
        end
        
        function saveDrawInsong_lineChart_2D_Cumulative(a)
            img1 = a.Insong_drawPitchHarmoLineChart;
            img2 = a.Insong_draw2DPitchVsHarmo;
            img3 = a.Insong_drawCumulativePitchDiff;
            img = horzcat(img1,img2,img3);
            imwrite(img,sprintf('汉Insong_lineChart_2D_Cumulative_%s.png',a.formated_name));
        end
        
        function BinaryThresholdDrawMeanFeaturesVsRespAsLineChart(a) % draw the distribution of mean features
            
            fraglist = a.judgeFragResponse;
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
            switch a.formated_name
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
            
            saveas(gcf,sprintf('New_BinaThres_LineChartMeanFeaturesVsResp-%s.png',a.formated_name));
            close(gcf);
        end
        
        function saveDrawPairwiseFragmentsMeanFeaturesDistribution(a)
            
            fraglist = a.judgeFragResponse;
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
            a.drawPCABasedOnAllFeatures;
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
            imwrite(IMG,sprintf('PairwiseFragmentsMeanFeaturesDistribution_%s.png',a.formated_name));
            
        end
        
        function saveDrawDTWSimilarityMatrix(a)
            
            fraglist = a.judgeFragResponse;
            if isempty(fraglist)
                return;
            end
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).maxvalue;
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
            imwrite(IMG,sprintf('DTWSimilarityMatrix_%s.png',a.formated_name));
            
        end
        
        function saveDrawDTWSimilarityMatrixBasedOnZscoredData(a)
            
            fraglist = a.judgeFragResponse;
            if isempty(fraglist)
                return;
            end
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).maxvalue;
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
            imwrite(IMG,sprintf('ZscoredFeatureSimilarity_%s.png',a.formated_name));
            
        end
        
        function saveDrawCoeffOfFeaturesLinearFit(a)
            
            fraglist = a.judgeFragResponse;
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
            imwrite(IMG,sprintf('N-order-coefficient_%s.png',a.formated_name));
            
            
            
        end
        
        function saveDrawMeanFeaturesInSongVsRespAsLineChart(a) % draw the distribution of mean features
            
            fraglist = a.calResponseToWithinSongFragsFromEleinf;
            if isempty(fraglist)
                return;
            end
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).maxvalue;
            end
            
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).amplitude = fraglist(u).mean_amplitude;
                toshow(u).pitch = fraglist(u).mean_pitch;
                toshow(u).AM = fraglist(u).mean_AM;
                toshow(u).mean_frequency = fraglist(u).mean_mean_frequency;
                toshow(u).FM = fraglist(u).mean_FM;
                toshow(u).peak_frequency = fraglist(u).mean_peak_frequency;
                toshow(u).entropy = fraglist(u).mean_entropy;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            
            toshow  = table2struct(sortrows(struct2table(toshow), 'resp','ascend'));
            
            
            num_toshow = [toshow.amplitude; toshow.pitch; toshow.AM; toshow.mean_frequency; toshow.FM; ...
                toshow.peak_frequency; toshow.entropy].';
            
            znum_toshow = zscore(num_toshow,0,1);
            %            figure
            %            plot(znum_toshow.');
            
            figure('Position',[1997 233 1388 658],'Color','w');
            hold on
            c = 1- rescale([toshow.resp].',0.1,1);
            for r = 1: length(znum_toshow)
                plot(znum_toshow(r,:),'Color',repmat(c(r),3,1));
            end
            
            colormap(flip(repmat(unique(c),1,3)))
            colorbar
            xlim([0,8]);
            xticks([0 1 2 3 4 5 6 7 8])
            xticklabels({'','amplitude','pitch','AM','mean_frequency','FM','peak-frequency','entropy',''});
            
            ylabel('Zscored Feature(averaged)');
            title(sprintf('Totally %u song elements',length(fraglist)));
            saveas(gcf,sprintf('V1-WithinSongsLineChartMeanFeaturesVsResp-%s.png',a.formated_name));
            close(gcf);
        end
        
        function saveDrawAlignedConsDegsReplasFrags(a)  % Align all together and draw
            dbstop if error
            tic
            
            % about songs been presented % Redudant code
            songlist = Neuron(a.neurons{a.song_id}).normlist;
            [~,postunique] = unique(cellstr(cellfun(@Convert.bid,{songlist.stimuliname}.','Uni',0)))
            songlist = songlist(postunique);
            
            
            % About Frag
            fragids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'frag|Frag|syl|Syl') ));
            if ~isempty(fragids)
                
                fraglist = a.list(fragids);
                %                 normlist = Neuron(a.neurons{a.song_id}).normlist;
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
            degids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'Deg|deg') ));
            if ~isempty(degids)
                
                deglist = a.list(degids);
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
                        % draw a figure for every four ote fragments
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
            
            a.drawFirstWaveform;
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
            
            imwrite(Iall,sprintf('Aligned_Normfrag_%s.png',a.formated_name));
            toc
            
            
        end
        
        function IMG = draw_pltthree_subfile(a,neuronid)
            % 画此神经元的subfile对应的three plots
            dbstop if error
            n = a.neurons{neuronid};
            
            for idx = 1: length(n.e)
                n.e{idx}.pltthree;
                frame = getframe(gcf);
                I{idx} = frame.cdata;
                close(gcf);
            end
            
            figure('Color','w');
            n.draw_waveform;     % draw waveform
            frame = getframe(gcf);
            I{length(I)+ 1} = frame.cdata;
            close(gcf);
            
            n.drawPC12;
            frame = getframe(gcf);
            I{length(I)+ 1} = frame.cdata;
            close(gcf);
            
            n.drawValleyVsTime;
            frame = getframe(gcf);
            I{length(I)+ 1} = frame.cdata;
            close(gcf);
            
            n.drawPC1VsTime;
            frame = getframe(gcf);
            I{length(I)+ 1} = frame.cdata;
            close(gcf);
            
            % draw blank white
            if length(I) > 30
                lieshu = 15;
            else
                lieshu = 4; % previously it was 6
            end
            hangshu = ceil(length(I)/lieshu);
            rest = lieshu*hangshu - length(I);
            white = uint8(255*ones(size(I{1})));
            
            if rest > 0
                for k = 1:rest
                    I = [I,white];
                end
            end
            
            reshapedI = reshape(I, lieshu,[])';
            %clear I
            IMG = cell2mat(reshapedI);
            imwrite(IMG,sprintf('RawThree_%s.png',a.formated_name));
            
            
        end
        
        
    end
end

