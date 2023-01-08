classdef Neuron_basicDrawings <handle
    
    %Axillary functions of Neuron
    
    
    methods
        function drawAllFigures(neu)
            
            img1 = A.Three;
            img2 = A.saveDrawAlignedConsDegs;
            img3 = A.saveDrawAlignedConsFrag;
            img4 = A.saveDrawAlignedConsReplas;
            img5 = A.saveDrawSortedRespToFrags;
            img234 = padFigures({img2,img3,img4});
            allimg = padFigures({img1,img234,img5});
            imwrite(allimg,sprintf('%s.png',neu.formated_name));
            
            % 首先生成一张大图
            % 之后可以生成一个报告
            %生成的A文件需要有自我更新的能力
            %             A.drawAllWaveform;
            %             A.drawFragScatter
            % probably draw it into one big figure
            
        end
        
        function Three(neu,keyword)
            % keyword可以是norm
            
            % draw three plot, using plt-data
            es = getAllEphysObject(neu);
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

            figure('Position',PM.size1,'color','w');
            neu.waveform.draw1st;     % draw 第一个 waveform
            frame = getframe(gcf);
            I{length(I)+ 1} = frame.cdata;
            close(gcf);

            figure('Position',PM.size1,'color','w');
            neu.waveform.drawAll;    % draw 所有的 waveform
            frame = getframe(gcf);
            I{length(I)+ 1} = frame.cdata;
            close(gcf);

            figure('Position',PM.size1,'color','w');
            neu.waveform.drawSeparated;   % draw 所有的 waveform
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
            imwrite(img,sprintf('齐Three_%s.png',neu.info.formated_name));
            
        end
           
        function Deprecated_drawFirstWaveform(neu)
            
            % temporialriy neu.neurons{1}
            waveforms = neu.experiments{1}.waves.waveform;

            hold on
            plot(waveforms.',':','Color',[.5,.5,.5]);
            plot(max(waveforms),'--','Color','blue');
            plot(min(waveforms),'--','Color','blue');
            plot(mean(waveforms),'Color','red');
            
        end
        
        function Deprecated_drawAllWaveform(neu)
            % this function is used to draw waveform
            % what this works for concatenating all neurons together???
            for k = 1: length(neu.neurons)
                waveforms{k} = neu.neurons{k}.waveform;
            end
            concat_waveforms = vertcat(waveforms{:});
            figure('Color','w');
            hold on
            plot(concat_waveforms.',':','Color',[.5,.5,.5]);
            plot(max(concat_waveforms),'--','Color','blue');
            plot(min(concat_waveforms),'--','Color','blue');
            plot(mean(concat_waveforms),'Color','red');
            
        end
        
        function Deprecated_drawSeparatedWaveform(neu,handle)
            
            fig = figure('Color','w');
            cmap = colormap(flip(hsv(5)));
            for k = 1:length(neu.experiments)
                local_waveform = neu.experiments{k}.info.waveform;
                hold on
                plot(local_waveform.',':','Color',[cmap(k,:),0.3]);
                plot(max(local_waveform),'--','Color',cmap(k,:)*0.7);
                plot(min(local_waveform),'--','Color',cmap(k,:)*0.7);
                plot(mean(local_waveform),'-*','Color',cmap(k,:)*0.5);
                
            end
            colorbar('southoutside')
            
            title(sprintf('%u plexon files',length(neu.experiments)));
            if exist('handle','var')
                saveas(gcf,sprintf('WaveformsSeparated-%s.png', neu.info.formated_name));
            end
        end
        

        function drawSelectedStimuliResp(neu,stimulinames,range)
            % range是从plty的起始为始，以秒计
            dbstop if error; tic
            if isa(stimulinames,'cell')
                selectedids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),strjoin(stimulinames,'|') )));
            else
                selectedids = stimulinames;
            end
            
            
            songids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'norm|spe') ));
            songids = intersect(songids,selectedids);
            songlist = neu.list(songids);
            %             [~,postunique] = unique(cellstr(cellfun(@Convert.bid,{songlist.stimuliname}.','Uni',0)))
            %             songlist = songlist(postunique);
            RONGYU = 0.5;
            range = range + RONGYU;
            for k = 1: length(songlist)
                songlist(k).pady = [zeros(RONGYU*songlist(k).fs,1);songlist(k).plty];
                songlist(k).padsptimes = cellfun( @(x) x + RONGYU, songlist(k).pltsptimes,'uni',0);
            end
            
            % About Frag
            fragids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'frag|Frag|syl|Syl') ));
            fragids = intersect(fragids,selectedids);
            if ~isempty(fragids)
                fraglist = neu.list(fragids);
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
            degids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'Deg|deg') ));
            degids = intersect(degids,selectedids);
            if ~isempty(degids)
                deglist = neu.list(degids);
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
            replaids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'Repla|repla|catego|Catego') ));
            replaids = intersect(replaids,selectedids);
            if ~isempty(replaids)
                replalist = neu.list(replaids);
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
        
        function drawStimuliOnly(neu,name_or_id,range)
            % range是从plty的起始为始，以秒计
            dbstop if error; tic
            if isnumeric(name_or_id)
                selectedid = name_or_id;
            else
                selectedid = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),name_or_id )));
            end
            
            k = selectedid;
            if exist('range','var')
                truncated_y = neu.list(k).plty(range(1)*neu.list(k).fs:range(2)*neu.list(k).fs);
                figure('Position',[-36 437 length(truncated_y)/18 528]);
                
                
            else
                figure('Position',[-36 437 length(neu.list(k).plty)/18 528]);
                Draw.spec(neu.list(k).plty,neu.list(k).fs);
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
        
        function deg_songnames = which_songs_are_targeted_for_further_modification(neu)
            
            degids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'deg|Deg')));
            
            deglist = neu.list(degids);
            
            deg_songnames = unique(cellstr(regexp([deglist.stimuliname].','[OGBYR]\d{3}','match')));
            
            
        end
        
 
        function IMG = draw_pltthree_subfile(neu,neuronid)
            % 画此神经元的subfile对应的three plots
            dbstop if error
            n = neu.neurons{neuronid};
            
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
            imwrite(IMG,sprintf('RawThree_%s.png',neu.formated_name));
            
            
        end
        
        
    end
end

