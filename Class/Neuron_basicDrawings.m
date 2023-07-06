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

