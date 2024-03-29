classdef SpikeWaveform
    %有关 Spike waveform的一切事宜，比如画spike waveform的图

    properties
        all % all the waveform
        separated % waveform of each individual experiment(不包括merged的情况）
        formated_name
        adfreq % 
        meanWL
        firstWL
    end

    methods
        function wv = SpikeWaveform(input,mode) % timesdurations是针对mergedfiles
            dbstop if error

            warning('@SpikeWaveform ：Not constructed yet');
            
            switch mode

                

                case 'single or merge'
                    wv.adfreq = input.adfreq;
                    waveforms = input.waveforms;
                    times = input.times;
                    timeedges = input.timeedges;
                    wv.all = waveforms;
                    
                    for k = 1:length(timeedges) % 这个计算有大问题!!
                        ids = find(times >=timeedges{k}(1)&times <timeedges{k}(2));
                        wv.separated{k} = waveforms(ids,:);
                    end

                    wv.meanWL = wv.calMeanWaveLength;
                case 'multiple'
                    wv.all = vertcat(input{:});
                    wv.separated = input;
                   % wv.adfreq = input{1}.adfreq;
            end

        

        end

        function draw1st(wv,handle)


            if exist('handle','var') && handle == 1
                figure('Color','w','Position',[2031 376 807 600]);
            end
            
            % temporialriy neu.neurons{1}
            waveforms = wv.separated{1};

            % 如果只取前0.8ms
            waveforms = waveforms(:,1:32);

            hold on
            plot(waveforms.','-','Color',[.5,.5,.5],'LineWidth',1);
            plot(max(waveforms),'--','Color','blue');
            plot(min(waveforms),'--','Color','blue');
            plot(mean(waveforms),'Color','red');
            xlabel('1st waveform')

            if exist('handle','var') && handle == 1
                xlabel(sprintf('1st_Waveform_%s', wv.formated_name));
                img = getframe(gcf).cdata;
                imwrite(img,sprintf('兖州FirstWaveform_%s.png', wv.formated_name)); 
            end
            
        end
        
        function drawAll(wv, handle)
            % this function is used to draw waveform
            % what this works for concatenating all neurons together???
            if exist('handle','var') && handle == 1
                figure('Color','w','Position',[2031 376 807 600]);
            end

            for k = 1: length(wv.separated)
                waveforms{k} = wv.separated{k};
                waveforms{k} = waveforms{k}(:,1:32);
            end
            concat_waveforms = vertcat(waveforms{:});
            %figure('Color','w');
            hold on
            plot(concat_waveforms.','-','Color',[.5,.5,.5],'LineWidth',1);
            plot(max(concat_waveforms),'--','Color','blue');
            plot(min(concat_waveforms),'--','Color','blue');
            plot(mean(concat_waveforms),'Color','red');
            xlabel('all waveforms')


            % 把需要的Text 信息标记出来


            % Set the text content
            textStr = sprintf('mWL:%f\nnumSpikes:%u',wv.calMeanWaveLength,size(concat_waveforms,1));

            % Get the current x and y axis limits
            xLimits = xlim(gca);
            yLimits = ylim(gca);

            % Calculate the position of the text
            textX = xLimits(1) + 0.05 * diff(xLimits);
            textY = yLimits(2) - 0.05 * diff(yLimits);

            % Add the text to the figure
            text(textX, textY, textStr,'FontSize', 15,'Color','b');


            
        end

        function drawSeparated(wv,handle)


            %fig = figure('Color','w');
            cmap = colormap(flip(hsv(7)));
            for k = 1:length(wv.separated)
                subplot(length(wv.separated),1,k)
                local_waveform = wv.separated{k};
                local_waveform = local_waveform(:,1:32);
                hold on
                plot(local_waveform.',':','Color',[cmap(k,:),0.3]);
                plot(max(local_waveform),'--','Color',cmap(k,:)*0.7);
                plot(min(local_waveform),'--','Color',cmap(k,:)*0.7);
                plot(mean(local_waveform),'-*','Color',cmap(k,:)*0.5);

            end
           % colorbar('southoutside')

            %title(sprintf('%u plexon files',length(neu.experiments)));
            if exist('handle','var')
                saveas(gcf,sprintf('WaveformsSeparated-%s.png', neu.info.formated_name));
            end
        end
        function completeDraw(wv)

            %把所有的waveform的信息都画出来

            
            figure('Color','w','Position',[2031 376 807 600]);

            
            for k = 1: length(wv.separated)
                waveforms{k} = wv.separated{k};
                waveforms{k} = waveforms{k}(:,1:32);
            end
            concat_waveforms = vertcat(waveforms{:});
            %figure('Color','w');

            x = [1:size(concat_waveforms.',1)]* (1/wv.adfreq)*1000;
            hold on
            plot(concat_waveforms.','-','Color',[.5,.5,.5],'LineWidth',1.2);
            plot(max(concat_waveforms),'--','Color','blue');
            plot(min(concat_waveforms),'--','Color','blue');
            plot(mean(concat_waveforms),'Color','red');
            xlabel('all waveforms (ms)')
            title(wv.formated_name,'Interpreter', 'none');

            % 把需要的Text 信息标记出来


            % Set the text content
            textStr = sprintf('mWL:%f\nnumSpikes:%u',wv.calMeanWaveLength,size(concat_waveforms,1));

            % Get the current x and y axis limits
            xLimits = xlim(gca);
            yLimits = ylim(gca);

            % Calculate the position of the text
            textX = xLimits(1) + 0.05 * diff(xLimits);
            textY = yLimits(2) - 0.05 * diff(yLimits);

            % Add the text to the figure
            text(textX, textY, textStr,'FontSize', 15,'Color','b');





            img = {};
            img{1} = getframe(gcf).cdata;
            close(gcf);

%             
% 
%             cmap = colormap(flip(hsv(7)));
%             for k = 1:length(wv.separated)
%                 figure('Color','w','Position',[2031 376 807 600]);
%                 local_waveform = wv.separated{k};
%                 local_waveform = local_waveform(:,1:32);
%                 hold on
% 
%                 if ~isempty(local_waveform) 
%                     plot(x,local_waveform.',':','Color',[cmap(k,:),0.3]);
%                     plot(x,max(local_waveform),'--','Color',cmap(k,:)*0.7);
%                     plot(x,min(local_waveform),'--','Color',cmap(k,:)*0.7);
%                     plot(x,mean(local_waveform),'-*','Color',cmap(k,:)*0.5);
%                     xlabel(sprintf('第%u个文件',k));
%                     img{1+k} = getframe(gcf).cdata;
%                     close(gcf);
%                 end
% 
%             end

            allimg = vertcat(img{:});

            imwrite(allimg,sprintf('黄州Waveform_%s.png', wv.formated_name));


        end



        function meanWL = calMeanWaveLength(wv)
            % 计算 meanWL，需要考虑到仪器的fs
            allwaveforms = vertcat(wv.separated{:});
            %waveforms =  n.waveform;
            [~,troughstime] = min(allwaveforms,[],2);
            wavlen_units = [];

            for k = 1: size(allwaveforms,1) % length is dangerous!!!!!
                this_wf = allwaveforms(k,:);
                [~,wavlen_units(k)] =  max(this_wf (troughstime(k):end));
            end



            meanWL =  mean(wavlen_units*(1/wv.adfreq)*1000); % ms


        end


        function  firstWL = calFirstWavelength(wv)

            % 只计算第一个文件的waveforms

            firstwaveforms = wv.separated{1};
            %waveforms =  n.waveform;
            [~,troughstime] = min(firstwaveforms,[],2);
            wavlen_units = [];

            for k = 1: size(firstwaveforms,1) % length is dangerous!!!!!
                this_wf = firstwaveforms(k,:);
                [~,wavlen_units(k)] =  max(this_wf (troughstime(k):end));
            end


            firstWL =  mean(wavlen_units*(1/wv.adfreq)*1000); % ms

        end
    end
end