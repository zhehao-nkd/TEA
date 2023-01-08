classdef SpikeWaveform
    %有关 Spike waveform的一切事宜，比如画spike waveform的图

    properties
        all % all the waveform
        separated % waveform of each individual experiment(不包括merged的情况）
    end

    methods
        function wv = SpikeWaveform(input,mode) % timesdurations是针对mergedfiles
            dbstop if error
            switch mode

                case 'single or merge'
                    waveforms = input.waveforms;
                    times = input.times;
                    timedurations = input.timedurations;
                    wv.all = waveforms;
                    for k = 1:length(timedurations)
                        if k == length(timedurations) % 如果是最后一段
                            ids = find(times >=timedurations(k));
                        else % 否则
                            ids = find(times >=timedurations(k)&times <timedurations(k+1));
                        end
                        wv.separated{k} = waveforms(ids,:);
                    end
                case 'multiple'
                    wv.all = vertcat(input{:});
                    wv.separated = input;
            end
           
        end

        function draw1st(wv)
            
            % temporialriy neu.neurons{1}
            waveforms = wv.separated{1};

            hold on
            plot(waveforms.',':','Color',[.5,.5,.5]);
            plot(max(waveforms),'--','Color','blue');
            plot(min(waveforms),'--','Color','blue');
            plot(mean(waveforms),'Color','red');
            xlabel('1st waveform')
            
        end
        
        function drawAll(wv)
            % this function is used to draw waveform
            % what this works for concatenating all neurons together???
            for k = 1: length(wv.separated)
                waveforms{k} = wv.separated{k};
            end
            concat_waveforms = vertcat(waveforms{:});
            %figure('Color','w');
            hold on
            plot(concat_waveforms.',':','Color',[.5,.5,.5]);
            plot(max(concat_waveforms),'--','Color','blue');
            plot(min(concat_waveforms),'--','Color','blue');
            plot(mean(concat_waveforms),'Color','red');
            xlabel('all waveforms')
            
        end

        function drawSeparated(wv,handle)

            %fig = figure('Color','w');
            cmap = colormap(flip(hsv(5)));
            for k = 1:length(wv.separated)
                subplot(length(wv.separated),1,k)
                local_waveform = wv.separated{k};
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

    end
end