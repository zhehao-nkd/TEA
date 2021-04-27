classdef draw
    
    methods(Static)
        
        function spec(y,fs)
            
            %             idx = find(y,1);
            %             zpt = idx/fs; % zp means Zero Point
            
            try
                spectrogram(y,hamming(512),512-round(fs/1e3),512,fs,'yaxis'); colorbar('off');
                
                axh = gca;                                             %add this
                axh.Toolbar = matlab.ui.controls.AxesToolbar();       %add this
                hold on;
            catch Err
                disp('Error!!')
            end
            
            ylabel('Hz','FontSize',18);
            xlabel('Time (s)');
            ax.XLim             = [0 length(y)/fs];
            
            ydata = get(gca,'Ylim');
            % line([idx/fs,idx/fs],[min(ydata),max(ydata)],'color','r');
            
        end
        
        function raster(sptimes,y,fs,color)
            %find initial of stimuli
            
            
            idx = find(y,1);
            
            %initial_timestamps = idx;
            % y = y - initial_timestamps;
            
            zpt = idx/fs;
            zpt = 0; %%%% close zpt compensation
            
            % calculation...
            hold on
            ylabel('Trials','FontSize',18);
            % calculation...
            for nUnits = 1:length(sptimes)
                
                spks            = sptimes{nUnits}';  % Get all spikes of respective trial
                
                xspikes         = repmat(spks,2,1);         % Replicate array
                yspikes      	= nan(size(xspikes));       % NaN array
                
                % if ~isempty(yspikes)
                yspikes(1,:) = nUnits-1;                % Y-offset for raster plot
                yspikes(2,:) = nUnits;
                % end
                if ~isempty(xspikes) % in case that xspikes is empty
                    plot(xspikes-zpt, yspikes, 'Color', color,'linewidth',1.3)
                end
                
            end
            
            %figure configuration
            
            
            ax.YLim             = [0 length(sptimes)];
            
            ax.XLabel.String  	= 'Time (s)';
            
            
            
            ydata = get(gca,'Ylim');
            line([0,0],[min(ydata),max(ydata)],'color','r');
            
            xlim([0-zpt length(y)/fs-zpt]);
        end
        
        
        function psth(sptimes)
            
            all = [];
            for iTrial = 1:length(sptimes)
                all             = [all; sptimes{iTrial}];               % Concatenate spikes of all trials
            end
            
            ax                  = subplot(4,1,2);
            nbins               = 100;
            h                   = histogram(all,nbins);
            h.FaceColor         = 'k';
            
            mVal                = max(h.Values)+round(max(h.Values)*.1);
            ax.XLim             = [-.2 .3];
            ax.YLim             = [0 mVal];
            ax.XTick            = [0 .2];
            ax.XLabel.String  	= 'Time [s]';
            ax.YLabel.String  	= 'Spikes/Bin';
            
            slength             = 500;                                  % Length of signal trace [ms]
            bdur                = slength/nbins;                        % Bin duration in [ms]
            nobins              = 1000/bdur;                            % No of bins/sec
            for iLab = 1:length(ax.YTickLabel)
                lab             = str2num(ax.YTickLabel{iLab});
                conv            = (lab / length(sptimes)) * nobins; 	% Convert to [Hz]: avg spike count * bins/sec
                newlabel{iLab}  = num2str(round(conv));                 % Change YLabel
            end
            ax.YTickLabel       = newlabel;
            ax.YLabel.String  	= 'Firing Rate [Hz]';
        end
        
        function sdf(y, fs, sptimes)
            resolution = 0.001;
            gausswidth = 0.02;
            msdf = cal.sdf(sptimes,y,fs,resolution,gausswidth);
            %figure; 
            x = linspace(0,length(msdf),length(msdf))*resolution;
            plot(x,msdf);
            xlim([0,max(x)])
        
        end
        
        function three(y,fs,sptimes)
            subplot(3,1,1)
            draw.spec(y,fs);
            subplot(3,1,2)
            draw.raster(sptimes, y, fs, 'k');
            subplot(3,1,3)
            draw.sdf(y,fs,sptimes);
        end
        
        function sixplot
        end
        
        
    end
    
end