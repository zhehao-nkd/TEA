classdef draw
    
    methods(Static)
        
        function sono(y,fs)
            plot(y);
        end
        
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
        
        function spec3(y,fs)
            [~,F,T,P] = spectrogram(y,hamming(1024),1024-round(fs/1e3),1024,32000,'yaxis');
            %imagesc(T, F, 10*log10(P+eps)) % add eps like pspectrogram does
            lightness = 10*log10(P+eps);
            %lightness = imadjust(rescale(lightness,0,1));
            lightness = histeq(rescale(lightness,0,1),1024);
            figure; imagesc(T,F,lightness);
            axis xy
            ylabel('Frequency (Hz)')
            xlabel('Time (s)')
            colorbar off
           % colormap('jet')
           % h = colorbar;
           % h.Label.String = 'Power/frequency (dB/Hz)';
        end
        
        function spec2(y,fs)
            [~,F,T,P] = spectrogram(y,hamming(512),512-round(fs/1e3),512,32000,'yaxis');
            imagesc(T, F, 10*log10(P+eps)) % add eps like pspectrogram does
            axis xy
            ylabel('Frequency (Hz)')
            xlabel('Time (s)') % force as seconds
            colorbar off
            colormap('jet')
           % h = colorbar;
           % h.Label.String = 'Power/frequency (dB/Hz)';
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
        
        function appspec(ax,y,fs)
            %myAxe = app.UIAxes;
            [S,F,T] = spectrogram(y,hamming(512),512-round(fs/1e3),512,fs);
            imagesc(ax, T, F, log(1+abs(S)) ); %plot the log spectrum
            set(ax,'YDir', 'normal');
            set(ax,'XLim',[min(T),max(T)]);
            set(ax,'YLim',[min(F),max(F)]);
            
            %title(myAxe,sprintf('%s-%u',app.sylinf(app.idx).sound,app.sylinf(app.idx+value).number));
        end
        
        function appraster(ax,sptimes,y,fs)
             
            
            idx = find(y,1);
            
            
            zpt = idx/fs;
            zpt = 0; %%%% close zpt compensation
            
            % calculation...
            
            % calculation...
            for nUnits = 1:length(sptimes)
                
                spks            = sptimes{nUnits}';  % Get all spikes of respective trial
                
                ypair = [nUnits-1,nUnits];
                
                if ~isempty(spks)
                    
                    for each = 1: length(spks)
                        xpair = repmat(spks(each),1,2);
                        plot(ax,xpair,ypair,'Color', 'k','linewidth',1.3);
                        hold(ax,'on')
%                         pause(0.5)
%                         drawnow
                    end
                    
                end
            
                
            end
            
            %figure configuration
            
            ax.YLabel.String = 'Trials';
            ax.YLim             = [0 max(length(sptimes),10)];
            ax.XLabel.String  	= 'Time (s)';
           % ydata = get(ax,'Ylim');
           % line([0,0],[min(ydata),max(ydata)],'color','r');
            ax.XLim = [0-zpt length(y)/fs-zpt];
             hold(ax,'off')
            
        end
        
        function appsdf(ax, y, fs, sptimes)
            resolution = 0.001;
            gausswidth = 0.02;
            msdf = cal.sdf(sptimes,y,fs,resolution,gausswidth);
            %figure; 
            x = linspace(0,length(msdf),length(msdf))*resolution;
            plot(ax,x,msdf);
            xlim([0,max(x)])
        end
        
        function raster2(sptimes,ylen,fs)
            %find initial of stimuli
            
            
            %idx = find(y,1);
            
            %initial_timestamps = idx;
            % y = y - initial_timestamps;
            
            %zpt = idx/fs;
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
                    plot(xspikes-zpt, yspikes, 'Color', 'k','linewidth',1.3)
                end
                
            end
            
            %figure configuration
            
            
            ax.YLim             = [0 length(sptimes)];
            
            ax.XLabel.String  	= 'Time (s)';
            
            
            
            ydata = get(gca,'Ylim');
            line([0,0],[min(ydata),max(ydata)],'color','r');
            
            xlim([0-zpt ylen/fs-zpt]);
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