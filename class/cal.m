classdef cal
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods(Static)
        function msdf = sdf(sptimes,y,fs,resolution,gausswidth)
            
            %idx = find(y,1);
            
            %initial_timestamps = idx;
            % y = y - initial_timestamps;
            
            %zpt = idx/fs; % zp means Zero Point
            if exist('resolution','var')
                tstep = resolution; % Resolution for SDF [s]
            else
                tstep     	= .001;
            end
            
            if  exist('gausswidth','var')
                sigma = gausswidth;
            else
                sigma = .005; % Width of gaussian/window [s]
            end
            
            time = 0:tstep:length(y)/fs;
            %time     	= tstep-.25:tstep:.25;
            % Time vector
            
            for iTrial = 1:length(sptimes)
                spks    = [];
                gauss   = [];
                spks   	= sptimes{iTrial}';          % Get all spikes of respective trial
                if isempty(spks)
                    sdf(iTrial,:)	= zeros(1,length(time));    % Add zero vector if no spikes
                else
                    % For every spike
                    for iSpk = 1:length(spks)
                        
                        % Center gaussian at spike time
                        mu              = spks(iSpk);
                        %tbase = tstep-.25:tstep:.25;
                        %time = tbase + mu;
                        % Calculate gaussian
                        p1              = -.5 * ((time - mu)/sigma) .^ 2;
                        p2              = (sigma * sqrt(2*pi));
                        gauss(iSpk,:)   = exp(p1) ./ p2;
                    end
                    % Sum over all distributions to get spike density function
                    sdf(iTrial,:)       = sum(gauss,1);
                end
            end
            % Single trial display
            %[len1,len2] = size(sdf);
            
            
            % Average response
            
            %plot((1:len2)/1000-zpt,mean(sdf), 'Color', 'k', 'LineWidth', 1);
            % mVal = max(mean(sdf)) + round(max(mean(sdf))*.1);
            
            
            
            msdf= mean(sdf);
        end
        
        
        function threshold = thres(sdf,range)  % calculate the threshold that above this threshold, percentage of samples is smaller than the range
            [N,edges] = histcounts(sdf,100);
            
            countN = 0;
            for  n = 1: length(N)
                reverse = length(N) +1 -n;
                countN = countN  + N(reverse);
                percentage = countN/sum(N);
                if percentage > range
                    out = reverse + 1;
                    break
                end
            end
            
            threshold = edges(out);
            
        end
        
        
        function image = img(y,fs)
            figure('visible','off')
            draw.spectrogram(y,fs);
            f = getframe(gcf);
            [rgb,~] = frame2im(f);
            image = rgb2gray(rgb);
            close(gcf);
        end
        
        
        function [STFT, f, t] = stft(x, win, hop, nfft, fs)
            
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %              Short-Time Fourier Transform            %
        %               with MATLAB Implementation             %
        %                                                      %
        % Author: Ph.D. Eng. Hristo Zhivomirov        12/21/13 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            % Input:
            % x - signal in the time domain
            % win - analysis window function
            % hop - hop size
            % nfft - number of FFT points
            % fs - sampling frequency, Hz
            %
            % Output:
            % STFT - STFT-matrix (only unique points, time
            %        across columns, frequency across rows)
            % f - frequency vector, Hz
            % t - time vector, s
            
            % representation of the signal as column-vector
            x = x(:);
            
            % determination of the signal length
            xlen = length(x);
            
            % determination of the window length
            wlen = length(win);
            
            % stft matrix size estimation and preallocation
            NUP = ceil((1+nfft)/2);     % calculate the number of unique fft points
            L = 1+fix((xlen-wlen)/hop); % calculate the number of signal frames
            STFT = zeros(NUP, L);       % preallocate the stft matrix
            
            % STFT (via time-localized FFT)
            for l = 0:L-1
                % windowing
                xw = x(1+l*hop : wlen+l*hop).*win;
                
                % FFT
                X = fft(xw, nfft);
                
                % update of the stft matrix
                STFT(:, 1+l) = X(1:NUP);
            end
            
            % calculation of the time and frequency vectors
            t = (wlen/2:hop:wlen/2+(L-1)*hop)/fs;
            f = (0:NUP-1)*fs/nfft;
            
        end
        
        
    end
    
end

