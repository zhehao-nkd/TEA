% functions to calculate sth.
classdef Cal
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods(Static)
        
        function diffs = allPairDiff(array)
            array = array(:)'; %强制转为行向量
            % Cal all pairwise diffs between members of an array
            part1 = nchoosek(array,2);
            part2 = part1;
            part2(:,[2,1]) = part2(:,[1,2]);
            concat = [part1;part2];
            diffs = concat(:,1)-concat(:,2);
            
        end
        function msdf = sdf(sptimes,y_or_length,fs,resolution,gausswidth)
            
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
            
            if length(y_or_length)==1
                thelength = y_or_length;
            else
                thelength = length(y_or_length);
            end
            time = 0:tstep:thelength/fs;
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
        
        function msdf = sdf_concated_sptimes(concated_sptimes,y,fs,resolution,gausswidth)
            
            sptimes{1} = concated_sptimes;
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
        
        
        function img = spec(y,fs) % calculate the image matrix of spectrogram
            y = 2*y;%myAxe = app.UIAxes;
            [S,F,T] = spectrogram(y,hamming(512),512-round(fs/1e3),512,fs);
            img = flip(log(1+abs(S)));
        end
        
        
        function img = Deprecated_spec1024(y,fs) % calculate the image matrix of spectrogram
            y = 2*y;%myAxe = app.UIAxes;
            [S,F,T] = spectrogram(y,hamming(1024),1024-round(fs/1e3),1024,fs);
            img = flip(log(1+abs(S)));
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
            Draw.spec(y,fs);
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
        
        function N = psth(sptimes,y,fs)
            
            all = [];
            for iTrial = 1:length(sptimes)
                all             = [all; sptimes{iTrial}];               % Concatenate spikes of all trials
            end
            
            len =  length(y)/fs
            edges = linspace(0,len,500) % 500 is the bin number!!!!!! wrong
            
            N = histcounts(all,edges);
        end
        
        function N = psth_SingleTrial(sptimes,ylen,fs) % calculate psth using length(y); Do not consider multiple trials
%             
%             all = [];
%             for iTrial = 1:length(sptimes)
%                 all             = [all; sptimes{iTrial}];               % Concatenate spikes of all trials
%             end
            
            sptimes = sptimes * 1000;
            len =  ylen/fs*1000; % 多少hao秒
            edges = linspace(0,len,len); %len is the bin number
            
            N = histcounts(sptimes,edges);
        end
        
        
        function [num,detail] = sylnum(dir) % calculate in total how many syllables in a song folder
            num = 0;
            files = Extract.filename (dir, '*.wav');
            for idx = 1: length(files)
                [yall,fs] = audioread(files{idx}); 
                y = yall(:,2); % yall means 2-channels y
                frag = segment(y,fs).seg5;
                num = length(frag) + num;
                detail(idx).name = files{idx};
                detail(idx).num = length(frag);
            end
            
        end 
        
        function N =  psth_frag(y,fs,sptimes) % specific for syllables/elements
            dur = 0.01; % 10 ms
            len = length(y)/fs;
            edges = linspace(0,len,ceil(len/dur));
            concat =  vertcat(sptimes{:});
            %figure; histogram(concat,edges);
            if isempty(concat)
                N = zeros(length(edges) - 1,1);
            else
                N = histcounts(concat,edges);
            end
            
        end
        
        
        function Num_Trigger_binaryCode = BinaryCodeTriggerNum(num_of_stimuli)
        end
        
        function firing_rates = eachTrialFiringRate(sptimes,duration) % duration in seconds
            
            for k = 1: length(sptimes)
               firing_rates(k) =  length(sptimes{k})/duration;
            end
            
        end
        
        
        function timeStimuliPresentation(dirpath)
            
            filenames = Extract.filename(dirpath);
        end
    end
    
end

