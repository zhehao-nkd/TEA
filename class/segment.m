% a collection of segment functions

% I will have to rewrite segment function to a much broader one which can
% use either eleinf,autoseg, simpelseg algorithms
classdef segment < handle
    
    properties
        y
        fs
    end

    methods
        
        function s = segment(y,fs) % constructor function
            s.y = y;
            s.fs = fs;
        end

        function fragment = seg1(s) % function for segmenting the sound
            
            measurey = abs(s.y);
            measurey = highpass(measurey,400,s.fs);
            
            
            
            [yup,~] = envelope(measurey,800,'rms');
            
            [lowpks,lowlcs] = findpeaks(-yup,'MinPeakHeight',-0.005,'MinPeakDistance',250);
            
            [highpks,highlcs] = findpeaks(yup,'MinPeakHeight',0.012,'MinPeakDistance',500);
            
            %             figure
            %             findpeaks(-yup,'MinPeakHeight',-0.005,'MinPeakDistance',250);
            %             figure
            %             findpeaks(yup,'MinPeakHeight',0.015,'MinPeakDistance',500);
            %
            idx = 0;
            
            for lcsN = 2: length(lowlcs)
                if  ~isempty(highlcs(lowlcs(lcsN-1)<highlcs&highlcs<lowlcs(lcsN)))
                    idx = idx + 1;
                    fragment(idx).initial = lowlcs(lcsN-1);
                    fragment(idx).terminal = lowlcs(lcsN);
                    fragment(idx).y = s.y(lowlcs(lcsN-1):lowlcs(lcsN));
                    
                end
            end
            
            
        end
      
        function fragment = seg2(s)
            
            
            
            
            config = configSegSyl('default');
            
            
            
            
            %% Extract envelope
            
            Hd = fil_hp(config.fil_freq); %high pass filter, this is the first parameter
            
            rms_all = rms(s.y); rms_ratio = 0.05/rms_all;
            
            y_normalized = s.y* rms_ratio;    % Normalize
            
            y_filtered = filter(Hd,y_normalized); %Pass filter with configuration as Hd
            
            y_filtered_abs = abs(y_filtered); %Absolute value
            
            [y_env, ~] = envelope(y_filtered_abs,120,'peak'); % Get Envelope, this is the second parameter
            
            y_env(y_env < 0) = 0;
            
            %% Transform signal samples to 0 and 1
            % Amplitue threshold, this is the third parameter
            y_env_add0 = [0; y_env;0];             %Add a 0 to the front to avoid k=1 or k = length(y)
            y_env_add0(y_env_add0 < config.amp_thre) = 0; %Convert singals  to 0 and 1
            y_env_add0(y_env_add0 >= config.amp_thre) = 1; % based on the threshold
            
            
            k = find(y_env_add0);   % Find all the 1
            
            
            initials = k(y_env_add0(k-1)==0)-1; % Find the initials of all the segements
            ends = k(y_env_add0(k+1)==0)-1;  % Find the ends of all the segements
            
            prior_initials = initials;
            prior_ends = ends;
            
            
            
            %% Remove small signals
            
            for n = 1:length(initials)
                if ends(n)-initials(n)< config.min_syl*s.fs
                    initials(n) = NaN;
                    ends(n)= NaN;
                end
            end
            
            initials (isnan(initials)) = [];
            ends (isnan(ends)) = [];
            
            
            %% Remove small gaps
            
            
            for n = 1:length(ends)-1
                if initials(n+1)-ends(n)< config.min_gap*s.fs
                    ends(n)= NaN;
                    initials(n+1) = NaN;
                end
            end
            
            initials (isnan(initials)) = [];
            ends (isnan(ends)) = [];
            
            %% write to sptimes
            
            fragment.prior_initials = prior_initials;
            fragment.prior_ends = prior_ends;
            fragment.initials = initials;
            fragment.ends = ends;
            fragment.y = s.y;
            fragment.fs = s.fs;
            fragment.y_filtered_abs = y_filtered_abs;
            fragment.y_env = y_env;
            
            
        end
        
        function frag = seg3(s)  % specific for dataset, will be merged into seg1 in the future
            
            %gpuy = gpuArray(y);
            measurey = abs(s.y);
            measurey = highpass(measurey,400,s.fs);
            
            
            [yup,~] = envelope(measurey,800,'rms');
            %     figure
            %     plot(yup)
            %     figure
            [lowpks,lowlcs] = findpeaks(-yup,'MinPeakHeight',-0.005,'MinPeakDistance',500);
            %     figure
            %     findpeaks(yup,'MinPeakHeight',0.01,'MinPeakDistance',500);
            [highpks,highlcs] = findpeaks(yup,'MinPeakHeight',0.015,'MinPeakDistance',500);
            
            frag = struct;
            frag(1).label = 0;
            for lcsN = 2: length(lowlcs)
                
                
                
                if  ~isempty(highlcs(lowlcs(lcsN-1)<highlcs&highlcs<lowlcs(lcsN)))
                    
                    frag(lcsN).label = 1;
                    frag(lcsN).initial = lowlcs(lcsN-1);
                    frag(lcsN).terminal = lowlcs(lcsN);
                    %syllable(lcsN).y = s.y(lowlcs(lcsN-1):lowlcs(lcsN));
                    
                else
                    frag(lcsN).label = 0;
                end
                
                %                 figure('Visible','off');
                %                 mySpectrogram(syllable(idx).y,fs,'Decide');
                %                 temp = getframe(gcf)
                %                 syllable(idx).I = temp.cdata;
                %                 close(gcf);
            end
            
            ind=find([frag.label]);  % remove fake data
            if ~isempty(ind)
                frag=frag(ind);
            else
                frag = [];
            end
            
            % frag = rmfield(frag,'label');
        end
        
        function frag = seg4(s)
            mingap = 0.001;
            
            allidx = find(abs(s.y)>0);
            
            initials = allidx(s.y(allidx-1)==0); % Find the initials of all the segements
            ends = allidx(s.y(allidx+1)==0);  % Find the ends of all the segements
            
            
            % Remove small gaps
            
            
            for n = 1:length(ends)-1
                if initials(n+1)-ends(n)< mingap*s.fs
                    ends(n)= NaN;
                    initials(n+1) = NaN;
                end
            end
            
            initials (isnan(initials)) = [];
            ends (isnan(ends)) = [];
            
            for idx = 1: length(initials)
                frag(idx).initial = initials(idx);
                frag(idx).terminal = ends(idx);
                frag(idx).y = s.y(initials(idx):ends(idx));
            end
            
        end
        
        function frag = seg5(s)
            % suitable for already normalized signal
%              mingap = 0.0000;
%             thres = 0.00056;
%             env =  medfilt1(abs(s.y),700,'truncate');
            
            mingap = 0.0001;
            thres = 0.00056;
            env =  medfilt1(abs(s.y),750,'truncate');
            env(env<=thres)= 0;
            env(env>thres)= 1;
            if env(1) == 1
                env = [0;env]; % if the pre-gap is zero  % dangerous code!!!!!
            end
            allidx = find(env>0);
            
            initials = allidx(env(allidx-1)==0); % Find the initials of all the segements
            ends = allidx(env(allidx+1)==0);  % Find the ends of all the segements
            
            
            % Remove small gaps
            
            
            for n = 1:length(ends)-1
                if initials(n+1)-ends(n)< mingap*s.fs
                    ends(n)= NaN;
                    initials(n+1) = NaN;
                end
            end

            initials (isnan(initials)) = [];
            ends (isnan(ends)) = [];
            
            
            % Remove small segments
            minseg = 0.005; % larger than 5 miliseconds
            
            for n = 1:length(ends)
                if ends(n)-initials(n)< minseg*s.fs
                    ends(n)= NaN;
                    initials(n) = NaN;
                end
            end
            
            initials (isnan(initials)) = [];
            ends (isnan(ends)) = [];
            
            
            % asemble calculated initials and ends back to struct frag
            for idx = 1: length(initials)
                frag(idx).initial = initials(idx);
                frag(idx).terminal = ends(idx);
                frag(idx).y = s.y(initials(idx):ends(idx));
            end
            
        end
        
    end
    
end



