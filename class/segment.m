% a collection of segment functions


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
        
        function fragment = seg3(s)  % segment already preprocessed song
            
        end
        
        
    end
    
end
    
