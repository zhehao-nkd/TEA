% [y,fs] = audioread("D:\R644\CON2.wav");
%
% audiowrite('CON2Single.wav',y(:,2), 32000)
%
%
%
%
% [y2,fs2] = audioread("D:\R644\CON2.wav");
% y2 = 0.26 *y2
%
% audiowrite('CON2realsmaller.wav',y2(:,2), 32000)

%global extend

classdef Sound < handle
    
    
    properties
        path
        raw
        name
        fs
        y
        corey
        fragment
        trigger
        initial
    end
    
    methods
        
        function s = Sound(path_wav) % constrcutor
            s.path = path_wav;
            
            [s.raw,s.fs] =  audioread(path_wav);
            
            [~,s.name,~] = fileparts(path_wav);
            
            if size(s.raw, 2)== 1 % if only a single channel
                s.y = s.raw;
            elseif size(s.raw,2)== 2
                s.y = s.raw(:,2);
                s.iniTrigger; % number of trigger pulse
            end
            
            s.segment;
            s.setcorey;
            s.initial = [s.fragment.initial].';
            
        end
        
        function s = segment(s) % function for segmenting the sound
            
            s.fragment = segment(s.y,s.fs).seg1;
            
        end
        
        function s = normalize(s) % function for normalize the sound
        end
        
        function s = preprocess(s) % preprocess the sound to remove the interval silence and normalize the syllables
        end
        
        function drawSegment(s)
            ini = [s.fragment.initial].'
            ter = [s.fragment.terminal].'
            
            % spectrogram
            figure
            
            mySpectrogram(s.y, s.fs, '')
            
            for idx = 1: length(ini)
                hold on
                line([ini(idx)/s.fs,ini(idx)/s.fs],[0,16],'color','b')
                hold on
                line([ter(idx)/s.fs,ter(idx)/s.fs],[0,16],'color','r')
            end
            
            % sonogram
            figure
            
            filteredY = highpass(abs(s.y),400,s.fs);
            plot(filteredY)
            
            for idx = 1: length(ini)
                hold on
                line([ini(idx),ini(idx)],[0,max(filteredY)],'color','b')
                hold on
                line([ter(idx),ter(idx)],[0,max(filteredY)],'color','r')
            end
            
        end
        
        function y = getY(s)
            y = s.y;
        end
        
        function s = iniTrigger(s)
            yT = s.raw(:,1);
            
            
            THRESHOLD = 0.5; %not necessary here
            
            yT (yT < THRESHOLD ) = 0; %Convert singals of pulse channel to 0 and 1
            yT (yT  >= THRESHOLD ) = 1; % based on the threshold
            
            yT = [0; yT];
            
            k = find(yT);                        % Find all the 1
            initials = k(yT(k-1)==0)-1;      % Find the initial of all the pulses
            
            s.trigger = length(initials);
        end
        
        function s = setcorey(s) % y with real sound signal, without silent duration
            
            global extend
            if ~isempty(extend)
                s.corey = s.y(extend*s.fs:end-extend*s.fs);
            else
                s.corey = s.y(s.fragment(1).initial:s.fragment(end).terminal);
            end
            
        end
        
        
        
    end
    
    
    methods(Static)
        function files= split(folder_wav)
            files = extract.filename(folder_wav,'*.wav');
        end
        
        function intercepted = intercept(y,fs,start,stop) % tpoint is the time in seconds
            dbstop if error
            
            if start <= 0
                frontpad = NaN(int32([-(start-1)*fs,1]));  % pad a little bit more
            else
                frontpad = [];
            end
            
            if stop >= length(y)/fs
                bottompad = NaN(int32([(stop+1-(length(y)/fs))*fs,1])); % pad a little bit more
            else
                bottompad = [];
            end
            
            padded = [frontpad;y;bottompad];
            
            lag = length(frontpad)/fs;
            
            intercepted = padded(int32((start+lag)*fs):int32((stop+lag)*fs));
        end
        
        
        
        
    end
end









