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
        wavinfo
        path
        raw
        name
        fs
        y
        corey
        fragment
        trigger
        initial
        terminal
        feature
        sapfragment
        sapinitial
        sapterminal
        sapexcel % an xls-fomratted excel table containing sap-segmented syllable/element
    end
    
    methods
        
        function s = Sound(path_wav) % constrcutor
            
            s.setexcel("C:\Users\Zhehao\Dropbox (OIST)\My_Luscinia\sap1276.xls");
            s.path = path_wav;
            
            
            
            [s.raw,s.fs] =  audioread(path_wav);
            
            [~,s.name,~] = fileparts(path_wav);
            
            if size(s.raw, 2)== 1 % if only a single channel
                s.y = s.raw;
            elseif size(s.raw,2)== 2
                s.y = s.raw(:,2);
                s.iniTrigger; % number of trigger pulse
            end
%            s.inifeature; % initialize feature   trmporaliy delete this
%            part

            wavinfo = audioinfo(path_wav);
            
            if isempty(wavinfo.Comment) % if segmentation information is not stored in the  wav file
                s.segment;
                s.setcorey;
                s.initial = [s.fragment.initial].';
                s.terminal = [s.fragment.terminal].';
            end
            % unfortunately I have not giot time to write about the another
            % possibility
            
          %  s.sapsegment;
          %  s.sapinitial = [s.sapfragment.initial].';
          %  s.sapterminal = [s.sapfragment.terminal].';
            disp('Sound一回！');
        end
        
        function s = setexcel(s, path)
            s.sapexcel = path;
        end
        function s = inifeature(s)
            temp = SAT_hijack(s.y,s.fs);
            s.feature = temp.features;
        end
        
        
        function thissap = getsap(s)
            dbstop if error
            if isempty(s.sapexcel)  % if there is no sapexcel
            [folder,b,c] = fileparts(s.path);
            xlsfile = extract.filename(folder,'*.xls');
            xlsfile = xlsfile{1};
            else
            [~,b,c] = fileparts(s.path);
            temp = strsplit(b,'_');   % this modification is dangerous  !!!!!!
            b = temp{1};
            xlsfile = s.sapexcel;
            end
            
            sap = table2struct(readtable(xlsfile));
            thisfile = sprintf('%s',b);
          
            thisidx = find(~cellfun(@isempty,regexp({sap(:).fileName},thisfile)));
            thissap = sap(thisidx);
            
        end
        
        function s = sapsegment(s)
            dbstop if error
            thissap = s.getsap;
            
            for idx = 1: length(thissap)
                s.sapfragment(idx).initial = thissap(idx).syllable_start/1000*s.fs;
                s.sapfragment(idx).terminal = (thissap(idx).syllable_start+ thissap(idx).syllable_duration )/1000*s.fs;
                s.sapfragment(idx).y = s.y(s.sapfragment(idx).initial:s.sapfragment(idx).terminal);
            end
            
        end  % segment function specific for sap
        
        function s = segment(s) % function for segmenting the sound
            
            s.fragment = segment(s.y,s.fs).seg5;
            warning('Now the segment method is seg5!!');
        end
        
        function rmsy = normalize(s) % function for normalize the sound
            filtery = highpass(s.y,400,s.fs);
            rmsy = filtery*rms(filtery)/0.5;  % when 0.5 is the goal
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
            
            difftype = length(unique(round(diff(initials),-1)));
            if difftype == 1||difftype == 0
            s.trigger = length(initials);
            elseif difftype == 2 % When the stimuli is using binary trigger code
                onevalue = max(unique(round(diff(initials),-1)));
                zerovalue = min(unique(round(diff(initials),-1)));
                
                temp = round(diff(initials),-1);
                temp(temp == onevalue ) = 1;
                temp(temp == zerovalue ) = 0;
                s.trigger = bi2de(temp(:).');
            end
                
                
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









