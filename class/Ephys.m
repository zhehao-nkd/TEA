classdef Ephys < handle
    %EPHYS2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        spike
        trigger
        sound
        sptimes
        fs
        y
        corey
        presptimes
    end
    
    methods
        function e = Ephys(spike, trigger, sound)
            
            
            e.spike = spike;
            e.trigger = trigger;
            e.sound = sound;
            e.y = e.sound.y;
            e.corey = e.sound.corey;
            
            % disp('MAGA!');
            e.fs = sound.fs;
            e.allocate;
            
            
        end
        
        function e = allocate(e) % here the response is a collection of response
            
            n = e.sound.trigger;
            initials = e.trigger.info(n).time;
            ylength = length(e.sound.corey)/e.sound.fs;
            
            for k = 1:length(initials)
                e.sptimes{k} = e.spike.time( initials(k)<e.spike.time& e.spike.time<initials(k) + ylength)...
                    - initials(k);
                e.presptimes{k} = e.spike.time( initials(k)- ylength<e.spike.time & e.spike.time<initials(k))- (initials(k)-ylength);
            end
            
        end
        
        function e = raster(e)% draw raster
            color = 'k';
            title = ' ';
            figure;
            draw.raster(e.sptimes,e.sound.corey,e.fs,color,title);
            %draw.raster(e.sound.corey, e.fs, e.sptimes);
        end
        
        
        function e = spec(e) % draw spectrogram
            
            figure;
            draw.spectrogram(e.sound.corey,e.fs);
            %draw.raster(e.sound.corey, e.fs, e.sptimes);
        end
        
        function e = sdf(e)    % draw sdf
            figure;
            draw.sdf(e.sound.corey,e.fs,e.sptimes);
        end
        
        function e = three(e)% draw three plots
            figure;
            draw.three(e.sound.corey,e.fs,e.sptimes);
        end
        
        
        function locs = peak(e)% time of response peak
            resolution = .001;
            gausswidth = resolution*20;
            
            minpeakprominence = 10;
            minpeakdistance = 90;
            range = 0.05;
            presdf = cal.sdf(e.presptimes,e.corey,e.fs,resolution,gausswidth);
            thres = cal.thres(presdf,range);
            
            sdf = cal.sdf(e.sptimes,e.corey,e.fs,resolution,gausswidth);
            %             figure
            %             findpeaks(sdf,"MinPeakHeight",thres,...
            %                 'MinPeakProminence',minpeakprominence,'MinPeakDistance',minpeakdistance);
            [~,locs] = findpeaks(sdf,"MinPeakHeight",thres,...
                'MinPeakProminence',minpeakprominence,'MinPeakDistance',minpeakdistance);
            locs = locs*0.001 ;% convert seconds to miliseconds
        end
        
        function sylidx = sig(e) % index of significant syllable
            dbstop if error
            initial = ( e.sound.initial -e.sound.initial(1) )/e.fs;  % this is for cvonverting initial to seconds relative to the first syllable
            locs = e.peak;
            if ~isempty(locs)
                for idx = 1: length(locs)
                    sylidx = length(find(~(locs(idx)<initial))); % return the idx of the significant syllables
                end
            else
                sylidx = [];
            end
            
            sylidx = unique(sylidx,'stable'); % remove repeated element
            
        end
        
        function siginf = siginf(e) % full info of sig syllable
            
            siginf = sylinf([e.sylinf.label] == 1); % derive from sylinf
            if isempty(e.sig)
                siginf = [];
            end
        end
        
        
        function syl = sylinf(e) % full info of all syllables
            dbstop if error
            number = e.sig;
            syl = struct;
            for n = 1:length(e.sound.initial)
                
                
                syl(n).sound = convertStringsToChars(e.sound.name);
                syl(n).number = n;
                syl(n).plx = convertStringsToChars(e.trigger.plxname);
                syl(n).channel = e.spike.channel;
                syl(n).unit = e.spike.unit;
                syl(n).y = e.sound.fragment(n).y;
                syl(n).image = cal.img(syl(n).y,e.fs); % store the image matrix
                temp = extract.feature(syl(n).y,e.fs);
                syl(n).goodness = temp.goodness;
                syl(n).meanf = temp.mean_frequency;
                syl(n).fm = temp.fm;
                syl(n).amplitude = temp.amplitude;
                syl(n).entropy = temp.entropy;
                syl(n).pitch = temp.pitch;
                syl(n).rawpitch = temp.rawpitch;
                syl(n).am = temp.AM;
                syl(n).initial = e.sound.fragment(n).initial;
                syl(n).terminal = e.sound.fragment(n).terminal;
                if n ~= 1
                    syl(n).pregap = (length(syl(n).initial) - length(syl(n-1).terminal))/e.fs;
                else
                    syl(n).pregap = inf; % pre-gap duration
                end
                syl(n).dur = length(syl(n).y)/e.fs;
                
                if ismember(n,number)
                    syl(n).label = 1; % significant
                else
                    syl(n).label = 0; % not significant
                end
                
                
            end
            
        end
        
        function syl = avgn(e) % full info of all syllables
            dbstop if error
            number = e.sig;
            syl = struct;
            for n = 1:length(e.sound.initial)
                syl(n).birdid = convert.bid(convertStringsToChars(e.sound.name));
                syl(n).filename = convertStringsToChars(e.sound.name);
                syl(n).number = n;
                syl(n).plx = convertStringsToChars(e.trigger.plxname);
                syl(n).channel = e.spike.channel;
                syl(n).unit = e.spike.unit;
                syl(n).y = e.sound.fragment(n).y;
                syl(n).hpy = highpass(syl(n).y,400,e.fs); % high passed y, threshold is 400
                syl(n).image = cal.img(syl(n).y,e.fs); % store the image matrix

                if ismember(n,number)
                    syl(n).label = 1; % significant
                else
                    syl(n).label = 0; % not significant
                end
            end

            
        end
        
        function response = resp(e) % mimic the old response function
            
            response.plxname = e.trigger.plxname;
            response.channel_name = e.spike.channel;
            response.unit_name = e.spike.unit;
            response.stimuli_name = convertStringsToChars(e.sound.name);
            response.wav_Y = e.sound.y;
            response.initial = e.sound.initial;
            response.terminal = [e.sound.fragment.terminal].';
            response.sigsylidx = e.sig;
            
        end
        
        function pre = preinf(e, duration)% info of pre-peak duration
            
            if ~exist('duration','var')
                dur = 0.200 ;% second , this is default
            else
                dur = duration;
            end
            
            locs =  e.peak;
            
            if ~isempty(locs)
                
                for n = 1: length(locs)
                    pre(n).y = Sound.intercept(e.corey,e.fs,locs(n)-dur,locs(n));
                    pre(n).fs = e.fs;
                    temp = extract.feature(pre(n).y,e.fs);
                    pre(n).goodness = temp.goodness;
                    pre(n).meanf = temp.mean_frequency;
                    pre(n).fm = temp.fm;
                    pre(n).amplitude = temp.amplitude;
                    pre(n).entropy = temp.entropy;
                    pre(n).pitch = temp.pitch;
                    pre(n).am = temp.AM;
                    pre(n).plx = e.trigger.plxname;
                    pre(n).channel = e.spike.channel{1};
                    pre(n).unit = e.spike.unit;
                    pre(n).sound = e.sound.name;
                    
                end
            else
                pre = [];
            end
            
        end
    end
end








