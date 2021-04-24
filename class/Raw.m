classdef Raw < handle & Audio
    
    properties
        
        path
        folder
        filename
        initial
        terminal
    end
    
    methods
        
        function r = Raw(filepath)
            % path = "Z:\Yazaki-SugiyamaU\Bird-song\Blue318\2016-06-22_20-15-43.WAV"
            [y,fs] = audioread(filepath);
            r@Audio(y,fs);
            r.path = filepath;
            [dir,r.filename,~] = fileparts(filepath);
            splited = split(dir,'\');
            r.folder = splited{end}; % shows in which folder the file locates
            seg = segment(r.y,r.fs).seg3;
            if ~isempty(seg)
                r.initial = [seg.initial]';
                r.terminal = [seg.terminal]';
            end
        end 
        
        
        function syl = avgn(r)
            dbstop if error
            syl = struct;
            
            if isempty(r.initial)
                syl = [];   % to avoid empty segment
                return
            end
            
            for n = 1:length(r.initial)
                syl(n).birdid = convertStringsToChars(convert.bid(r.folder));
                syl(n).filename = convertStringsToChars(r.filename);
                syl(n).number = n;
                syl(n).plx = nan; % these are only available for ephys object
                syl(n).channel = nan; 
                syl(n).unit = nan;
                syl(n).y = r.y(r.initial(n):r.terminal(n));
                syl(n).hpy = highpass(syl(n).y,pa.hp,r.fs); % high passed y, threshold is 400
               %syl(n).stft = cal.stft(syl(n).y,400,48,400,r.fs);
                syl(n).image = cal.img(syl(n).hpy,r.fs); % store the image matrix 
                syl(n).label = 0; % not significant
            end
        end
        
        
    end
        
    
    
    
end