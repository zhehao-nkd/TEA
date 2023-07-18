%%% constitute : syllable/element constitute


classdef constitute < handle
    
    properties
        y
        fs
        normalized
    end
    
    methods
        
        function c = constitute(y,fs)
            c.y = y;
            c.fs = fs;
            c.norm(0.05);
        end
        
        function c = norm(c,tarrms) % generate normalized syllable
            
            rawrms = rms(c.y); rms_ratio = tarrms/rawrms;
            c.normalized = c.y*rms_ratio; 
           % disp('maga');
        end
        
        function mirrored = mirr(c)  % generate mirrored syllable
            mirrored = flip(c.normalized);
        end
       % some further methods to measure the feature of each syllables
    end
    
end