classdef Audio < handle
    
    properties
        y
        fs
    end
    
    methods
        function a = Audio(y,fs)
            a.y = y;
            a.fs = fs;
        end
    end
    
end
