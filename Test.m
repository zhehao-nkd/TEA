classdef Test < handle & EphysAnalysis
    properties
        value = 0
    end
    methods
        function t = Test
           
        end
        
        function add(t)
            t.value = t.value + 1;
        end
        
    end
    
    
end