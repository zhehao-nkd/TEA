classdef Sort
    %SORT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods(Static)
        function output = songname(info)
            %SORT Construct an instance of this class
            %   Detailed explanation goes here
            tinfo = struct2table(info);
            output = table2struct(sortrows(tinfo,{'songname','fragid'}));
        end
        
        function output = stiminfo(info)
             output = table2struct(sortrows(struct2table(info),{'songname','fragid'}));
             
        end
        
       
    end
end

