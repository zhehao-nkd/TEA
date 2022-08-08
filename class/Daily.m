classdef Daily
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods(Static)
        function d20220705
            dbstop if error
            % update BS_NS_Neurons
            sul = Sultan("E:\BS_NS_Neurons");
            sul.update_Analysis_Files("E:\Updated_BS_NS_Neurons");
            
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

