classdef utl
    % store some utilis functions that I have not yet decided where to put
    
  
    methods(Static)
        
        function padded = pad(y,targetlen)
            % pad zero to data in both side
            % here length(y) must < targetlen
            y = y(:); % force y to be column vector
            if mod(targetlen - length(y),2) == 0 % if even number
                frontside = (targetlen - length(y))/2;
                backside = (targetlen - length(y))/2; % length in one side
            else
                frontside = (targetlen - length(y)-1)/2;
                backside = targetlen - length(y)-frontside;
            end
            
           padded = [zeros(frontside,1);y;zeros(backside,1)];
           
            
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

