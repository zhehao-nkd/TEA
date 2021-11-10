
% Core of stimuli generating function 
classdef stimcore
    %STIMCORE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
          inf% inf struct or eleinf
          fs
    end
    
    methods
        
        function s = stimcore(inf)
            %STIMCORE Construct an instance of this class
            %   Detailed explanation goes here
            s.inf = inf;
            s.fs = unique([inf.fs].');
           
        end

        
        function r = replace(inf)
        end
        
        function new = normalize2(s,tarrms)
            % every syllable have rms as 0.05
            
            old = s.inf;
            new = old;
            
            for k = 1: length(s.inf)
                rawrms = rms(old(k).y); rms_ratio = tarrms/rawrms;
                new(k).y = old(k).y*rms_ratio; 
            end

        end
        
%         function new = normalize2(s,tarrms)
%              %normalize to make the song have the same rms
%       
%             old = s.inf;
%             new = old;
%             ys = {old(:).y}.';
%             sumy = vertcat(ys{:});
%             
%             sumrms = rms(sumy);
%             
%             for k = 1: length(old)
%                
%                 ratio = tarrms/sumrms;
%                 new(k).y = ratio*old(k).y;
%             end
%  
%         end

   
        
        function summed = assemble(s) % a function to assemble constitute together
            
            
            old = s.inf;
            ys = {old.y}.';
            gaps = [old(:).pregap].';
            gaps(1) = 0; % replace the first pregap to 0,which is a inf
            
            
            summed = [];
            
            for haris = 1: length(ys)
                summed = [summed;zeros(int64(gaps(haris)*s.fs),1);ys{haris}];
            end
            
            
        end  % ?? 这个函数貌似没用了
        
        

 

        
        
    end
    
    methods(Static)
        
       function new = normalize(fraginf,tarrms)
            % normalize to make the song have the same rms
           if ~exist('tarrms','var')
               tarrms = 0.05; % default high-pass filter
           end
           
            old = fraginf;
            new = old;
            ys = {old(:).y}.';
            sumy = vertcat(ys{:});
            
            sumrms = rms(sumy);
            
            for k = 1: length(old)
               
                ratio = tarrms/sumrms;
                new(k).y = ratio*old(k).y;
            end
 
       end
       
       function new = highpass(fraginf,hpf) % hpf is the high pass frequency
           % this normalize all element/syllable to have the same rms
           if ~exist('hpf','var')
               hpf = 420; % default high-pass filter
           end
           
           old = fraginf;
           new = old;
           
           for k = 1: length(old)
               new(k).y = highpass(old(k).y, hpf,old(k).fs);
           end
            
        end
        
    end
        
   
    
    
end

