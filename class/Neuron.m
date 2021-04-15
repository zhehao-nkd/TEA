classdef Neuron < handle
    %ANALYSIS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        e
    end
    
    methods
        
        function n = Neuron(neuron, path_plx, folder_wav)
            
             
        %             for m = 1: length(neurons) % for each neuron
        %                 sp = Spike(neurons{m});
        %
        %                 temp = struct;
        %                 parfor n = 1: length(songs)
        %                     so = Sound(songs{n});
        %                     temp(n).e = Ephys2(sp,t,so);
        %                 end
        %                 temp2 = [temp.e];
        %             end
        
            songs = Sound.split(folder_wav);
            t = Trigger(path_plx);
            
            
            % for this neuron
            sp = Spike(neuron);
            
            
            parfor k = 1: length(songs)
                so = Sound(songs{k});
                e{k} = Ephys(sp,t,so);
            end
            n.e = e;
        end
          
        
        function sigonly = siginf(n)
            
            collect = cellfun(@(obj) obj.siginf, n.e,'UniformOutput',false);
            sigonly =  vertcat(collect{:});
            
        end
        
        
        function response = sylinf(n)
  
            collect = cellfun(@(obj) obj.sylinf, n.e,'UniformOutput',false);
            response =  horzcat(collect{:});
            
        end
        
        
        function response = resp(n) % mimic old response
            
            collect = cellfun(@(obj) obj.resp, n.e,'UniformOutput',false);
            response =  horzcat(collect{:}).';
            
        end
        
        
        function pre = preinf(n)
            pre = cellfun(@(obj) obj.preinf, n.e,'UniformOutput',false);
        end
        %METHOD1 Summary of this method goes here
        %   Detailed explanation goes here
        
    end
end


