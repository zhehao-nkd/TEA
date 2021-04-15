classdef batch < handle
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        n
    end
    
    methods
        
        function b = batch(varargin)
            dbstop if error
            
            info = batch.input(varargin);
            
            for k = 1: length(info)
                n{k} = Neuron(info(k).spike,info(k).path_plx,info(k).path_folder);
            end
            b.n = n;
        end
        
    end
    
    methods(Static)
        
        function info = input(cells)
            if length(cells) == 1 % the path of the csv folder
                path = table2struct(readtable(cells{1}{1}));
            elseif length(cells) == 3
                path.path_txt = cells{1};
                path.path_plx = cells{2};
                path.path_folder = cells{3};
            end
            
            info = struct;
            idx = 0;
            for k = 1:size(path,1)
                spikes = Spike.split(path(k).path_txt);
                for m = 1: length(spikes)
                    idx = idx + 1;
                    info(idx).spike = spikes{m};
                    info(idx).path_plx = path(k).path_plx;
                    info(idx).path_folder = path(k).path_folder;
                end
            end
            

            
        end
    end
end



