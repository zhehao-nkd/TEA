classdef batch < handle
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        input
        path % just a list of all path
        neu % these are the inputs for class Neuron
        plx % 
        wavfolder%
        sneu  % selected
        splx
        swavfolder
        nlist  % to show the list of neuron name
    end
    
    methods
        
        function b = batch(varargin)
            b.input = varargin;
            b.split;
            b.nlist = b.nlist';
        end
        
        
        function b = split(b) % generate input
            if length(b.input) == 1 % the path of the csv folder
                b.path = table2struct(readtable(b.input{1}{1}));
            elseif length(b.input) == 3
                b.path.path_txt = b.input{1};
                b.path.path_plx = b.input{2};
                b.path.path_folder = b.input{3};
            end
            
            idx = 0;
            for k = 1:size(b.path,1)
                spikes = Spike.split(b.path(k).path_txt);
                for m = 1: length(spikes)
                    idx = idx + 1;
                    b.neu{idx} = spikes{m};
                    b.plx{idx} = b.path(k).path_plx;
                    b.wavfolder{idx} = b.path(k).path_folder;
                    [~,plxname,~] = fileparts(b.plx{idx});
                    channelname = unique(b.neu{idx}.channelname);
                    channelname = channelname{1};
                    unitname = unique(b.neu{idx}.unit);
                     b.nlist(idx).idx = idx;
                    b.nlist(idx).neuronname = sprintf('%s_%s_%u',plxname,channelname,unitname);
                    
                end
            end  
        end
        
        
        function avgn(b) % write avgn mat files
            outdir = 'batch_avgn';
            mkdir(outdir);
            
            for idx = 1: length(b.neu)
                syllables = Neuron(b.neu{idx},b.plx{idx},b.wavfolder{idx}).avgn;
                syllables = syllables';
                
                [~,rawid,~] = fileparts(b.plx{idx});
                
                save(sprintf('%s\\%s.mat',outdir,rawid),'syllables');
                
            end
            
            
        end
        
        function select(b,idx)
            if ~exist('idx','var')
                b.sneu = b.neu;
                b.splx = b.splx;
                b.swavfolder = b.wavfolder;
            else
                b.sneu = b.neu(idx);
                b.splx = b.plx(idx);
                b.swavfolder = b.wavfolder(idx);
            end
        end
        
        function neurons = getn(b) % initiatialize neurons

            for idx = 1: length(b.sneu)
                neurons{idx} = Neuron(b.sneu{idx},b.splx{idx},b.swavfolder{idx});
            end
        end
      
        
        
    end
    
   
end



