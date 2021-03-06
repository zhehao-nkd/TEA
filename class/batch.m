% batch for data analysis
classdef Batch < handle
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
        
        function b = Batch(varargin)
            b.input = varargin;
            b.split;
            b.nlist = b.nlist';
        end
        
        function spikeinf = manspike(b)
            for ii = 1: length(b.neu)
                b.select(ii);
                temp = b.getn;
                thisn = temp{1};
                thisn.manspike;
            end
        end
        
        function b = split(b) % split a recording file to neurons with different channel and unit names
            
            if length(b.input) == 1 % the path of the csv folder
                b.path = table2struct(readtable(b.input{1}{1}));
            elseif length(b.input) == 3
                b.path.path_txt = b.input{1};
                b.path.path_plx = b.input{2};
                b.path.path_folder = b.input{3};
            end
            
            idx = 0;
            for k = 1:size(b.path,1)
                spikes = Spike.split(b.path(k).path_txt); % in this step, the unsorted spikes has been removed
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
        
%         function b = add_same_channel_spikes(b) % For each neuron object from the bacth, write same_channel_spikes info into it
%             
%         end
        
        function avgn(b) % write avgn mat files
            dbstop if error
            outdir = 'batch_avgn';
            mkdir(outdir);
            % 48 is jumped out
            for idx = 49: length(b.neu) %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                temp = Neuron(b.neu{idx},b.plx{idx},b.wavfolder{idx});
                syllables = temp.avgn;
                %syllables = syllables';
                
                %[~,rawid,~] = fileparts(b.plx{idx});
                fullid = temp.neuronname; %%%%%%%%%%%%%%%%%%
                
                save(sprintf('%s\\%s.mat',outdir,fullid),'syllables');
                
            end
            
            
        end
        
        function select(b,idx)
            if ~exist('idx','var')
                b.sneu = b.neu;
                b.splx = b.plx;
                b.swavfolder = b.wavfolder;
            else
                b.sneu = b.neu(idx);
                b.splx = b.plx(idx);
                b.swavfolder = b.wavfolder(idx);
            end
        end
        
        function neurons = getn(b,mergeIdx) % initiatialize neurons
            
            if exist('mergeIdx','var')
                for idx = 1: length(b.sneu)
                    neurons{idx} = Neuron(b.sneu{idx},b.splx{idx},b.swavfolder{idx},mergeIdx);
                end
            else
                for idx = 1: length(b.sneu)
                    neurons{idx} = Neuron(b.sneu{idx},b.splx{idx},b.swavfolder{idx});
                end
            end
            
            for w = 1: length(neurons) % for each neuron, write the same_channel_spikes info
                
                
               same_channel_spikes = Spike.extract_specific_channel(b.path.path_txt,neurons{w}.channelname);
               neurons{w}.sameChannelSpikes = same_channel_spikes;
                
            end
            
        end
      
        function featuretsne(b)
            for idx = 11: length(b.neu) %%%%%%%% modified here
                Neuron(b.neu{idx},b.plx{idx},b.wavfolder{idx}).featuretsne;
            end
        end
        
        function three(b)
            
            
            for idx = 1: length(b.neu)
                b.select(idx);
                neurons = b.getn;
                neurons{1}.three;
                disp('HAhahahahahaha');
            end
        end
        
        function sapscatter(b)
            
            for idx = 1: length(b.neu)
                
                b.select(idx);
                neuron = b.getn;
                neuron = neuron{1};
                try
                    neuron.sapscatter;
                catch Error
                end
            end
        end
        
        
    end
    
   
end



