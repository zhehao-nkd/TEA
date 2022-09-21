
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
        
        function b = Batch(path_txt,path_plx,path_folder)
            
            if nargin== 0
                return
            end
            
            b.path.path_txt = path_txt;
            b.path.path_plx = path_plx;
            b.path.path_folder = path_folder;
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
            
            idx = 0;
            for k = 1:size(b.path,1)
                spikes = Spike.split(b.path(k).path_txt); % in this step, the unsorted spikes has been removed
                for m = 1: length(spikes)
                    idx = idx + 1;
                    b.neu{idx} = spikes{m};
                    b.plx{idx} = b.path(k).path_plx;
                    b.wavfolder{idx} = b.path(k).path_folder;
                    if isa(b.plx{idx},'Trigger')
                        [~,plxname,~] = fileparts(b.plx{idx}.inputpath);
                    elseif isa(b.plx{idx},'string')||isa(b.plx{idx},'char')
                        [~,plxname,~] = fileparts(b.plx{idx});
                    end
                    channelname = unique(b.neu{idx}.channelname);
                    channelname = channelname{1};
                    unitname = unique(b.neu{idx}.unit);
                    b.nlist(idx).idx = idx;
                    b.nlist(idx).neuronname = sprintf('%s_%s_%u',plxname,channelname,unitname);
                    
                end
            end
        end
        
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

        function neurons = getn_shift(b,shift_value)
           
            parfor idx = 1: length(b.sneu) % heer should be parfor , i edited here just tio chelc the bug
                NN = Neuron;
                neurons{idx} = NN.shiftNeuron(b.sneu{idx},b.splx{idx},b.swavfolder{idx},shift_value);
            end
            
            parfor w = 1: length(neurons) % for each neuron, write the same_channel_spikes info
                
                
                same_channel_spikes = Spike.extract_specific_channel(b.path.path_txt,neurons{w}.channelname);
                neurons{w}.sameChannelSpikes = same_channel_spikes;
                
            end
             
        end
        
        function neurons = getn(b) % initiatialize neurons
            
            neurons = {};
            parfor idx = 1: length(b.sneu) % Here should be parfor
                neurons{idx} = Neuron(b.sneu{idx},b.splx{idx},b.swavfolder{idx});
                
            end

            
            parfor w = 1: length(neurons) % for each neuron, write the same_channel_spikes info
                
                
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
        
        function simuThree(b)
            
            IMG = {};
            
            neuronlist = b.getn;
            
            for k = 1:length(neuronlist)
                IMG{k} = neuronlist{k}.OneRowThree;
            end
            
            final_IMG = vertcat(IMG{:});
            [~,plxname,~] = fileparts(b.input{2});
            imwrite(final_IMG,sprintf('Three_%s.png',plxname));
        end
        
    end
    
   methods(Static)
       
       function neuronlist = pipline(path_txt,path_plx,path_folder)
           dbstop if error
           addpath(genpath("C:\Users\Zhehao\Dropbox (OIST)\My_Matlab\TEA"))
           b = Batch(path_txt,path_plx,path_folder);
           b.select;
           neuronlist = b.getn;
           
           for k = 1: length(neuronlist)
               thisn = neuronlist{k};
               
              % A = Analysis(thisn);
               %save(thisn.neuronname,'A','-v7.3');
               %thisn.three;
               thisn.pltthree(1);
               %thisn.ResponseBasedOrderedThreePlots;
           end
       end

   end
end



