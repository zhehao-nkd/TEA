classdef Analysis < handle
    %ANALYZE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        neurons
        lists
        list_sum
        group
       
        source    % read xlxs into struct
        list
       
        fragnames % a cell of unique frag names
        degnames
        replanames
        
        normlist % list of norm-songs
        all_eleinf
        conspe_eleinf
        
        targets % unique(targetnames)
        
        neuinfo
        uniqueid
        birdid
        unique_neuronname
        
        % ids of neuros corresponding to different stimuli
        song_id_redundant
        song_id
        song_only_id % the neuron which has only song stimuli
        frag_id
        deg_id
        repla_id
        plx_data_fs
        
    end
    
    methods
        
        
        function a = Analysis(neurons)
            
           if exist('neurons','var')
               input = neurons;
               
               if isa(input,'Neuron')
                   a.neurons{1} = input;
               else
                   a.neurons = input;
               end
               
               a.setNeuronInfo;
               a.setStimuliCorrespondingNeuronId;
               
               % if norm stimuli exist in multiple Neuron files, then
               % delect the first one which has only norm stimuli, as it is
               % far away from others in time
               to_remove_id = intersect(a.song_only_id,setdiff(a.song_id_redundant,a.song_id));
               
               neuron_count = 0;
               for k = 1: length(a.neurons)
                   if ~ismember(k,to_remove_id)
                       neuron_count = neuron_count + 1;
                       lists{k} = a.neurons{neuron_count}.todisplay;
                   end
               end
               a.list = horzcat(lists{:}); % The list are merging individual lists of all member Neurons
               % sort
               a.sort  % !!!!!!!!!!!!
               
               temp = regexp(a.neurons{1}.plxname,'[RBOYRG]\d{3}','match');
               
               % set birdid uniqueid and unique_neuronname
               a.birdid = temp{1};
               a.uniqueid = a.neurons{1}.uniqueid;
               a.unique_neuronname = sprintf('%s_%u',a.birdid,a.uniqueid);
           end   
            
        end

        function e_objects = getAllEphysObject(a)
           collects = {};
            for k = 1: length(a.neurons)  
                collects{k} = a.neurons{k}.e(:);
            end
            e_objects = vertcat(collects{:});
        end
        
         function IMG = three(a) % draw three plot
            
             es = getAllEphysObject(a);
            for idx = 1: length(es)
                es{idx}.three_with_sig_resp; % newer version of threeplot drawing method
                frame = getframe(gcf);
                I{idx} = frame.cdata;
                close(gcf);
            end
            
            a.draw_waveform;     % draw waveform
            frame = getframe(gcf);
            I{length(I)+ 1} = frame.cdata;
            close(gcf);
            
            % draw blank white
            lieshu = 12;
            hangshu = ceil(length(I)/lieshu);
            rest = lieshu*hangshu - length(I);
            white = uint8(255*ones(size(I{1})));
            
            if rest > 0
                for k = 1:rest
                    I = [I,white];
%                     ax = gcf;
%                     ax.Position(3) = 560;
%                     ax.Position(4) = 420;
                end
            end
            
            reshapedI = reshape(I, lieshu,[])';     
            clear I     
            IMG = cell2mat(reshapedI);
            imwrite(IMG,sprintf('Three_%s.png',a.unique_neuronname));
            
         end
        
        function a =setStimuliCorrespondingNeuronId(a)
            
            a.song_id_redundant = [];
            a.frag_id = [];
            a.deg_id = [];
            a.repla_id = [];
            a.song_only_id = [];
            for k = 1: length(a.neuinfo)
                if ismember("song",a.neuinfo(k).keywords)
                    a.song_id_redundant = [a.song_id_redundant,k];
                end
                
                if ismember("frag",a.neuinfo(k).keywords)
                    a.frag_id = [a.frag_id,k];
                end
                
                 if ismember("deg",a.neuinfo(k).keywords)
                    a.deg_id = [a.deg_id,k];
                 end
                
                 if ismember("repla",a.neuinfo(k).keywords)
                    a.repla_id = [a.repla_id,k];
                 end
                 if strcmp("song",a.neuinfo(k).keywords)
                     a.song_only_id = [a.song_only_id,k];
                 end
                 
                
            end
            
            if length(a.song_id_redundant) == 1
                a.song_id = a.song_id_redundant;
            elseif length(a.song_id_redundant)~=0
                [~,a.song_id] = max(cellfun(@length,{a.neuinfo(a.song_id_redundant).keywords}));
                % find the id which have max length of keywords, if the
                % result is multipl ids, then select the first one
                % But maybe the last one will be much proper!
            else
                a.song_id = 1; % Very bad meaningless code
            end
            
        end
        
        function unique_neuronname = neuronname(a)
            unique_neuronname = sprintf('%s_%u',a.birdid,a.uniqueid);
        end
        
        function drawWaveform(a)
            
            % temporialriy a.neurons{1}
            figure
            a.neurons{1}.draw_waveform;
            
        end
        
        function meanWL = calMeanWaveLength(a)
            
            wavecollect = {};
            for s = 1: length(a.neurons)
                wavecollect{s} = a.neurons{s}.waveform;
            end
            waveforms = vertcat(wavecollect{:});
            %waveforms =  n.waveform;
            [~,troughstime] = min(waveforms,[],2);
            wavlen_units = [];
            
            for k = 1: size(waveforms,1) % length is dangerous!!!!!
                this_wf = waveforms(k,:);
               [~,wavlen_units(k)] =  max(this_wf (troughstime(k):end));
            end
            
            
            a.plx_data_fs = 30000; %hard code !!!!!! Dangerous
            meanWL =  mean(wavlen_units*(1/a.plx_data_fs)*1000); % ms
            
            
        end
        
        function a = setNeuronInfo(a)
            %  to generate info of neurons of an Analysis object, used for
            %  Augustus
            
            a.neuinfo = struct;
            
            for k = 1: length(a.neurons)
                a.neuinfo(k).uniqueid = a.neurons{k}.uniqueid;
                a.neuinfo(k).neuronname = a.neurons{k}.neuronname;
                a.neuinfo(k).keywords = [];
                
                if length(find(~cellfun(@isempty,regexp([a.neurons{k}.slist.name].','norm|Norm|song|Song')))) > 16
                    a.neuinfo(k).keywords = [a.neuinfo(k).keywords,"song"];
                end
                
                if ~isempty(find(~cellfun(@isempty,regexp([a.neurons{k}.slist.name].','syl|Syl|Ele|ele|frag|Frag'))))
                    a.neuinfo(k).keywords = [a.neuinfo(k).keywords,"frag"];
                end
                
                if ~isempty(find(~cellfun(@isempty,regexp([a.neurons{k}.slist.name].','deg|Deg'))))
                    a.neuinfo(k).keywords = [a.neuinfo(k).keywords,"deg"];
                end
                
                
                if ~isempty(find(~cellfun(@isempty,regexp([a.neurons{k}.slist.name].','repla|Repla|catego|Catego'))))
                    a.neuinfo(k).keywords = [a.neuinfo(k).keywords,"repla"];
                end
                
            end
        end
        
        function a = set_eleinf(a,eleinf) %set eleinf,eleinf is the info-data used to construct all the stimuli
            if isa(eleinf,'struct')
                a.all_eleinf = eleinf;
            elseif isa(eleinf,'string')|| isa(eleinf,'char')
                loaded = load(eleinf);
                a.all_eleinf = loaded.all_eleinf;
            end
            
             conspe_ids = find( ~cellfun(@isempty, regexp([a.all_eleinf.songname].','CON|SPE')) );
             a.conspe_eleinf = a.all_eleinf(conspe_ids);
        end
        
        function a = update(a)  % update stimuli-response list after adding extra plx-txt-folder files
            for k = 1: length(a.neurons)
                lists{k} = a.neurons{k}.todisplay;
            end
            a.list = vertcat(lists{:}); 
            a.sort;
        end
        
        function a = add_neuron(new_neuron) % add neuron directly
            a.neurons = {a.neurons,new_neuron};
            a.update;
        end
        
        function a = add_from_path(txt,plx,folder,channel,unit) % add neuron from path
        
        end
        
        function a = sort(a) % sort is the function to split stimuli-response pairs with different types (frags, replas,norms)
            % frag
            fragidx = find(~cellfun(@isempty, regexp([a.list(:).stimuliname].','Frag')));
            fraglist = a.list(fragidx);
            for k = 1: length(fraglist)
                
                temp = strsplit(fraglist(k).stimuliname,'-');
                
                fragnames{k} = sprintf('%s-%s-%s',temp{1},temp{2},temp{3});
            end
            
            if exist('fragnames','var')
                a.fragnames = unique(fragnames);
            end
            
            
            % deg
            degidx = find(~cellfun(@isempty, regexp([a.list.stimuliname].','deg')));
            deglist = a.list(degidx);
            for k = 1: length(deglist)
                
                temp = strsplit(deglist(k).stimuliname,'-');
                
                degnames{k} = sprintf('%s-%s-%s',temp{1},temp{2},temp{3});
            end
            if exist('degnames','var')
                a.degnames = unique(degnames);
            end
           

            % repla 
            replaidx = find(~cellfun(@isempty, regexp([a.list.stimuliname].','Repla')));
            replalist = a.list(replaidx);
            for k = 1: length(replalist)
                
                temp = strsplit(replalist(k).stimuliname,'-');
                
                replanames{k} = sprintf('%s-%s-%s',temp{1},temp{2},temp{3});
            end
            if exist('replanames','var')
                a.replanames = unique(replanames);
            end
            
             % norm
            normidx = find(~cellfun(@isempty, regexp([a.list.stimuliname].','norm')));
            a.normlist = a.list(normidx);
            
            % target
            targetidx = find(~cellfun(@isempty, regexp([a.list.stimuliname].','Repla')));
            targetlist = a.list(targetidx);
            for k = 1: length(targetlist)
                
                temp = strsplit(targetlist(k).stimuliname,'-');
                
                targetnames{k} = sprintf('%s-%s-%s',temp{6},temp{7},temp{8});
            end
            
            if exist('targetnames','var')
                a.targets = unique(targetnames);
            end
            
        end
        
        function threePlotsWithPitch(a)
            
            e_objects = getAllEphysObject(a);
            
            for idx = 1: length(e_objects)
                e_objects{idx}.threeWithFreqRelatedFeatures; % newer version of threeplot drawing method
                frame = getframe(gcf);
                I{idx} = frame.cdata;
                close(gcf);
            end
            
            figure;
            a.neurons{a.song_id}.draw_waveform;     % draw waveform
            % a.draw_waveform; % this is the original draw function
            frame = getframe(gcf);
            I{length(I)+ 1} = frame.cdata;
            close(gcf);
            
            % draw blank white
            lieshu = 10;
            hangshu = ceil(length(I)/lieshu);
            rest = lieshu*hangshu - length(I);
            white = uint8(255*ones(size(I{1})));
            
            if rest > 0
                for k = 1:rest
                    I = [I,white];
                    %                     ax = gcf;
                    %                     ax.Position(3) = 560;
                    %                     ax.Position(4) = 420;
                end
            end
            
            reshapedI = reshape(I, lieshu,[])';
            clear I
            IMG = cell2mat(reshapedI);
            imwrite(IMG,sprintf('Three_%s.png',a.unique_neuronname));
         
        end
        
        function drawselect(a,ids,range)
            d = Display(a.list);
            d.showselect(ids,range);
        end
        
        function drawfrag(a,keyword,rangeratio,ids)
            
            mergedeleinf = a.all_eleinf;
            d = Display(a.list);
            
            if exist('keyword','var')
                
                if exist('rangeratio','var')
                    if exist('ids','var')
                        d.showfrag(keyword,mergedeleinf,rangeratio,ids);
                    else
                        d.showfrag(keyword,mergedeleinf,rangeratio);
                    end
                    
                else
                    d.showfrag(keyword,mergedeleinf);
                end  
            else
                for k = 1: length(a.fragnames)
                    d.showfrag(a.fragnames{k},mergedeleinf);
                end
            end
            
        end
        
        function drawdeg(a)
            %load("C:\Users\Zhehao\Dropbox (OIST)\My_Stimuli\Y661@06282021\mergedeleinf.mat");
            d = Display(a.list);
            for k = 1: length(a.degnames)
                d.showdeg(a.degnames{k});
            end
        end
        
        function drawrepla(a)
       
            d = Display(a.list);
            
            for k = 1: length(a.replanames)
                d.showrepla(a.replanames{k});
            end
            
        end
        
        function drawpartrepla(a,initial ,terminal) % show part of all
       
            d = Display(a.list);
            repla = d.findrepla;
            
            repla = repla(initial: terminal);
            
            d = Display(repla);
            d.showrepla;
            
        end  
        
        function drawallrepla(a)
            d = Display(a.list);
            d.showrepla;
            
        end
            
        function drawsinglenorm(a,keyword)
            d = Display(a.list);
            d.shownorm(keyword);
            
        end
        
        function drawsimplefrag(a,keyword)
            
            d = Display(a.list);
            
            if exist('keyword','var')
                d.showsimplefrag(a.eleinf,keyword);
            else
                d.showsimplefrag(a.eleinf);
            end
        end
        
        function drawtransform(a,keyword)
        end
        
        function fraglist = judgeFragResponse(a) %%% To judge whether the neuron response to a frag or not
            
           ids = find(~cellfun(@isempty, regexp({a.list.stimuliname}.','Frag|syl'))); % find all frags
           
           % ’syl'可以兼容旧的stimuli命名规则
           
           fraglist = a.list(ids);
           
           for n = 1: length(fraglist)
            tempsum = cal.psth_frag(fraglist(n).plty,fraglist(n).fs,fraglist(n).pltsptimes);
            halfsum = sum(tempsum(end/2:end));
            fullsum = sum(tempsum);
            maxvalue = max(cal.psth_frag(fraglist(n).plty,fraglist(n).fs,fraglist(n).pltsptimes));
            fraglist(n).maxvalue = maxvalue;
            fraglist(n).halfsum = halfsum;
            fraglist(n).fullsum = fullsum;
            if maxvalue > 6 % here the threshold is very important % originally set as 8
                fraglist(n).label = 1;
            else
                fraglist(n).label = 0;
            end    
           end
           
        end
        
        function Conlist = evaluateConResponse(a) % Eveluate Conspecific song response
            
            ids = find(~cellfun(@isempty, regexp({a.list.stimuliname}.','norm'))); % find all frags
            
            % ’syl'可以兼容旧的stimuli命名规则
            
            Conlist = a.list(ids);
            
            for n = 1: length(Conlist)
                sdf = cal.sdf(Conlist(n).pltsptimes,Conlist(n).plty,Conlist(n).fs,0.001,0.004);
                [maxsdf,maxidx] = max(sdf);
                
                percentage_max = maxidx/length(sdf);
                time_max = length(Conlist(n).plty)/Conlist(n).fs*percentage_max;
                % check whether the surroding are has spikes in most of the
                % trials
                extracted_sptimes = extract.sptimes(Conlist(n).pltsptimes,time_max - 0.15, time_max + 0.15); % 前后 100ms
                
                num_of_not_empty_trials = length(find(~cellfun(@isempty, extracted_sptimes)));
                
                mean_maxsdf = maxsdf/length(Conlist(n).pltsptimes);
                tempsum = cal.psth_frag(Conlist(n).plty,Conlist(n).fs,Conlist(n).pltsptimes);
                halfsum = sum(tempsum(end/2:end));
                fullsum = sum(tempsum);
                maxvalue = max(cal.sdf(Conlist(n).pltsptimes,Conlist(n).plty,Conlist(n).fs)); % I changed the value here from cal.psth_frag to cal.sdf
                Conlist(n).sdf = cal.sdf(Conlist(n).pltsptimes,Conlist(n).plty,Conlist(n).fs);
                Conlist(n).maxvalue = maxvalue;
                Conlist(n).mean_maxsdf = mean_maxsdf;
                Conlist(n).halfsum = halfsum;
                Conlist(n).fullsum = fullsum;
                Conlist(n).maxsdf = maxsdf;
                %if maxvalue > 20 % here the threshold is very important originally 24
                if (mean_maxsdf > 2.22) 
                    Conlist(n).label = 1;
                else
                    Conlist(n).label = 0;
                end
                
                if fullsum > 40 % if not 60% of the trails are not empty
                    Conlist(n).label = 1;
                end
                
                
                if (num_of_not_empty_trials/length(Conlist(n).pltsptimes)<0.5)||(num_of_not_empty_trials<5) % if not 60% of the trails are not empty
                    Conlist(n).label = 0;
                end
                
                if fullsum > 48 % if not 60% of the trails are not empty
                    Conlist(n).label = 1;
                end
                
                % A rescue
                slim_extracted_sptimes = extract.sptimes(Conlist(n).pltsptimes,time_max - 0.8, time_max + 0.8);
                slim_trials = length(find(~cellfun(@isempty,slim_extracted_sptimes)));
                if slim_trials/length(Conlist(n).pltsptimes)> 0.5 % if not 60% of the trails are not empty
                    Conlist(n).label = 1;
                end
                
                
                if (maxsdf) > 19.5&& strcmp(Conlist(n).stimuliname,'norm-Y515A-21Pulses') && strcmp(Conlist(n).plxname,'Y661_Z17')
                    disp('Incredible bug in pltsptimes of function Analysis.evaluateConResponse !!');
                    Conlist(n).label = 1;
                end
                    
            end
           
        end
        
        function fraglist = to2ndAcousticSpace(a) %%%
            ids = find(~cellfun(@isempty, regexp({a.list.stimuliname}.','repla'))); % find all frags
           
           fraglist = a.list(ids);
           
           for n = 1: length(fraglist)
            tempsum = cal.psth_frag(fraglist(n).rawy,fraglist(n).fs,fraglist(n).rawsptimes);
            range = 1 % the very fisrt 1 second
            beginmax = max(tempsum(1: ceil(length(tempsum)*range/(length(fraglist(n).rawy)/fraglist(n).fs)) ));% the maximum value of the begining 0.5 second
            halfsum = sum(tempsum(end/2:end));
            fullsum = sum(tempsum);
            maxvalue = max(cal.psth_frag(fraglist(n).rawy,fraglist(n).fs,fraglist(n).rawsptimes));
            fraglist(n).maxvalue = maxvalue;
            fraglist(n).halfsum = halfsum;
            fraglist(n).fullsum = fullsum;
             fraglist(n).beginmax = beginmax;
            if maxvalue > 8 % here the threshold is very important
                fraglist(n).label = 1;
            else
                fraglist(n).label = 0;
            end    
           end
           
        end
        
        function draw_frag_scatter(a,not_tested_handle)  % not_tested_handle = 1 means draw
            
            if ~exist('not_tested_handle','var') % to judge whether to draw not tested elements or not
                not_tested_handle = 0;
            end
            
            global_eleinf = a.all_eleinf;
            fraglist = a.judgeFragResponse;
            
            
            for k = 1: length(global_eleinf)
                uniqueid = sprintf('-%s-%u-',global_eleinf(k).songname, global_eleinf(k).fragid);
                
               if ~isempty (find(~cellfun(@isempty, regexp({fraglist.stimuliname}.',uniqueid)), 1)) % if this element was tested
                     id_in_fraglist = find(~cellfun(@isempty, regexp({fraglist.stimuliname}.',uniqueid)));
                   if fraglist(id_in_fraglist).label == 0 % if the tested ele does not trigger reponse
                       
                   global_eleinf(k).scatter = 0;
                   elseif fraglist(id_in_fraglist).label == 1 % if the tested ele trigger response
                       global_eleinf(k).scatter = 1;
                   end
               else
                   global_eleinf(k).scatter = -1;
               end
            end
            
            % section for drawing
            figure
            hold on
            for k = 1: length(global_eleinf)
                if global_eleinf(k).scatter == -1
                    if not_tested_handle == 1
                        scatter(global_eleinf(k).coor_1,global_eleinf(k).coor_2,[],'k','filled'); % black for not tested
                    elseif not_tested_handle == 0 % if not draw test handle, do nothing
                    end
                elseif global_eleinf(k).scatter == 0
                    scatter(global_eleinf(k).coor_1,global_eleinf(k).coor_2,[],'g','filled');   % green for tested but not response-eliciting
                   % text(double(global_eleinf(k).coor_1),double(global_eleinf(k).coor_2),sprintf('%s-%u',global_eleinf(k).songname, global_eleinf(k).fragid) )
                elseif global_eleinf(k).scatter == 1
                    scatter(global_eleinf(k).coor_1,global_eleinf(k).coor_2,[],'r','filled');
                    % text(double(global_eleinf(k).coor_1),double(global_eleinf(k).coor_2),sprintf('%s-%u',global_eleinf(k).songname, global_eleinf(k).fragid) )
                end
                
                
            end
            
            
            % label the targets
            global_names = [global_eleinf.songname].';
            global_fragids = [global_eleinf.fragid].';
            
            global_merged = {};
           for w = 1: length(global_names)
               global_merged{w} = sprintf('%s-%u',global_names(w),global_fragids(w));
           end
           global_merged = global_merged.';
            
            
            for u = 1: length(a.targets)
                [~,beta] = ismember( a.targets{u}, global_merged );
                scatter( global_eleinf(beta).coor_1,global_eleinf(beta).coor_2,[],'k','h');
            end
            hold off
            
            
        end
        
        function sort_frags_by_response_strength_and_then_draw(a)  % good explanation % Use maxvalue in judgeResponse
            
            fraglist =  judgeFragResponse(a);
            if isempty(fraglist)
                return
            end
            sorted_fraglist = table2struct(sortrows( struct2table(fraglist) ,'maxvalue','descend'));
            
          
            I = {}; % collection of frag-response-three images
            for k = 1: length(sorted_fraglist)
                h =  figure('Position',[681 403 523 696],'Color','w');
                %h.WindowState = 'maximized';
                draw.two(sorted_fraglist(k).plty,sorted_fraglist(k).fs,sorted_fraglist(k).pltsptimes);
                xlabel(sorted_fraglist(k).stimuliname);
                temp = getframe(gcf);
                I{k} = temp.cdata;
                
                
                close(h)
            end
            
%             figure;
%             n.draw_waveform;     % draw waveform
%             frame = getframe(gcf);
%             I{length(I)+ 1} = frame.cdata;
%             close(gcf);
            
            % draw blank white
            lieshu = 10;
            hangshu = ceil(length(I)/lieshu);
            rest = lieshu*hangshu - length(I);
            white = uint8(255*ones(size(I{1})));
            
            if rest > 0
                for k = 1:rest
                    I = [I,white];
                    %                     ax = gcf;
                    %                     ax.Position(3) = 560;
                    %                     ax.Position(4) = 420;
                end
            end
            
            reshapedI = reshape(I, lieshu,[])';
            clear I
            IMG = cell2mat(reshapedI);
            imwrite(IMG,sprintf('NeuralReponseToFragments_%s.png',a.unique_neuronname));
            
        end
        
        function draw_waveform(a) % this function is used to draw waveform
            % what this works for concatenating all neurons together???
            for k = 1: length(a.neurons)
                waveforms{k} = a.neurons{k}.waveform;
            end
            
            concat_waveforms = vertcat(waveforms{:});
            
            figure('Color','w');
            hold on
            plot(concat_waveforms.',':','Color',[.5,.5,.5]);
            plot(max(concat_waveforms),'--','Color','blue');
            plot(min(concat_waveforms),'--','Color','blue');
            plot(mean(concat_waveforms),'Color','red');
            
           
            
        end
        
        function draw_waveform_each(a) % this function is used to draw waveform
            % what this works for concatenating all neurons together???
            for k = 1: length(a.neurons)
                
                this_waveform =  a.neurons{k}.waveform;
                
                figure;
                hold on
                plot(this_waveform.','Color',[.5,.5,.5]);
                plot(max(this_waveform),'--','Color','blue');
                plot(min(this_waveform),'--','Color','blue');
                plot(mean(this_waveform),'Color','red');
            
            end
            
            
            
           
            
        end
        
        function reordered_ids = reorder_three_plots(a)
            
            fragids = find(~cellfun(@isempty, regexp({a.list.stimuliname}.','Frag|syl|Syl'))); % find all frags,兼容 syl|Syl
            if ~isempty( fragids)
                fraglist = a.list(fragids);
                for n = 1: length(fraglist)
                    tempsum = cal.psth_frag(fraglist(n).rawy,fraglist(n).fs,fraglist(n).rawsptimes);
                    halfsum = sum(tempsum(end/2:end));
                    fullsum = sum(tempsum);
                    maxvalue = max(cal.psth_frag(fraglist(n).rawy,fraglist(n).fs,fraglist(n).rawsptimes));
                    fraglist(n).maxvalue = maxvalue;
                    fraglist(n).halfsum = halfsum;
                    fraglist(n).fullsum = fullsum;
                    if maxvalue > 6 % here the threshold is very important % originally set as 8
                        fraglist(n).label = 1;
                    else
                        fraglist(n).label = 0;
                    end
                end
                [~,fragindex] = sortrows( struct2table(fraglist) ,'maxvalue','descend');
                reordered_fragids = fragids(fragindex);
            else
                reordered_fragids = [];
            end
            
            replaids = find(~cellfun(@isempty, regexp({a.list.stimuliname}.','Repla|catego'))); % find all frags,兼容catego
            if ~isempty( replaids)
                replalist = a.list(replaids);
                for n = 1: length(replalist)
                    tempsum = cal.psth_frag(replalist(n).rawy,replalist(n).fs,replalist(n).rawsptimes);
                    halfsum = sum(tempsum(end/2:end));
                    fullsum = sum(tempsum);
                    maxvalue = max(cal.psth_frag(replalist(n).rawy,replalist(n).fs,replalist(n).rawsptimes));
                    replalist(n).maxvalue = maxvalue;
                    replalist(n).halfsum = halfsum;
                    replalist(n).fullsum = fullsum;
                    if maxvalue > 6 % here the threshold is very important % originally set as 8
                        replalist(n).label = 1;
                    else
                        replalist(n).label = 0;
                    end
                end
                [~,replaindex] = sortrows( struct2table(replalist) ,'maxvalue','descend');
                reordered_replaids = replaids(replaindex);
            else
                reordered_replaids = [];
            end
            
            normids = find(~cellfun(@isempty, regexp({a.list.stimuliname}.','norm'))); % find all frags
            if ~isempty( normids )
                normlist = a.list( normids);
                for n = 1: length( normlist)
                    tempsum = cal.psth_frag( normlist(n).rawy, normlist(n).fs, normlist(n).rawsptimes);
                    halfsum = sum(tempsum(end/2:end));
                    fullsum = sum(tempsum);
                    maxvalue = max(cal.psth_frag( normlist(n).rawy, normlist(n).fs, normlist(n).rawsptimes));
                    normlist(n).maxvalue = maxvalue;
                    normlist(n).halfsum = halfsum;
                    normlist(n).fullsum = fullsum;
                    if maxvalue > 6 % here the threshold is very important % originally set as 8
                        normlist(n).label = 1;
                    else
                        normlist(n).label = 0;
                    end
                end
                
                [~, normindex] = sortrows( struct2table( normlist) ,'maxvalue','descend');
                reordered_normids =  normids( normindex);
            else
                reordered_normids = [];
            end
            
            
            
            otherids = find(cellfun(@isempty, regexp({a.list.stimuliname}.','Frag|Repla|norm|catego|Syl|syl'))); % find all frags
            if ~isempty(otherids)
                otherlist = a.list(otherids);
                for n = 1: length(otherlist)
                    tempsum = cal.psth_frag( otherlist(n).rawy,otherlist(n).fs,otherlist(n).rawsptimes);
                    halfsum = sum(tempsum(end/2:end));
                    fullsum = sum(tempsum);
                    maxvalue = max(cal.psth_frag(otherlist(n).rawy,otherlist(n).fs,otherlist(n).rawsptimes));
                    otherlist(n).maxvalue = maxvalue;
                    otherlist(n).halfsum = halfsum;
                    otherlist(n).fullsum = fullsum;
                    if maxvalue > 6 % here the threshold is very important % originally set as 8
                        otherlist(n).label = 1;
                    else
                        otherlist(n).label = 0;
                    end
                end
                
                
                [~, otherindex] = sortrows( struct2table(otherlist) ,'maxvalue','descend');
                reordered_otherids =  otherids( otherindex);
                
            else
                reordered_otherids = [];
            end
            
            reordered_ids = [reordered_fragids;reordered_replaids;reordered_otherids;reordered_normids];  
                 
        end
        
        function drawSSIMSimlarityMatrix(a)
            % This function firstly order the frags by response strength
            % then measure the pairwise similarity between elements
            fraglist = a.judgeFragResponse;
            sorted_fraglist = table2struct(sortrows( struct2table(fraglist) ,'halfsum','descend'));  %此处用 maxvalue或许不太对
            
            for k = 1: length(sorted_fraglist)
                fiy = bandpass(sorted_fraglist(k).y,[900 6000],sorted_fraglist(k).fs); % preiously 5000
                
                envy = rescale(smooth(abs(fiy),150)); % amplitude envelope of y
                %powery = downsample(fiy.^2/length(fiy),fs/1000);
                %downy = downsample(abs(fiy),fs/1000);
                I = cal.spec(fiy,sorted_fraglist(k).fs); % I is the image of the whole song
                sorted_fraglist(k).normalized_img = imresize(I,[257,50]); 
            end
            
            for m = 1:length(sorted_fraglist)
                parfor p = 1:length(sorted_fraglist)
                    sim(m,p) = ssim(sorted_fraglist(m).normalized_img,sorted_fraglist(p).normalized_img);
                end
            end
            
            figure;
            imagesc(sim);
            
            saveas(gcf,sprintf('SSIMSimlarityMatrix-%s.png',a.unique_neuronname));
            close(gcf);
        end
        
        function alignFragsWithSongThenDraw(a)
            dbstop if error
            
            normlist = a.normlist;
            
            ids = find(~cellfun(@isempty, regexp({a.list.stimuliname}.','Frag'))); % find all frags 
            fraglist = a.list(ids);
            
            
            % generate addlagy and addlagsptimes for plotting
            for k = 1:length(normlist)
                
                splited = split( normlist(k).stimuliname,'-');
                songname = splited{3}; % name of the song
                this_fragids = find( ~cellfun(@isempty, regexp({fraglist.stimuliname}.',songname)));
                lensong = length(normlist(k).plty);
                
                normlist(k).fragids = this_fragids;
                
                for f = 1: length(this_fragids)
                    timelag = finddelay(highpass(fraglist(this_fragids(f)).y,450,32000),normlist(k).plty); % Dangerous!!!!!!
                    disp('Dangerous code exist!!!')
                    fraglist(this_fragids(f)).addlagy = [zeros(timelag,1);fraglist(this_fragids(f)).y; zeros((lensong- timelag - length(fraglist(this_fragids(f)).y)),1)]; 
                    fraglist(this_fragids(f)).addlagsptimes = cellfun(@(x) x+timelag/32000-fraglist(this_fragids(f)).pltext, fraglist(this_fragids(f)).pltsptimes,'un',0); % every sptimes will be added a timelag time
                    fraglist(this_fragids(f)).timelag = timelag;
                end

            end
            
            
            % draw Two plots
            
            I = {};
            for k = 1: length(normlist)
                
                figure('Position',[282 759 1444 272],'Color','w');
                draw.two(normlist(k).plty,normlist(k).fs,normlist(k).pltsptimes);
                frame = getframe(gcf);
                jihe{1} = frame.cdata;
                close(gcf);
                
                for v = 1:length(normlist(k).fragids)
                    localids = normlist(k).fragids(v);
                    figure('Position',[282 759 1444 272],'Color','w');
                    draw.two(fraglist(localids ).addlagy,fraglist(localids ).fs,fraglist(localids ).addlagsptimes);
                    frame = getframe(gcf);
                    jihe{1 + v} = frame.cdata;
                    close(gcf);
                end
                
                
                I{k} = vertcat(jihe{:});
                
            end
            
             for w = 1: length(I)
                 imwrite(I{w},sprintf('Three-%u.png',w));
             end
            
%             n.draw_waveform;     % draw waveform
%             frame = getframe(gcf);
%             I{length(I)+ 1} = frame.cdata;
%             close(gcf);
%             
            % draw blank white
%             lieshu = 9;
%             hangshu = ceil(length(I)/lieshu);
%             rest = lieshu*hangshu - length(I);
%             white = uint8(255*ones(size(I{1})));
%             
%             if rest > 0
%                 for k = 1:rest
%                     I = [I,white];
%                     %                     ax = gcf;
%                     %                     ax.Position(3) = 560;
%                     %                     ax.Position(4) = 420;
%                 end
%             end
%             
%             reshapedI = reshape(I, lieshu,[])';
%             clear I
%             IMG = cell2mat(reshapedI);
%             imwrite(IMG,sprintf('Three.png'));
%             
            
            
            
        end
        
        function drawMeanFeatureVsResp(a) % draw the distribution of mean features
            
%             fragids = find(~cellfun(@isempty, regexp({a.list.stimuliname}.','frag|Frag|syl|Syl')));
%             
%             if isempty(fragids)
%                 return
%             end
%             fraglist = a.list(fragids);

            fraglist = a.judgeFragResponse;
            if isempty(fraglist)
                return;
            end
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).maxvalue;
            end
            fraglist = table2struct(sortrows(struct2table(fraglist), 'responseMeasure'));
            
           figure;
           subplot(5,2,1)
           toshow = struct;
           for u = 1: length(fraglist)
               toshow(u).value = fraglist(u).meanfeatures.amplitude;
               toshow(u).resp = fraglist(u).responseMeasure;
           end
           toshow  = table2struct(sortrows(struct2table(toshow), 'value'));  
           plot([toshow.value].',[toshow.resp].');    
           xlabel('Amplitude')
           ylabel('Maximum FR');
           
           subplot(5,2,2)
           toshow = struct;
           for u = 1: length(fraglist)
               toshow(u).value = fraglist(u).meanfeatures.pitch;
               toshow(u).resp = fraglist(u).responseMeasure;
           end
           toshow  = table2struct(sortrows(struct2table(toshow), 'value'));  
           plot([toshow.value].',[toshow.resp].');    
           xlabel('Pitch')
           ylabel('Maximum FR');
           
           subplot(5,2,3)
           toshow = struct;
           for u = 1: length(fraglist)
               toshow(u).value = fraglist(u).meanfeatures.mean_frequency;
               toshow(u).resp = fraglist(u).responseMeasure;
           end
           toshow  = table2struct(sortrows(struct2table(toshow), 'value'));  
           plot([toshow.value].',[toshow.resp].');    
           xlabel('Mean Frequency')
           ylabel('Maximum FR');
           
           subplot(5,2,4)
           toshow = struct;
           for u = 1: length(fraglist)
               toshow(u).value = fraglist(u).meanfeatures.FM;
               toshow(u).resp = fraglist(u).responseMeasure;
           end
           toshow  = table2struct(sortrows(struct2table(toshow), 'value'));  
           plot([toshow.value].',[toshow.resp].');    
           xlabel('FM')
           ylabel('Maximum FR');
           
           subplot(5,2,5)
           toshow = struct;
           for u = 1: length(fraglist)
               toshow(u).value = fraglist(u).meanfeatures.AM;
               toshow(u).resp = fraglist(u).responseMeasure;
           end
           toshow  = table2struct(sortrows(struct2table(toshow), 'value'));  
           plot([toshow.value].',[toshow.resp].');    
           xlabel('AM')
           ylabel('Maximum FR');
           
           subplot(5,2,6)
           toshow = struct;
           for u = 1: length(fraglist)
               toshow(u).value = fraglist(u).meanfeatures.goodness;
               toshow(u).resp = fraglist(u).responseMeasure;
           end
           toshow  = table2struct(sortrows(struct2table(toshow), 'value'));  
           plot([toshow.value].',[toshow.resp].');    
           xlabel('Goodness')
           ylabel('Maximum FR');
           
           subplot(5,2,7)
           toshow = struct;
           for u = 1: length(fraglist)
               toshow(u).value = fraglist(u).meanfeatures.entropy;
               toshow(u).resp = fraglist(u).responseMeasure;
           end
           toshow  = table2struct(sortrows(struct2table(toshow), 'value'));  
           plot([toshow.value].',[toshow.resp].');    
           xlabel('Entropy')
           ylabel('Maximum FR');
           
           subplot(5,2,8)
           toshow = struct;
           for u = 1: length(fraglist)
               toshow(u).value = fraglist(u).meanfeatures.peak_frequency;
               toshow(u).resp = fraglist(u).responseMeasure;
           end
           toshow  = table2struct(sortrows(struct2table(toshow), 'value'));  
           plot([toshow.value].',[toshow.resp].');    
           xlabel('Peak Frequency')
           ylabel('Maximum FR');
           
           subplot(5,2,9)
           toshow = struct;
           for u = 1: length(fraglist)
               toshow(u).value = fraglist(u).meanfeatures.continuity_t;
               toshow(u).resp = fraglist(u).responseMeasure;
           end
           toshow  = table2struct(sortrows(struct2table(toshow), 'value'));  
           plot([toshow.value].',[toshow.resp].');    
           xlabel('Time continuity')
           ylabel('Maximum FR');
           
           subplot(5,2,10)
           toshow = struct;
           for u = 1: length(fraglist)
               toshow(u).value = fraglist(u).meanfeatures.continuity_f;
               toshow(u).resp = fraglist(u).responseMeasure;
           end
           toshow  = table2struct(sortrows(struct2table(toshow), 'value'));  
           plot([toshow.value].',[toshow.resp].');    
           xlabel('frequency continuity')
           ylabel('Maximum FR');
           
        
           saveas(gcf, sprintf('SeparatedFeatureVsResponse_%s.png',a.unique_neuronname));
           close(gcf);
           
           
        end
        
        function drawMeanFeaturesVsRespAsLineChart(a) % draw the distribution of mean features
           
            fraglist = a.judgeFragResponse;
            if isempty(fraglist)
                return;
            end
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).maxvalue;
            end
            
            toshow = struct;
           for u = 1: length(fraglist)
               toshow(u).amplitude = fraglist(u).meanfeatures.amplitude;
               toshow(u).pitch = fraglist(u).meanfeatures.pitch;
               toshow(u).AM = fraglist(u).meanfeatures.AM;
               toshow(u).mean_frequency = fraglist(u).meanfeatures.mean_frequency;
               toshow(u).FM = fraglist(u).meanfeatures.FM;
               toshow(u).peak_frequency = fraglist(u).meanfeatures.peak_frequency;
               toshow(u).entropy = fraglist(u).meanfeatures.entropy;
               toshow(u).resp = fraglist(u).responseMeasure;
           end
           
           toshow  = table2struct(sortrows(struct2table(toshow), 'resp','ascend'));  
           
           
           num_toshow = [toshow.amplitude; toshow.pitch; toshow.AM; toshow.mean_frequency; toshow.FM; ...
               toshow.peak_frequency; toshow.entropy].';
           
           znum_toshow = zscore(num_toshow,0,1);
%             
%            figure
%            plot(znum_toshow.');
           
           figure('Position',[1997 233 1388 658],'Color','w');
           hold on
           c = 1- rescale([toshow.resp].',0.1,1);
           for r = 1: size(znum_toshow,1)
               plot(znum_toshow(r,:),'Color',repmat(c(r),3,1));
               %drawnow
              % pause(0.5)
           end
           colormap(flip(repmat(unique(c),1,3)))
           colorbar
           xlim([0,8]);
           xticks([0 1 2 3 4 5 6 7 8])
           xticklabels({'','amplitude','pitch','AM','mean_frequency','FM','peak-frequency','entropy',''});
           title(sprintf('Totally %u song elements',length(fraglist)));
           ylabel('Zscored Feature(averaged)');
           
           saveas(gcf,sprintf('LineChartMeanFeaturesVsResp-%s.png',a.unique_neuronname));
           close(gcf);         
        end
        
         function BinaryThresholdDrawMeanFeaturesVsRespAsLineChart(a) % draw the distribution of mean features
           
            fraglist = a.judgeFragResponse;
            if isempty(fraglist)
                return;
            end
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).maxvalue;
            end
            
            toshow = struct;
           for u = 1: length(fraglist)
               toshow(u).amplitude = fraglist(u).meanfeatures.amplitude;
               toshow(u).pitch = fraglist(u).meanfeatures.pitch;
               toshow(u).AM = fraglist(u).meanfeatures.AM;
               toshow(u).mean_frequency = fraglist(u).meanfeatures.mean_frequency;
               toshow(u).FM = fraglist(u).meanfeatures.FM;
               toshow(u).peak_frequency = fraglist(u).meanfeatures.peak_frequency;
               toshow(u).entropy = fraglist(u).meanfeatures.entropy;
               toshow(u).resp = fraglist(u).responseMeasure;
           end
           
           toshow  = table2struct(sortrows(struct2table(toshow), 'resp','ascend'));  
           
           
           num_toshow = [toshow.amplitude; toshow.pitch; toshow.AM; toshow.mean_frequency; toshow.FM; ...
               toshow.peak_frequency; toshow.entropy].';
           
           znum_toshow = zscore(num_toshow,0,1);
%             
%            figure
%            plot(znum_toshow.');
           
           figure('Position',[1997 233 1388 658],'Color','w');
           hold on
          % c = 1- rescale([toshow.resp].',0.1,1);
          
          
          % very artifical bad code
          switch a.unique_neuronname
              case 'O686_34'
                  thres = 2;
              case 'R677_55'
                  thres = 6;
              case 'Y661_5'
                  thres = 6;
              case 'Y661_8'
                  thres = 2;
              case 'Y675_22'
                  thres = 3;
              case 'R677_59'
                  thres = 3;
              otherwise
                  thres = 5;
          end
          
          kemal = 0;
          for r = 1: size(znum_toshow,1)
              if toshow(r).resp > thres
                  kemal = kemal + 1;
                  defaultLW = get(gca,'LineWidth');
                  plot(znum_toshow(r,:),'Color',[1 0 0],'LineWidth',defaultLW*1.8);
              else
                  plot(znum_toshow(r,:),'Color',[0.6 0.6 0.6]);
              end
          end
        
           xlim([0,8]);
           xticks([0 1 2 3 4 5 6 7 8])
           xticklabels({'','amplitude','pitch','AM','mean_frequency','FM','peak-frequency','entropy',''});
           title(sprintf('%u of  %u song elements',kemal,length(fraglist)));
           ylabel('Zscored Feature(averaged)');
           
           saveas(gcf,sprintf('New_BinaThres_LineChartMeanFeaturesVsResp-%s.png',a.unique_neuronname));
           close(gcf);         
        end
        
        function drawPairwiseFragmentsMeanFeaturesDistribution(a)
            
            fraglist = a.judgeFragResponse;
            if isempty(fraglist)
                return
            end 
            if isempty(fraglist(1).meanfeatures) %%% Too bad !!
                return
            end
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).maxvalue;
            end
            
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).amplitude = fraglist(u).meanfeatures.amplitude;
                toshow(u).pitch = fraglist(u).meanfeatures.pitch;
                toshow(u).AM = fraglist(u).meanfeatures.AM;
                toshow(u).mean_frequency = fraglist(u).meanfeatures.mean_frequency;
                toshow(u).FM = fraglist(u).meanfeatures.FM;
                toshow(u).peak_frequency = fraglist(u).meanfeatures.peak_frequency;
                toshow(u).entropy = fraglist(u).meanfeatures.entropy;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
           
            name = {'amplitude','pitch','AM','mean_frequency','FM','peak_frequency','entropy' };
            ncomb = nchoosek(name,2);

            I = {};
            parfor idx = 1:size(ncomb,1)
                
              %subplot(hangshu,lieshu,idx);
              figure('Color','w');
              hold on
              cmap = 1 - rescale([toshow.resp].');   
              for  md = 1: length(toshow)  
              scatter(toshow(md).(ncomb{idx,1})', toshow(md).(ncomb{idx,2}),[],repmat(cmap(md,:),1,3),'filled');
              end
              xlabel(replace(ncomb{idx,1},'_','-'));
              ylabel(replace(ncomb{idx,2},'_','-'));
              
              frame = getframe(gcf);
              I{idx} = frame.cdata;
              close(gcf);
                
            end
            
            % PCA plot
            a.drawPCABasedOnAllFeatures;
            frame = getframe(gcf);
            I{length(I)+ 1} = frame.cdata;
            close(gcf);
            
            
            lieshu = 7;
            hangshu = ceil(length(I)/lieshu);
            
            rest = lieshu*hangshu - length(I);
            white = uint8(255*ones(size(I{1})));
            
            if rest > 0
                for k = 1:rest
                    I = [I,white];
                end
            end
            
            reshapedI = reshape(I, lieshu,[])';     
            clear I     
            IMG = cell2mat(reshapedI);
            imwrite(IMG,sprintf('PairwiseFragmentsMeanFeaturesDistribution_%s.png',a.unique_neuronname));       

        end
        
        function drawPCABasedOnAllFeatures(a)
            
            fraglist = a.judgeFragResponse;
            if isempty(fraglist)
                return;
            end
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).maxvalue;
            end
            
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).amplitude = fraglist(u).meanfeatures.amplitude;
                toshow(u).pitch = fraglist(u).meanfeatures.pitch;
                toshow(u).AM = fraglist(u).meanfeatures.AM;
                toshow(u).mean_frequency = fraglist(u).meanfeatures.mean_frequency;
                toshow(u).FM = fraglist(u).meanfeatures.FM;
                toshow(u).peak_frequency = fraglist(u).meanfeatures.peak_frequency;
                toshow(u).entropy = fraglist(u).meanfeatures.entropy;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            
            toshow  = table2struct(sortrows(struct2table(toshow), 'resp','ascend'));
            
            
            num_toshow = [toshow.amplitude; toshow.pitch; toshow.AM; toshow.mean_frequency; toshow.FM; ...
                toshow.peak_frequency; toshow.entropy].';
            
            [coe,score,~] = pca(num_toshow);
            
           cmap = 1- rescale( [toshow.resp].');
            % 1- to plot the data
            figure('Color','white');
            hold on
            scatter(score(:,1),score(:,2),[],repmat(cmap,1,3),'filled');
            %colorbar;
            xlabel('PCA Dim-1');
            ylabel('PCA Dim-2');
            %title('PCA data');
            
        end   
        
        function drawDTWSimilarityMatrix(a)
            
            fraglist = a.judgeFragResponse;
             if isempty(fraglist)
                return;
            end
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).maxvalue;
            end
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).amplitude = fraglist(u).features.amplitude;
                toshow(u).pitch = fraglist(u).features.pitch;
                toshow(u).AM = fraglist(u).features.AM;
                toshow(u).mean_frequency = fraglist(u).features.mean_frequency;
                toshow(u).FM = fraglist(u).features.FM;
                toshow(u).peak_frequency = fraglist(u).features.peak_frequency;
                toshow(u).entropy = fraglist(u).features.entropy;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            toshow  = table2struct(sortrows(struct2table(toshow), 'resp','ascend'));
            
            
            names = {'amplitude','pitch','AM','FM','mean_frequency','peak_frequency','entropy'};
            
            I = {}; % inatialize the Image collection
            f = waitbar(0,'Start');

            for na = 1: length(names)
                waitbar(na/length(names),f,replace(sprintf('Generating %s similarity matrix ...',names{na}),'_','-'));
                figure('Color','w');
                sim = [];
                for s = 1: length(toshow)
                    parfor b = 1: length(toshow)
                        sim(s,b) = 1 - dtw(toshow(s).(names{na}),toshow(b).(names{na}));
                    end
                end
                imagesc(sim);
                xlabel(sprintf('%s-Similarity',names{na})); 
                frame = getframe(gcf);
                I{na} = frame.cdata;
                close(gcf);
            end
            close(f);
          
            
            lieshu = 4;
            hangshu = ceil(length(I)/lieshu);
            rest = lieshu*hangshu - length(I);
            white = uint8(255*ones(size(I{1})));
            
            if rest > 0
                for k = 1:rest
                    I = [I,white];
                end
            end
            
            reshapedI = reshape(I, lieshu,[])';
            clear I
            IMG = cell2mat(reshapedI);
            imwrite(IMG,sprintf('DTWSimilarityMatrix_%s.png',a.unique_neuronname));
                
        end
            
        function drawDTWSimilarityMatrixBasedOnZscoredData(a)
            
            fraglist = a.judgeFragResponse;
            if isempty(fraglist)
                return;
            end
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).maxvalue;
            end
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).amplitude = zscore(fraglist(u).features.amplitude);
                toshow(u).pitch = zscore(fraglist(u).features.pitch);
                toshow(u).AM = zscore(fraglist(u).features.AM);
                toshow(u).mean_frequency = zscore(fraglist(u).features.mean_frequency);
                toshow(u).FM = zscore(fraglist(u).features.FM);
                toshow(u).peak_frequency = zscore(fraglist(u).features.peak_frequency);
                toshow(u).entropy = zscore(fraglist(u).features.entropy);
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            toshow  = table2struct(sortrows(struct2table(toshow), 'resp','ascend'));
            
            
            names = {'amplitude','pitch','AM','FM','mean_frequency','peak_frequency','entropy'};
            
            I = {}; % inatialize the Image collection
            f = waitbar(0,'Start');

            for na = 1: length(names)
                waitbar(na/length(names),f,replace(sprintf('Generating %s similarity matrix ...',names{na}),'_','-'));
                figure('Color','w');
                sim = [];
                for s = 1: length(toshow)
                    parfor b = 1: length(toshow)
                        sim(s,b) = 1 - dtw(toshow(s).(names{na}),toshow(b).(names{na}));
                    end
                end
                imagesc(sim);
                xlabel(sprintf('Zscored-DTWSimilarityMatrix-%s',names{na})); 
                frame = getframe(gcf);
                I{na} = frame.cdata;
                close(gcf);
            end
            close(f);
          
            
            lieshu = 4;
            hangshu = ceil(length(I)/lieshu);
            rest = lieshu*hangshu - length(I);
            white = uint8(255*ones(size(I{1})));
            
            if rest > 0
                for k = 1:rest
                    I = [I,white];
                end
            end
            
            reshapedI = reshape(I, lieshu,[])';
            clear I
            IMG = cell2mat(reshapedI);
            imwrite(IMG,sprintf('ZscoredFeatureSimilarity_%s.png',a.unique_neuronname));
                
        end

        function drawCoeffOfFeaturesLinearFit(a)
            
                
            fraglist = a.judgeFragResponse;
            if isempty(fraglist)
                return
            end
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).maxvalue;
            end
            toshow = struct;
            for u = 1: length(fraglist)
                % for polyfits
                commonx = 1: length(fraglist(u).features.amplitude);
                toshow(u).amplitude_coef = polyfit(commonx, fraglist(u).features.amplitude,1);
                toshow(u).pitch_coef = polyfit(commonx, fraglist(u).features.pitch,1);
                toshow(u).AM_coef = polyfit(commonx, fraglist(u).features.AM,1);
                toshow(u).mean_frequency_coef = polyfit(commonx, fraglist(u).features.mean_frequency,1);
                toshow(u).FM_coef = polyfit(commonx, fraglist(u).features.FM,1);
                toshow(u).peak_frequency_coef = polyfit(commonx, fraglist(u).features.peak_frequency,1);
                toshow(u).entropy_coef = polyfit(commonx, fraglist(u).features.entropy,1);
                toshow(u).resp = fraglist(u).responseMeasure;
  
            end
            toshow  = table2struct(sortrows(struct2table(toshow), 'resp','ascend'));
            
            
            names = {'amplitude_coef','pitch_coef','AM_coef','FM_coef','mean_frequency_coef','peak_frequency_coef','entropy_coef'};
            I = {}; % inatialize the Image collection
            f = waitbar(0,'Start');
            
            for na = 1: length(names)
                waitbar(na/length(names),f,replace(sprintf('Creating %s scatter plot ...',names{na}),'_','-'));
                figure('Color','w');
                hold on
                feature_coefs = {toshow.(names{na})}.';
                feature_coefs = vertcat(feature_coefs{:});
                resp = [toshow.resp].';
                
                scatter(feature_coefs(:,1),feature_coefs(:,2),[],resp,'filled');
                xline(0,'--r');
                
                xlabel('1-th order coefficient','Interpreter','none');
                ylabel('0-th order coefficient','Interpreter','none');
                
                title(replace(sprintf('%s-1 and 0 coefficient',names{na}),'-','_'));
                frame = getframe(gcf);
                I{na} = frame.cdata;
                close(gcf);
            end
            close(f);
          
            
            lieshu = 4;
            hangshu = ceil(length(I)/lieshu);
            rest = lieshu*hangshu - length(I);
            white = uint8(255*ones(size(I{1})));
            
            if rest > 0
                for k = 1:rest
                    I = [I,white];
                end
            end
            
            reshapedI = reshape(I, lieshu,[])';
            clear I
            IMG = cell2mat(reshapedI);
            imwrite(IMG,sprintf('N-order-coefficient_%s.png',a.unique_neuronname));
            
            
            
        end
        
        function drawResponseToRepla(a)
            
            % About the classification of song elements, should I classify
            % elements by 2nd order modulation???
            repla_ids = find( ~cellfun(@isempty, regexp({a.list.stimuliname}.','catego|repla|Repla')));
            
            repla_list = a.list(repla_ids);
            % paused here
            regexp('catego-1-B346A-1-before-Y616X-16-0.0187-61Pulses','(?<=before-)[OYBGR]\d{3}','match')
        end
        
        function split_fraginf = calResponseToWithinSongFragsFromEleinf(a)
            
            latency = 50*0.001; % 50 ms
            
            split_fraginf = struct;
            
            counts = 0;
            if isempty(a.normlist)
                 split_fraginf = struct([]);
                 return
            end
            
            for k = 1: length(a.normlist)
               %songname = regexp(a.normlist(k).stimuliname,'(CON|SPE|norm)-([BRGOY]\d{3}|Fcall|Mcall|WNS|HET|Het)','match');
               songname = regexp(a.normlist(k).stimuliname,'([BRGOY]\d{3}|Fcall|Mcall|WNS|HET|Het)','match');
               
               % use autoseg to segment the elements
%                [rawy,fiy,I,syledge,eleedge] = main(fiy,fs,birdid,CONFIG)
%                
               eleids = find(~cellfun(@isempty, regexp([a.conspe_eleinf.songname].',songname)));
               local_eleinf = a.conspe_eleinf(eleids);
               
%                for w = 1: length(local_eleinf)
%                    local_eleinf(w).yini = local_eleinf(w).initial - local_eleinf(1).initial;
%                    local_eleinf(w).yter = local_eleinf(w).terminal - local_eleinf(1).initial;
%                end

% 为了方便命名和统一形式
                for w = 1: length(local_eleinf)
                    local_eleinf(w).yini = local_eleinf(w).initial - 0;
                    local_eleinf(w).yter = local_eleinf(w).terminal - 0;
                end

               
               for m = 1:length(local_eleinf)
                   counts = counts + 1;
                   split_fraginf(counts).y = local_eleinf(m).y;
                    split_fraginf(counts).fs = local_eleinf(m).fs;
                   split_fraginf(counts).padded_y = [local_eleinf(m).y;zeros(latency*local_eleinf(m).fs,1)]; % pad zeors with latency length
                   
                   split_fraginf(counts).initial = local_eleinf(m).yini/local_eleinf(m).fs;
                   split_fraginf(counts).terminal = local_eleinf(m).yter/local_eleinf(m).fs;
               
                   split_fraginf(counts).padded_sptimes = extract.sptimes(a.normlist(k).rawsptimes, split_fraginf(counts).initial...
                       ,split_fraginf(counts).terminal + latency);
                   
                   front_percentage = local_eleinf(m).yini/length(a.normlist(k).rawy);
                   back_percentage = local_eleinf(m).yter/length(a.normlist(k).rawy);
                   
                   featurenames = fieldnames(a.normlist(k).rawfeatures); % iterate for each feature
                   featurenames = setdiff(featurenames,{'file_index','file_name'});
                   for omega = 1: length(featurenames)
                       lenfeature = length(a.normlist(k).rawfeatures.(featurenames{omega}));
                       if  round(front_percentage*lenfeature)== 0
                           
                           split_fraginf(counts).(featurenames{omega}) = a.normlist(k).rawfeatures.(featurenames{omega})(...
                               1:round(back_percentage*lenfeature));
                       else
                           split_fraginf(counts).(featurenames{omega}) = a.normlist(k).rawfeatures.(featurenames{omega})(...
                               round(front_percentage*lenfeature):round(back_percentage*lenfeature));
                       end
                       split_fraginf(counts).(sprintf('mean_%s',featurenames{omega})) = mean(split_fraginf(counts).(featurenames{omega}));
                       
                   end
                   
               end
                
            end
            
            for n = 1: length(split_fraginf)
                tempsum = cal.psth_frag(split_fraginf(n).padded_y,split_fraginf(n).fs,split_fraginf(n).padded_sptimes);
                if isempty(tempsum)
                    halfsum = 0;
                else
                    halfsum = sum(tempsum(end/2:end));
                end
                fullsum = sum(tempsum);
                maxvalue = max(cal.psth_frag(split_fraginf(n).padded_y,split_fraginf(n).fs,split_fraginf(n).padded_sptimes));
                split_fraginf(n).maxvalue = maxvalue;
                split_fraginf(n).halfsum = halfsum;
                split_fraginf(n).fullsum = fullsum;
                if maxvalue > 6 % here the threshold is very important % originally set as 8
                    split_fraginf(n).label = 1;
                else
                    split_fraginf(n).label = 0;
                end
            end
            
            
        end
            
        function V1drawMeanFeaturesInSongVsRespAsLineChart(a) % draw the distribution of mean features
           
            fraglist = a.calResponseToWithinSongFragsFromEleinf;
            if isempty(fraglist)
                return;
            end
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).maxvalue;
            end
            
            toshow = struct;
           for u = 1: length(fraglist)
               toshow(u).amplitude = fraglist(u).mean_amplitude;
               toshow(u).pitch = fraglist(u).mean_pitch;
               toshow(u).AM = fraglist(u).mean_AM;
               toshow(u).mean_frequency = fraglist(u).mean_mean_frequency;
               toshow(u).FM = fraglist(u).mean_FM;
               toshow(u).peak_frequency = fraglist(u).mean_peak_frequency;
               toshow(u).entropy = fraglist(u).mean_entropy;
               toshow(u).resp = fraglist(u).responseMeasure;
           end
           
           toshow  = table2struct(sortrows(struct2table(toshow), 'resp','ascend'));  
           
           
           num_toshow = [toshow.amplitude; toshow.pitch; toshow.AM; toshow.mean_frequency; toshow.FM; ...
               toshow.peak_frequency; toshow.entropy].';
           
           znum_toshow = zscore(num_toshow,0,1);
%             
%            figure
%            plot(znum_toshow.');
           
           figure('Position',[1997 233 1388 658],'Color','w');
           hold on
           c = 1- rescale([toshow.resp].',0.1,1);
           for r = 1: length(znum_toshow)
               plot(znum_toshow(r,:),'Color',repmat(c(r),3,1));
           end
           
           colormap(flip(repmat(unique(c),1,3)))
           colorbar
           xlim([0,8]);
           xticks([0 1 2 3 4 5 6 7 8])
           xticklabels({'','amplitude','pitch','AM','mean_frequency','FM','peak-frequency','entropy',''});
           
           ylabel('Zscored Feature(averaged)');
           title(sprintf('Totally %u song elements',length(fraglist)));
           saveas(gcf,sprintf('V1-WithinSongsLineChartMeanFeaturesVsResp-%s.png',a.unique_neuronname));
           close(gcf);         
        end
        
        function V2drawMeanFeaturesInSongVsRespAsLineChart(a)
            songneuron = a.neurons{a.song_id};
            songneuron.drawInSongFragsMeanFeaturesVsRespAsLineChart
        end
        
        function new_draw_repla(a)
            
            d = Display(a.list);
            repla = d.findrepla();
            
            if isempty(repla)
                return
            end
            
            normlist = d.findnorm();
            
            replanames = {};
            for k = 1:length(repla)
                temp = strsplit(repla(k).stimuliname,'before')
                temp = regexp(temp{2},'[BGYOR]\d{3}','match');
                replanames{k} = temp{1};
            end
            replanames = replanames.';
            unique_replanames = unique(replanames);
            
            norm_ids = [];
            for k = 1: length(unique_replanames)
                norm_ids(k) = find(~cellfun(@isempty, regexp({normlist.stimuliname}.',unique_replanames{k})));
                repla(length(repla) + 1) = normlist(norm_ids(k));  
            end
            repla = Display.f0(repla);

            
            I = {};
            for idx = 1: length(repla)
                figure('Color','w','Position',[2031 536 732 281]);
                draw.three(repla(idx).f0y,repla(idx).fs,repla(idx).f0sptimes); % newer version of threeplot drawing method
                frame = getframe(gcf);
                I{idx} = frame.cdata;
                close(gcf);
            end
            
            
            lieshu = 3;
            hangshu = ceil(length(I)/lieshu);
            rest = lieshu*hangshu - length(I);
            white = uint8(255*ones(size(I{1})));
            
            if rest > 0
                for k = 1:rest
                    I = [I,white];
                    %                     ax = gcf;
                    %                     ax.Position(3) = 560;
                    %                     ax.Position(4) = 420;
                end
            end
            
            reshapedI = reshape(I, lieshu,[])';
            clear I
            IMG = cell2mat(reshapedI);
            
            imwrite(IMG,sprintf('ReplaThree-%s.png',a.unique_neuronname));
            
            %deg = Display.descend(deg);
            
              %replalist = d.info(find(~cellfun(@isempty, regexp({d.info.stimuliname}.','Repla|catego'))));
            
        end
        
        
        function new_draw_deg(a)

            dbstop if error
            
            d = Display(a.list);
            deg = d.finddeg();
            
            if isempty(deg)
                return
            end
            
            normlist = d.findnorm();
            
            deganames = {};
            for k = 1:length(deg)
                temp = regexp(deg(k).stimuliname,'[BGYOR]\d{3}|Fcall|Mcall|Het','match');
                deganames{k} = temp{1};
            end
            deganames = deganames.';
            unique_replanames = unique(deganames);
            
            norm_ids = [];
            for k = 1: length(unique_replanames)
                temp = find(~cellfun(@isempty, regexp({normlist.stimuliname}.',unique_replanames{k})));
                norm_ids(k) = temp(end); % if there are multiple values in temp, using the last one
                deg(length(deg) + 1) = normlist(norm_ids(k));  
            end
            deg = Display.f0(deg);

            
            I = {};
            for idx = 1: length(deg)
                figure('Color','w','Position',[2031 536 732 281]);
                draw.three(deg(idx).f0y,deg(idx).fs,deg(idx).f0sptimes); % newer version of threeplot drawing method
                frame = getframe(gcf);
                I{idx} = frame.cdata;
                close(gcf);
            end
            
            
            lieshu = 3;
            hangshu = ceil(length(I)/lieshu);
            rest = lieshu*hangshu - length(I);
            white = uint8(255*ones(size(I{1})));
            
            if rest > 0
                for k = 1:rest
                    I = [I,white];
                    %                     ax = gcf;
                    %                     ax.Position(3) = 560;
                    %                     ax.Position(4) = 420;
                end
            end
            
            reshapedI = reshape(I,[],hangshu)';
            clear I
            IMG = cell2mat(reshapedI);
            
            imwrite(IMG,sprintf('DegThree-%s.png',a.unique_neuronname));
            
            %deg = Display.descend(deg);
            
              %replalist = d.info(find(~cellfun(@isempty, regexp({d.info.stimuliname}.','Repla|catego'))));
            
        end

        function drawAlignedNormFragTwoPlots(a)
            dbstop if error
            tic
            
            % about Norms % Redudant code
            normlist = Analysis(a.neurons{a.song_id}).normlist;
            [~,postunique] = unique(cellfun(@convert.bid,[normlist.stimuliname].','Uni',0))
            normlist = normlist(postunique);


            % About Frag
            fragids = find(~cellfun(@isempty, regexp([a.list.stimuliname].','frag|Frag|syl|Syl') ));
            if ~isempty(fragids)
                
                fraglist = a.list(fragids);
%                 normlist = Analysis(a.neurons{a.song_id}).normlist;
%                 [~,postunique] = unique(cellfun(@convert.bid,[normlist.stimuliname].','Uni',0))
%                 normlist = normlist(postunique);
                
                for m = 1: length(fraglist)
                    
                    birdid = convert.bid(fraglist(m).stimuliname);
                    ids_norm = find(~cellfun(@isempty, regexp([normlist.stimuliname].',birdid) ) );
                    
                    if ~isempty(ids_norm)& length(ids_norm) == 1
                        fraglist(m).sylIni = Analysis.findIni(normlist(ids_norm).plty,fraglist(m).y);
                        fraglist(m).pady = [zeros([fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext,1]);fraglist(m).plty;zeros(length(normlist(ids_norm).plty)...
                            - (fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext) - length(fraglist(m).plty),1)];
                        fraglist(m).padsptimes = cellfun( @(x) x + (fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext)/fraglist(m).fs, fraglist(m).pltsptimes,'uni',0);
                    end
                end
            end
            
            % About Deg
            degids = find(~cellfun(@isempty, regexp([a.list.stimuliname].','Deg|deg') ));
            if ~isempty(degids)
 
                deglist = a.list(degids);
                for m = 1: length(deglist)
                    birdid = convert.bid(deglist(m).stimuliname);
                    ids_norm = find(~cellfun(@isempty, regexp([normlist.stimuliname].',birdid) ) );
                    if ~isempty(ids_norm)& length(ids_norm) == 1
                        [deglist(m).sylIni,trump_diffvalue] = Analysis.findIni(normlist(ids_norm).plty,deglist(m).y);
                        fprintf('%s has the DIFFVALUE as: %u and INI as: %f \n ',deglist(m).stimuliname,trump_diffvalue,deglist(m).sylIni/32000);
                        deglist(m).pady = [zeros([deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext,1]);deglist(m).plty;zeros(length(normlist(ids_norm).plty)...
                            - (deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext) - length(deglist(m).plty),1)];
                        deglist(m).padsptimes = cellfun( @(x) x + (deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext)/deglist(m).fs, deglist(m).pltsptimes,'uni',0);
                    end
                end
            end
            
            % merge the new fraglist and the deglist with the normlist
            I_of_each_column = {};
            for w = 1: length(normlist)

                if ~isempty(degids)
                    birdid = convert.bid(normlist(w).stimuliname);
                    ids_indeg = find(~cellfun(@isempty, regexp([deglist.stimuliname].',birdid) ) );
                    selected_deglist = deglist(ids_indeg);
                    [~,temp_index] = sortrows([selected_deglist.sylIni].');
                    selected_deglist = selected_deglist(temp_index);
                end

                if ~isempty(fragids)
                    birdid = convert.bid(normlist(w).stimuliname);
                    ids_infrag = find(~cellfun(@isempty, regexp([fraglist.stimuliname].',birdid) ) );
                    selected_fraglist = fraglist(ids_infrag);
                    [~,temp_index] = sortrows([selected_fraglist.sylIni].');
                    selected_fraglist = selected_fraglist(temp_index);
                end
                
                
                % draw the basic figure
                Icollect = {};
                figure('Color','w');

                draw.two(normlist(w).plty,normlist(w).fs,normlist(w).pltsptimes);
                xlabel(normlist(w).stimuliname);
                frame = getframe(gcf);
                Icollect{1} = frame.cdata;
                close(gcf)
                
                if ~isempty(degids)
                    for hh = 1: length(selected_deglist)
                        
                        figure('Color','w');
                        draw.two(selected_deglist(hh).pady,selected_deglist(hh).fs,selected_deglist(hh).padsptimes);
                        xlabel(selected_deglist(hh).stimuliname);
                        frame = getframe(gcf);
                        Icollect{1 + hh} = frame.cdata;
                        close(gcf);
                        
                    end
                end
                
                frozen_Icollect_len = length(Icollect);
                
                if ~isempty(fragids)
                    for bb = 1: length(selected_fraglist)
                        
                        figure('Color','w');
                        draw.two(selected_fraglist(bb).pady,selected_fraglist(bb).fs,selected_fraglist(bb).padsptimes);
                        xlabel(selected_fraglist(bb).stimuliname);
                        frame = getframe(gcf);
                        Icollect{frozen_Icollect_len + bb} = frame.cdata;
                        close(gcf);
                        
                    end
                end
                
                I_of_each_column{w} = vertcat(Icollect{:});
            end
            
            % this part draw individual elements from OTE songs
            if exist('fraglist','var')
                oteids = find(~cellfun(@isempty, regexp([fraglist.stimuliname].','OTE|ote|Ote') ));
                if ~isempty(oteids)
                    otelist = fraglist(oteids);
                    ote_Icollect = {};
                    for k = 1: ceil(length(otelist)/4)
                        % draw a figure for every four ote fragments
                        figure('Color','w');
                        for kk = flip([0:3])
                            if 4*k - kk <length(otelist)
                                subplot(2,4,4-kk);

                                draw.spec(otelist(4*k - kk).plty,otelist(4*k - kk).fs);
                                subplot(2,4,4-kk + 4 );
                                draw.raster(otelist(4*k - kk).pltsptimes,otelist(4*k - kk).plty,otelist(4*k - kk).fs);
                                ylabel('')
                            end
                        end
                        frame = getframe(gcf);
                        ote_Icollect{k} = frame.cdata;
                        close(gcf)

                    end
                    ote_img = vertcat( ote_Icollect{:});
                    I_of_each_column{length(I_of_each_column)+ 1} = ote_img;

                end
            end
            
            a.compareWaveforms;
            temp = getframe(gcf);
            w_img = temp.cdata;
            I_of_each_column{length(I_of_each_column)+ 1} = w_img;

            % padding each I based on the maximum size of local I
            size1 = [];

             for oo = 1: length(I_of_each_column)
                size1(oo) = size(I_of_each_column{oo},1);
            end

            [max_size1,max_oo] = max(size1);

            Ipad = {};
            for oo = 1: length(I_of_each_column)
                localI = I_of_each_column{oo};
                Ibase= uint8(256*ones(size(I_of_each_column{max_oo})));
                Ibase(1:size(localI,1),1:size(localI,2),1:size(localI,3)) = localI;
                Ipad{oo} = Ibase;
            end

            Iall = horzcat(Ipad{:});

            imwrite(Iall,sprintf('Aligned_Normfrag_%s.png',a.unique_neuronname));
            toc

        end

        function compareWaveforms(a)

            fig = figure('Color','w');
            cmap = colormap(flip(hsv(5)));
            for k = 1:length(a.neurons)
                local_waveform = a.neurons{k}.waveform;
                
              
                hold on
                plot(local_waveform.',':','Color',[cmap(k,:),0.45]);

                plot(max(local_waveform),'--','Color',cmap(k,:)*0.7);
                plot(min(local_waveform),'--','Color',cmap(k,:)*0.7);
                plot(mean(local_waveform),'-*','Color',cmap(k,:)*0.5);

            end

            title(sprintf('%u plexon files',length(a.neurons)));

           

        end

        function drawSeparatedWaveforms(a)
             a.compareWaveforms;
             saveas(gcf,sprintf('WaveformsSeparated-%s.png', a.unique_neuronname));
             close(fig);

        end
    end
    
    methods(Static)
        
         function check_detail_folder(dirpath,global_eleinf)
            % This function draw distribution of song elements in spectral
            % space ( scatter / small spectrogram-rasterPlot)
            
            % extract .wav file names in the target folder
            filenames = extract.filename(dirpath,'*.wav');
            nameonly = {};
            for k = 1: length(filenames)
                [~,nameonly{k},~] = fileparts(filenames{k});
            end
            nameonly = nameonly.';
            nameonly = [nameonly{:}].';
            
              % concatenate each songname and fragid for comparing
             global_merge = {};
            for w = 1: length(global_eleinf)
                global_merge{w} = sprintf('%s-%02u',global_eleinf(w).songname,global_eleinf(w).fragid);
            end
            
            
            % extract the target song-element from the name of the target
            % directory
            temp = split(dirpath,'\');
            temp = temp{end};
            temp = split(temp,'-');
            targetname = sprintf('%s-%s-%s',temp{2},temp{3},temp{4});
            [~,targetid] = ismember( targetname,global_merge);
            
            % extract element data of song element Fragments/Replacements
            fragids = find(~cellfun(@isempty, regexp(nameonly,'Frag')))
            fragnames = nameonly(fragids);
            fragname_remove_frag = split(fragnames,'Frag-');
            fragname_remove_frag = fragname_remove_frag(:,2);
            
            % label global_eleinf by frag data
            for k = 1: length(global_eleinf)
                global_eleinf(k).whether_frag_is_tested = 0;
            end
            
            for u = 1: length(fragname_remove_frag)
                [~,beta] = ismember( fragname_remove_frag{u},global_merge);
                global_eleinf(beta).whether_frag_is_tested = 1;
            end
            
            frag_1_eleinf = global_eleinf(find([global_eleinf.whether_frag_is_tested].' == 1));
            frag_0_eleinf = global_eleinf(find([global_eleinf.whether_frag_is_tested].' == 0));
            %--% draw
            figure; 
            hold on
            scatter([frag_0_eleinf.coor_1].', [frag_0_eleinf.coor_2].','k','filled');
            scatter([frag_1_eleinf.coor_1].', [frag_1_eleinf.coor_2].','r','filled');
            scatter(global_eleinf(targetid).coor_1,global_eleinf(targetid ).coor_2,[],'g','h');
            title(sprintf('Frag-%s',targetname));
            
            
             % extract data of replaced song
             replaids = find(~cellfun(@isempty, regexp(nameonly,'Repla')));
            replanames = nameonly(replaids);
            for e = 1: length(replanames)
                temp = split(replanames{e},'-before-');
                replanames{e} = temp{1};
            end
            replaname_remove_repla = split(replanames,'Repla-');
            if isempty(replaname_remove_repla)
                return % if there is no repla, the just return
            end
            replaname_remove_repla = replaname_remove_repla(:,2);
            
            
              % label global_eleinf by Repla data
            for k = 1: length(global_eleinf)
                global_eleinf(k).whether_repla_is_tested = 0;
            end
            
            for u = 1: length(replaname_remove_repla)
                [~,beta] = ismember( replaname_remove_repla{u},global_merge);
                global_eleinf(beta).whether_repla_is_tested = 1;
            end
            
            repla_1_eleinf = global_eleinf(find([global_eleinf.whether_repla_is_tested].' == 1));
            repla_0_eleinf = global_eleinf(find([global_eleinf.whether_repla_is_tested].' == 0));
            %--% draw
            figure; 
            hold on
            scatter([repla_0_eleinf.coor_1].', [repla_0_eleinf.coor_2].','k','filled');
            scatter([repla_1_eleinf.coor_1].', [repla_1_eleinf.coor_2].','r','filled');
            
            if global_eleinf(targetid).fragid ~=1
             scatter(global_eleinf(targetid).coor_1,global_eleinf(targetid-1).coor_2,[],'g','h');
            end
            title(sprintf('Repla-%s',targetname));
            
            
         end

         function [SynIni,diff_value] = findIni(y, yfrag) % find the correspoding initial timestamps by aliging two time series

             %y = B521; yfrag = B521_7;

%              counts = 0;
             diff_info = struct;
             parfor k = 1: length(y) - length(yfrag)
                 totest = y(k:k+ length(yfrag) -1);

                 if sum(totest) ~= 0

                     diff = sum(abs(totest - yfrag));

                     diff_info(k).diff = diff;
                     diff_info(k).kvalue = k;
                     %                      end
                 else

                     diff_info(k).diff = Inf;
                     diff_info(k).kvalue = k;
                 end


             end

             [~,min_ids] = min([diff_info.diff].');
             SynIni = diff_info(min_ids).kvalue;
             diff_value = diff_info(min_ids).diff;

         end

    end
    
end

