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
        
        targets % unique(targetnames)
        
    end
    
    methods
        
        
        function a = Analysis(neurons)
            
           if exist('neurons','var')
               a.neurons = neurons;
               for k = 1: length(a.neurons)
                   lists{k} = a.neurons{k}.todisplay;
               end
               a.list = vertcat(lists{:});
               % sort
               a.sort  % !!!!!!!!!!!!

           end   
            
        end
        
        function a = set_eleinf(a,eleinf) %set eleinf,eleinf is the info-data used to construct all the stimuli
            if isa(eleinf,'struct')
                a.all_eleinf = eleinf;
            elseif isa(eleinf,'string')|| isa(eleinf,'char')
                loaded = load(eleinf);
                a.all_eleinf = loaded.all_eleinf;
            end
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
            
            a.targets = unique(targetnames);
            
        end
        
        function primthree(a)
            
            
            for k = 1: length(neurons)
                selected{k}.three;
            end
         
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
            
           ids = find(~cellfun(@isempty, regexp([a.list.stimuliname].','Frag'))); % find all frags
           
           fraglist = a.list(ids);
           
           for n = 1: length(fraglist)
            tempsum = cal.psth_frag(fraglist(n).rawy,fraglist(n).fs,fraglist(n).rawsptimes);
            halfsum = sum(tempsum(end/2:end));
            fullsum = sum(tempsum);
            maxvalue = max(cal.psth_frag(fraglist(n).rawy,fraglist(n).fs,fraglist(n).rawsptimes))
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
            
            if ~exist('not_tested_handle','var') % to judge whether to dra not tested elements or not
                not_tested_handle = 0;
            end
            
            global_eleinf = a.all_eleinf;
            fraglist = a.judgeFragResponse;
            
            
            for k = 1: length(global_eleinf)
                uniqueid = sprintf('-%s-%u-',global_eleinf(k).songname, global_eleinf(k).fragid);
                
               if ~isempty (find(~cellfun(@isempty, regexp([fraglist.stimuliname].',uniqueid)), 1)) % if this element was tested
                     id_in_fraglist = find(~cellfun(@isempty, regexp([fraglist.stimuliname].',uniqueid)));
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
        
        function sort_frags_by_response_strength_and_then_draw(a)
            
            fraglist =  judgeFragResponse(a);
            sorted_fraglist = table2struct(sortrows( struct2table(fraglist) ,'halfsum','ascend'));
            
            % use montage
            I = {}; % collection of frag-response-three images
            for k = 1: length(sorted_fraglist)
                h =  figure;
                h.WindowState = 'maximized';
                draw.three(sorted_fraglist(k).plty,sorted_fraglist(k).fs,sorted_fraglist(k).sptimes);
                temp = getframe(gcf);
                I{k} = temp.cdata;
                
                close(h)
            end
            
            figure;
            montage(I);
            
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
    end
    
end

