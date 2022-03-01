classdef Sultan
    %Sultan
    
    % So Sultan controls everything, more powerful than Consul
    properties
        anas % a cell collection of neurons
        ana_filelist
        anainfo
        con_anas
    end
    
    methods
        function s = Sultan(analysis_filelist)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            f = waitbar(0,'Loading...');
            s.ana_filelist = analysis_filelist;
            for k = 1: length(analysis_filelist)
                
                load(analysis_filelist{k});
                s.anas{k} = A;
                waitbar(k/length(analysis_filelist),f,sprintf(' Now loading %u // %u files',k,length(analysis_filelist)));
            end
            
            close(f);
        end
        
        function s = selectConNeurons(s)
            
            selected_ids = [];
            for k = 1: length(s.neurons)
                
                if length(cellfun(@isempty, regexp('norm',{s.neurons{k}.slist.name}.')))> 12
                    selected_ids = [selected_ids,k];
                end
            end
            s.con_neurons = s.neurons{selected_ids};
        end % 需要修改
        
        function drawCONResponse(s)
            % extract binarized neuron's responses to CONs
            dbstop if error
            con_info = struct;
            counts = 0;
            
            common_cons = {};
            for k = 1: length(s.neurons)
                Conlist = Analysis(s.neurons{k}).evaluateConResponse;
                if length(Conlist) >= 18 % hard code here
                    counts = counts + 1;
                    con_match = regexp({Conlist.stimuliname}.','[BGYRO]\d{3}','match');
                    con_ids = find(~cellfun(@isempty,con_match));
                    spe_match = regexp({Conlist.stimuliname}.','TUT|BOS|Fcall|Mcall|WNS','match');
                    spe_ids = find(~cellfun(@isempty,spe_match));
                    % remove TUT and BOS
                    con_ids = setdiff(con_ids,spe_ids);
                    con_match = [con_match{con_ids}].';
                    % extract con_label
                    con_info(counts).neuronname = s.neurons{counts}.neuronname;
                    con_info(counts).con_match = con_match;
                    con_info(counts).con_resp = [Conlist(con_ids).label].';
                    % find out Common Conspecific songs
                    if counts == 1
                        common_cons = con_info(counts).con_match;
                    else
                        temp = intersect(common_cons,con_info(counts).con_match);
                        
                        common_cons = temp; % weird bugs happened when temp was not used
                        
                        if length(common_cons)< 18
                            pause
                        end
                        
                    end 
                end
            end
            
                 
            % extract binary resp map
            respmap = [];
            for m = 1: length(con_info)
                [~,loc] = ismember (common_cons,con_info(m).con_match);
                
               % new_con_info(m).name = {con_info(m).con_match{loc}}.';
                respmap(m,:) = con_info(m).con_resp(loc);
            end
            
            % Sort the matrix
            [~,new_rowids] = sort(sum(respmap,2),'descend');
            [~,new_columnids] = sort(sum(respmap,1),2,'descend');
            tempo =  respmap(new_rowids,:);
            new_respmap = tempo(:,new_columnids);
            new_cnames = {common_cons{new_columnids}}; %#ok<CCAT1>
            tempo = {con_info.neuronname};
            new_rnames = {tempo{new_rowids}};
            
%             figure
%             subplot(1,3,1)
%             imagesc(respmap)
%             subplot(1,3,2)
%             imagesc( respmap(new_rowids,:));
%             subplot(1,3,3)
%             tempo =  respmap(new_rowids,:);
%             imagesc( tempo(:,new_columnids));


            respT = array2table(respmap);
            respT.Properties.VariableNames = new_cnames;
            respT.Properties.RowNames = new_rnames;
            
            
        end
        
        function wl_info = calWavelength(s)
            wl_info = struct;
            
            for k = 1: length(s.anas)
                wl_info(k).wl = s.anas{k}.calMeanWaveLength;
                wl_info(k).neuronname = s.anas{k}.neuronname;
            end
            
        end
        
        function frags_exist_names = markNeuronsWithFragsAsStimuli(s)
            
            frags_exist_names = {};
            counts = 0;
            for k = 1: length(s.anas)
                if ~isempty(~arrayfun(@isempty, regexp([s.anas{k}.neuinfo.keywords],"frag")) )
                    counts = counts + 1;
                   frags_exist_names{counts} = s.anas{k}.unique_neuronname;
                end
            end
        end
        
        function histWaveLength(s)
            dbstop if error       
%             figure
%             imagesc(new_respmap);
%             xticks(1:length(new_cnames));
%             xticklabels(new_cnames);
%             yticks(1:length(new_rnames));
%             yticklabels(new_rnames);
%             set(gca,'TickLabelInterpreter','none');
           
            wl_info = struct;
            
            for k = 1: length(s.anas)
                 wl_info(k).wl = s.anas{k}.calMeanWaveLength;
                 wl_info(k).neuronname = s.anas{k}.unique_neuronname;
            end
            
            figure('Color','w');
            histogram([wl_info.wl].',30); % originally 15
            
            xlabel('Spike width (ms)');
            ylabel('Neurons per bin');
            
        end
        
        function drawCONSPEResponse(s)
            
             % extract binarized neuron's responses to CONs
            dbstop if error
            con_info = struct;
            counts = 0;
            common_cons = {};
            common_spes = {};
            
            for k = 1: length(s.anas)
                this_neuron = s.anas{k}.neurons{s.anas{k}.song_id};
                Conlist = Analysis(this_neuron).evaluateConResponse;
                if length(Conlist) >= 18 % hard code here
                    counts = counts + 1;
                    con_match = regexp({Conlist.stimuliname}.','[BGYRO]\d{3}','match');
                    con_ids = find(~cellfun(@isempty,con_match));
                    spe_match = regexp({Conlist.stimuliname}.','TUT|BOS|Fcall|Mcall|WNS','match');
                    spe_ids = find(~cellfun(@isempty,spe_match));
                    spe_match = [spe_match{spe_ids}].';
                    % remove TUT and BOS
                    con_ids = setdiff(con_ids,spe_ids);
                    con_match = [con_match{con_ids}].';
                    
                   
                    % add TUT-BOS-Fcall-Mcall-WNS (Special)
                    
                    % extract con_label
                    con_info(counts).wav_len = s.anas{k}.calMeanWaveLength;
                    con_info(counts).neuronname = s.anas{k}.neuronname;
                    con_info(counts).con_match = con_match;
                    con_info(counts).con_resp = [Conlist(con_ids).label].';
                    con_info(counts).spe_match = spe_match;
                    con_info(counts).spe_resp = [Conlist(spe_ids).label].';
                    % find out Common Conspecific songs
                    if counts == 1
                        common_cons = con_info(counts).con_match;
                        %common_spes = con_info(counts).spe_match;
                    else
                        temp = intersect(common_cons,con_info(counts).con_match);
                        
                        common_cons = temp; % weird bugs happened when temp was not used
                        
%                         
%                         temp = intersect(common_spes,con_info(counts).spe_match);
%                         common_spes = temp;
                        
                        if length(common_cons)< 18
                            pause
                        end
                        
                    end 
                end
            end
            
                 
            % extract binary con-resp map
            respmap = [];
            for m = 1: length(con_info)
                [~,loc] = ismember (common_cons,con_info(m).con_match);
                
               % new_con_info(m).name = {con_info(m).con_match{loc}}.';
                respmap(m,:) = con_info(m).con_resp(loc);
            end
            
            
            % extract spe-resp map
            common_spes = {'TUT','BOS','Fcall','Mcall','WNS'};
            spemap = [];
            for m = 1: length(con_info)
                
                if ~isempty(con_info(m).spe_match)
                    [~,loc] = ismember (common_spes,con_info(m).spe_match);
                else
                    [~,loc] = ismember (common_spes,{});
                end
                
                % new_con_info(m).name = {con_info(m).con_match{loc}}.';
                for mm = 1: length(loc)
                    if loc(mm) == 0
                        spemap(m,mm) = nan;
                    else
                        spemap(m,mm) = con_info(m).spe_resp(loc(mm));
                    end
                end
            end
            
            
            % Sort the matrix
            [~,new_rowids] = sort(sum(respmap,2),'descend');
            [~,new_columnids] = sort(sum(respmap,1),2,'descend');
            tempo =  respmap(new_rowids,:);
            new_respmap = tempo(:,new_columnids);
            new_cnames = {common_cons{new_columnids}}; %#ok<CCAT1>
            tempo = {con_info.neuronname};
            new_rnames = {tempo{new_rowids}};

            respT = array2table(respmap);
            respT.Properties.VariableNames = new_cnames;
            respT.Properties.RowNames = new_rnames;
             
            % sort spe-map
            new_spemap = spemap(new_rowids,:);
            
            
            % reordered waveforms length
            
            WLs = [con_info.wav_len].';
            WLs = WLs(new_rowids);
            
            
            for k = 1: size(new_respmap,1)
                value_is_one_ids = find(new_respmap(k,:));
                new_respmap(k,value_is_one_ids) = WLs(k);    
            end
            
            figure;
            colormap('summer');
            imagesc(new_respmap);
            
            colorbar;
            
            xticks(1:size(new_respmap,2));
            xticklabels(new_cnames);
            
            
            yticks(1:size(new_respmap,1));
            yticklabels(new_rnames);
            set(gca,'TickLabelInterpreter','none');
            set(gca,'TickLength',[0.001, 0.001]);
            
        end
        
        function drawCONSPEResponse_NSBS(s)
            
             % extract binarized neuron's responses to CONs
            dbstop if error
            con_info = struct;
            counts = 0;
            common_cons = {};
            common_spes = {};
            
            for k = 1: length(s.anas)
                this_neuron = s.anas{k}.neurons{s.anas{k}.song_id};
                Conlist = Analysis(this_neuron).evaluateConResponse;
                if length(Conlist) >= 18 % hard code here
                    counts = counts + 1;
                    con_match = regexp({Conlist.stimuliname}.','[BGYRO]\d{3}','match');
                    con_ids = find(~cellfun(@isempty,con_match));
                    spe_match = regexp({Conlist.stimuliname}.','TUT|BOS|Fcall|Mcall|WNS','match');
                    spe_ids = find(~cellfun(@isempty,spe_match));
                    spe_match = [spe_match{spe_ids}].';
                    % remove TUT and BOS
                    con_ids = setdiff(con_ids,spe_ids);
                    con_match = [con_match{con_ids}].';
                    
                   
                    % add TUT-BOS-Fcall-Mcall-WNS (Special)
                    
                    % extract con_label
                    con_info(counts).wav_len = this_neuron.calMeanWaveLength;
                    con_info(counts).neuronname = this_neuron.unique_neuronname;
                    con_info(counts).con_match = con_match;
                    con_info(counts).con_resp = [Conlist(con_ids).label].';
                    con_info(counts).spe_match = spe_match;
                    con_info(counts).spe_resp = [Conlist(spe_ids).label].';
                    % find out Common Conspecific songs
                    if counts == 1
                        common_cons = con_info(counts).con_match;
                        %common_spes = con_info(counts).spe_match;
                    else
                        temp = intersect(common_cons,con_info(counts).con_match);
                        
                        common_cons = temp; % weird bugs happened when temp was not used                               
%                         temp = intersect(common_spes,con_info(counts).spe_match);
%                         common_spes = temp;
                        
                        if length(common_cons)< 18
                            pause
                        end
                    end 
                    
                end
            end
            
                 
            % extract binary con-resp map
            respmap = [];
            for m = 1: length(con_info)
                [~,loc] = ismember (common_cons,con_info(m).con_match);
                
               % new_con_info(m).name = {con_info(m).con_match{loc}}.';
                respmap(m,:) = con_info(m).con_resp(loc);
            end
            
            
            % extract spe-resp map
            common_spes = {'TUT','BOS','Fcall','Mcall','WNS'};
            spemap = [];
            for m = 1: length(con_info)
                
                if ~isempty(con_info(m).spe_match)
                    [~,loc] = ismember (common_spes,con_info(m).spe_match);
                else
                    [~,loc] = ismember (common_spes,{});
                end
                
                % new_con_info(m).name = {con_info(m).con_match{loc}}.';
                for mm = 1: length(loc)
                    if loc(mm) == 0
                        spemap(m,mm) = nan;
                    else
                        spemap(m,mm) = con_info(m).spe_resp(loc(mm));
                    end
                end
            end
            
            
            % Sort the matrix
            [~,new_rowids] = sort(sum(respmap,2),'descend');
            [~,new_columnids] = sort(sum(respmap,1),2,'descend');
            tempo =  respmap(new_rowids,:);
            new_respmap = tempo(:,new_columnids);
            new_cnames = {common_cons{new_columnids}}; %#ok<CCAT1>
            tempo = {con_info.neuronname};
            new_rnames = {tempo{new_rowids}};

            respT = array2table(respmap);
            respT.Properties.VariableNames = new_cnames;
            respT.Properties.RowNames = new_rnames;
             
            % sort spe-map
            new_spemap = spemap(new_rowids,:);
            
            
            % reordered waveforms length
            
            WLs = [con_info.wav_len].';
            WLs = WLs(new_rowids);
            
            wlfr_info = plotWavelengthVsFiringRate(s);
            
            names_in_wlfr = {wlfr_info.neuronname}.';
            [~,ids_for_reorder] = ismember(new_rnames,names_in_wlfr);
            new_wlfr_info = wlfr_info(ids_for_reorder);
            
            for k = 1: size(new_respmap,1)
                value_is_one_ids = find(new_respmap(k,:));
                if new_wlfr_info(k).isBS == 0
                    new_respmap(k,value_is_one_ids) = 0.5;
                elseif new_wlfr_info(k).isBS == 1
                    new_respmap(k,value_is_one_ids) = 1;
                end  
            end
            
            fig = figure;
            % colormap('summer');
             imagesc(new_respmap);
            
            map = [ 1 1 1
                0 0.4470 0.7410
                0.8500 0.3250 0.0980];
            
            colormap(fig,map);
            
            %colorbar;
            
            xticks(1:size(new_respmap,2));
            xticklabels(new_cnames);
            
            
            yticks(1:size(new_respmap,1));
            yticklabels(new_rnames);
            set(gca,'TickLabelInterpreter','none');
            set(gca,'TickLength',[0.001, 0.001]);
            
        end        
        
        function drawCONSPEResponse_NSBS_markFrag(s)
            
             % extract binarized neuron's responses to CONs
            dbstop if error
            con_info = struct;
            counts = 0;
            common_cons = {};
            common_spes = {};
            
            for k = 1: length(s.anas)
                this_neuron = s.anas{k}.neurons{s.anas{k}.song_id};
                Conlist = Analysis(this_neuron).evaluateConResponse;
                if length(Conlist) >= 18 % hard code here
                    counts = counts + 1;
                    con_match = regexp({Conlist.stimuliname}.','[BGYRO]\d{3}','match');
                    con_ids = find(~cellfun(@isempty,con_match));
                    spe_match = regexp({Conlist.stimuliname}.','TUT|BOS|Fcall|Mcall|WNS','match');
                    spe_ids = find(~cellfun(@isempty,spe_match));
                    spe_match = [spe_match{spe_ids}].';
                    % remove TUT and BOS
                    con_ids = setdiff(con_ids,spe_ids);
                    con_match = [con_match{con_ids}].';
                    
                   
                    % add TUT-BOS-Fcall-Mcall-WNS (Special)
                    
                    % extract con_label
                    con_info(counts).wav_len = this_neuron.calMeanWaveLength;
                    con_info(counts).neuronname = this_neuron.unique_neuronname;
                    con_info(counts).con_match = con_match;
                    con_info(counts).con_resp = [Conlist(con_ids).label].';
                    con_info(counts).spe_match = spe_match;
                    con_info(counts).spe_resp = [Conlist(spe_ids).label].';
                    % find out Common Conspecific songs
                    if counts == 1
                        common_cons = con_info(counts).con_match;
                        %common_spes = con_info(counts).spe_match;
                    else
                        temp = intersect(common_cons,con_info(counts).con_match);
                        
                        common_cons = temp; % weird bugs happened when temp was not used                               
%                         temp = intersect(common_spes,con_info(counts).spe_match);
%                         common_spes = temp;
                        
                        if length(common_cons)< 18
                            pause
                        end
                    end 
                    
                end
            end
            
                 
            % extract binary con-resp map
            respmap = [];
            for m = 1: length(con_info)
                [~,loc] = ismember (common_cons,con_info(m).con_match);
                
               % new_con_info(m).name = {con_info(m).con_match{loc}}.';
                respmap(m,:) = con_info(m).con_resp(loc);
            end
            
            % Sort the matrix
            [~,new_rowids] = sort(sum(respmap,2),'descend');
            [~,new_columnids] = sort(sum(respmap,1),2,'descend');
            tempo =  respmap(new_rowids,:);
            new_respmap = tempo(:,new_columnids);
            new_cnames = {common_cons{new_columnids}}; %#ok<CCAT1>
            tempo = {con_info.neuronname};
            new_rnames = {tempo{new_rowids}};

            respT = array2table(respmap);
            respT.Properties.VariableNames = new_cnames;
            respT.Properties.RowNames = new_rnames;
             
            % sort spe-map
           
            
            
            % reordered waveforms length
            
            WLs = [con_info.wav_len].';
            WLs = WLs(new_rowids);
            
            wlfr_info = plotWavelengthVsFiringRate(s);
            
            names_in_wlfr = {wlfr_info.neuronname}.';
            [~,ids_for_reorder] = ismember(new_rnames,names_in_wlfr);
            new_wlfr_info = wlfr_info(ids_for_reorder);
            
            for k = 1: size(new_respmap,1)
                value_is_one_ids = find(new_respmap(k,:));
                if new_wlfr_info(k).isBS == 0
                    new_respmap(k,value_is_one_ids) = 1;
                elseif new_wlfr_info(k).isBS == 1
                    new_respmap(k,value_is_one_ids) = 3;
                end  
            end
            
            frags_exist_names = s.markNeuronsWithFragsAsStimuli;
            
            for k = 1: size(new_respmap,1)
                value_is_one_ids = find(new_respmap(k,:));
                if ismember(new_wlfr_info(k).neuronname,frags_exist_names) && new_wlfr_info(k).isBS == 0
                    new_respmap(k,value_is_one_ids) = 2;
                elseif ismember(new_wlfr_info(k).neuronname,frags_exist_names) && new_wlfr_info(k).isBS == 1
                    new_respmap(k,value_is_one_ids) = 4;
                end  
            end
            
            
            
            fig = figure;
            % colormap('summer');
             imagesc(new_respmap);
             
             map = [ 1 1 1  % white for 0
                 0 0.4470 0.7410 % shallow blue for 1
                 0 0 1          % bright blue for 2
                 0.8500 0.3250 0.0980  % shallow red for 3
                 1 0 0];        % bright red for 4
           % figure; viscircles([0,0],10,'Color',[0 0 1]);
            colormap(fig,map);
            caxis manual
            caxis([0 4]);
           % colorbar
            
            %colorbar;
            
            xticks(1:size(new_respmap,2));
            xticklabels(new_cnames);
            
            
            yticks(1:size(new_respmap,1));
            yticklabels(new_rnames);
            set(gca,'TickLabelInterpreter','none');
            set(gca,'TickLength',[0.001, 0.001]);
            
        end
        
        function drawCONSPEResponse_NSBS_FragOnly(s)
            
             % extract binarized neuron's responses to CONs
            dbstop if error
            con_info = struct;
            counts = 0;
            common_cons = {};
            common_spes = {};
            
            for k = 1: length(s.anas)
                this_neuron = s.anas{k}.neurons{s.anas{k}.song_id};
                Conlist = Analysis(this_neuron).evaluateConResponse;
                if length(Conlist) >= 18 % hard code here
                    counts = counts + 1;
                    con_match = regexp({Conlist.stimuliname}.','[BGYRO]\d{3}','match');
                    con_ids = find(~cellfun(@isempty,con_match));
                    spe_match = regexp({Conlist.stimuliname}.','TUT|BOS|Fcall|Mcall|WNS','match');
                    spe_ids = find(~cellfun(@isempty,spe_match));
                    spe_match = [spe_match{spe_ids}].';
                    % remove TUT and BOS
                    con_ids = setdiff(con_ids,spe_ids);
                    con_match = [con_match{con_ids}].';
                    
                   
                    % add TUT-BOS-Fcall-Mcall-WNS (Special)
                    
                    % extract con_label
                    con_info(counts).wav_len = this_neuron.calMeanWaveLength;
                    con_info(counts).neuronname = this_neuron.unique_neuronname;
                    con_info(counts).con_match = con_match;
                    con_info(counts).con_resp = [Conlist(con_ids).label].';
                    con_info(counts).spe_match = spe_match;
                    con_info(counts).spe_resp = [Conlist(spe_ids).label].';
                    % find out Common Conspecific songs
                    if counts == 1
                        common_cons = con_info(counts).con_match;
                        %common_spes = con_info(counts).spe_match;
                    else
                        temp = intersect(common_cons,con_info(counts).con_match);
                        
                        common_cons = temp; % weird bugs happened when temp was not used                               
%                         temp = intersect(common_spes,con_info(counts).spe_match);
%                         common_spes = temp;
                        
                        if length(common_cons)< 18
                            pause
                        end
                    end 
                    
                end
            end
            
                 
            % extract binary con-resp map
            respmap = [];
            for m = 1: length(con_info)
                [~,loc] = ismember (common_cons,con_info(m).con_match);
                
               % new_con_info(m).name = {con_info(m).con_match{loc}}.';
                respmap(m,:) = con_info(m).con_resp(loc);
            end
            
            % Sort the matrix
            [~,new_rowids] = sort(sum(respmap,2),'descend');
            [~,new_columnids] = sort(sum(respmap,1),2,'descend');
            tempo =  respmap(new_rowids,:);
            new_respmap = tempo(:,new_columnids);
            new_cnames = {common_cons{new_columnids}}; %#ok<CCAT1>
            tempo = {con_info.neuronname};
            new_rnames = {tempo{new_rowids}};

            respT = array2table(respmap);
            respT.Properties.VariableNames = new_cnames;
            respT.Properties.RowNames = new_rnames;
             
           
            
            % reordered waveforms length
            
            WLs = [con_info.wav_len].';
            WLs = WLs(new_rowids);
            
            wlfr_info = plotWavelengthVsFiringRate(s);
            
            names_in_wlfr = {wlfr_info.neuronname}.';
            [~,ids_for_reorder] = ismember(new_rnames,names_in_wlfr);
            new_wlfr_info = wlfr_info(ids_for_reorder);
            
            for k = 1: size(new_respmap,1)
                value_is_one_ids = find(new_respmap(k,:));
                if new_wlfr_info(k).isBS == 0
                    new_respmap(k,value_is_one_ids) = 1;
                elseif new_wlfr_info(k).isBS == 1
                    new_respmap(k,value_is_one_ids) = 3;
                end  
            end
            
            frags_exist_names = s.markNeuronsWithFragsAsStimuli;
            
            for k = 1: size(new_respmap,1)
                value_is_one_ids = find(new_respmap(k,:));
                if ismember(new_wlfr_info(k).neuronname,frags_exist_names) && new_wlfr_info(k).isBS == 0
                    new_respmap(k,value_is_one_ids) = 2;
                elseif ismember(new_wlfr_info(k).neuronname,frags_exist_names) && new_wlfr_info(k).isBS == 1
                    new_respmap(k,value_is_one_ids) = 4;
                end  
            end
            
            labled_tokeep = [];
            for k = 1: size(new_respmap,1)
                if ismember(4,new_respmap(k,:))
                    labled_tokeep = [labled_tokeep,k];
                end
            end
            
            labeled_respmap = new_respmap(labled_tokeep,:);
            
            
            fig = figure;
            % colormap('summer');
             imagesc(labeled_respmap);
             
             map = [ 1 1 1  % white for 0
%                  0 0.4470 0.7410 % shallow blue for 1
%                  0 0 1          % bright blue for 2
%                  0.8500 0.3250 0.0980  % shallow red for 3
                 1 0 0];        % bright red for 4
           % figure; viscircles([0,0],10,'Color',[0 0 1]);
            colormap(fig,map);
            caxis manual
            caxis([0 4]);
           % colorbar
            
            %colorbar;
            
            xticks(1:size(labeled_respmap,2));
            xticklabels(new_cnames);
            
            
            yticks(1:size(new_respmap,1));
            yticklabels({new_rnames{labled_tokeep}});
            set(gca,'TickLabelInterpreter','none');
            set(gca,'TickLength',[0.001, 0.001]);
            
            
            figure;
            
            summed_resps = sum(labeled_respmap,2)/4;
            
            histogram(summed_resps);
            xlabel('Number of response-eliciting songs');
            ylabel ('Number of BS neurons');
            
        end
          
        function neuronStimuliInfo = getStimuliInfo(s)
            
            neuronStimuliInfo = struct;
            
            for k = 1: length(s.anas)
                these_names = {s.anas{k}.slist.name}.';
                
                % initialize
                neuronStimuliInfo(k).neuronname = s.neurons{k}.neuronname;
                neuronStimuliInfo(k).norm = 0;
                neuronStimuliInfo(k).frag = 0;
                neuronStimuliInfo(k).repla = 0;
                
                if ~isempty(find(~cellfun(@isempty, regexp(these_names,'norm'))))
                    neuronStimuliInfo(k).norm = 1; 
                end
                
                if ~isempty(find(~cellfun(@isempty, regexp(these_names,'frag|Frag|syl|Syl'))))
                    neuronStimuliInfo(k).frag = 1; 
                end
                
                if ~isempty(find(~cellfun(@isempty, regexp(these_names,'Repla|repla|catego'))))
                    neuronStimuliInfo(k).repla = 1; 
                end
                
            end
            
            sumnorm = sum([neuronStimuliInfo.norm].');
            sumfrag = sum([neuronStimuliInfo.frag].');
            sumrepla = sum([neuronStimuliInfo.repla].');
            
            %????????????????????????????????? something is wrong
            
        end
        
        function fr_info = multiRepeatsFiringRate(s)
            fr_info = struct;
             
            for k = 1: length(s.anas)
                % for each neuron
               
                thisA = s.anas{k};
                fr_info(k).neuronname = s.anas{k}.neuronname;
                
                sum_prelen = 0; % summed prey length 
                concat_presptimes = []; % concatenated prey sptimes
                
                sum_pltlen = 0; %summed prey( stimuli y, not plty or rawy) length 
                concat_pltsptimes = []; %  % concatenated y sptimes
                
                all_es = thisA.getAllEphysObject;
                for m = 1: length(all_es)
                    
                    
                    % for prey 
                    fr_info(k).presptimes{m} = all_es{m}.presptimes
                    fr_info(k).preylen{m} = length(all_es{m}.y)/all_es{m}.fs;
                    fr_info(k).repnum{m} = size(all_es{m}.presptimes,2);
                    temp = all_es{m}.presptimes.';
                    concat_presptimes = [concat_presptimes;vertcat(vertcat(temp{:}))+ sum_prelen];
                    sum_prelen = sum_prelen +  fr_info(k).preylen{m};
                    
                    % for plty
                    fr_info(k).pltsptimes{m} = all_es{m}.pltsptimes
                    fr_info(k).pltlen{m} = length(all_es{m}.plty)/all_es{m}.fs;
                    temp = all_es{m}.pltsptimes.';
                    concat_pltsptimes = [concat_pltsptimes;vertcat(vertcat(temp{:}))+ sum_pltlen];
                    sum_pltlen = sum_pltlen +  fr_info(k).pltlen{m};
                    
                end
                
                % for pre_y
                fr_info(k).concat_pre_sptimes = concat_presptimes;
                fr_info(k).concat_pre_len = sum_prelen;
                fr_info(k).mean_pre_fr = length(concat_presptimes)/sum_prelen;
                
                % for plt_y
                fr_info(k).concat_plt_sptimes = concat_pltsptimes;
                fr_info(k).concat_plt_len = sum_pltlen;
                fr_info(k).mean_plt_fr = length(concat_pltsptimes)/sum_pltlen;
                
                
            end
            
        end
        
        function wlfr_info = plotWavelengthVsFiringRate(s)
            
            wl_info = calWavelength(s);
            fr_info = multiRepeatsFiringRate(s);

            ids_in_fr = [];
            for k = 1: length(wl_info)
                this_nname = wl_info(k).neuronname;
                ids_in_fr(k) = find(strcmp(this_nname,{fr_info.neuronname}.'));
            end
            
            fr_info = fr_info(ids_in_fr); % re-order
            
            fr_info = rmfield( fr_info,'neuronname');

            wlfr_info = table2struct([struct2table(wl_info),struct2table(fr_info)]);

            fig = figure;
            h = scatter([wlfr_info.wl],[wlfr_info.mean_plt_fr]);
            xlabel('Mean Spike width');
            ylabel('Spontaneous Firing Rate');
            title('Brush BS neurons!');
            brush on
            pause(10) 
            brush off
            title('');
            brushdata = logical(get(h, 'BrushData'));
            brushed_ids = find(brushdata);
            not_brushed_ids = find(brushdata == 0);
            [wlfr_info(:).isBS] = deal(0);
            for w = 1: length(brushed_ids)
                wlfr_info(brushed_ids(w)).isBS = 1;
            end
            close(fig);
            
            
            
            figure;
            hold on
            scatter([wlfr_info(brushed_ids).wl],[wlfr_info(brushed_ids).mean_plt_fr],[], [0.8500 0.3250 0.0980],'filled');
            scatter([wlfr_info(not_brushed_ids).wl],[wlfr_info(not_brushed_ids).mean_plt_fr],[],[0 0.4470 0.7410],'filled');
            xlabel('Mean Spike width');
            ylabel('Spontaneous Firing Rate');
             
            %disp(brushed);

%             figure;
%             scatter([wlfr_info.wl],[wlfr_info.mean_plt_fr]);
%             xlabel('Mean Spike width');
%             ylabel('Evoked Firing Rate');
        end
        
        function arrangeThreePlotsByRespMapAndDraw(s)
            
            % extract binarized neuron's responses to CONs
            dbstop if error
            con_info = struct;
            counts = 0;
            common_cons = {};
            common_spes = {};
            
            for k = 1: length(s.anas)
                this_neuron = s.anas{k}.neurons{s.anas{k}.song_id};
                Conlist = Analysis(this_neuron).evaluateConResponse;
                if length(Conlist) >= 18 % hard code here
                    counts = counts + 1;
                    con_match = regexp({Conlist.stimuliname}.','[BGYRO]\d{3}','match');
                    con_ids = find(~cellfun(@isempty,con_match));
                    spe_match = regexp({Conlist.stimuliname}.','TUT|BOS|Fcall|Mcall|WNS','match');
                    spe_ids = find(~cellfun(@isempty,spe_match));
                    spe_match = [spe_match{spe_ids}].';
                    % remove TUT and BOS
                    con_ids = setdiff(con_ids,spe_ids);
                    con_match = [con_match{con_ids}].';
                    
                   
                    % add TUT-BOS-Fcall-Mcall-WNS (Special)
                    
                    % extract con_label
                    con_info(counts).wav_len = this_neuron.calMeanWaveLength;
                    con_info(counts).neuronname = this_neuron.unique_neuronname;
                    con_info(counts).con_match = con_match;
                    con_info(counts).con_resp = [Conlist(con_ids).label].'; % resp is just label
                    con_info(counts).spe_match = spe_match;
                    con_info(counts).spe_resp = [Conlist(spe_ids).label].';
                    % find out Common Conspecific songs
                    if counts == 1
                        common_cons = con_info(counts).con_match;
                        %common_spes = con_info(counts).spe_match;
                    else
                        temp = intersect(common_cons,con_info(counts).con_match);
                        
                        common_cons = temp; % weird bugs happened when temp was not used                               
%                         temp = intersect(common_spes,con_info(counts).spe_match);
%                         common_spes = temp;
                        
                        if length(common_cons)< 18
                            pause
                        end
                    end 
                    
                end
            end
                 
            % extract binary con-resp map
            respmap = [];
            for m = 1: length(con_info)
                [~,loc] = ismember (common_cons,con_info(m).con_match);
                
               % new_con_info(m).name = {con_info(m).con_match{loc}}.';
                respmap(m,:) = con_info(m).con_resp(loc);
            end
            
            
            % extract spe-resp map
            common_spes = {'TUT','BOS','Fcall','Mcall','WNS'};
            spemap = [];
            for m = 1: length(con_info)
                
                if ~isempty(con_info(m).spe_match)
                    [~,loc] = ismember (common_spes,con_info(m).spe_match);
                else
                    [~,loc] = ismember (common_spes,{});
                end
                
                % new_con_info(m).name = {con_info(m).con_match{loc}}.';
                for mm = 1: length(loc)
                    if loc(mm) == 0
                        spemap(m,mm) = nan;
                    else
                        spemap(m,mm) = con_info(m).spe_resp(loc(mm));
                    end
                end
            end
            
            % Sort the matrix
            [~,new_rowids] = sort(sum(respmap,2),'descend');
            [~,new_columnids] = sort(sum(respmap,1),2,'descend');
            tempo =  respmap(new_rowids,:);
            new_respmap = tempo(:,new_columnids);
            new_cnames = {common_cons{new_columnids}}; %#ok<CCAT1>
            tempo = {con_info.neuronname};
            new_rnames = {tempo{new_rowids}};

            respT = array2table(respmap);
            respT.Properties.VariableNames = new_cnames;
            respT.Properties.RowNames = new_rnames;
             
            % sort spe-map
            new_spemap = spemap(new_rowids,:);
            
            % reordered waveforms length
            
            WLs = [con_info.wav_len].';
            WLs = WLs(new_rowids);
            
            wlfr_info = plotWavelengthVsFiringRate(s);
            
            names_in_wlfr = {wlfr_info.neuronname}.';
            [~,ids_for_reorder] = ismember(new_rnames,names_in_wlfr);
            new_wlfr_info = wlfr_info(ids_for_reorder);
            
            for k = 1: size(new_respmap,1)
                value_is_one_ids = find(new_respmap(k,:));
                if new_wlfr_info(k).isBS == 0
                    new_respmap(k,value_is_one_ids) = 1;
                elseif new_wlfr_info(k).isBS == 1
                    new_respmap(k,value_is_one_ids) = 3;
                end  
            end
            
            frags_exist_names = s.markNeuronsWithFragsAsStimuli;
            
            for k = 1: size(new_respmap,1)
                value_is_one_ids = find(new_respmap(k,:));
                if ismember(new_wlfr_info(k).neuronname,frags_exist_names) && new_wlfr_info(k).isBS == 0
                    new_respmap(k,value_is_one_ids) = 2;
                elseif ismember(new_wlfr_info(k).neuronname,frags_exist_names) && new_wlfr_info(k).isBS == 1
                    new_respmap(k,value_is_one_ids) = 4;
                end  
            end
               
            % draw the three plots
            Icollect = {}; % to collect figure frames for each pairwise three plots
            
            for z = 1: length(new_rowids)
               new_N_id =  new_rowids(z);  % reorder the neurons ranked by numnber of eliciting songs
               
               [e_songs,slist_songs] = s.anas{new_N_id}.neurons{s.anas{new_N_id}.song_id}.onlyExtractEphysOfSongs;
               converted_names = cellfun(@convert.bid,{slist_songs.name}.','UniformOutput',0);
               [~,order_for_three_plots] = ismember(new_cnames, converted_names);
               
               
               for hh = 1: length(order_for_three_plots)
                   e_songs{order_for_three_plots(hh)}.three;
                   subplot(3,1,3);
                   xlabel(sprintf('Neuron:%s---Stimuli:%s',s.anas{new_N_id}.unique_neuronname,...
                       convert.bid(e_songs{order_for_three_plots(hh)}.sound.name)),'Interpreter', 'none');
                   
                   switch(new_respmap(z,hh)) % change figure color based on resp and Neuron type
                       case 0
                           set(gcf,'Color','[1 1 1 ]')
                       case 1
                           set(gcf,'Color','[0 0.4470 0.7410]')
                       case 2
                           set(gcf,'Color','[0 0 1 ]')
                       case 3
                           set(gcf,'Color','[0.8500 0.3250 0.0980]')
                       case 4
                           set(gcf,'Color','[1 0 0]')
                   end
                   
                   
                
                frame = getframe(gcf);
                Icollect{z,hh} = frame.cdata;
                close(gcf);
                
                
                
               end
            
            end
            
            %column_shuliang = size(Icollect,2);
            
            wfIcollect = {};
            for g = 1: length(new_rowids)
                s.anas{new_rowids(g)}.neurons{s.anas{new_rowids(g)}.song_id}.draw_waveform; % draw waveform plot
                %s.anas{new_rowids(g)}.draw_waveform; % draw waveform plot
                frame = getframe(gcf);
                wfIcollect{g,1} = frame.cdata;
                close(gcf)
            end
            Icollect = horzcat(Icollect,wfIcollect);
            
            
             % draw neuron ids column
              nameIcollect = {};
              for g = 1: length(new_rowids)
                  figure('menubar','none','Color','w') ;
                  ah = gca ;
                  th = text(1,1,s.anas{new_rowids(g)}.unique_neuronname,'Interpreter','none','FontSize',57);
                  set(ah,'visible','off','xlim',[0 2],'ylim',[0 2],'Position',[0 0 1 1]) ;
                  set(th,'visible','on','HorizontalAlignment','center','VerticalAlignment','middle');
                  frame = getframe(gcf);
                  nameIcollect{g,1} = frame.cdata;
                  close(gcf)
              end
              
              finalIcollect = horzcat(nameIcollect,Icollect);
              finalI = cell2mat(finalIcollect);
              
              imwrite(finalI,sprintf('AllNeurons_ColoredThree.png'));

        end
        
        function arrangeThreePlotsByRespMapAndDraw_NoWaveform(s)
            
            % extract binarized neuron's responses to CONs
            dbstop if error
            con_info = struct;
            counts = 0;
            common_cons = {};
            common_spes = {};
            
            for k = 1: length(s.anas)
                this_neuron = s.anas{k}.neurons{s.anas{k}.song_id};
                Conlist = Analysis(this_neuron).evaluateConResponse;
                if length(Conlist) >= 18 % hard code here
                    counts = counts + 1;
                    con_match = regexp({Conlist.stimuliname}.','[BGYRO]\d{3}','match');
                    con_ids = find(~cellfun(@isempty,con_match));
                    spe_match = regexp({Conlist.stimuliname}.','TUT|BOS|Fcall|Mcall|WNS','match');
                    spe_ids = find(~cellfun(@isempty,spe_match));
                    spe_match = [spe_match{spe_ids}].';
                    % remove TUT and BOS
                    con_ids = setdiff(con_ids,spe_ids);
                    con_match = [con_match{con_ids}].';
                    
                    % extract con_label
                    con_info(counts).wav_len = this_neuron.calMeanWaveLength;
                    con_info(counts).neuronname = this_neuron.unique_neuronname;
                    con_info(counts).con_match = con_match;
                    con_info(counts).con_resp = [Conlist(con_ids).label].'; % resp is just label
                    con_info(counts).spe_match = spe_match;
                    con_info(counts).spe_resp = [Conlist(spe_ids).label].';
                    % find out Common Conspecific songs
                    if counts == 1
                        common_cons = con_info(counts).con_match;
                        %common_spes = con_info(counts).spe_match;
                    else
                        temp = intersect(common_cons,con_info(counts).con_match);
                        
                        common_cons = temp; % weird bugs happened when temp was not used                               
%                         temp = intersect(common_spes,con_info(counts).spe_match);
%                         common_spes = temp;
                        
                        if length(common_cons)< 18
                            pause
                        end
                    end 
                    
                end
            end
                 
            % extract binary con-resp map
            respmap = [];
            for m = 1: length(con_info)
                [~,loc] = ismember (common_cons,con_info(m).con_match);
                
               % new_con_info(m).name = {con_info(m).con_match{loc}}.';
                respmap(m,:) = con_info(m).con_resp(loc);
            end
            
            
            % extract spe-resp map
            common_spes = {'TUT','BOS','Fcall','Mcall','WNS'};
            spemap = [];
            for m = 1: length(con_info)
                
                if ~isempty(con_info(m).spe_match)
                    [~,loc] = ismember (common_spes,con_info(m).spe_match);
                else
                    [~,loc] = ismember (common_spes,{});
                end
                
                % new_con_info(m).name = {con_info(m).con_match{loc}}.';
                for mm = 1: length(loc)
                    if loc(mm) == 0
                        spemap(m,mm) = nan;
                    else
                        spemap(m,mm) = con_info(m).spe_resp(loc(mm));
                    end
                end
            end
            
            % Sort the matrix
            [~,new_rowids] = sort(sum(respmap,2),'descend');
            [~,new_columnids] = sort(sum(respmap,1),2,'descend');
            tempo =  respmap(new_rowids,:);
            new_respmap = tempo(:,new_columnids);
            %frozen_new_respmap = new_respmap;
            new_cnames = {common_cons{new_columnids}}; %#ok<CCAT1>
            tempo = {con_info.neuronname};
            new_rnames = {tempo{new_rowids}};

            respT = array2table(respmap);
            respT.Properties.VariableNames = new_cnames;
            respT.Properties.RowNames = new_rnames;
             
            % sort spe-map
            new_spemap = spemap(new_rowids,:);
            
            % reordered waveforms length
            
            WLs = [con_info.wav_len].';
            WLs = WLs(new_rowids);
            
            wlfr_info = plotWavelengthVsFiringRate(s);
            
            names_in_wlfr = {wlfr_info.neuronname}.';
            [~,ids_for_reorder] = ismember(new_rnames,names_in_wlfr);
            new_wlfr_info = wlfr_info(ids_for_reorder);
            
            for k = 1: size(new_respmap,1)
                value_is_one_ids = find(new_respmap(k,:));
                if new_wlfr_info(k).isBS == 0
                    new_respmap(k,value_is_one_ids) = 1;
                elseif new_wlfr_info(k).isBS == 1
                    new_respmap(k,value_is_one_ids) = 3;
                end  
            end
            
            frags_exist_names = s.markNeuronsWithFragsAsStimuli;
            
            for k = 1: size(new_respmap,1)
                value_is_one_ids = find(new_respmap(k,:));
                if ismember(new_wlfr_info(k).neuronname,frags_exist_names) && new_wlfr_info(k).isBS == 0
                    new_respmap(k,value_is_one_ids) = 2;
                elseif ismember(new_wlfr_info(k).neuronname,frags_exist_names) && new_wlfr_info(k).isBS == 1
                    new_respmap(k,value_is_one_ids) = 4;
                end  
            end
               
            % draw the three plots
            Icollect = {}; % to collect figure frames for each pairwise three plots
            
            for z = 1: length(new_rowids)
               new_N_id =  new_rowids(z);  % reorder the neurons ranked by numnber of eliciting songs
               
               [e_songs,slist_songs] = s.anas{new_N_id}.neurons{s.anas{new_N_id}.song_id}.onlyExtractEphysOfSongs;
               converted_names = cellfun(@convert.bid,{slist_songs.name}.','UniformOutput',0);
               [~,order_for_three_plots] = ismember(new_cnames, converted_names);
               
               
               for hh = 1: length(order_for_three_plots)
                   e_songs{order_for_three_plots(hh)}.three;
                   subplot(3,1,3);
                   xlabel(sprintf('Neuron:%s---Stimuli:%s',s.anas{new_N_id}.unique_neuronname,...
                       convert.bid(e_songs{order_for_three_plots(hh)}.sound.name)),'Interpreter', 'none');
                   
                   switch(new_respmap(z,hh)) % change figure color based on resp and Neuron type
                       case 0
                           set(gcf,'Color','[1 1 1 ]')
                       case 1
                           set(gcf,'Color','[0 0.4470 0.7410]')
                       case 2
                           set(gcf,'Color','[0 0 1 ]')
                       case 3
                           set(gcf,'Color','[0.8500 0.3250 0.0980]')
                       case 4
                           set(gcf,'Color','[1 0 0]')
                   end  
                frame = getframe(gcf);
                Icollect{z,hh} = frame.cdata;
                close(gcf);   
               end
            
            end
            
            
             % draw neuron ids column
              nameIcollect = {};
              for g = 1: length(new_rowids)
                  figure('menubar','none','Color','w') ;
                  ah = gca ;
                  th = text(1,1,s.anas{new_rowids(g)}.unique_neuronname,'Interpreter','none','FontSize',57);
                  set(ah,'visible','off','xlim',[0 2],'ylim',[0 2],'Position',[0 0 1 1]) ;
                  set(th,'visible','on','HorizontalAlignment','center','VerticalAlignment','middle');
                  frame = getframe(gcf);
                  nameIcollect{g,1} = frame.cdata;
                  close(gcf)
              end
              
              finalIcollect = horzcat(nameIcollect,Icollect);
              finalI = cell2mat(finalIcollect);
              
              imwrite(finalI,sprintf('NoWaveform_AllNeurons_ColoredThree.png'));

        end
        
        
        function drawSongOnlyWaveformVsAllWaveform(s)
            
            dbstop if error
            con_info = struct;
            counts = 0;
            common_cons = {};
            common_spes = {};
            
            for k = 1: length(s.anas)
                this_neuron = s.anas{k}.neurons{s.anas{k}.song_id};
                Conlist = Analysis(this_neuron).evaluateConResponse;
                if length(Conlist) >= 18 % hard code here
                    counts = counts + 1;
                    con_match = regexp({Conlist.stimuliname}.','[BGYRO]\d{3}','match');
                    con_ids = find(~cellfun(@isempty,con_match));
                    spe_match = regexp({Conlist.stimuliname}.','TUT|BOS|Fcall|Mcall|WNS','match');
                    spe_ids = find(~cellfun(@isempty,spe_match));
                    spe_match = [spe_match{spe_ids}].';
                    % remove TUT and BOS
                    con_ids = setdiff(con_ids,spe_ids);
                    con_match = [con_match{con_ids}].';
                    % add TUT-BOS-Fcall-Mcall-WNS (Special)
                    % extract con_label
                    con_info(counts).wav_len = this_neuron.calMeanWaveLength;
                    con_info(counts).neuronname = this_neuron.unique_neuronname;
                    con_info(counts).con_match = con_match;
                    con_info(counts).con_resp = [Conlist(con_ids).label].';
                    con_info(counts).spe_match = spe_match;
                    con_info(counts).spe_resp = [Conlist(spe_ids).label].';
                    % find out Common Conspecific songs
                    if counts == 1
                        common_cons = con_info(counts).con_match;
                    else
                        temp = intersect(common_cons,con_info(counts).con_match);
                        common_cons = temp; % weird bugs happened when temp was not used                               
                        if length(common_cons)< 18
                            pause
                        end
                    end 
                    
                end
            end
            % extract binary con-resp map
            respmap = [];
            for m = 1: length(con_info)
                [~,loc] = ismember (common_cons,con_info(m).con_match);
                
               % new_con_info(m).name = {con_info(m).con_match{loc}}.';
                respmap(m,:) = con_info(m).con_resp(loc);
            end
            % Sort the matrix
            [~,new_rowids] = sort(sum(respmap,2),'descend');
            
            
            Icollect = {};
            
            column_shuliang = size(Icollect,2);
            
            % first column song only waveform
            
            for g = 1: length(new_rowids)
                s.anas{new_rowids(g)}.neurons{s.anas{new_rowids(g)}.song_id}.draw_waveform; % draw waveform plot
                frame = getframe(gcf);
                Icollect{g, column_shuliang + 1} = frame.cdata;
                close(gcf)
            end
            
             column_shuliang = size(Icollect,2);
             
            % second column all waveforms
              for g = 1: length(new_rowids)
                s.anas{new_rowids(g)}.draw_waveform; % draw waveform plot
                frame = getframe(gcf);
                Icollect{g, column_shuliang + 1} = frame.cdata;
                close(gcf)
              end
              
              
              % neuron ids
              nameIcollect = {};
              parfor g = 1: length(new_rowids)
                  figure('menubar','none','Color','w') ;
                  ah = gca ;
                  th = text(1,1,s.anas{new_rowids(g)}.unique_neuronname,'Interpreter','none','FontSize',57);
                  set(ah,'visible','off','xlim',[0 2],'ylim',[0 2],'Position',[0 0 1 1]) ;
                  set(th,'visible','on','HorizontalAlignment','center','VerticalAlignment','middle');
                  frame = getframe(gcf);
                  nameIcollect{g,1} = frame.cdata;
                  close(gcf)
                  
              end
              
              Icollect = horzcat(nameIcollect,Icollect);
            
              
          finalI = cell2mat(Icollect);
           
          
          imwrite(finalI,sprintf('SongOnlyWaveforms_Vs_AllWaveforms.png'));
          
          
          
            
        end
        
        function info = get_NumResp_WL_FR_Info(s)
            
             dbstop if error
            con_info = struct;
            counts = 0;
            common_cons = {};
            common_spes = {};
            
            for k = 1: length(s.anas)
                this_neuron = s.anas{k}.neurons{s.anas{k}.song_id};
                Conlist = Analysis(this_neuron).evaluateConResponse;
                if length(Conlist) >= 18 % hard code here
                    counts = counts + 1;
                    con_match = regexp({Conlist.stimuliname}.','[BGYRO]\d{3}','match');
                    con_ids = find(~cellfun(@isempty,con_match));
                    spe_match = regexp({Conlist.stimuliname}.','TUT|BOS|Fcall|Mcall|WNS','match');
                    spe_ids = find(~cellfun(@isempty,spe_match));
                    spe_match = [spe_match{spe_ids}].';
                    % remove TUT and BOS
                    con_ids = setdiff(con_ids,spe_ids);
                    con_match = [con_match{con_ids}].';
                    
                   
                    % add TUT-BOS-Fcall-Mcall-WNS (Special)
                    
                    % extract con_label
                    con_info(counts).wav_len = this_neuron.calMeanWaveLength;
                    con_info(counts).neuronname = this_neuron.unique_neuronname;
                    con_info(counts).con_match = con_match;
                    con_info(counts).con_resp = [Conlist(con_ids).label].';
                    con_info(counts).spe_match = spe_match;
                    con_info(counts).spe_resp = [Conlist(spe_ids).label].';
                    % find out Common Conspecific songs
                    if counts == 1
                        common_cons = con_info(counts).con_match;
                        %common_spes = con_info(counts).spe_match;
                    else
                        temp = intersect(common_cons,con_info(counts).con_match);
                        
                        common_cons = temp; % weird bugs happened when temp was not used                               
%                         temp = intersect(common_spes,con_info(counts).spe_match);
%                         common_spes = temp;
                        
                        if length(common_cons)< 18
                            pause
                        end
                    end 
                    
                end
            end
            
                 
            % extract binary con-resp map
            respmap = [];
            for m = 1: length(con_info)
                [~,loc] = ismember (common_cons,con_info(m).con_match);
                
               % new_con_info(m).name = {con_info(m).con_match{loc}}.';
                respmap(m,:) = con_info(m).con_resp(loc);
            end
            
            % Sort the matrix
            [~,new_rowids] = sort(sum(respmap,2),'descend');
            [~,new_columnids] = sort(sum(respmap,1),2,'descend');
            tempo =  respmap(new_rowids,:);
            new_respmap = tempo(:,new_columnids);
            new_cnames = {common_cons{new_columnids}}; %#ok<CCAT1>
            tempo = {con_info.neuronname};
            new_rnames = {tempo{new_rowids}};

            respT = array2table(respmap);
            respT.Properties.VariableNames = new_cnames;
            respT.Properties.RowNames = new_rnames;
             
           
            
            % reordered waveforms length
            
            WLs = [con_info.wav_len].';
            WLs = WLs(new_rowids);
            
            wlfr_info = plotWavelengthVsFiringRate(s);
            
            names_in_wlfr = {wlfr_info.neuronname}.';
            [~,ids_for_reorder] = ismember(new_rnames,names_in_wlfr);
            new_wlfr_info = wlfr_info(ids_for_reorder);
            
            for k = 1: size(new_respmap,1)
                value_is_one_ids = find(new_respmap(k,:));
                if new_wlfr_info(k).isBS == 0
                    new_respmap(k,value_is_one_ids) = 1;
                elseif new_wlfr_info(k).isBS == 1
                    new_respmap(k,value_is_one_ids) = 3;
                end  
            end
            
            frags_exist_names = s.markNeuronsWithFragsAsStimuli;
            
            for k = 1: size(new_respmap,1)
                value_is_one_ids = find(new_respmap(k,:));
                if ismember(new_wlfr_info(k).neuronname,frags_exist_names) && new_wlfr_info(k).isBS == 0
                    new_respmap(k,value_is_one_ids) = 2;
                elseif ismember(new_wlfr_info(k).neuronname,frags_exist_names) && new_wlfr_info(k).isBS == 1
                    new_respmap(k,value_is_one_ids) = 4;
                end  
            end
            
            labled_tokeep = [];
            for k = 1: size(new_respmap,1)
                if ismember(4,new_respmap(k,:))
                    labled_tokeep = [labled_tokeep,k];
                end
            end
            
            
            % add number of responsive songs into info
            for k = 1: length(new_wlfr_info)
                new_wlfr_info(k).numResp = length(find(new_respmap(k,:)));
                if ismember(new_wlfr_info(k).neuronname,frags_exist_names)
                    new_wlfr_info(k).fragexist = 1;
                else
                    new_wlfr_info(k).fragexist = 0;
                end
            end
        
            info = new_wlfr_info;
   
        end
        
        function WLvsNumberOfRespnsiveSongs(s)
            
            info = s.get_NumResp_WL_FR_Info;
            
            BSinfo = info(find([info.isBS]==1));
            
            figure('Color','w');
            subplot(1,3,1)
            scatter([BSinfo.numResp],[BSinfo.wl],[],[0.8500 0.3250 0.0980],'filled');
            xL=xlim;
            yL=ylim;
            ylim([yL(1),yL(2)]);
            [cc,p] = corrcoef([BSinfo.numResp],[BSinfo.wl]);
            text(0.99*xL(2),0.99*yL(2),sprintf('corr-coef: %.3f \n p-value: %.3f',cc(1,2),p(1,2)),'HorizontalAlignment','right','VerticalAlignment','top');
            xlabel('Number of response-eliciting songs');
            ylabel('Mean wavelength');
            
            subplot(1,3,2)
            scatter([BSinfo.numResp],[BSinfo.mean_plt_fr],[],[0.8500 0.3250 0.0980],'filled');
            xL=xlim;
            yL=ylim;
            [cc,p] = corrcoef([BSinfo.numResp],[BSinfo.mean_plt_fr]);
            text(0.99*xL(2),0.99*yL(2),sprintf('corr-coef: %.3f \n p-value: %.3f',cc(1,2),p(1,2)),'HorizontalAlignment','right','VerticalAlignment','top');
            xlabel('Number of response-eliciting songs');
            ylabel('Evoked firing rate');
            
            subplot(1,3,3)
            scatter([BSinfo.numResp],[BSinfo.mean_pre_fr],[],[0.8500 0.3250 0.0980],'filled');
            xL=xlim;
            yL=ylim;
            [cc,p] = corrcoef([BSinfo.numResp],[BSinfo.mean_pre_fr]);
            text(0.99*xL(2),0.99*yL(2),sprintf('corr-coef: %.3f \n p-value: %.3f',cc(1,2),p(1,2)),'HorizontalAlignment','right','VerticalAlignment','top');
            xlabel('Number of response-eliciting songs');
            ylabel('Spontaneous firing rate');
            
        end
            
    end
    
    methods(Static)
        
        
        function runAllAnalysis(ana_path)
            dbstop if error

            ANAfiles = extract.filename(ana_path,'*.mat');
            
            wuhu = waitbar(0,'Start processing');
            
            for m = 1 : length(ANAfiles)
                
                load(ANAfiles{m});
                A.set_eleinf("C:\Users\Zhehao\Dropbox (OIST)\SaveAllMatXlsData\My_eleinf\all_eleinf.mat");
                A.V1drawMeanFeaturesInSongVsRespAsLineChart;
                A.drawDTWSimilarityMatrixBasedOnZscoredData;
                A.drawDTWSimilarityMatrix;
                A.drawCoeffOfFeaturesLinearFit;
                A.drawMeanFeatureVsResp;
                A.drawMeanFeaturesVsRespAsLineChart;
                A.drawPairwiseFragmentsMeanFeaturesDistribution;
                A.threePlotsWithPitch;
                A.V2drawMeanFeaturesInSongVsRespAsLineChart
                waitbar(m/length(ANAfiles),wuhu,sprintf('%u of totally %u files',m,length(ANAfiles)));
            end
            
            close(wuhu);
        end
        
        function writeAnalysisObjectForMergedPlexonFiles(tablepath)

            %tablepath = "C:\Users\Zhehao\Desktop\AllInOne (27).xlsx";
            
            
            T = table2struct(readtable(tablepath));
            
            IDs = [T.UniqueID].';
            num_ids = rmmissing(unique(IDs(IDs~=0)));
            
            %kbad = [28,31,32,41,42,44,64,65];
            
            
            wb = waitbar(0,'Creating Neuron analysis objects');
            for k = 11: length(num_ids)
                waitbar(k/length(num_ids),wb,sprintf('%u of %u Neuron',k,length(num_ids)));
                ids_in_T = find(num_ids(k) == IDs);
                
                Ns = {};
                for i = 1:length(ids_in_T)
                    b = Batch(T(ids_in_T(i)).MergedTxtPath,T(ids_in_T(i)).MergedPlxPath,T(ids_in_T(i)).StimuliPath);
                    
                    this_channel = T(ids_in_T(i)).ChannelName;
                    this_unit = T(ids_in_T(i)).UnitName;
                    neu_list = {b.nlist.neuronname}.';
                    
                    channel_ids = find(~cellfun(@isempty,regexp(neu_list,this_channel)));
                    unit_ids = find(~cellfun(@isempty,regexp(neu_list,sprintf('_%u',this_unit))));
                    neuron_ids = intersect(channel_ids,unit_ids);
                    
                    if T(ids_in_T(i)).MergedIndex == 0 || isempty(T(ids_in_T(i)).MergedIndex ) % When the raw data is not a merged object (Or not treated as)
                        b.select(neuron_ids);
                        N = b.getn{1};         
                    else
                        b.select(neuron_ids);
                        tempN = b.getn(T(ids_in_T(i)).MergedIndex);
                        N = tempN{1};  
                    end
                
                    %N.mergeIdx = T(ids_in_T(i)).mergedIdx;
                    N.signalGoodness = T(ids_in_T(i)).Goodness;
                    N.set_uniqueid(T(ids_in_T(i)).UniqueID);
                    % allocate sap-based feature information to each neuron
                    sorted_data = Neuron.extractFeaturesFromSapRawDataJson(T(ids_in_T(i)).FeatureData, T(ids_in_T(i)).FeatureInfo);
                    N.setEachStimuliSapFeatures(sorted_data);
                    N.calMeanFeatures;
                    Ns{i} = N;
                end
                
                A = Analysis(Ns);
                A.uniqueid = num_ids(k);
                
                
                save(sprintf('%s_%u',A.birdid,A.uniqueid),'A','-v7.3');
                
                
            end
            
            close(wb);




        end
        
    end
    
end

