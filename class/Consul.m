
classdef Consul < handle
    % has been deprecated by Sultan
    % A consul held the highest elected political office of the Roman Republic 
    % So Consul do analysis for the whole data
    
    properties
        neurons % a cell collection of neurons
        neuronlist
        neuroninfo
        con_neurons
        brushed
    end
    
    methods
        function c = Consul(neuronlist)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            f = waitbar(0,'Loading...');
            c.neuronlist = neuronlist;
            for k = 1: length(neuronlist)
                
                load(neuronlist{k});
                c.neurons{k} = N;
                waitbar(k/length(neuronlist),f,sprintf(' Now loading %u // %u files',k,length(neuronlist)));
            end
            
            close(f);
        end
        
        function c = selectConNeurons(c)
            
            selected_ids = [];
            for k = 1: length(c.neurons)
                
                if length(cellfun(@isempty, regexp('norm',{c.neurons{k}.slist.name}.')))> 12
                    selected_ids = [selected_ids,k];
                end
            end
            c.con_neurons = c.neurons{selected_ids};
        end
        
        function drawCONResponse(c)
            % Extract binarized neuron's responses to CONs
            dbstop if error
            con_info = struct;
            counts = 0;
            
            common_cons = {};
            for k = 1: length(c.neurons)
                Conlist = Analysis(c.neurons{k}).evaluateConResponse;
                if length(Conlist) >= 18 % hard code here
                    counts = counts + 1;
                    con_match = regexp({Conlist.stimuliname}.','[BGYRO]\d{3}','match');
                    con_ids = find(~cellfun(@isempty,con_match));
                    spe_match = regexp({Conlist.stimuliname}.','TUT|BOS|Fcall|Mcall|WNS','match');
                    spe_ids = find(~cellfun(@isempty,spe_match));
                    % remove TUT and BOS
                    con_ids = setdiff(con_ids,spe_ids);
                    con_match = [con_match{con_ids}].';
                    % Extract con_label
                    con_info(counts).neuronname = c.neurons{counts}.neuronname;
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
            
                 
            % Extract binary resp map
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
        
        function wl_info = calWavelength(c)
              wl_info = struct;
            
            for k = 1: length(c.neurons)
                 wl_info(k).wl = c.neurons{k}.calMeanWaveLength;
                 wl_info(k).neuronname = c.neurons{k}.neuronname;
            end
            
        end
        
        function histWaveLength(c)
            dbstop if error       
%             figure
%             imagesc(new_respmap);
%             xticks(1:length(new_cnames));
%             xticklabels(new_cnames);
%             yticks(1:length(new_rnames));
%             yticklabels(new_rnames);
%             set(gca,'TickLabelInterpreter','none');
           
            wl_info = struct;
            
            for k = 1: length(c.neurons)
                 wl_info(k).wl = c.neurons{k}.calMeanWaveLength;
                 wl_info(k).neuronname = c.neurons{k}.neuronname;
            end
            
            figure
            histogram([wl_info.wl].',30); % originally 15
            
            xlabel('Spike width (ms)');
            ylabel('Neurons per bin');
            
        end
        
        function drawCONSPEResponse(c)
            
             % Extract binarized neuron's responses to CONs
            dbstop if error
            con_info = struct;
            counts = 0;
            common_cons = {};
            common_spes = {};
            
            for k = 1: length(c.neurons)
                Conlist = Analysis(c.neurons{k}).evaluateConResponse;
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
                    
                    % Extract con_label
                    con_info(counts).wav_len = c.neurons{k}.calMeanWaveLength;
                    con_info(counts).neuronname = c.neurons{k}.neuronname;
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
            
                 
            % Extract binary con-resp map
            respmap = [];
            for m = 1: length(con_info)
                [~,loc] = ismember (common_cons,con_info(m).con_match);
                
               % new_con_info(m).name = {con_info(m).con_match{loc}}.';
                respmap(m,:) = con_info(m).con_resp(loc);
            end
            
            
            % Extract spe-resp map
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
        
        function drawCONSPEResponse_NSBS(c)
            
             % Extract binarized neuron's responses to CONs
            dbstop if error
            con_info = struct;
            counts = 0;
            common_cons = {};
            common_spes = {};
            
            for k = 1: length(c.neurons)
                Conlist = Analysis(c.neurons{k}).evaluateConResponse;
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
                    
                    % Extract con_label
                    con_info(counts).wav_len = c.neurons{k}.calMeanWaveLength;
                    con_info(counts).neuronname = c.neurons{k}.neuronname;
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
            
                 
            % Extract binary con-resp map
            respmap = [];
            for m = 1: length(con_info)
                [~,loc] = ismember (common_cons,con_info(m).con_match);
                
               % new_con_info(m).name = {con_info(m).con_match{loc}}.';
                respmap(m,:) = con_info(m).con_resp(loc);
            end
            
            
            % Extract spe-resp map
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
            
            wlfr_info = plotWavelengthVsFiringRate(c);
            
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
        
        function neuronStimuliInfo = getStimuliInfo(c)
            
            neuronStimuliInfo = struct;
            
            for k = 1: length(c.neurons)
                these_names = {c.neurons{k}.slist.name}.';
                
                % initialize
                neuronStimuliInfo(k).neuronname = c.neurons{k}.neuronname;
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
        
        function fr_info = multiRepeatsFiringRate(c)
            fr_info = struct;
             
            for k = 1: length(c.neurons)
                % for each neuron
               
                thisn = c.neurons{k};
                fr_info(k).neuronname = c.neurons{k}.neuronname;
                
                sum_prelen = 0; % summed prey length 
                concat_presptimes = []; % concatenated prey sptimes
                
                sum_pltlen = 0; %summed prey( stimuli y, not plty or rawy) length 
                concat_pltsptimes = []; %  % concatenated y sptimes
                
                for m = 1: length(thisn.e)
                    
                    
                    % for prey 
                    fr_info(k).presptimes{m} = thisn.e{m}.presptimes
                    fr_info(k).preylen{m} = length(thisn.e{m}.y)/thisn.e{m}.fs;
                    fr_info(k).repnum{m} = size(thisn.e{m}.presptimes,2);
                    temp = thisn.e{m}.presptimes.';
                    concat_presptimes = [concat_presptimes;vertcat(vertcat(temp{:}))+ sum_prelen];
                    sum_prelen = sum_prelen +  fr_info(k).preylen{m};
                    
                    % for plty
                    fr_info(k).pltsptimes{m} = thisn.e{m}.pltsptimes
                    fr_info(k).pltlen{m} = length(thisn.e{m}.plty)/thisn.e{m}.fs;
                    temp = thisn.e{m}.pltsptimes.';
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
        
        function wlfr_info = plotWavelengthVsFiringRate(c)
            
            wl_info = calWavelength(c);
            fr_info = multiRepeatsFiringRate(c);

            ids_in_fr = [];
            for k = 1: length(wl_info)
                this_nname = wl_info(k).neuronname;
                ids_in_fr(k) = find(strcmp(this_nname,{fr_info.neuronname}.'));
            end
            
            fr_info = fr_info(ids_in_fr); % re-order
            
            fr_info = rmfield( fr_info,'neuronname');

            wlfr_info = table2struct([struct2table(wl_info),struct2table(fr_info)]);

            fig = figure;
            h = scatter([wlfr_info.wl],[wlfr_info.mean_pre_fr]);
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
            scatter([wlfr_info(brushed_ids).wl],[wlfr_info(brushed_ids).mean_pre_fr],[], [0.8500 0.3250 0.0980],'filled');
            scatter([wlfr_info(not_brushed_ids).wl],[wlfr_info(not_brushed_ids).mean_pre_fr],[],[0 0.4470 0.7410],'filled');
            xlabel('Mean Spike width');
            ylabel('Spontaneous Firing Rate');
             
            %disp(brushed);

%             figure;
%             scatter([wlfr_info.wl],[wlfr_info.mean_plt_fr]);
%             xlabel('Mean Spike width');
%             ylabel('Evoked Firing Rate');


    
            
            
        end
        
        
        
    end
    


end

