classdef simultaneousNeurons < Sultan
    % inherit Sultan, process simutaneously recorded neuron files
  
    
    properties
        simugroups
        groups
    end
    
    methods
        
        function s = simultaneousNeurons(dirs_of_analysis)
            
            dbstop if error
            s@Sultan(dirs_of_analysis);
            
            grouping_info = {};
            for k = 1: length(s.anas)
                grouping_info{k} = regexp(s.anas{k},'[OGBRY]\d{3}_[ZP]\d{2}','match');   
            end
            
            unqs = unique(cellstr(grouping_info.'));
            
            s.groups = {};
            for k = 1:length(unqs)
                
                ids = find(~cellfun(@isempty, regexp(cellstr(grouping_info),unqs{k})));
                
                s.groups{k} = s.anas(ids);
            end
            
        end
        
        function  How_Do_NCM_Neurons_respond_to_Songs(s)
            
            dbstop if error
            simu_neuron_larger_than1_ids = find (cellfun(@(x) length(x), s.groups) > 1);
            
            wb = waitbar(0,'Start processing');
            num_files = length(simu_neuron_larger_than1_ids );
            Utl.UpdateParforWaitbar(num_files, wb);
            D = parallel.pool.DataQueue;
            afterEach(D, @ult.UpdateParforWaitbar);   
            
            for k = 1:length(simu_neuron_larger_than1_ids)
                simultaneousNeurons.How_Do_simu_Neurons_respond_to_Songs(s.groups{simu_neuron_larger_than1_ids(k)});
                
                send(D, 1);
            end
            
            % Extract binarized neuron's responses to CONs
            
        end
        
    end
    
    methods (Static)
        
         function conallneurons = How_Do_simu_Neurons_respond_to_Songs(ana_pathes)
            
            % Extract binarized neuron's responses to CONs
            dbstop if error
            
            conkeywords = {'B346','B512','B521','B554','B606','G429','G506','G518','G548','G573',...
                'G578','O331','O507','O509','O540','Y515','Y606','Y616'};
            
            spekeywords = {'BOS','TUT','Fcall','Mcall','Het','WNS'};
            
            wb = waitbar(0,'Start processing');
            num_files = length(ana_pathes);
            % Dummy call to nUpdateWaitbar to initialise
            Utl.UpdateParforWaitbar(num_files, wb);
            % Go back to simply calling nUpdateWaitbar with the data
            D = parallel.pool.DataQueue;
            afterEach(D, @ Utl.UpdateParforWaitbar);
            
            conallneurons = struct;
            parfor k = 1: length(ana_pathes) % should be par-for
                loaded = load(ana_pathes{k});
                A = loaded.A;
                A.judgeConResp_FR;
                %A.judgeConResp; % update con resp labels
                FRINFO  = A.multiRepeatsFiringRate;
                conallneurons(k).wl =  A.neurons{A.song_id}.calMeanWaveLength;
                conallneurons(k).neuronname = A.formated_name;
                conallneurons(k).mean_used_fr = FRINFO.mean_pre_fr;
                
                for kk = 1: length(conkeywords) % for Cons
                    
                    thisid = find(~cellfun(@isempty,regexp(cellstr({A.list.stimuliname}.'),['norm\S+(?!(TUT|BOS))',conkeywords{kk},'(?!(TUT|BOS))'] )));
                    if length(thisid) == 1 % 当stimuli sets 有且只有一个对应song时
                        conallneurons(k).figcon{kk} = A.list(thisid).image;
                        conallneurons(k).con_biresp(kk) = A.list(thisid).label;
                    elseif isempty(thisid) % 没有对应song时
                        conallneurons(k).figcon{kk} = uint8(255*ones(size(A.list(1).image)));
                        conallneurons(k).con_biresp(kk) = 0;
                    else
                        conallneurons(k).figcon{kk} = A.list(thisid(1)).image; % 有多个对应时
                        conallneurons(k).con_biresp(kk) = A.list(thisid(1)).label;
                    end
                end
                
                for kk = 1:length(spekeywords) %for Spes
                    
                    thisid = find(~cellfun(@isempty,regexp(cellstr({A.list.stimuliname}.'),['norm\S+',spekeywords{kk}] )));
                    if length(thisid) == 1
                        conallneurons(k).figspe{kk} = A.list(thisid).image;
                        conallneurons(k).spe_biresp(kk) = A.list(thisid).label;
                    elseif isempty(thisid)
                        conallneurons(k).figspe{kk} = uint8(255*ones(size(A.list(1).image)));
                        conallneurons(k).spe_biresp(kk) = 0;
                    else
                        conallneurons(k).figspe{kk} = A.list(thisid(1)).image;
                        conallneurons(k).spe_biresp(kk) = A.list(thisid(1)).label;
                    end
                end
                
                figure('Position',[2108 544 895 672],'color','w','Visible','off')
                A.neurons{A.song_id}.draw_waveform; % draw waveform plot
                frame = getframe(gcf);close(gcf);
                conallneurons(k).figwf = frame.cdata; % waveform figure
                %waitbar(k/length(s.anas),wb,sprintf('%u of %u files',k,length(s.anas)));
                
                % Note we send only an "increment" for the waitbar.
                send(D, 1);
                
            end
            
            close(wb);
            
            con_respmap = cell2mat({conallneurons(:).con_biresp}.'); % old
            
            % Sort the matrix
            [~,new_rowids] = sort(sum(con_respmap,2),'descend');
            [~,new_columnids] = sort(sum(con_respmap,1),'descend');
            
            % edit（reorder） the struct based on matrix  sorting
            conallneurons = conallneurons(new_rowids);
            for k = 1: length(conallneurons)
                conallneurons(k).figcon = conallneurons(k).figcon(new_columnids);
                conallneurons(k).con_biresp = conallneurons(k).con_biresp(new_columnids);
            end
            
            new_con_respmap = cell2mat({conallneurons.con_biresp}.'); % new
            new_spe_respmap = cell2mat({conallneurons.spe_biresp}.'); % new
            new_cnames = {conkeywords{new_columnids}};
%             wlfr_info = Sultan.plotWLvsFR(conallneurons); % To judge whether it is BS or NS
%             
%             for k = 1: size(new_con_respmap,1)
%                 is1 = find(new_con_respmap(k,:));
%                 if wlfr_info(k).isBS == 0
%                     new_con_respmap(k,is1) = 1;
%                 elseif wlfr_info(k).isBS == 1
%                     new_con_respmap(k,is1) = 2;
%                 end
%             end
%             
%             for k = 1: size(new_spe_respmap,1)
%                 is1 = find(new_spe_respmap(k,:));
%                 if wlfr_info(k).isBS == 0
%                     new_spe_respmap(k,is1) = 1;
%                 elseif wlfr_info(k).isBS == 1
%                     new_spe_respmap(k,is1) = 2;
%                 end
%             end
            
%             % draw binary response-map
%             new_rnames = {conallneurons.neuronname}.';
%             %concat_respmap = horzcat(new_con_respmap,new_spe_respmap);
%             concat_respmap = new_con_respmap;
%             %             ysize = size(concat_respmap,1);
%             %             xsize = size(concat_respmap,2);
%             %             yvector = linspace(1,ysize,ysize);
%             %             xvector = linspace(1,xsize,xsize);
%             pcolormap = flip(concat_respmap);
%             pcolormap = vertcat(zeros(1,size(pcolormap,2)),pcolormap);
%             pcolormap = horzcat(pcolormap,zeros(size(pcolormap,1),1));
%             bfig = figure; s = pcolor(pcolormap);
%             s.EdgeColor = 'k';
%             s.LineWidth = 1.6;
%             %contour(concat_respmap,'LineColor','k');
%             map = [ 1 1 1
%                 0 0.4470 0.7410
%                 0.8500 0.3250 0.0980];
%             
%             colormap(bfig,map);
%             %colorbar;
%             xticks(1:size(concat_respmap,2));
%             xticklabels(horzcat(new_cnames,spekeywords));
%             yticks(1:size(concat_respmap,1));
%             yticklabels(new_rnames);
%             set(gca,'TickLabelInterpreter','none');
%             set(gca,'TickLength',[0.001, 0.001]);
%             xtickangle(45)
%             savefig(bfig,'Binary_Resp_Map.fig');
            
            
            % draw NS only map
            % ns_rows = find([wlfr_info(k).isBS] == 0)
            % draw the three plots
            Icollect = {}; % to collect figure frames for each pairwise three plots
            for k = 1: length(conallneurons)
                for kk = 1: length(conallneurons(k).figcon)
                    Icollect{k,kk} = conallneurons(k).figcon{kk};
                    if new_con_respmap(k,kk) ==1
                        Icollect{k,kk} = Convert.colorEdge(Icollect{k,kk},'r'); %NS neurons-red
                    elseif new_con_respmap(k,kk) ==2
                        Icollect{k,kk} = Convert.colorEdge(Icollect{k,kk},'b');
                    end
                end
            end
            
            specollect = {}; % to collect figure frames for each pairwise three plots
            for k = 1: length(conallneurons)
                for kk = 1: length(conallneurons(k).figspe)
                    specollect{k,kk} = conallneurons(k).figspe{kk};
                    if new_spe_respmap(k,kk) ==1
                        specollect{k,kk} = Convert.colorEdge(specollect{k,kk},'r'); %NS neurons-red
                    elseif new_spe_respmap(k,kk) ==2
                        specollect{k,kk} = Convert.colorEdge(specollect{k,kk},'b');
                    end
                end
            end
            
            wfIcollect = {};
            for g = 1: length(new_rowids)
                wfIcollect{g} = conallneurons(g).figwf;
            end
            Icollect = horzcat(Icollect,specollect,wfIcollect.');
            
            
            % draw neuron ids column
            nameIcollect = {};
            for g = 1: length(conallneurons)
                figure('Position',[2108 544 895 672],'color','w','menubar','none','Visible','off')
                ah = gca ;
                th = text(1,1,conallneurons(g).neuronname,'Interpreter','none','FontSize',51);
                set(ah,'visible','off','xlim',[0 2],'ylim',[0 2],'Position',[0 0 1 1]) ;
                set(th,'visible','on','HorizontalAlignment','center','VerticalAlignment','middle');
                frame = getframe(gcf);
                nameIcollect{g,1} = frame.cdata;
                close(gcf)
            end
            
            finalIcollect = horzcat(Icollect,nameIcollect);
            %finalI = cell2mat(finalIcollect);
            
            numperfig = 50;
            numsections = ceil(length(ana_pathes)/numperfig);
            
            for k = 1: numsections
                
                if k < numsections
                    Icollect_this_section = finalIcollect((k-1)*numperfig+1:numperfig*k,:);
                elseif k == numsections
                    Icollect_this_section = finalIcollect((k-1)*numperfig+1:end,:);
                end
                
                 grouping_info = regexp(ana_pathes{1},'[OGBRY]\d{3}_[ZP]\d{2}','match')
                I_this_section = cell2mat(Icollect_this_section);
                t = Tiff(sprintf('SimuNeurons_%s_Part%u.tiff',grouping_info,k),'w8');
                setTag(t,'ImageLength',size(I_this_section,1));
                setTag(t,'ImageWidth',size(I_this_section,2));
                setTag(t,'Photometric',Tiff.Photometric.RGB);
                setTag(t,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
                setTag(t,'BitsPerSample',8);
                setTag(t,'SamplesPerPixel',3);
                % write data
                write(t,I_this_section);
                close(t);
                
            end
            
            end
        
        
    end
end

