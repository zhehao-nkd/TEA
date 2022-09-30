classdef Deg
    
    properties
        deglist
        list
    end
    methods
        function dg = Deg(list)
            dg.list = list;
            dg.deglist = list( find(~cellfun(@isempty, regexp(cellstr({dg.list.stimuliname}.'),'Deg|deg'))));
        end
        
        
          function Iall = saveDrawAlignedConsDegs(dg,songnames)
            dbstop if error
            
            %只有一个模式： 只针对二次播放里包含的norm songs进行degs的对齐
            tic
            
            degids = find(~cellfun(@isempty, regexp(cellstr({dg.list.stimuliname}.'),'deg|Deg')));
            
            deglist = dg.list(degids);
            deg_Fid = unique({deglist.Fid}.');
            % normlist = Neuron(neu.neurons{neu.song_id}).normlist;
            
            
            
            subfile_deg_ids = find(~cellfun(@isempty,regexp(cellstr({dg.list.Fid}.'),strjoin(deg_Fid,'|'))));
            hard_to_name_ids = subfile_deg_ids;
            if exist('songnames','var')
                songnameids = find(~cellfun(@isempty, regexp(cellstr({dg.list.stimuliname}.'),strjoin(songnames,'|'))));
                hard_to_name_ids = intersect(subfile_deg_ids,songnameids);
            end
            
            fucklist = dg.list(hard_to_name_ids);
            
            normlist = fucklist(find(~cellfun(@isempty, regexp(cellstr({fucklist.stimuliname}.'),'norm'))));
            %             [~,postunique] = unique(cellfun(@Convert.bid,cellstr({fucklist.stimuliname}.'),'Uni',0));
            %             normlist = normlist(postunique);
            
            % About Deg
            degids = find(~cellfun(@isempty, regexp(cellstr({dg.list.stimuliname}.'),'Deg|deg') ));
            if ~isempty(degids)
                
                deglist = dg.list(degids);
                for m = 1: length(deglist)
                    birdid = Convert.bid(deglist(m).stimuliname);
                    ids_norm = find(~cellfun(@isempty, regexp(cellstr({normlist.stimuliname}.'),birdid) ) );
                    if ~isempty(ids_norm)& length(ids_norm) == 1
                        [deglist(m).sylIni,trump_diffvalue] = Neuron.findIni(normlist(ids_norm).plty,deglist(m).y);
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
                    birdid = Convert.bid(normlist(w).stimuliname);
                    ids_indeg = find(~cellfun(@isempty, regexp(cellstr({deglist.stimuliname}.'),birdid) ) );
                    selected_deglist = deglist(ids_indeg);
                    [~,temp_index] = sortrows([selected_deglist.sylIni].');
                    selected_deglist = selected_deglist(temp_index);
                end
                
                
                % draw the basic figure
                Icollect = {};
                figure('Color','w','Position',PM.size_wide);
                
                Draw.two(normlist(w).plty,normlist(w).fs,normlist(w).pltsptimes);
                xlabel(normlist(w).stimuliname);
                frame = getframe(gcf);
                Icollect{1} = frame.cdata;
                close(gcf)
                
                
                
                if ~isempty(degids)
                    for hh = 1: length(selected_deglist)
                        figure('Color','w','Position',PM.size_wide);
                        Draw.two(selected_deglist(hh).pady,selected_deglist(hh).fs,selected_deglist(hh).padsptimes);
                        xlabel(selected_deglist(hh).stimuliname);
                        frame = getframe(gcf);
                        Icollect{1 + hh} = frame.cdata;
                        close(gcf);
                    end
                end
                
                frozen_Icollect_len = length(Icollect);
                
                
                I_of_each_column{w} = vertcat(Icollect{:});
            end
            
            %             neu.drawFirstWaveform;
            %             temp = getframe(gcf);
            %             w_img = temp.cdata;
            %             I_of_each_column{length(I_of_each_column)+ 1} = w_img;
            
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
            
            % imwrite(Iall,sprintf('Aligned_ConsDegs_%s.png',neu.neurons{1}.neuronname));
            imwrite(Iall,sprintf('Aligned_ConsDegs_%s.png',dg.formated_imagename));
            toc
            
        end
       
    
    

    
    end
    
end