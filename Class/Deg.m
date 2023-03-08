classdef Deg < handle
    
    properties
        deglist % it is separated for each birdid
        list
        corresp_normlist
        formated_name
    end
    
    methods

        function dg = Deg(list)
            dg.list = list;
            degids = find(~cellfun(@isempty, regexp(cellstr({dg.list.stimuliname}.'),'Deg|deg') ));
            list_contain_degs = dg.list(degids);
            deg_birdids = cellfun(@Convert.bid, cellstr({list_contain_degs.stimuliname}.'),'Uni',0);
            unique_deg_birdids = unique(deg_birdids);
            unpack = @(x) x{1};
            summer = {};
            for k = 1:length(unique_deg_birdids)
                correspids = strcmp(unique_deg_birdids{k},deg_birdids); % 找到有这个birdid的所有deg的ids,这ids是关于list_contain_degs的
                dg.deglist{k} = list_contain_degs(correspids);

                corresp_fid = unpack(unique({dg.deglist{k}.Fid}.')); % 找到deg stimuli对应的fid

                corresp_normids = mintersect( find(~cellfun(@isempty, regexp(cellstr({dg.list.stimuliname}.'),unique_deg_birdids{k}))),...
                    find(~cellfun(@isempty, regexp(cellstr({dg.list.stimuliname}.'),'norm'))),...
                    find(~cellfun(@isempty, regexp(cellstr({dg.list.Fid}.'),corresp_fid))) );

%                 if isempty(corresp_normids)
                    temp = intersect( find(~cellfun(@isempty, regexp(cellstr({dg.list.stimuliname}.'),unique_deg_birdids{k}))),...
                    find(~cellfun(@isempty, regexp(cellstr({dg.list.stimuliname}.'),'norm'))));

                    summer{k} = dg.list(temp(1)); % 取第一位


%                 end


                dg.deglist{k} = horzcat(dg.list(corresp_normids),dg.deglist{k});

            end

            dg.corresp_normlist = vertcat(summer{:});

            dg.list = []; % 为了节省存储，在最后一步把此项设为零

        end
        
        
        function Iall = saveDrawAlignedConsDegs(dg)
            dbstop if error

            %只有一个模式： 只针对二次播放里包含的norm songs进行degs的对齐
            tic
            for k = 1:length(dg.deglist)

                sublist = dg.deglist{k};

                degids = find(~cellfun(@isempty, regexp(cellstr({sublist.stimuliname}.'),'deg|Deg')));

                deglist = sublist(degids);
                deg_Fid = unique({deglist.Fid}.');
                % normlist = Neuron(neu.neurons{neu.song_id}).normlist;



                subfile_deg_ids = find(~cellfun(@isempty,regexp(cellstr({sublist.Fid}.'),strjoin(deg_Fid,'|'))));
                hard_to_name_ids = subfile_deg_ids;
%                 if exist('songnames','var')
%                     songnameids = find(~cellfun(@isempty, regexp(cellstr({sublist.stimuliname}.'),strjoin(songnames,'|'))));
%                     hard_to_name_ids = intersect(subfile_deg_ids,songnameids);
%                 end

                fucklist = sublist(hard_to_name_ids);

                normlist = fucklist(find(~cellfun(@isempty, regexp(cellstr({fucklist.stimuliname}.'),'norm'))));

                if isempty(normlist)
                    normlist = dg.corresp_normlist(k);
                end
                %             [~,postunique] = unique(cellfun(@Convert.bid,cellstr({fucklist.stimuliname}.'),'Uni',0));
                %             normlist = normlist(postunique);

                % About Deg
                degids = find(~cellfun(@isempty, regexp(cellstr({sublist.stimuliname}.'),'Deg|deg') ));
                if ~isempty(degids)

                    deglist = sublist(degids);
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
%                         [~,temp_index] = sortrows([selected_deglist.sylIni].');
%                         selected_deglist = selected_deglist(temp_index);
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
                imwrite(Iall,sprintf('对齐的Degs_%s_stimuli是%s.png',dg.formated_name,Convert.bid(normlist.stimuliname)));
                toc
                %degnames{k} = Convert.bid(dg.deglist{k}(1).stimuliname);
            end



        end


        function dg = judgeDegResp(dg)
            dbstop if error
            % evaluate the responsiveness of song degressive deletion

            %ids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'deg'))); % find all norms

            for k = 1:length(dg.deglist) % 当有多个degsong target的时候
                
                locallist = dg.deglist{k};


                ids_norm = find(~cellfun(@isempty, regexp(cellstr({locallist.stimuliname}.'),'norm')));
                % ids_norm 有很多时取第一个吧？？ 2023.01.06 无法推进，暂停在此处

                for m = 1: length(locallist) %setdiff([1: length(locallist)],ids_norm)

                    [locallist(m).sylIni,trump_diffvalue] = Neuron.findIni(locallist(ids_norm).plty,locallist(m).y);%此处目的是找到deg与norm的分歧点，所以具体是哪个E的norm不重要
                    fprintf('%s has the DIFFVALUE as: %u and INI as: %f \n ',locallist(m).stimuliname,trump_diffvalue,locallist(m).sylIni/32000);
                    locallist(m).pady = [zeros([locallist(m).sylIni-1-locallist(m).fs*locallist(m).pltext,1]);locallist(m).plty;zeros(length(locallist(ids_norm).plty)...
                        - (locallist(m).sylIni-1-locallist(m).fs*locallist(m).pltext) - length(locallist(m).plty),1)];
                    locallist(m).padsptimes = cellfun( @(x) x + (locallist(m).sylIni-1-locallist(m).fs*locallist(m).pltext)/locallist(m).fs, locallist(m).pltsptimes,'uni',0);
                    padsdf = Cal.sdf(locallist(m).padsptimes,locallist(m).pady,locallist(m).fs,0.001,0.02); % 0.001,0.004
                    %                     padsdfcollect{m} = padsdf;

                    [locallist(m).pks,locallist(m).locs] = findpeaks(padsdf,'MinPeakHeight',10);
                    if isempty(locallist(m).pks)
                        locallist(m).pks = nan;
                        locallist(m).locs = nan;

                    end


                end

                for v = 1:length(locallist(1).locs) % v means value
                    %                     find(rmmissing([locallist.locs].') -locallist(1).locs(v) <=3)
                    consistent_ids = find( arrayfun(@(x) min(abs(x - locallist(1).locs(v))), [locallist.locs].')<=3);
                    % 比如 结果是 1 2 4 5 也就是说对应序号的stimuli里存在对应的神经元反应
                    for vv = consistent_ids.'
                        if ismember(vv,consistent_ids)
                            locallist(vv).label = 1;
                        end


                    end
                   

                end

                % 暂且放这里，后面继续改
                %figure;
                img = {};

                for m = 1: length(locallist)%setdiff([1: length(locallist)],ids_norm)
                    figure('Position',[293 772 1232 319]);
                    Draw.specSdf(locallist(m).pady,locallist(m).fs,locallist(m).padsptimes);
                    img{m} = getframe(gcf).cdata;
                    close(gcf)
                end

                summer = vertcat(img{:});
                figure('Position',[1938 -295 1014 1560]);
                imshow(summer)
                imwrite(summer,'陈留_DegFigure.png');
                disp('Done')

            end
            % based on the name of replaced song, find out the corresponding % norm song
        end

    
    end
    
end