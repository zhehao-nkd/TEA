classdef ResponseList
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Property1
    end

    methods(Static)

        function AlignedFragsDegs(inputlist, formated_name,outdir)

            if ~exist('outdir','var')
                outdir = '.';
            end

            dbstop if error
            tic

            % 首先，定义songlist
            songlist = inputlist(find(~cellfun(@isempty, regexp(cellstr({inputlist.stimuliname}.'),'norm'))));
            degids = find(~cellfun(@isempty, regexp(cellstr({inputlist.stimuliname}.'),'Deg|deg') ));
            if isempty(degids); return; end
            deglist = inputlist(degids);
            deg_bids = unique(cellfun(@Convert.bid,cellstr({deglist.stimuliname}.'),'Uni',0));

            I_song = {};
            for w = 1: length(deg_bids)
                Icollect = {};
                degexist_norm_list  = songlist(find(~cellfun(@isempty, regexp(cellstr({songlist.stimuliname}.'),deg_bids{w}))));
                for k = 1; length(degexist_norm_list)
                    smallfig = figure('Color','w','Position',[29 525 1953 420]);

                    newy = degexist_norm_list(k).plty(int64(1.8*degexist_norm_list(k).fs):...
                        int64(length(degexist_norm_list(k).plty) -1.8*degexist_norm_list(k).fs)  );

                    newsptimes = Extract.sptimes_resetSP(degexist_norm_list(k).pltsptimes,1.8,...
                        length(degexist_norm_list(k).plty)/degexist_norm_list(k).fs-1.8);

                    Draw.two(newy,degexist_norm_list(k).fs,newsptimes);

                    %Draw.two(degexist_norm_list(k).plty,degexist_norm_list(k).fs,degexist_norm_list(k).pltsptimes);
                    xlabel(degexist_norm_list(k).stimuliname);
                    frame = getframe(smallfig); Icollect{1} = frame.cdata;close(gcf)
                end
                I_song{w} = vertcat(Icollect{:});
            end
            %  [~,postunique] = unique(cellstr(cellfun(@Convert.bid,{songlist.stimuliname}.','Uni',0)));
            % songlist = songlist(postunique);

            % 其二，找到deressive songs对应的birdid
            degids = find(~cellfun(@isempty, regexp(cellstr({inputlist.stimuliname}.'),'Deg|deg') ));
            if isempty(degids); return; end
            deglist = inputlist(degids);
            for m = 1: length(deglist)
                birdid_collect{m} = Convert.bid(deglist(m).stimuliname);
                ids_norm = find(~cellfun(@isempty, regexp(cellstr({songlist.stimuliname}.'),birdid_collect{m}) ) );
                ids_norm = ids_norm(1); % ids_norm 可能有多个
                if ~isempty(ids_norm)& length(ids_norm) == 1
                    [deglist(m).sylIni,trump_diffvalue] = Neuron.findIni(songlist(ids_norm).plty,deglist(m).y);
                    fprintf('%s has the DIFFVALUE as: %u and INI as: %f \n ',deglist(m).stimuliname,trump_diffvalue,deglist(m).sylIni/32000);
                    deglist(m).pady = [zeros([deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext,1]);deglist(m).plty;zeros(length(songlist(ids_norm).plty)...
                        - (deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext) - length(deglist(m).plty),1)];
                    deglist(m).padsptimes = cellfun( @(x) x + (deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext)/deglist(m).fs, deglist(m).pltsptimes,'uni',0);
                end
            end

            unique_bid = unique(cellstr(birdid_collect));
            deg_fids = unique(cellstr({deglist.Fid}.'));
            deglist(1).Fid
            normlist = songlist(find(~cellfun(@isempty, regexp(cellstr({songlist.stimuliname}.'),strjoin(unique_bid,'|')) )));

            normlist = normlist(find(~cellfun(@isempty, regexp(cellstr({normlist.Fid}.'),strjoin(deg_fids,'|')))));

            I_Deg = {}; % 最初定义
            for w = 1: length(normlist)

                birdid = Convert.bid(normlist(w).stimuliname);
                ids_indeg = find(~cellfun(@isempty, regexp(cellstr({deglist.stimuliname}.'),birdid) ) );
                selected_deglist = deglist(ids_indeg);
                selected_deglist = selected_deglist(find(~cellfun(@isempty, regexp(cellstr({selected_deglist.Fid}.'),normlist(w).Fid))));
                [~,temp_index] = sortrows([selected_deglist.sylIni].');
                selected_deglist = selected_deglist(temp_index);

                % draw the norm figure
                Icollect = {};figure('Color','w','Position',[406 675 1378 420]);


                newy = normlist(w).plty(int64(1.8*normlist(w).fs):...
                    int64(length(normlist(w).plty) -1.8*normlist(w).fs)  );

                newsptimes = Extract.sptimes_resetSP(normlist(w).pltsptimes,1.8,...
                    length(normlist(w).plty)/normlist(w).fs-1.8);

                Draw.two(newy,normlist(w).fs,newsptimes);


                %Draw.two(normlist(w).plty,normlist(w).fs,normlist(w).pltsptimes);
                xlabel(sprintf('%s-%s',normlist(w).Fid,normlist(w).stimuliname));
                frame = getframe(gcf); Icollect{1} = frame.cdata;close(gcf);

                % draw the deg figure
                for hh = 1: length(selected_deglist)
                    figure('Color','w','Position',[406 675 1378 420]);

                    newy = selected_deglist(hh).pady(int64(1.8*selected_deglist(hh).fs):...
                        int64(length(selected_deglist(hh).pady) -1.8*selected_deglist(hh).fs)  );

                    newsptimes = Extract.sptimes_resetSP(selected_deglist(hh).padsptimes,1.8,...
                        length(selected_deglist(hh).pady)/selected_deglist(hh).fs-1.8);

                    Draw.two(newy,selected_deglist(hh).fs,newsptimes);

                    %Draw.two(selected_deglist(hh).pady,selected_deglist(hh).fs,selected_deglist(hh).padsptimes);
                    xlabel(sprintf('%s-%s',selected_deglist(hh).Fid,selected_deglist(hh).stimuliname));
                    frame = getframe(gcf); Icollect{1 + hh} = frame.cdata;close(gcf);
                end

                I_Deg{w} = vertcat(Icollect{:});
            end

            % 其三，找到对应birdid的frags
            fragids1 = find(~cellfun(@isempty, regexp(cellstr({inputlist.stimuliname}.'),'frag|Frag|syl|Syl') ));
            unique_bid = unique(cellstr(birdid_collect));
            fragids2 = find(~cellfun(@isempty, regexp(cellstr({inputlist.stimuliname}.'),strjoin(unique_bid,'|')) ));
            fragids = intersect(fragids1,fragids2);
            if ~isempty(fragids)
                fraglist = inputlist(fragids);
                for m = 1: length(fraglist)

                    birdid = Convert.bid(fraglist(m).stimuliname);
                    ids_norm = find(~cellfun(@isempty, regexp(cellstr({songlist.stimuliname}.'),birdid) ) );


                    if ~isempty(ids_norm)
                        ids_norm = ids_norm(1);
                        fraglist(m).sylIni = Neuron.findIni(songlist(ids_norm).plty,fraglist(m).y);
                        fraglist(m).pady = [zeros([fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext,1]);fraglist(m).plty;zeros(length(songlist(ids_norm).plty)...
                            - (fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext) - length(fraglist(m).plty),1)];
                        fraglist(m).padsptimes = cellfun( @(x) x + (fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext)/fraglist(m).fs, fraglist(m).pltsptimes,'uni',0);
                    end
                end
            end

            I_Frag = {};
            for w = 1: length(normlist)

                if ~isempty(fragids)
                    birdid = Convert.bid(normlist(w).stimuliname);
                    ids_infrag = ~cellfun(@isempty, regexp(cellstr({fraglist.stimuliname}.'),birdid) ) ;
                    selected_fraglist = fraglist(ids_infrag);
                    selected_fraglist = selected_fraglist(find(~cellfun(@isempty, regexp(cellstr({selected_fraglist.Fid}.'),normlist(w).Fid))));

       
                    [~,temp_index] = sortrows([selected_fraglist.sylIni].');
                    selected_fraglist = selected_fraglist(temp_index);
                end

                Icollect = {};figure('Color','w','Position',[406 675 1378 420]);

                newy = normlist(w).plty(int64(1.8*normlist(w).fs):...
                    int64(length(normlist(w).plty) -1.8*normlist(w).fs)  );

                newsptimes = Extract.sptimes_resetSP(normlist(w).pltsptimes,1.8,...
                    length(normlist(w).plty)/normlist(w).fs-1.8);

                Draw.two(newy,normlist(w).fs,newsptimes);


                % Draw.two(normlist(w).plty,songlist(w).fs,normlist(w).pltsptimes); xlabel(sprintf('%s-%s',normlist(w).Fid,normlist(w).stimuliname));
                frame = getframe(gcf);   Icollect{1} = frame.cdata; close(gcf)

                if ~isempty(fragids)
                    for bb = 1: length(selected_fraglist)
                        figure('Color','w','Position',[406 675 1378 420]);
                        newy = selected_fraglist(bb).pady(int64(1.8*selected_fraglist(bb).fs):...
                            int64(length(selected_fraglist(bb).pady) -1.8*selected_fraglist(bb).fs)  );

                        newsptimes = Extract.sptimes_resetSP(selected_fraglist(bb).padsptimes,1.8,...
                            length(selected_fraglist(bb).pady)/selected_fraglist(bb).fs-1.8);

                        Draw.two(newy,selected_fraglist(bb).fs,newsptimes);

                        %Draw.two(selected_fraglist(bb).pady,selected_fraglist(bb).fs,selected_fraglist(bb).padsptimes);
                        xlabel(sprintf('%s-%s',selected_fraglist(bb).Fid,selected_fraglist(bb).stimuliname));
                        frame = getframe(gcf);Icollect{1 + bb} = frame.cdata;close(gcf);
                    end
                end
                I_Frag{w} = vertcat(Icollect{:});
            end


            % 其末 draw and save
            %             figure;
            %             neu.waveform.draw1st;
            %             temp = getframe(gcf);close(gcf);
            %             w_img = temp.cdata;
            %             I_WF{1} = w_img;

            %             try
            I_of_each_column = horzcat(I_song,I_Deg,I_Frag);
            %             catch
            % %                 I_of_each_column = horzcat(I_song,I_Deg,I_Frag,I_WF);
            %             end
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

            imwrite(Iall,sprintf('%s\\对齐AlignedFragsDegs_%s.png',outdir,formated_name));
            toc


        end


        function result = RatioBufenZhengti(list)

            % Context effect 是多少？计算 syllable分开后的和值/Syllables as sequence的值
            raw_samesongfrags = list(find(~cellfun(@isempty,regexp( cellstr({list.stimuliname}.'),'samesong'))));

            if length(raw_samesongfrags) == 0
                result.ratio_bufen_zhengti = nan;  % 这个Neuron是否对single syllable之和的反应是否高于对整首song的反应？
                result.normalized_ratio = nan;
                return
            end

            [counts, groupnames] = groupcounts({raw_samesongfrags.Fid}.');
            [maxvalue, maxid] = max(counts);
            if length(counts) > 1
                warning('Dangerous!!!!!!');
                %pause
            end
            selected_fid = groupnames(maxid);
            samesongfrags = raw_samesongfrags(find(~cellfun(@isempty,regexp( cellstr({raw_samesongfrags.Fid}.'),selected_fid))));

            songname = unique(cellstr({samesongfrags.stimuliname}.'));

            songid = {};
            for h = 1:length(samesongfrags)

                songid{h} = regexp(convertCharsToStrings(samesongfrags(h).stimuliname),'[OGBYR]\d+','match');

            end

            replalist = list(find(~cellfun(@isempty,regexp( cellstr({list.stimuliname}.'),'Repla'))));
            
            replasongid = {};
            for h = 1:length(replalist)

                parts = strsplit(replalist(h).stimuliname,'before');
                replasongid{h} = regexp(convertCharsToStrings(parts{2}),'[OGBYR]\d+','match');
            end

            if ~isempty(replasongid )
                targetsongid = intersect(cellstr(songid.'),cellstr(replasongid.'));

            else
                targetsongid = songid;
            end


            if length(unique(cellstr(targetsongid.'))) == 1

                corresp_songid = targetsongid{1};

                criteria1 = find(~cellfun(@isempty,regexp( cellstr({list.stimuliname}.'),corresp_songid)));
                criteria2 = find(~cellfun(@isempty,regexp( cellstr({list.stimuliname}.'),'norm')));
                criteria3 = find(~cellfun(@isempty,regexp( cellstr({list.Fid}.'),selected_fid)));
                hitted_song = list(intersect(criteria1,intersect(criteria2,criteria3)));

                if length(hitted_song) == 0
                    result.ratio_bufen_zhengti = nan;  % 这个Neuron是否对single syllable之和的反应是否高于对整首song的反应？
                    result.normalized_ratio = nan;
                    return
                end

                song_prespikes = hitted_song.presptimes;
                song_duringspikes = hitted_song.sptimes;
                song_numspikes = length(vertcat(song_duringspikes{:})) - length(vertcat(song_prespikes{:}));


            end

            [~,index1] = sortrows(cellstr({samesongfrags.stimuliname}.')); samesongfrags = samesongfrags(index1); clear index1


            summer_frags_numspikes = [];
            for h = 1:length(samesongfrags)

                thisrow = samesongfrags(h);
                frag_zpt = thisrow.zpt;


                if h ~= length(samesongfrags)

                    frag_sptimes = Extract.sptimes_resetSP( thisrow.rawsptimes, frag_zpt, frag_zpt + 0.5); % 取500ms的间距
                    frag_presptimes = Extract.sptimes_resetSP( thisrow.rawsptimes, frag_zpt -0.5, frag_zpt ); % 取500ms的间距
                else
                    frag_sptimes = thisrow.sptimes; % 取500ms的间距
                    frag_presptimes =thisrow.presptimes; % 取500ms的间距

                end


                summer_frags_numspikes(h) = length(vertcat(frag_sptimes{:})) - length(vertcat(frag_presptimes{:}));

            end


            frags_numspikes = sum(summer_frags_numspikes);

            result.ratio_bufen_zhengti = frags_numspikes/song_numspikes;  % 这个Neuron是否对single syllable之和的反应是否高于对整首song的反应？
            result.normalized_ratio = (frags_numspikes - song_numspikes)/(abs(frags_numspikes)+abs(song_numspikes));


        end




        function drawSelectedStimuliResp(list,stimulinames,range)
            % range是从plty的起始为始，以秒计
            dbstop if error; tic
            if isa(stimulinames,'cell')
                selectedids = find(~cellfun(@isempty, regexp(cellstr({list.stimuliname}.'),strjoin(stimulinames,'|') )));
            else
                selectedids = stimulinames;
            end
            
            
            songids = find(~cellfun(@isempty, regexp(cellstr({list.stimuliname}.'),'norm|spe') ));
            songids = intersect(songids,selectedids);
            songlist = list(songids);
            %             [~,postunique] = unique(cellstr(cellfun(@Convert.bid,{songlist.stimuliname}.','Uni',0)))
            %             songlist = songlist(postunique);
            RONGYU = 0.5;
            range = range + RONGYU;
            for k = 1: length(songlist)
                songlist(k).pady = [zeros(RONGYU*songlist(k).fs,1);songlist(k).plty];
                songlist(k).padsptimes = cellfun( @(x) x + RONGYU, songlist(k).pltsptimes,'uni',0);
            end
            
            % About Frag
            fragids = find(~cellfun(@isempty, regexp(cellstr({list.stimuliname}.'),'frag|Frag|syl|Syl') ));
            fragids = intersect(fragids,selectedids);
            if ~isempty(fragids)
                fraglist = list(fragids);
                for m = 1: length(fraglist)
                    birdid = Convert.bid(fraglist(m).stimuliname);
                    ids_norm = find(~cellfun(@isempty, regexp(cellstr({songlist.stimuliname}.'),birdid) ) );
                    if ~isempty(ids_norm)& length(ids_norm) == 1
                        fraglist(m).sylIni = Neuron.findIni(songlist(ids_norm).pady,fraglist(m).y);
                        fraglist(m).pady = [zeros([fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext,1]);fraglist(m).plty;zeros(length(songlist(ids_norm).plty)...
                            - (fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext) - length(fraglist(m).plty),1)];
                        fraglist(m).padsptimes = cellfun( @(x) x + (fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext)/fraglist(m).fs, fraglist(m).pltsptimes,'uni',0);
                    end
                end
            end
            
            % About Deg
            degids = find(~cellfun(@isempty, regexp(cellstr({list.stimuliname}.'),'Deg|deg') ));
            degids = intersect(degids,selectedids);
            if ~isempty(degids)
                deglist = list(degids);
                for m = 1: length(deglist)
                    birdid = Convert.bid(deglist(m).stimuliname);
                    ids_norm = find(~cellfun(@isempty, regexp(cellstr({songlist.stimuliname}.'),birdid) ) );
                    if ~isempty(ids_norm)& length(ids_norm) == 1
                        [deglist(m).sylIni,trump_diffvalue] = Neuron.findIni(songlist(ids_norm).pady,deglist(m).y);
                        fprintf('%s has the DIFFVALUE as: %u and INI as: %f \n ',deglist(m).stimuliname,trump_diffvalue,deglist(m).sylIni/32000);
                        deglist(m).pady = [zeros([deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext,1]);deglist(m).plty;zeros(length(songlist(ids_norm).plty)...
                            - (deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext) - length(deglist(m).plty),1)];
                        deglist(m).padsptimes = cellfun( @(x) x + (deglist(m).sylIni-1-deglist(m).fs*deglist(m).pltext)/deglist(m).fs, deglist(m).pltsptimes,'uni',0);
                    end
                end
            end
            
            % About Repla
            replaids = find(~cellfun(@isempty, regexp(cellstr({list.stimuliname}.'),'Repla|repla|catego|Catego') ));
            replaids = intersect(replaids,selectedids);
            if ~isempty(replaids)
                replalist = list(replaids);
                ids_norm_collecting_box = [];
                for m = 1: length(replalist)
                    
                    afterBefore = regexp(replalist(m).stimuliname,'(?<=before-)\S*','match');
                    afterBefore = afterBefore{1};
                    birdid = Convert.bid(afterBefore);
                    ids_norm = find(~cellfun(@isempty, regexp(cellstr({songlist.stimuliname}.'),birdid) ) );
                    
                    if ~isempty(ids_norm)%& length(ids_norm) == 1
                        ids_norm_1st = ids_norm(1);
                        ids_norm_collecting_box = [ids_norm_collecting_box, ids_norm];
                        afterpad_length = length(songlist(ids_norm_1st).plty) + RONGYU*songlist(ids_norm_1st).fs;  % +0.5s
                        
                        replalist(m).pady = [zeros(afterpad_length- length(replalist(m).plty),1);replalist(m).plty];
                        replalist(m).padsptimes = cellfun( @(x) x + length(zeros(afterpad_length- length(replalist(m).plty),1))/fraglist(m).fs,replalist(m).pltsptimes,'uni',0);
                    end
                end
            end
            
            % merge the Cons,Degs,Frags,Replas lists together
            for w = 1: length(songlist)
                
                if ~isempty(degids)
                    birdid = Convert.bid(songlist(w).stimuliname);
                    ids_indeg = find(~cellfun(@isempty, regexp(cellstr({deglist.stimuliname}.'),birdid) ) );
                    selected_deglist = deglist(ids_indeg);
                    [~,temp_index] = sortrows([selected_deglist.sylIni].');
                    selected_deglist = selected_deglist(temp_index);
                else
                    selected_deglist = [];
                end
                
                if ~isempty(fragids)
                    birdid = Convert.bid(songlist(w).stimuliname);
                    ids_infrag = find(~cellfun(@isempty, regexp(cellstr({fraglist.stimuliname}.'),birdid) ) );
                    selected_fraglist = fraglist(ids_infrag);
                    [~,temp_index] = sortrows([selected_fraglist.sylIni].');
                    selected_fraglist = selected_fraglist(temp_index);
                else
                    selected_fraglist = [];
                end
                
                if ~isempty(replaids)
                    birdid = Convert.bid(songlist(w).stimuliname);
                    ids_inrepla = find(~cellfun(@isempty, regexp(cellstr({replalist.stimuliname}.'),birdid) ) );
                    selected_replalist = replalist(ids_inrepla);
                else
                    selected_replalist = [];
                end
                
                % draw figures
                if ~isempty(selected_deglist)
                    selected_deglist = rmfield(selected_deglist,'sylIni');
                end
                
                if ~isempty(selected_fraglist)
                    selected_fraglist = rmfield(selected_fraglist,'sylIni');
                end
                alllist = horzcat(songlist(w),selected_deglist,selected_fraglist,selected_replalist);
                len = length(alllist);
                figure('Position',[1935 -207 520 1458/9*len],'Color','none');
                
                ax = tight_subplot(2*len, 1, 0.002, 0.02, 0);
                for k = 1: len
                    
                    axes(ax(2*(k-1)+ 1)); % Draw.two(,,);
                    %ax(2*(k-1)+ 1).Position(4) =  ax(2*(k-1)+ 1).Position(4);
                    if exist('range','var')
                        truncated_y = alllist(k).pady(range(1)*alllist(k).fs:range(2)*alllist(k).fs);
                        Draw.spec(truncated_y,alllist(k).fs);
                    else
                        Draw.spec(alllist(k).pady,alllist(k).fs);
                    end
                    xlabel('')
                    ylabel('')
                    
                    set(gca,'TickLength',[0 .01])
                    set(gca,'Yticklabel',[])
                    set(gca,'Xticklabel',[])
                    
                    axes(ax(2*k));
                    
                    if exist('range','var')
                        truncated_sptimes = Extract.sptimes_resetSP(alllist(k).padsptimes, range(1), range(2));
                        truncated_y = alllist(k).pady(range(1)*alllist(k).fs:range(2)*alllist(k).fs);
                        %Draw.raster(truncated_sptimes,truncated_y,alllist(k).fs);
                        Draw.raster(truncated_sptimes,truncated_y,alllist(k).fs,2.8,'k');
                    else
                        %Draw.raster(alllist(k).padsptimes,alllist(k).pady,alllist(k).fs);
                        Draw.raster(alllist(k).padsptimes,alllist(k).pady,alllist(k).fs,2.8,'k');
                    end
                    
                    ylabel('')
                    set(gca,'TickLength',[0 .01])
                    
                    
                    
                    %                         set(gca,'Yticklabel',[])
                    %                         set(gca,'Xticklabel',[])
                    %                     else
                    %                         xlabel(alllist(k).stimuliname);
                    %                     end
                    
                end
                
                for k = 1:len
                    ax(2*k-1).Position(4) =  ax(2*k-1).Position(4) -0.008;
                    ax(2*k-1).Position(2) =  ax(2*k-1).Position(2)+0.006;
                    box(ax(2*k-1),'off');
                    set(ax(2*k-1),'XColor','none','YColor','none')
                    %set(ax(2*k-1), 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
                end
                
                for k = 1:len
                    ax(2*k).Position(2) =  ax(2*k).Position(2)+0.000;
                    ax(2*k).Position(4) =  ax(2*k).Position(4)+0.002;
                    box(ax(2*k),'off');
                    set(ax(2*k),'XColor','none','YColor','none')
                    %set(ax(2*k), 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
                end
                disp('pause here')
                
                
                
                
            end
        end
     
    
    end


end