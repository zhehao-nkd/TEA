% manFindEli
function  [responsiveinfo,neweleinf] = manPickEli(neuron,mergedeleinf,mergedsylinf)
    mode = 2; % mode 2 means elemode
    dbstop if error
%     b = Batch(xlsxdir);
%     b.select(neuronnum);
% 
%     temp = b.getn;
%     neuron = temp{1};

    mergedeleinf = stimuli.adduniqueid(mergedeleinf);
    mergedsylinf = stimuli.adduniqueid(mergedsylinf);
    % using object stimuli

    sylstim = stimuli(mergedsylinf);
    splited = sylstim.splited;
    elestim = stimuli(mergedeleinf);
    elesplited = elestim.splited;

    %%%%%% sum of response-eliciting syllable/element id
    elicitingid = [];
    molihuacha = 0; % idx controller

    %thesamplesongs = unique(cellstr({frag.songname}.'))

    for k = 1: length(neuron.e)   % 对 每个 three plot

        fs = 32000;
        thiss = splited{k};
        ys = {thiss.y}.';     % generate summed y
        gaps = [thiss(:).pregap].';
        gaps(1) = 0; % replace the first pregap to 0,which is a inf
        summed = [];

        for haris = 1: length(ys)
            summed = [summed;zeros(int64(gaps(haris)*fs),1);ys{haris}];
        end
        % get all initials and all terminals
        if mode == 1
            initials = [thiss.initial].';
            % normalize initials
            initials = initials - initials(1);
            initials =initials/fs;
            terminals = initials + [thiss.dur].';
        elseif mode == 2 % ele mode
            thise = elesplited{k};
            initials = [thise.initial].';
            % normalize initials
            initials = initials - initials(1);
            initials =initials/fs;
            lenys = cellfun(@length,{thise.y}.');
            durs = lenys/fs;
            terminals = initials + durs;
        end

        %%%%%%%%%%%%%%%%%%%%% draw subplot 1
        figure('Position',[172 81 1715 1008]);
        subplot(2,1,1);
        %%%%% there are something wrong with the eleinf's terminal calculation

        draw.spec6(summed,fs);
        justatemp = unique(cellstr({thiss.songname}.'));
        xlabel(justatemp{1});
        hold on
        for a = 1: length(initials)
            thei = initials(a); % this initial
            thet = terminals(a); % this terminal
            line([thei, thei],[0,16000],'Color','r');
            hold on
            line([thet, thet],[0,16000],'Color','b');
            hold on
        end

        hold off

        %%%%%%%%%%%%%%%%%%%% draw subplot 2

        subplot(2,1,2);
        color = 'k';
        title = ' ';
        e = neuron.e{k};
        draw.raster(e.sptimes,e.y,e.fs);  % this draw the raster based on Neuron based data
        justanothertemp = e.sound.name;
        xlabel(justanothertemp);




        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ginput, then mark them
        subplot(2,1,1)
        hold on
        [xs,ys]=ginput;
        for m = 1: length(xs)
            x = xs(m);
            line([x,x],[0,16],'Color','g');
            hold on
        end

        % calculate which element the ginput locates
        close(gcf)
        for j = 1: length(xs)

            thisx = xs(j);
            initialleft = find(initials<thisx);
            terminalright = find(terminals>thisx);
            rank = intersect(initialleft,terminalright);
            
            if isempty(rank)
                break  % in case you click some weird location
            end

            molihuacha = molihuacha + 1;

            if mode == 1
                elicitingid(molihuacha) = thiss(rank).uniqueid;
            elseif mode == 2
                elicitingid(molihuacha) = thise(rank).uniqueid;
            end
        end
    end



    % based on the eliciting id, get the info;
    neweleinf = mergedeleinf;

    for d = 1: length(mergedeleinf)
        if ismember(d,elicitingid)
            neweleinf(d).respinsong = 1; % neuron response is responsive
        else
            neweleinf(d).respinsong = 0;
        end


    end

    % get the new mergedeleinf 
    
    

    %%%%% then plot all responsive syllables

    responsiveinfo = neweleinf(elicitingid);
    figure;

    for v = 1: length(responsiveinfo)
        subplot(1,length(responsiveinfo),v);
        draw.spec(responsiveinfo(v).y,responsiveinfo(v).fs);
        titletext =  sprintf('%s-%u',responsiveinfo(v).songname,responsiveinfo(v).fragid) ;
        xlabel(titletext);
        ylabel('')
    end
end
