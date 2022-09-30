
% using object batch and ephys
mode = 2; % mode 2 means elemode
dbstop if error
b = Chorus("C:\Users\Zhehao\Dropbox (OIST)\My_EphysInput\input2021.xlsx");

b.select(1);

temp = b.getn;
neuron = temp{1};

% using object stimuli
load("C:\Users\Zhehao\Dropbox (OIST)\My_Stimuli\B659@06242021\mergedsylinf.mat");  % syl
load("C:\Users\Zhehao\Dropbox (OIST)\My_Stimuli\B659@06242021\mergedeleinf.mat");  % ele

frag = stimuli.polish(mergedsylinf);
stim = stimuli(frag);
splited = stim.splited;

frag2 = stimuli.polish(eleinf);
elestim = stimuli(frag2);
elesplited = elestim.splited;

%%%%%% sum of response-eliciting syllable/element id 
elicitingid = [];
molihuacha = 0; % idx controller

%thesamplesongs = unique(cellstr({frag.songname}.'))


for k = 1: length(neuron.e)

    fs = 32000;
   
    
    % generate summed y
    thiss = splited{k};
    ys = {thiss.y}.';
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
        terminals = initials + [thise.dur].';
    end
   
        
    % draw subplot 1
    figure('Position',[172 81 1715 1008]); 
    subplot(2,1,1);
    %%%%% there are something wrong with the eleinf's terminal calculation
    
    Draw.spec(summed,fs);
    justatemp = unique(cellstr({thiss.songname}.'));
    xlabel(justatemp{1});
    hold on
    for a = 1: length(initials)
        thei = initials(a); % this initial
        thet = terminals(a); % this terminal
        line([thei, thei],[0,16],'Color','r');
        hold on
        line([thet, thet],[0,16],'Color','b');
        hold on
    end
    
    hold off
    
    % draw subplot 2
    
    subplot(2,1,2);
    
    color = 'k';
    title = ' ';
    
    e = neuron.e{k};
    Draw.raster(e.sptimes,e.y,e.fs,color);  % this draw the raster based on Experiment based data
    justanothertemp = e.sound.name;
    xlabel(justanothertemp);
    % ginput, then mark them
    subplot(2,1,1)
    hold on
    [xs,ys]=ginput;
    for m = 1: length(xs)
        x = xs(m);
        line([x,x],[0,16],'Color','g');
        hold on
    end
    
    close(gcf)
    % calculate which element the ginput locates
    
    for j = 1: length(xs)
        
        thisx = xs(j);
        initialleft = find(initials<thisx);
        terminalright = find(terminals>thisx);
        rank = intersect(initialleft,terminalright);
        
        molihuacha = molihuacha + 1;
        
        if mode == 1
         elicitingid(molihuacha) = thiss(rank).uniqueid;
        elseif mode == 2
         elicitingid(molihuacha) = thise(rank).uniqueid;
        end
        
        
    end
        
        
    
    
end

% based on the eliciting id, get the info;

elicitinginfo = frag(elicitingid);


