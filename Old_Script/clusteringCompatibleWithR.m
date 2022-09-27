
meanf = readtable("C:\Users\Zhehao\Dropbox (OIST)\My_table\featurelist.xlsx",'Sheet', 2);
meanf = table2array(meanf);


ff = readtable("C:\Users\Zhehao\Dropbox (OIST)\My_table\featurelist.xlsx",'Sheet', 3);
ff = table2array(ff);

ffc = readtable("C:\Users\Zhehao\Dropbox (OIST)\My_table\featurelist.xlsx",'Sheet', 4);
ffc = table2array(ffc);

harmo = readtable("C:\Users\Zhehao\Dropbox (OIST)\My_table\featurelist.xlsx",'Sheet', 5);
harmo = table2array(harmo);


fea = {'meanf','ff','ffc','harmo'};
weight = [1,1,1,1];

distcollect = {};
for trump = 1: length(fea)
    
    eval(['feature = ',fea{trump},';']);
    
    
    dist = [];
    for a = 1: length(feature)
        for b = 1: length(feature)
            alpha = rmmissing(feature(a,:));
            beta = rmmissing(feature(b,:));
            dist(a,b) = dtw(alpha,beta);
        end
    end
    
    dist01 = rescale(dist,0,1);
    distcollect{trump} = dist01;

end

distsum = zeros(size(distcollect{1}));
for biden = 1: length(distcollect)
    distsum = distsum + distcollect{biden}*weight(biden);
end



[dims,stress] = mdscale(distsum,10);  % dimensions in the mds non-paramteric\

disp(['Stress is',stress]);

dimscell = num2cell(dims,2); % Convert dims matrix to cell


%%%%%%%%%% another section
% read exported song
wavs = Extract.filename("C:\Users\Zhehao\Dropbox (OIST)\My_Luscinia\ExportedSong",'*.wav');

for w = 1: length(wavs)
    [song(w).y,song(w).fs] = audioread(wavs{w});
    [~,song(w).name,~] = fileparts(wavs{w});
end

% element inf
rawinf = table2struct(readtable("C:\Users\Zhehao\Dropbox (OIST)\My_Luscinia\eleinf.csv"));

% construct the element inf structure
eleinf = struct;
for e = 1: length(rawinf)
    
    eleinf(e).sound = rawinf(e).SongName;
    eleinf(e).number = rawinf(e).ElementNumber; % element number
    songidx = find(~cellfun(@isempty, regexp([song(:).name].', eleinf(e).sound)));
    songy = song(songidx).y;
    eleinf(e).initial = rawinf(e).StartTime*0.001*song(songidx).fs;
    eleinf(e).terminal = (rawinf(e).StartTime+rawinf(e).Length)*0.001*song(songidx).fs;
    eleinf(e).y = songy( eleinf(e).initial: eleinf(e).terminal);
    
end
% Convert sylinf to new sylinf
[eleinf.dims] = dimscell{:};


% this is for clusteirng syllables using kmedoids clustering method
nclu = 4; % cluster to 10 clusters
cluidx = x;
%cluidx = kmedoids(dims,nclu);

cluidxcell = num2cell(cluidx,2); % Convert dims matrix to cell

% Convert sylinf to new sylinf
[sylinf.cluidx] = cluidxcell{:};

figure('Color','w');
for i = 1: nclu
    thisclu = find(cluidx==i);
    scatter(dims(thisclu,1),dims(thisclu,2),'filled');
    hold on
end
xlabel('Dim1'); ylabel('Dim2')
