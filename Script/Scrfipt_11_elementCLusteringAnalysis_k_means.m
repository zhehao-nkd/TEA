
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

dimscell = num2cell(dims,2); % convert dims matrix to cell

% convert sylinf to new sylinf
[sylinf.dims] = dimscell{:};


% this is for clusteirng syllables using kmedoids clustering method
nclu = 10; % cluster to 10 clusters
cluidx = kmedoids(dims,nclu);

cluidxcell = num2cell(cluidx,2); % convert dims matrix to cell

% convert sylinf to new sylinf
[sylinf.cluidx] = cluidxcell{:};

figure('Color','w');
for i = 1: nclu
    thisclu = find(cluidx==i);
    scatter(dims(thisclu,1),dims(thisclu,2),'filled');
    hold on
end
xlabel('Dim1'); ylabel('Dim2')
