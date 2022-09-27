%%%% Convert sylinf to new sylinf with cluster information

% three example features

meanf = {sylinf.meanf}.';   % mean frequency 
rawpitch = {sylinf.rawpitch}.';  % fundamental frequency
fm = {sylinf.fm}.';            % fundamental frequency change
entro = {sylinf.entropy}.';    % will be changed to harmonicity
% change their name 

fea = {'meanf','rawpitch','fm','entro'};
weight = [1,1,1,1];

distcollect = {};
for trump = 1: length(fea)
    
    eval(['feature = ',fea{trump},';']);
    
    
    dist = [];
    for a = 1: length(feature)
        for b = 1: length(feature)
            feature{a}(isnan( feature{a}))=0; % Convert nan to zero
            feature{b}(isnan( feature{b}))=0;% Convert nan to zero
            dist(a,b) = dtw(feature{a},feature{b});
        end
    end
    
    dist01 = rescale(dist,0,1);
    distcollect{trump} = dist01;

end


distsum = zeros(size(distcollect{1}));
for biden = 1: length(distcollect)
    distsum = distsum + distcollect{biden}*weight(biden);
end



dims = mdscale(distsum,10);  % dimensions in the mds non-paramteric

dimscell = num2cell(dims,2); % Convert dims matrix to cell

% Convert sylinf to new sylinf
[sylinf.dims] = dimscell{:};


% this is for clusteirng syllables using kmedoids clustering method
nclu = 10; % cluster to 10 clusters
cluidx = kmedoids(dims,nclu);

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






% figure; scatter(Y(:,1),Y(:,2));
