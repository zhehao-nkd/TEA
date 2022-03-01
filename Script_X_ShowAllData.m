% This script read data tables downloaded from google sheet
dbstop if error

birdid = {'B659','Y661','Y675','Y687','O686','G656','R677'};
rangeid = {'B_659','Y_661','Y_675','Y_687','O_686','G_656','R_677'};
tablepath = "C:\Users\Zhehao\Desktop\AllInOne (23).xlsx";

Tcollect = {}; % to collect tables
for bid = 1: length(birdid)
    
    
    Tcollect{bid} = table2struct(readtable(tablepath,'Sheet',birdid{bid},'Range',rangeid{bid}));
end


XLSRangeNames('C:\Users\Zhehao\Desktop\AllInOne (23).xlsx');



Tall = vertcat(Tcollect{:});

candidates = find([Tall.Goodness].' == 1|[Tall.Goodness].' == 2|[Tall.Goodness].' == 3|[Tall.Goodness].' == 4);

%candidates = find([Tall.Goodness].' == 4.5);

%candidates = find([Tall.Goodness].' == 1|[Tall.Goodness].' == 2|[Tall.Goodness].' == 3);


% candidates = find([Tall.Goodness].' == 4);
Tcan = Tall(candidates);

% k = 7 ,8 , 24 are dead!!!!!!!!!! for rank 1 2 3

% k = 21 27 39 40 47 , 60, 61 are dead for rank 4
for k = 1: length(Tcan)
    
    b = Batch(Tcan(k).TxtPath,Tcan(k).PlxPath,Tcan(k).StimuliPath);
    
    this_channel = Tcan(k).ChannelName;
    this_unit = Tcan(k).UnitName;
    neu_list = {b.nlist.neuronname}.';
    
    channel_ids = find(~cellfun(@isempty,regexp(neu_list,this_channel)));
    unit_ids = find(~cellfun(@isempty,regexp(neu_list,sprintf('_%u',this_unit))));
    neuron_ids = intersect(channel_ids,unit_ids);
    
    b.select(neuron_ids);
    N = b.getn{1};  
    N.signalGoodness = Tcan(k).Goodness;
    
    % allocate sap-based feature information to each neuron
    sorted_data = Neuron.extractFeaturesFromSapRawDataJson(Tcan(k).FeatureData, Tcan(k).FeatureInfo);
    N.setEachStimuliSapFeatures(sorted_data);
    N.calMeanFeatures;
    
    save(N.neuronname,'N');
    
    %A = Analysis(N);
    
end

neuronfiles = extract.filename('./','*.mat');

% m = 39 is wrong
for m =40 : length(neuronfiles)
    
    load(neuronfiles{m});
    
    A = Analysis(N);
   
    N.threePlotsWithFreqRelatedFeatures;
    A.drawFeatureDistrutionInDifferentDimensions;
end



neuronfiles = extract.filename('./','*.mat');

% m = 39 is wrong
for m =40 : length(neuronfiles)
    
    load(neuronfiles{m});
    A = Analysis(N);
    A.linearfitFeature;
    N.calSTA;
end


% for all neurons

% m = 39 is wrong!!!
for m =40 : length(neuronfiles)
    
    load(neuronfiles{m});
    
    
    N.threePlotsWithFreqRelatedFeatures;  % ThreePlots for all the stimuli
    N.calSTA; % estimate STA
    A = Analysis(N);
    A.drawDTWSimilarityMatrixBasedOnZscoredData;
    A.drawDTWSimilarityMatrix;
    A.drawCoeffOfFeaturesLinearFit;
    A.drawMeanFeatureVsResp;
    A.drawMeanFeaturesVsRespAsLineChart;

end


neuronfiles = extract.filename('./','*.mat');

% m = 39 is wrong
f = waitbar(0,'Starting');
for m = 40 : length(neuronfiles)
    load(neuronfiles{m});
    N.calSTA;
    A = Analysis(N);
    A.drawPairwiseFragmentsMeanFeaturesDistribution
    
    waitbar(m/length(neuronfiles),f,sprintf( 'Number %u',m/length(neuronfiles)));
end

close(f);

neuronfiles = extract.filename('./','*.mat');
for k = 1: length(neuronfiles)
    load(neuronfiles{k});
    A = Analysis(N);
    A.set_eleinf("C:\Users\Zhehao\Dropbox (OIST)\SaveAllMatXlsData\My_eleinf\all_eleinf.mat");
    A.drawMeanFeaturesInSongVsRespAsLineChart;
end


