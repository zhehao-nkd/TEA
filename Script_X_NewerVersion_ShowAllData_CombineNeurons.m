% This script read data tables downloaded from google sheet
dbstop if error

birdid = {'B659','Y661','Y675','Y687','O686','G656','R677'};
rangeid = {'B_659','Y_661','Y_675','Y_687','O_686','G_656','R_677'};
tablepath = "C:\Users\Zhehao\Desktop\AllInOne (27).xlsx";

Tcollect = {}; % to collect tables
for bid = 1: length(birdid)

    Tcollect{bid} = table2struct(readtable(tablepath,'Sheet',birdid{bid},'Range',rangeid{bid}));
    
end


T = vertcat(Tcollect{:});
IDs = [T.UniqueID].';
num_ids = rmmissing(unique(IDs(IDs~=0)));

%kbad = [28,31,32,41,42,44,64,65];


wb = waitbar(0,'Creating Neuron analysis objects');
for k = 1: length(num_ids)
    waitbar(k/length(num_ids),wb,sprintf('%u of %u Neuron',k,length(num_ids)));
    ids_in_T = find(num_ids(k) == IDs);
    
    Ns = {};
    for i = 1:length(ids_in_T)
        b = Batch(T(ids_in_T(i)).TxtPath,T(ids_in_T(i)).PlxPath,T(ids_in_T(i)).StimuliPath);
        
        this_channel = T(ids_in_T(i)).ChannelName;
        this_unit = T(ids_in_T(i)).UnitName;
        neu_list = {b.nlist.neuronname}.';
        
        channel_ids = find(~cellfun(@isempty,regexp(neu_list,this_channel)));
        unit_ids = find(~cellfun(@isempty,regexp(neu_list,sprintf('_%u',this_unit))));
        neuron_ids = intersect(channel_ids,unit_ids);
        
        b.select(neuron_ids);
        N = b.getn{1};
        N.signalGoodness = T(ids_in_T(i)).Goodness;
        N.set_uniqueid(T(ids_in_T(i)).UniqueID);
        % allocate sap-based feature information to each neuron
        sorted_data = Neuron.extractFeaturesFromSapRawDataJson(T(ids_in_T(i)).FeatureData, T(ids_in_T(i)).FeatureInfo);
        N.setEachStimuliSapFeatures(sorted_data);
        N.calMeanFeatures;
        Ns{i} = N;
    end
    
    A = Analysis(Ns);
    A.uniqueid = num_ids(k);
    
    
    save(sprintf('%s_%u',A.birdid,A.uniqueid),'A','-v7.3');
    
    
end

close(wb);


ANAfiles = extract.filename('./','*.mat');


for m = 6 : length(ANAfiles)
    
    load(neuronfiles{m});
    
    A.set_eleinf("C:\Users\Zhehao\Dropbox (OIST)\SaveAllMatXlsData\My_eleinf\all_eleinf.mat");
    A.V1drawMeanFeaturesInSongVsRespAsLineChart;
    %     N.threePlotsWithFreqRelatedFeatures;  % ThreePlots for all the stimuli
    %     N.calSTA; % estimate STA
    %     A = Analysis(N);
    A.drawDTWSimilarityMatrixBasedOnZscoredData;
    A.drawDTWSimilarityMatrix;
    A.drawCoeffOfFeaturesLinearFit;
    A.drawMeanFeatureVsResp;
    A.drawMeanFeaturesVsRespAsLineChart;
    A.drawPairwiseFragmentsMeanFeaturesDistribution;
    A.threePlotsWithPitch;
    A.V2drawMeanFeaturesInSongVsRespAsLineChart
end


%drawAlignedNormFragTwoPlots

wb = waitbar(0,'Reading Analysis objects');
ANAfiles = extract.filename('./','*.mat');
for k = 5: length(ANAfiles)
    load(ANAfiles{k});
    
    A.drawAlignedNormFragTwoPlots;

    waitbar(k/length(ANAfiles),wb,sprintf('%u of %u Neuron',k,length(ANAfiles)));
end

close(wb)


% Neuron.cal_latency is
% wrong!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
