dbstop if error

ANAfiles = extract.filename('./','*.mat');

wuhu = waitbar(0,'Start processing');

for m =  [11,43]

    load(ANAfiles{m});
    A.set_eleinf("C:\Users\Zhehao\Dropbox (OIST)\SaveAllMatXlsData\My_eleinf\all_eleinf.mat");
    A.V1drawMeanFeaturesInSongVsRespAsLineChart;
    A.drawDTWSimilarityMatrixBasedOnZscoredData;
    A.drawDTWSimilarityMatrix;
    A.drawCoeffOfFeaturesLinearFit;
    A.drawMeanFeatureVsResp;
    A.drawMeanFeaturesVsRespAsLineChart;
    A.drawPairwiseFragmentsMeanFeaturesDistribution;
    A.threePlotsWithPitch;
    A.V2drawMeanFeaturesInSongVsRespAsLineChart
    waitbar(m/length(ANAfiles),wuhu,sprintf('%u of totally %u files',m,length(ANAfiles)));
    %A.sort_frags_by_response_strength_and_then_draw
end

close(wuhu);