%drawAlignedNormFragTwoPlots
addpath(genpath("C:\Users\Zhehao\Dropbox (OIST)\My_Matlab\TEA"))
wb = waitbar(0,'Reading Analysis objects');
ANAfiles = extract.filename('./','*.mat');
for k = 21: length(ANAfiles)
    load(ANAfiles{k});
    
    A.drawAlignedNormFragTwoPlots;

    %A.compareWaveforms;

    waitbar(k/length(ANAfiles),wb,sprintf('%u of %u Neuron',k,length(ANAfiles)));
end

close(wb)
