% revise lots of segmentation file as a series 


dirpath = "E:\Stimuli_Source\R690_twoMotif\SegData"; %"E:\Stimuli_Source\senatusTwoMotif\SegData";
segfiles = Extract.filename(dirpath,'*.mat');



% dirpath = "E:\Stimuli_Source\allBirdsSong";
% segfiles = Extract.filesAllLevel(dirpath,'*.mat');



for k = 1: length(segfiles)
    
    reviseSeg(segfiles{k});
    uiwait(gcf);
end