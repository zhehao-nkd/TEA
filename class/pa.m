classdef pa
   % ps is a paramter class
    properties (Constant)
        wav = "C:\Users\Zhehao\Desktop\B464A.wav";
        folder = "Z:\Yazaki-SugiyamaU\Zhehao\My_Stimuli\G606\G606-1st";
        txt = "E:\consongPlx\G606_Z01SO.txt";
        plx  = "E:\consongPlx\G606_Z01SO.plx";
        T = "E:\DiskDBackUp\ComprehensiveAnalysis\pathTable.xlsx"
        a = "C:\Users\Remote\Documents\G606_Z01SO.txt"
        b = "C:\Users\Remote\Documents\G606_Z01SO.plx"
        c = "C:\Users\Remote\Documents\G606-1st"
        bucket = "Z:\Yazaki-SugiyamaU\Bird-song"
        hp = 450;  % high pass threshold
        SYLLEN = 12800 % threshold for screnning the syllables
        wavfolders = unique(readtable(pa.T).path_folder);
% this class store the commonly used paramters( ect. pathes) for testing
% the function
    end
 

end