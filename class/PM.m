classdef PM % parameter
   % ps is a paramter class
    properties (Constant)
        wav = "C:\Users\Zhehao\Desktop\B464A.wav";
        folder = "Z:\Yazaki-SugiyamaU\Zhehao\My_Stimuli\G606\G606-1st";
        txt = "E:\consongPlx\G606_Z01SO.txt";
        plx  = "E:\consongPlx\G606_Z01SO.plx";
        %T = "E:\DiskDBackUp\ComprehensiveAnalysis\pathTable.xlsx"
        a = "C:\Users\Remote\Documents\G606_Z01SO.txt"
        b = "C:\Users\Remote\Documents\G606_Z01SO.plx"
        c = "C:\Users\Remote\Documents\G606-1st"
        bucket = "Z:\Yazaki-SugiyamaU\Bird-song"
        hp = 450;  % high pass threshold
        SYLLEN = 12800 % threshold for screnning the syllables
        %wavfolders = unique(readtable(pa.T).path_folder);
        t2021 = "C:\Users\Zhehao\Dropbox (OIST)\My_EphysInput\input2021.xlsx"
        TEAOutputFolder = "C:\Users\Zhehao\Dropbox (OIST)\My_Matlab\TEA\Output"
        selected = "C:\Users\Zhehao\Dropbox (OIST)\My_EphysInput\selected.xlsx"
        % this class store the commonly used paramters( ect. pathes) for testing
        % the function
        size3plots = [1986 97 895 672]
        size_wide = [1201 584 1732 407];
        
        
    end
 

end