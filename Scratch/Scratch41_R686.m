
%conspe_eleinf = autoseg("C:\Users\Zhehao\Dropbox (OIST)\My_Stimuli\R686@11022021").singleSegmenter;




ss = stimuli(conspe_eleinf);
ss.setoutdir("E:\Total_stimuli_R686_New");

ss.writenorm;


% generate stimuli degressives
allsongnames = [conspe_eleinf.songname].';
songnames = unique(allsongnames);

for v = 1: length(songnames)
    this_song_ids = find(~cellfun(@isempty, regexp(songnames{v},allsongnames)));
    this_song_inf = conspe_eleinf(this_song_ids);
    degS = stimuli(this_song_inf);
    degS.setoutdir("E:\Total_stimuli_R686_New");
    degS.writedegressive; % write degressive songs
end
