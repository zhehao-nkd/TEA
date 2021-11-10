con_eleinf = autoseg("Z:\Zhehao\Dropbox (OIST)\My_Stimuli\Y687@10252021").singleSegmenter;


all_eleinf = con_eleinf;


for p = 1: length(all_eleinf)
    all_eleinf(p).uniqueid = p;
end

preprocessVrae(all_eleinf);


for k = 1: length(all_eleinf)
    all_eleinf(k).len = length(all_eleinf(k).y);
end