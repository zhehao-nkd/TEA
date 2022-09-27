sap = table2struct(readtable("C:\Users\Zhehao\Dropbox (OIST)\My_Luscinia\sap1276.xls"));

wavfolder = 'Z:\Yazaki-SugiyamaU\Zhehao\My_Stimuli\G606\G606-1st';
wavs = Extract.filename(wavfolder,'*.wav');

nn = 0;

for mm = 1: length(wavs)
    thiswav = wavs{mm};
    [y ,fs] = audioread(thiswav);
    y = y(:,2);
    [~,b,c] = fileparts(thiswav);
    temp = strsplit(b,'_');   % this modification is dangerous  !!!!!!
    b = temp{1};
    
   
    
    thisfile = sprintf('%s',b);
    
    thisidx = find(~cellfun(@isempty,regexp({sap(:).fileName},thisfile)));
    thissap = sap(thisidx);
    
    
    for idx = 1: length(thissap)
        nn = nn + 1;
        frag(nn).initial = thissap(idx).syllable_start/1000*fs;
        frag(nn).terminal = (thissap(idx).syllable_start+ thissap(idx).syllable_duration )/1000*fs;
        frag(nn).y = y(frag(nn).initial:frag(nn).terminal);
        frag(nn).wholey = y;
        frag(nn).songid = b;
        frag(nn).sylnum = idx;
 
    end
    
end


% for ee = 1: length(frag)
%     judgeclass(frag{ee});
% end


