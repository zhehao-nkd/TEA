addpath(genpath('C:\Users\Zhehao\Dropbox (OIST)\My_Matlab\TEA'))

segment_workstate = 0; % 如果不启用新的CONs set就设为0
classify_workstate = 0;

% 第一部分

if segment_workstate == 1
    
    app = manSeg("C:\Users\Zhehao\Dropbox (OIST)\My_Stimuli\Y661@06282021");  % manSeg special songs
    waitfor(app,'outputstate',1);
    %output
    disp('Output acquired')
    elepart2 = output{3};
    sylpart2 = output{2}; 
    
    load("C:\Users\Zhehao\Dropbox (OIST)\My_Stimuli\SourceSongSet1@06232021\eleinf.mat");
    load("C:\Users\Zhehao\Dropbox (OIST)\My_Stimuli\SourceSongSet1@06232021\sylinf.mat");
    
    mergedeleinf = horzcat(eleinf, elepart2);
    mergedsylinf = horzcat(sylinf,sylpart2);
    
    
    disp('stupid section1 has been finished')
  
    %manClus(eleinf);
end


load("C:\Users\Zhehao\Dropbox (OIST)\My_Stimuli\Y661@06282021\mergedeleinf.mat");
load("C:\Users\Zhehao\Dropbox (OIST)\My_Stimuli\Y661@06282021\mergedsylinf.mat");
mergedeleinf = stimuli.adduniqueid(mergedeleinf); % add unique id for each syllable/note
mergedsylinf = stimuli.adduniqueid(mergedsylinf);

stimuli(mergedsylinf).writenorm2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlsxdir = "C:\Users\Zhehao\Dropbox (OIST)\My_EphysInput\input2021.xlsx";
responsiveinfo = manPickEli(xlsxdir,23,mergedeleinf,mergedsylinf);
%%% after processing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 人工选择一个 response eliciting syllables 的 idx


responsiveidx = 88


mansim = manFindSim(mergedeleinf, responsiveidx);
waitfor(mansim,'outputstate');
    %output
disp('Output acquired')

picked

% write picked ele into frags
stimuli(picked).writeintofrags;

%%% write single syllables of the whole songs from picked into .wav files
pickedsongs = cellstr({picked.songname}.');
eleidxfrompickedsongs = ismember(cellstr({mergedeleinf(:).songname}.'),pickedsongs);
pickedwholeeles = mergedeleinf(eleidxfrompickedsongs);

stimuli(pickedwholeeles2).writeintofrags;



% write replaced songs
responsiveidx = 88;
notresponsiveidx = 316;

stimuli(mergedeleinf).writereplace(responsiveidx,notresponsiveidx);
pause(2)

stimuli(mergedeleinf).writereplace(notresponsiveidx,responsiveidx); % 换位再来一下


% degressive 

stimuli(pickedwholeeles2).writedegressive1(8,19)
stimuli(pickedwholeeles3).writedegressive1(6,20)
