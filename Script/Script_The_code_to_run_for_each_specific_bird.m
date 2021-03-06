% The code to run for each specific bird id

% %4------calculate TUT-BOS de eleinf
tutbos_dir = "E:\Stimuli_Source\R690_twoMotif";
% as = autoseg(tutbos_dir);
% as.standard;
tutbos_eleinf = getInf.Eleinf(tutbos_dir,1,'SegData');
tutbos_eleinf  = categoFrags(tutbos_eleinf).eachsong;
% 
% TargetDir_specific_bird = "E:\R690_Special";
% 
% S_of_this_specific_bird = stimuli(tutbos_eleinf);
% S_of_this_specific_bird.setoutdir(TargetDir_specific_bird);
% S_of_this_specific_bird.writenorm;
% S_of_this_specific_bird.writeEachSongFrag;
% S_of_this_specific_bird.writeEachSongFragInOneFolder;
% 
% % generate stimuli degressives
% tutbossongnames = [tutbos_eleinf.songname].';
% songnames = unique(tutbossongnames);
% 
% for v = 1: length(songnames)
%     this_song_ids = find(~cellfun(@isempty, regexp(songnames{v},tutbossongnames)));
%     this_song_inf = tutbos_eleinf(this_song_ids);
%     specific_degS = stimuli(this_song_inf);
%     specific_degS.setoutdir(TargetDir_specific_bird);
%     specific_degS.writedegressive; % write degressive songs
% end
% 
% % generate stimuli details
% 
% 
% save(sprintf('tutbos_dir\%s.mat', convert.bid(tutbos_dir)),tutbos_dir);
% 
% % write fromat for coordinate analysis
% 
% for p = 1: length(tutbos_eleinf)
%     tutbos_eleinf(p).uniqueid = p;
% end
% preprocessVrae(tutbos_eleinf); % this code will write labeled_eleConOnly.mat
% 
% 
% 
% 
% 
% 
% 
% % Here I should load the coor_Z trained by the neural network
% for p = 1: length(tutbos_eleinf)
%     tutbos_eleinf(p).coor_1 = coor_Z(p,1);
%     tutbos_eleinf(p).coor_2 = coor_Z(p,2);
% end
% 
% all_eleinf = horzcat(all_eleinf,tutbos_eleinf);
% 
% for p = 1: length(all_eleinf)
%     all_eleinf(p).uniqueid = p;
% end



% generate detailed stimuli set
tic

TARGET_DIR = "E:\R690_Special"
tutbos_ids = find(~cellfun(@isempty,regexp([all_eleinf(:).songname].','TUT|BOS')));

ss = stimuli(all_eleinf);
disp('王勃')
ss.setoutdir(TARGET_DIR);
ss.set_far(45);
% ss.writeFrag_far_from_senatus;
disp('佩玉鸣鸾罢歌舞');
first_ids = find([ss.prepro.fragid].' == 1); % 意思是每首歌的起始的index

for m = tutbos_ids(:).' % for each element of CONs, generate the corresponding stimuli
    
    disp(m);
    
    ss.writeFrag_samesong(m); %
    
%     ss.set_near(15);
%     ss.writeFrag_near_from_senatus(m);
%     disp('画栋朝飞南浦云');
    
    %num_all_near = ss.from_senatus_near_find_all_near(ss.num_of_near,m);
    ss.set_near(20);
    ss.writeFrag_near_from_all(m);
    disp('珠帘暮卷西山雨');
       
    ss.set_far(30);
    ss.writeRepla_far_from_senatus(m);
    disp('物换星移几度秋');
    
    if ~ismember(m,first_ids)
        
%         ss.set_near(12);
%         ss.writeRepla_near_from_senatus(m);
%         disp('阁中帝子今何在');
        
        ss.set_near(15);
       % num_all_near = ss.from_senatus_near_find_all_near(ss.num_of_near,m-1);
        ss.writeRepla_near_from_all(m);
        disp('槛外长江空自流');
        
    end
    
    disp('滕王阁');
end


