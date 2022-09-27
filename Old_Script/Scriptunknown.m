%1----- new script to generate shared senatus data
dbstop if error

% % % Generate con & spe eleinf
% % senatus_dir = "E:\Stimuli_Source\senatusTwoMotif";
% % as = autoseg(senatus_dir);
% % as.standard;
% % raw_con_spe_eleinf = getInf.Eleinf(senatus_dir ,1,'SegData');
% con_spe_eleinf = categoFrags(raw_con_spe_eleinf).eachsong;
% 
% 
% % 2----- Generate OTE eleinf
% % 
% % CONSPE_dir = "E:\Stimuli_Source\senatusTwoMotif";
% % % 下面的部分直接生成stimuli set
% % TARGET_DIR = "E:\G656_ELE";
% % mkdir(TARGET_DIR);ConDirs = fromWavFindParentFolder(CONSPE_dir,"E:\Stimuli_Source\allBirdsSong");
% %  wav_dir = "E:\Stimuli_Source\allBirdsSong";
% %  subdirs = Extract.folder(wav_dir);
% %  [~,locb] = ismember([ConDirs{:}].',[subdirs{:}].');
% % noncon_dir_whole_collection = subdirs(setdiff([1: length(subdirs)],locb(locb>0)).').';
% 
% 
% % noncon_ote_eleinf = getInf.Eleinf(noncon_dir_whole_collection ,300,'SegData'); % 300 means remove folderes with bird id smaller than 300
% % tic
% % ote_eleinf = categoFrags(noncon_ote_eleinf).eachsong;
% % toc
% % unique_ote_eleinf = stimuli.unique_it( ote_eleinf );
% % remove too long or too short elements
% wbar = waitbar(0,'Removing...');
% for k = 1: length(unique_ote_eleinf )
%     unique_ote_eleinf(k).len = length(unique_ote_eleinf (k).y);
%     waitbar(k/length(unique_ote_eleinf ),wbar,'Removing');
% end
% close(wbar);
% unique_ote_eleinf = unique_ote_eleinf(find(8000 >[unique_ote_eleinf.len].'& [unique_ote_eleinf.len].'> 800));
% unique_ote_eleinf = rmfield(unique_ote_eleinf, 'len');
% %split ote_eleinf to non-con and con eleinf
% 
% 
% %4------calculate TUT-BOS de eleinf
% % tutbos_dir = "E:\Stimuli_Source\R677_twoMotif";
% % as = autoseg(tutbos_dir);
% % as.standard;
% tutbos_eleinf = getInf.Eleinf(tutbos_dir,1,'SegData');
% tutbos_eleinf  = categoFrags(tutbos_eleinf).eachsong;
% 
% 
% %3------merge the Cons and Non COns together
% temp = horzcat(con_spe_eleinf,unique_ote_eleinf);
% all_eleinf = horzcat(temp,tutbos_eleinf);
% for p = 1: length(all_eleinf)
%     all_eleinf(p).uniqueid = p;
% end
% preprocessVrae(all_eleinf); % this code will write labeled_eleConOnly.mat
% 
% 
% 
% 
% % Here I should load the coor_Z trained by the neural network
% for p = 1: length(all_eleinf)
%     all_eleinf(p).coor_1 = coor_Z(p,1);
%     all_eleinf(p).coor_2 = coor_Z(p,2);
% end
% 
% %

TARGET_DIR = "E:\R677_ELE_realDonaldTrump";
conS = stimuli(con_spe_eleinf);
conS.setoutdir(TARGET_DIR);
conS.writenorm;
conS.writeEachSongFrag;
conS.writeEachSongFragInOneFolder;

% generate stimuli degressives
allsongnames = [con_spe_eleinf.songname].';
songnames = unique(allsongnames);

for v = 1: length(songnames)
    this_song_ids = find(~cellfun(@isempty, regexp(songnames{v},allsongnames)));
    this_song_inf = con_spe_eleinf(this_song_ids);
    degS = stimuli(this_song_inf);
    degS.setoutdir(TARGET_DIR);
    degS.writedegressive; % write degressive songs
end


% % generate similar-but-different-motif fragments
% Trump = stimuli(con_merged_eleinf);
% Trump.setoutdir(TARGET_DIR);
% Trump.writeSimilarFragsButDifferentMotif;

% generate detailed stimuli set
tic

con_ids = find(~cellfun(@isempty,regexp([all_eleinf(:).songname].','CON|SPE')));

ss = stimuli(all_eleinf);
disp('王勃')
ss.setoutdir(TARGET_DIR);
ss.set_far(35);
ss.writeFrag_far_from_senatus;
disp('佩玉鸣鸾罢歌舞');
first_ids = find([ss.prepro.fragid].' == 1); % 意思是每首歌的起始的index

for m = con_ids(:).' % for each element of CONs, generate the corresponding stimuli
    
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

toc
send_mail_message('zhehao.cheng@oist.jp','Done','finished')
send_mail_message('379380788@qq.com','Done','finished')
send_mail_message('chengzh.nkd@gmail.com','Done','finished')