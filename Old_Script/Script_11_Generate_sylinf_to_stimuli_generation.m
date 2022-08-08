
ConDirs = fromWavFindParentFolder("C:\Users\Zhehao\Dropbox (OIST)\My_Stimuli\O684@test@20211019","E:\WavsCollection");
 wav_dir = "E:\WavsCollection";
 subdirs = extract.folder(wav_dir);
 [Lia,locb] = ismember([ConDirs{:}].',[subdirs{:}].');


%(((^~^))) Convert segdata files into eleinf, and then convert eleinf into
%the format suitable for the python script

noncon_dir_whole_collection = subdirs(setdiff([1: length(subdirs)],locb(locb>0)).');
noncon_ote_eleinf = buildFraginf.assemble_as_sylinf(noncon_dir_whole_collection ,300,'SegData'); % 300 means remove folderes with bird id smaller than 300
ote_eleinf = categoFrags(noncon_ote_eleinf).eachsong;
unique_ote_eleinf = stimuli.unique_it( ote_eleinf );
% remove too long or too short elements
for k = 1: length(unique_ote_eleinf )
    unique_ote_eleinf(k).len = length(unique_ote_eleinf (k).y);
end
%T_unique_ote_eleinf = struct2table(unique_ote_eleinf);
rm_eleinf = unique_ote_eleinf(find(8000 >[unique_ote_eleinf.len].'& [unique_ote_eleinf.len].'> 800));
%split ote_eleinf to non-con and con eleinf


% 这里的目的是为了让con-ELEINF的运算流程和ote-eleinf的运算流程一样
as = autoseg("E:\Selected_wavs");
as.standard;
con_spe_eleinf = getInf.Sylinf("E:\Selected_wavs" ,1,'SylData');

parent_dir_con = subdirs(locb(locb>0));
con_ote_eleinf = getInf.Sylinf(parent_dir_con ,1,'SegData');

con_merged_eleinf = categoFrags(con_ote_eleinf).include_con(con_spe_eleinf); % add catego info into the eleinf


% merge the Cons and Non COns together
all_eleinf = horzcat(con_spe_eleinf,rm_eleinf);
for p = 1: length(all_eleinf)
    all_eleinf(p).uniqueid = p;
end
preprocessVrae(all_eleinf); % this code will write labeled_eleConOnly.mat


% Here I should load the coor_Z trained by the neural network
for p = 1: length(all_eleinf)
    all_eleinf(p).coor_1 = coor_Z(p,1);
    all_eleinf(p).coor_2 = coor_Z(p,2);
end



% 下面的部分直接生成stimuli set
TARGET_DIR = "E:\O686_REAL_TMD_NEW";
mkdir(TARGET_DIR);

conS = stimuli(con_spe_eleinf);
conS.setoutdir(TARGET_DIR);
conS.writenorm;

% generate stimuli degressives
allsongnames = [conspe_eleinf.songname].';
songnames = unique(allsongnames);

for v = 1: length(songnames)
    this_song_ids = find(~cellfun(@isempty, regexp(songnames{v},allsongnames)));
    this_song_inf = conspe_eleinf(this_song_ids);
    degS = stimuli(this_song_inf);
    degS.setoutdir(TARGET_DIR);
    degS.writedegressive; % write degressive songs
end


% generate similar-but-different-motif fragments
Trump = stimuli(con_merged_eleinf);
Trump.writeSimilarFragsButDifferentMotif

% generate detailed stimuli set
tic

con_ids = find(~cellfun(@isempty,regexp([all_eleinf(:).songname].','CON|SPE')));

ss = stimuli(all_eleinf);
disp('岳阳楼')
ss.setoutdir(TARGET_DIR);
ss.set_far(45);
ss.writeFrag_far_from_senatus;
disp('佩玉鸣鸾罢歌舞');
first_ids = find([ss.prepro.fragid].' == 1); % 意思是每首歌的起始的index

for m = con_ids(:).' % for each element of CONs, generate the corresponding stimuli
    
    disp(m);
    
    ss.writeFrag_samesong(m); %
    
    ss.set_near(15);
    ss.writeFrag_near_from_senatus(m);
    disp('画栋朝飞南浦云');
    
    num_all_near = ss.from_senatus_near_find_all_near(ss.num_of_near,m);
    ss.set_near(20);
    ss.with_sampling_writeFrag_near_from_all(m,num_all_near);
    disp('珠帘暮卷西山雨');
       
    ss.set_far(36);
    ss.writeRepla_far_from_senatus(m);
    disp('物换星移几度秋');
    
    if ~ismember(m,first_ids)
        
        ss.set_near(12);
        ss.writeRepla_near_from_senatus(m);
        disp('阁中帝子今何在');
        
        ss.set_near(12);
        num_all_near = ss.from_senatus_near_find_all_near(ss.num_of_near,m-1);
        ss.with_sampling_writeRepla_near_from_all(m,num_all_near);
        disp('槛外长江空自流');
        
    end
    
    disp('滕王阁');
end

toc
send_mail_message('zhehao.cheng@oist.jp','Done','finished')
send_mail_message('379380788@qq.com','Done','finished')
send_mail_message('chengzh.nkd@gmail.com','Done','finished')



