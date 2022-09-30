% Special Stimuli Generation code for G649
addpath(genpath("C:\Users\Zhehao\Dropbox (OIST)\My_Matlab\TEA"))


TargetDir_specific_bird = "E:\Stimuli_Synthesis_Source\TheSource\B733_twoMotif";
TARGET_DIR = "E:\SpecialDepot\B733_Special"


%第一步  自动分割
tutbos_dir = "E:\StimuliSource\O737_twoMotif";
as = autoseg(tutbos_dir);
as.standard;

%第二部 检视 手动分割
dirpath = "E:\StimuliSource\O737_twoMotif\SegData";%"E:\Stimuli_Source\senatusTwoMotif\SegData";
segfiles = Extract.filename(dirpath,'*.mat');
for k = 1: length(segfiles)   
    reviseSeg(segfiles{k});
    uiwait(gcf);
end


% 第三步 读取为 eleinf
tutbos_eleinf = MetaStimuli.Eleinf(tutbos_dir,1,'SegData');
tutbos_eleinf  = categoFrags(tutbos_eleinf).eachsong;


% 第四步 生成 norm eachsongfrag
TargetDir_specific_bird = "E:\SpecialDepot\O737_Special";
mkdir(TargetDir_specific_bird);

S_of_this_specific_bird = stimuli(tutbos_eleinf);
S_of_this_specific_bird.setoutdir(TargetDir_specific_bird);
S_of_this_specific_bird.writenorm;
S_of_this_specific_bird.writeEachSongFrag;
S_of_this_specific_bird.writeEachSongFragInOneFolder;
S_of_this_specific_bird.writereverses;
S_of_this_specific_bird.writemirrors;

%第五步 生成 degressive songs

tutbossongnames = [tutbos_eleinf.songname].';
songnames = unique(tutbossongnames);

for v = 1: length(songnames)
    this_song_ids = find(~cellfun(@isempty, regexp(songnames{v},tutbossongnames)));
    this_song_inf = tutbos_eleinf(this_song_ids);
    specific_degS = stimuli(this_song_inf);
    specific_degS.setoutdir(TargetDir_specific_bird);
    specific_degS.writedegressive; % write degressive songs
end

% 第六步 生成 detailed stimulli

% save(sprintf('tutbos_dir\%s.mat', Convert.bid(tutbos_dir)),tutbos_dir);
% % write fromat for coordinate analysis
for p = 1: length(tutbos_eleinf)
    tutbos_eleinf(p).uniqueid = p;
end
preprocessVrae(tutbos_eleinf); % this code will write labeled_eleConOnly.mat
% 第七步 用python生成坐标 读取 Coordinate_ELE
for p = 1: length(tutbos_eleinf)
    tutbos_eleinf(p).coor_1 = coor_Z(p,1);
    tutbos_eleinf(p).coor_2 = coor_Z(p,2);
end

% 第八步 读取 all_eleinf 并且和 tutbos—eleinf 合并 还要删掉 tutbos-eleinf的uniqueid 这个 field
all_eleinf = horzcat(all_eleinf,tutbos_eleinf);
for p = 1: length(all_eleinf)
    all_eleinf(p).uniqueid = p;
end

% 第九步 即最后一步
% generate detailed stimuli set
tic


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

