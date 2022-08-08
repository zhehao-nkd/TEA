% new stimuli generatuion code


dbstop if error

% conspe_eleinf = autoseg("C:\Users\Zhehao\Dropbox (OIST)\My_Stimuli\Y687@10252021").singleSegmenter;
conS = stimuli(conspe_eleinf);
conS.setoutdir("E:\Total_stimuli_Y687_New");
conS.writenorm;

% generate stimuli degressives
allsongnames = [conspe_eleinf.songname].';
songnames = unique(allsongnames);

for v = 1: length(songnames)
    this_song_ids = find(~cellfun(@isempty, regexp(songnames{v},allsongnames)));
    this_song_inf = conspe_eleinf(this_song_ids);
    degS = stimuli(this_song_inf);
    degS.setoutdir("E:\Total_stimuli_Y687_New");
    degS.writedegressive; % write degressive songs
end

tic

con_ids = find(~cellfun(@isempty,regexp([all_eleinf(:).songname].','CON|SPE')));
% ss = stimuli(all_eleinf);
% disp('岳阳楼')
% ss.setoutdir("E:\Total_stimuli_Y687");
% ss.set_far(20);
% ss.writeFrag_far_from_all;
% disp('滕王高阁临江渚');
ss.setoutdir("E:\Total_stimuli_Y687_New");
ss.set_far(45);
ss.writeFrag_far_from_senatus;
disp('佩玉鸣鸾罢歌舞');


first_ids = find([ss.prepro.fragid].' == 1);
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
    
    %
    %     ss.set_far(15);
    %     ss.writeRepla_far_from_all(m);
    %     disp('闲云潭影日悠悠');
    
    ss.set_far(36); % Repla FarfromSenatus
    ss.writeRepla_far_from_senatus(m);
    disp('物换星移几度秋');
    
    if ~ismember(m,first_ids)
        
        ss.set_near(12); % Repla NearfromSenatus
        ss.writeRepla_near_from_senatus(m);
        disp('阁中帝子今何在');
        
        ss.set_near(12); % Repla Near from All
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