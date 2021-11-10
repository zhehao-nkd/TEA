dbstop if error

%conspe_eleinf = autoseg("C:\Users\Zhehao\Dropbox (OIST)\My_Stimuli\Y687@10252021").singleSegmenter;
%conS = stimuli(conspe_eleinf);
% conS.setoutdir("E:\Total_stimuli_Y687");
% conS.writenorm;
% 
% % generate stimuli degressives
% allsongnames = [conspe_eleinf.songname].';
% songnames = unique(allsongnames);
% 
% for v = 1: length(songnames)
%     this_song_ids = find(~cellfun(@isempty, regexp(songnames{v},allsongnames)));
%     this_song_inf = conspe_eleinf(this_song_ids);
%     degS = stimuli(this_song_inf);
%     degS.setoutdir("E:\Total_stimuli_Y687");
%     degS.writedegressive; % write degressive songs
% end


% con_ids = find(~cellfun(@isempty,regexp([all_eleinf(:).songname].','CON|SPE')));
% ss = stimuli(all_eleinf);
% disp('岳阳楼')
% ss.setoutdir("E:\Total_stimuli_Y687");
% ss.set_far(20);
% ss.writeFrag_far_from_all;
% disp('滕王高阁临江渚');
ss.set_far(45);
ss.writeFrag_far_from_senatus;
disp('佩玉鸣鸾罢歌舞');


first_ids = find([ss.prepro.fragid].' == 1);
parfor m = con_ids(:).' % for each element of CONs, generate the corresponding stimuli
    
    disp(m);
    ss.set_near(15);
    ss.writeFrag_near_from_all(m);
    disp('画栋朝飞南浦云');
    
    ss.set_near(25);
    ss.writeFrag_near_from_senatus(m);
    disp('珠帘暮卷西山雨');
%     
%     ss.set_far(15);
%     ss.writeRepla_far_from_all(m);
%     disp('闲云潭影日悠悠');
    
    ss.set_far(25);
    ss.writeRepla_far_from_senatus(m);
    disp('物换星移几度秋');
    
    if ~ismember(m,first_ids)
        
        ss.set_near(15);
        ss.writeRepla_near_from_all(m);
        disp('阁中帝子今何在');
        
        ss.set_near(25);
        ss.writeRepla_near_from_senatus(m);
        disp('槛外长江空自流');
    end
    
    disp('滕王阁');
end

send_mail_message('zhehao.cheng@oist.jp','Done','finished')
send_mail_message('379380788@qq.com','Done','finished')
send_mail_message('chengzh.nkd@gmail.com','Done','finished')