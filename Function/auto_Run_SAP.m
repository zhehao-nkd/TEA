dbstop if error
%[neuroster,malfunctions,cdrfroster,not_cdrfroster,nofeature_cdrfroster] = Archon.extractAnalysisInfo('D:\');
% noscdroster = neuroster([neuroster.SCDexist].' == 0);
% Archon.batch_createSinChFolder(noscdroster); % generate single channel folder for each neurons in roster
% 

% nf_roster = neuroster([neuroster.featureexist].' == 0);
% [unqiebirdneuronpair,ia,~] = unique([nf_roster.birdneuron].');
% uni_nf_roster = nf_roster(ia);
tic

for k = 1: length(uni_nf_roster) % k was 227
    
    pause(0.3)
    !C:/Program Files (x86)/SAP2011/SAP2011.exe&
    au = AutoGui;
    pause(0.8);
    au.click(304,426,0.2);
    au.click(188,108,0.1);
    
    
    au.click(565,100,0.1);
    au.typewrite(char(uni_nf_roster(k).birdneuron));
    au.click(513,197,0.1); % 注意输入法必须是英文
    au.click(514,261,0.1);
    au.click(243,95,0.1);
    au.click(112,114,0.1);
    au.click(520,433,0.1);
    au.click(505,121,0.1);
    
    au.click(735,85,0.1); % click Drive icon
    au.click(684,110,0.1);
    au.doubleclick(714,117,0.1);
    
    sharedX = 674;
    bY = 129;
    zpY = 145;
    scY = 160;
    gap = 16;
    
    [b_order,zp_order,sc_order] = AutoSap.getDirOrder(uni_nf_roster(k).birdname,uni_nf_roster(k).neuronid);
    
    dynamic_by = AutoSap.getcoory(b_order,  bY, gap);
    au.move(sharedX,dynamic_by);
    au.doubleclick(sharedX,dynamic_by,0.1); % double-click 鸟文件夹
    
    dynamic_zpy = AutoSap.getcoory(zp_order,  zpY, gap);
    au.move(sharedX,dynamic_zpy);
    au.doubleclick(sharedX,dynamic_zpy,0.1); % double-click ZP文件夹
    
    dynamic_scy = AutoSap.getcoory(sc_order,  scY, gap);
    au.move(sharedX,dynamic_scy);
    au.doubleclick(sharedX,dynamic_scy,0.1); % double-click ZP文件夹
    
    au.click(500,192,0.1)% click next
    au.click(853,111,0.1)%click the sound file
    au.click(500,445,0.1) % click next
    
    au.drag([84,544],[84,648],0.1); % drag threshold
    pause(0.1)
    au.doubleclick(889,110,0.1);
    au.click(516,448,0.1);
    
    au.click(488,435,0.1); % run the batch
    
    
    image = au.capture([82,75],[92,85]);
    
    %判断颜色
    while true
        
        pause(4)
        image = au.capture([82,75],[92,85]);
        mcolor = au.meancolor(image);
        fprintf('G通道值是：%u\n',mcolor(2));
        if mcolor(2) > 100 % 当颜色变绿
            break
        end
        pause(0.2);
    end
    
    
    !Taskkill/IM SAP2011.exe >nul&
    !Taskkill/IM cmd.exe >nul&

     fprintf('此次循环数是%u,当前用时：%f\n',k,toc);
    
end


