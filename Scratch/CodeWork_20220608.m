% 20220608 Codework
addpath(genpath("C:\Users\Zhehao\Dropbox (OIST)\My_Matlab\TEA"))


%[neuroster,malfunctions,cdrfroster,not_cdrfroster,nofeature_neuroster] = Archon.extractAnalysisInfo('D:\');
%Archon.batch_createSinChFolder(nofeature_neuroster)
% auto-run_sap
dbstop if error
% 
%unilist = autoSap.script_autorun_scp(neuroster);
% % 
% wb = waitbar(0,'Starting');
% for k = 1: length(unilist)
%     waitbar(k/length(unilist),wb,sprintf('%u of %u',k,length(unilist)));
%     disp(k)
%     imwrite(unilist(k).testify_img,sprintf('%s.png',unilist(k).birdneuron));
% end
% 
% % info/data files
% 
% autoSap.export_infofiles(0);
 %autoSap.export_datafiles(0);


% run Analysis for all neuroster


error_neuroster = Archon.batch_genAnalysis(neuroster);
%sul = Sultan("C:\Users\Zhehao\Downloads");
sul = Sultan("C:\Users\Zhehao\Downloads\New0616");
%sul = Sultan("C:\Users\Zhehao\Downloads\BadSignalNeurons");
sul.allSongallNeurons_FolderOrder;

Sultan.applySameFunctionForAll("C:\Users\Zhehao\Downloads\OtherNeurons");



% new model of autogui

% au = autogui;
% 
% draggap = (210 -125)/(26-18);
% 
% for k = 1:8
%     
%     au.click(741,85,0.2)
%     au.click(712,103,0.2)
%     au.click(741,85,0.2)
%     au.click(708,113,0.2) % 点击到birdid level
%     
%     au.doubleclick(692,274,0.5)
%     au.drag([747,123],[747,123+k*draggap],0.5)
%     au.click(668,418,0.8)
%     
% end


% coorfirst = [637,145];
% coorend = [637,417];
% endnum = 18; totalnum = 26;
% dragfirst = [745,131];
% dragrange = 85;
% 
% targetnum = 25;
% 
% 
% [draglength,clicky] = au.dragClickInBox(coorfirst(2),coorend(2),endnum, totalnum,dragrange,targetnum);
%             % firsty: 选择框第一项的y坐标      endy：选择框最末项的y坐标      endnum：选择框最末项的次序
%             % totalnum: 所有可选项的数目      dragrange：拖动条从起始到使得最末项位于endy位置时所行使的距离
%             % targetnum; 目标项的次序数
% au.drag(dragfirst,[dragfirst(1),dragfirst(2)- 20],0.1); % 先把进度条拽到最初位置
% 
% au.drag(dragfirst,[dragfirst(1),dragfirst(2)+draglength]);
% au.move(coorfirst(1),clicky)