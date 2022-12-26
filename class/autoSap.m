classdef AutoSap
    % 自动运行SAP2011 Automatically run sap2011
    % based on AutoGui
    
    properties
        Property1
    end
    
    methods(Static)

        
        function [birdorder,zporder,sinchorder,bmax,zpmax,sinmax] = getDirOrder(birdname,ZPid)
%             arch = Archon('D:/');
            % bird level
            alldirnames = cellstr(Extract.folder("F:\").');
            % 不知道为什么会出现RECYCLE.BIN 和 System Volume Information
            % 这两个怪folder，不得不去掉
            disp('删除隐藏folders @AutoSap.getDirOrder');
            bad1 = find(~cellfun(@isempty, regexp(alldirnames,'System Volume Information')));
            bad2 = find(~cellfun(@isempty, regexp(alldirnames,'\$')));

            notbad = alldirnames(setdiff([1:length(alldirnames)],[bad1,bad2]));
            birddirs = sortrows(notbad);
            bmax = length(birddirs);
            birdorder = find(~cellfun(@isempty,regexp(birddirs,birdname)));

            
            
            % zplevel
            hitbird_dir = birddirs{birdorder};
            zpdirs = sortrows(cellstr(Extract.folder(hitbird_dir).'));
            zpmax =length(zpdirs);
            zporder = find(~cellfun(@isempty,regexp(zpdirs,ZPid)));
            
            %
            hitzp_dir = zpdirs{zporder};
            withindirs = cellstr(Extract.folder(hitzp_dir).');
            [~,idxwithindirs] = sortrows(lower(withindirs)); % 不区分大小写的排序
            withindirs = withindirs(idxwithindirs);
            sinmax= length(withindirs);
            sinchorder = find(~cellfun(@isempty,regexp(withindirs,'SingleChannel_Merged_Stimuli')));
           
            
            
        end
        
        function coory = getcoory(order, starty, gap)
            coory = starty + (order-1)*gap;
        end
        
        function uni_roster = script_autorun_scp(all_zpdirs)
            %里的每一个神经元，用自动点击器运行计算提取features,提取后的features会被储存在sap2011 mysql
            %数据库中
            % 已知但未解决的bug, 录入的数据名称和实际点开的文件夹有偏差，这种错误存在于Y675
            % Z20,Z21,Z22,Z25----2022.12.10

            dbstop if error
            tic
            au = AutoGui;

            %切换输入法为ENG
            au.move(1783,1177);
            au.click(1783,1177,0.2);
            au.click(1748,1072,0.1);

%             all_possible_dirs = Extract.foldersAllLevel('D:');
%             all_zpdirs = all_possible_dirs(find(~cellfun(@isempty, regexp(all_possible_dirs,'[ZP]\d+$'))));

            for k = 1: length(all_zpdirs) % k was 227
                birdname = regexp(all_zpdirs{k},'[OGBYR]\d{3}','match');
                birdname = birdname{1};
                zpids = regexp(all_zpdirs{k},'[ZP]\d+','match');
                zpids = zpids{1};
                birdname_zpids = [birdname,zpids];

                pause(0.2)
                !C:/Program Files (x86)/SAP2011/SAP2011.exe&
              
                pause(0.9);
                if k == 1
                    pause(2);
                end
                au.click(304,426,0.2);
                au.click(188,108,0.1);
                pause(0.1)
                
                
                au.click(565,100,0.1); % 输入batch name
                au.typewrite(char(birdname_zpids));% 注意输入法必须是英文
                pause(0.1)

                au.click(513,197,0.1); 
                au.click(514,261,0.1); %点击进入 SAP settings 页面
                au.click(243,95,0.1);
                au.click(112,114,0.1);
                au.click(520,433,0.1);
                au.click(505,121,0.1);%点击进入 sound files 页面
                au.click(243,95,0.1);
                
                au.click(735,85,0.1); % click Drive icon
                au.click(703,127,0.1); % click disk F
                au.doubleclick(714,117,0.1);
                
                sharedX = 674;
                bY = 129;
                zpY = 145;
                scY = 160;
                gap = 16;
                dragstart = [745, 136];
                coorstart = [680,128];
                  
                [b_order,zp_order,sc_order,bmax,zpmax,sinmax]= AutoSap.getDirOrder(birdname,zpids);
                [drag1,click1] = au.dragClickInBox(128,416,19,22,40,b_order);
                % firsty: 选择框第一项的y坐标      endy：选择框最末项的y坐标      endnum：选择框最末项的次序
                % totalnum: 所有可选项的数目      dragrange：拖动条从起始到使得最末项位于endy位置时所行使的距离
                % targetnum; 目标项的次序数
                %au.drag(dragfirst,[dragfirst(1),dragfirst(2)- 10],0.1); % 先把进度条拽到最初位置
                au.drag(dragstart,[dragstart(1),dragstart(2)+drag1],0.1); % 拖动到目标位置
                au.doubleclick(coorstart(1),click1,0.1) % 拖动后点击

                     
                
                % click ZP folder
                coorfirst = [637,145];
                coorend = [637,417];
                endnum = 18; totalnum = zpmax;
                dragfirst = [745,131];
%                 dragrange = 105;
                targetnum = zp_order;
                pause(0.15); %截图前停一下，使图像稳定
                scroll_bar = rgb2gray(au.capture([739,124],[752,420]));
                %figure; imshow(scroll_bar);
                linear_fig = mean(scroll_bar.');
                dragrange = (1-find(linear_fig > 220,1)/length( linear_fig))*(420-124);
                [draglength,clicky] = au.dragClickInBox(coorfirst(2),coorend(2),endnum, totalnum,dragrange,targetnum);
                % firsty: 选择框第一项的y坐标      endy：选择框最末项的y坐标      endnum：选择框最末项的次序
                % totalnum: 所有可选项的数目      dragrange：拖动条从起始到使得最末项位于endy位置时所行使的距离
                % targetnum; 目标项的次序数
                %au.drag(dragfirst,[dragfirst(1),dragfirst(2)- 10],0.1); % 先把进度条拽到最初位置
                au.drag(dragfirst,[dragfirst(1),dragfirst(2)+draglength],0.1); % 拖动到目标位置
                au.doubleclick(coorfirst(1),clicky,0.1) % 拖动后点击

                disp([birdname,zpids]);
                pause(0.2);
                dynamic_scy = AutoSap.getcoory(sc_order,  scY, gap);
                au.move(sharedX,dynamic_scy);
                au.doubleclick(sharedX,dynamic_scy,0.1); % double-click sc文件夹
                
                au.click(500,192,0.1)% click next
                au.click(853,111,0.1)%click the sound file
                au.click(500,445,0.1) % click next
                
                au.drag([84,544],[84,648],0.1); % drag threshold
                pause(0.1);
                au.doubleclick(889,110,0.1);
                au.click(516,448,0.1); %点击进入review页面
                
                au.click(488,435,0.1); % run the batch
                
               
                %判断颜色
                while true
                    
                    pause(4)
                    image = au.capture([82,75],[92,85]);
                    mcolor = au.meancolor(image);
                    fprintf('G通道值是：%u\n',mcolor(2));
                    if mcolor(2) > 100; break;end; % 当颜色变绿
                    pause(0.2);
                end
                uni_roster(k).birdname = birdname;
                uni_roster(k).zpids = zpids;
                uni_roster(k).testify_img = au.capture([19,61],[1100,740]); %截一张runbatch的图
                pause(0.1);

                !Taskkill/IM SAP2011.exe >nul&
                !Taskkill/IM cmd.exe >nul&
                
                fprintf('此次循环数是%u,当前用时：%f\n',k,toc);
                
            end
            
            
        end
        
        function export_infofiles(overwrite)
            all_possible_dirs = Extract.foldersAllLevel('F:');
           all_zpdirs = all_possible_dirs(find(~cellfun(@isempty, regexp(all_possible_dirs,'[ZP]\d+$'))));
            %把数据库中储存的feature写成feature文件，存在对应的subfolder中
            if overwrite == 1
                disp('覆写模式');
            else
                disp('跳过模式');
            end
            sourcelist =  listDataSources();
            
            datasource = 'sap';
            username = 'root';
            password = 'sap2011';
            conn = database(datasource,username,password);
            
            tablelist = table2struct(sqlfind(conn,'','Catalog','sap'));
            alltablenames = {tablelist.Table}.';
            
            infonames = alltablenames(~cellfun(@isempty,regexp(alltablenames, 'file_table_[ogbyr]\d{3}[zp]\d+')));
            
            wb1 = waitbar(0,'开始写入Info文件'); % loop the query to export txt
            for k = 1:length(infonames) % info 文件
                
                
                bid = upper(regexp(convertCharsToStrings(infonames{k}),'[ogbyr]\d{3}','match'));
                zp = upper(regexp(convertCharsToStrings(infonames{k}),'[zp]\d+','match'));
                hitbname = find(~cellfun(@isempty, regexp(all_zpdirs,bid)));
                hitzpid = find(~cellfun(@isempty, regexp(all_zpdirs,zp)));
                hitboth = intersect(hitbname, hitzpid);
                zpdir = all_zpdirs{hitboth};
       
                featuredir = fullfile(zpdir,'Feature');
                mkdir(featuredir); % 建立Feature文件夹
                infofilename = sprintf('%s-%s-Info.txt',bid,zp);
                outpath = Convert.path(fullfile(featuredir,infofilename),'unix');
                sqlquery = sprintf(['(SELECT ''bird_ID'',''file_index'',''file_name'',''file_age'',''bird_age'')\n',...
                    'UNION\n',...
                    '(SELECT bird_ID,file_index,file_name,file_age,bird_age\n',...
                    'FROM sap.%s\n',...
                    'INTO OUTFILE ''%s''\n',...
                    'FIELDS ENCLOSED BY ''"'' TERMINATED BY '','' ESCAPED BY ''"''\n',...
                    'LINES TERMINATED BY '',\\r'');']...
                    ,infonames{k},outpath);
                
                if ~isfile(outpath)
                    execute(conn,sqlquery);
                else
                    if overwrite == 1
                        delete(outpath);
                        execute(conn,sqlquery);
                        disp('原文件被覆盖');
                    else
                        disp('此文件已经生成,未被覆写');
                    end
                end
                waitbar(k/length(infonames),wb1,sprintf('已完成全部%u中的第%u个Info文件',length(infonames),k));
            end
            
            
        end
        
        function export_datafiles(overwrite)
              %把数据库中储存的feature写成feature文件，存在对应的subfolder中
              all_possible_dirs = Extract.foldersAllLevel('F:');
              all_zpdirs = all_possible_dirs(find(~cellfun(@isempty, regexp(all_possible_dirs,'[ZP]\d+$'))));
              if overwrite == 1
                  disp('覆写模式');
              else
                  disp('跳过模式');
              end

              sourcelist =  listDataSources();

              datasource = 'sap';
              username = 'root';
              password = 'sap2011';
              conn = database(datasource,username,password);

              tablelist = table2struct(sqlfind(conn,'','Catalog','sap'));
              alltablenames = {tablelist.Table}.';


              datanames = alltablenames(~cellfun(@isempty,regexp(alltablenames, 'raw_[ogbyr]\d{3}[zp]\d+')));
              wb2 = waitbar(0,'开始写入Data文件');
              for k = 1:length(datanames)% data 文件
                  bid = upper(regexp(convertCharsToStrings(datanames{k}),'[ogbyr]\d{3}','match'));
                  zp = upper(regexp(convertCharsToStrings(datanames{k}),'[zp]\d+','match'));
                  hitbname = find(~cellfun(@isempty, regexp(all_zpdirs,bid)));
                  hitzpid = find(~cellfun(@isempty, regexp(all_zpdirs,zp)));
                  hitboth = intersect(hitbname, hitzpid);
                  zpdir = all_zpdirs{hitboth};
   
                
                featuredir = fullfile(zpdir,'Feature');
                mkdir(featuredir); % 建立Feature文件夹
                
                infofilename = sprintf('%s-%s-Data.txt',bid,zp);
                outpath = Convert.path(fullfile(featuredir,infofilename),'unix');
                sqlquery = sprintf(...
                    ['(SELECT ''time'',''file_index'',''amplitude'',''mean_frequency_amp'',''pitch'',''mean_frequency'',',...
                    '''FM'',''am'',''goodness'',''entropy'',''peak_frequency'',''DAS'',''continuity_t'',''continuity_f'')\n',...
                    'UNION\n',...
                    '(SELECT time,file_index,amplitude,mean_frequency_amp,pitch,mean_frequency,',...
                    'FM,am,goodness,entropy,peak_frequency,DAS,continuity_t,continuity_f\n',...
                    'FROM sap.%s',...
                    ' INTO OUTFILE ''%s''\n',...
                    'FIELDS ENCLOSED BY ''"'' TERMINATED BY '','' ESCAPED BY ''"''\n',...
                    'LINES TERMINATED BY '',\\r'');']...
                    ,datanames{k},outpath);
                
                if ~isfile(outpath)
                    execute(conn,sqlquery);
                else
                    
                    if overwrite == 1
                        delete(outpath);
                        execute(conn,sqlquery);
                        disp('原文件被覆盖');
                    else
                        disp('此文件已经生成,未被覆写');
                    end
                    
                end
                waitbar(k/length(datanames),wb2,sprintf('已完成全部%u中的第%u个Data文件',length(datanames),k));
            end
            
            
        end
        
        function write_sapscreenshot(unilist)
            
            %unilist = AutoSap.script_autorun_scp(neuroster);
            
            wb = waitbar(0,'Starting');
            for k = 1: length(unilist)
                waitbar(k/length(unilist),wb,sprintf('%u of %u',k,length(unilist)));
                disp(k)
                imwrite(unilist(k).testify_img,sprintf('%s%s.png',unilist(k).birdname,unilist(k).zpids));
            end
            
        end
        
    end
end

