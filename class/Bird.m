
classdef Bird < handle
    % A class for   (1) extracting bird songs from the bucket surver, birdsong folders
    %               (2) finding candidate birds for surgery.

    properties
        birdlist %读取的buckect里叫做birdlist的表格
        songdirs % Birdsong directories

        folders
        adultfolders % 有正常的adult song的folders
    end
    % from the targeted folder to Extract syllables
    methods

        function bd = Bird(~)

            bd.birdlist = Bird.readBirdlist;

            fdir = "Z:\Yazaki-SugiyamaU\Bird-song";
            bd.songdirs = cellstr(Extract.folder(fdir).').';


        end
 
        function answer = findSons(bd,fathername,type)
            %找到父鸟所有的子鸟
            % configs for reading the table
            %birdlist = bd.readBirdlist;


            %如果鸟的名字是这种形式  "r-97 (O218)"，那么这名字应该被转化为O218这种格式

            converted_fathername = regexp(convertCharsToStrings(fathername),'[OGBYR]\d{3}','match');


            children_index = find(~cellfun(@isempty, regexp(bd.birdlist.father,converted_fathername)));
            male_index = find(~cellfun(@isempty, regexp(bd.birdlist.Gender,"♂")));
            sons_index = intersect(children_index, male_index);
            sons_names = bd.birdlist.BirdID(sons_index);
            hatch_dates =bd.birdlist.Hatchdate(sons_index);


            summer = {};
            for k = 1: length(sons_names)
                summer{k} = bd.getRecordingSate(sons_names{k}); % recording state被分为[-1,0,0.5,1,1.5]这五种情况，0.5意思是只有juvenile song
                fprintf('His song is:  %s, hatching date: %s, recording state: %.1f\n',sons_names{k},hatch_dates{k},summer{k}.identifier);
                newline;
            end

            answer = vertcat(summer{:});%default结果

            %如果只是找recorded songs,覆盖default结果
            if exist('type','var') && ~isempty(answer)
                if strcmp(type, 'recorded')
                    answer = answer(find([answer.identifier].' >=1));
                end
            end


        end


        function answer = getRecordingSate(bd,birdname)
            %判断输入的鸟ID是否有对应的"记录好的Adult Bird Song"
            % Output: answer.birdname, answer.identifier, answer.dirpath

            answer.birdname = birdname;
            folders = bd.songdirs;
            foldernames = cellfun(@Utl.fileparts,folders,'Uni',0);

            color = regexp(birdname,'[OGBYR]','match');
            number = regexp(birdname,'\d+','match');
            colorids = find(~cellfun(@isempty, regexp(foldernames,color{1})));
            numberids = find(~cellfun(@isempty, regexp(foldernames,number{1})));
            hitindex = intersect(colorids,numberids); %如果压根没这个鸟对应的folder，recording state就是0

            if isempty(hitindex)
                answer.identifier = 0; %没有记录
                answer.dirpath = [];
                return
            elseif length(hitindex) >1 %如果有多个hitfolder,取最晚创建的folder
%                 date_created = {};
                datetimes = NaT(length(hitindex),1);
                for k = 1:length(hitindex)
                    datetimes(k) = datetime(py.os.path.getctime(folders{hitindex(k)}),...
                        'ConvertFrom','epochtime','TicksPerSecond',1,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
                end

                [~,maxindex] = max(datetimes);
                hitindex = hitindex(maxindex);
            end


            foldername = foldernames{hitindex}; % 对应于birdname的foldername，绝对路径
            hitfolder = folders{hitindex};
            pure_foldername = Utl.fileparts(foldername); % 只取文件夹名

            bname = Convert.bid(pure_foldername); %鸟的ID编号

            corresp_index = find(~cellfun(@isempty, regexp(bd.birdlist.BirdID,bname)));
            if length(corresp_index) > 1
                corresp_index = corresp_index(1); % 有时列表里有重复的行，这种时候只取第一个
            end
            hatchdate = regexp(convertCharsToStrings(bd.birdlist.Hatchdate(corresp_index)),'\d+','match');
            if isempty(hatchdate) %如果没有hatchdate，多半是购买的，可被认定为adult song
                answer.identifier = 1;
                answer.dirpath = hitfolder;
                return
            end

            rawformat_gender = bd.birdlist.Gender(corresp_index);
            father = bd.birdlist.father(corresp_index);
            if any(ismissing(father))||isempty(father)
                father = '0'; % 0 means unknown
            elseif length(father)>1
                father = father{1};
            end

            %path_log = Extract.filename(folders{hitfolder},'*.log'); %
            %corresp_log = readtable(path_log{1}); %一般情况下，一个鸟歌文件夹只会有一个log文件
            path_log = fullfile(folders{hitindex},strcat(foldernames{hitindex},'.log'));
            try
                corresp_log = readtable(path_log);  %这个方法很慢，所以要直接合成path_log    
            catch %如果有异常情况，认定为什么歌都没有记录

                answer.identifier = -1;
                answer.dirpath = hitfolder;
                return
            end
                   
            try
                trick = corresp_log.Var3([1,2,3,height(corresp_log)-2,height(corresp_log)-1]); % trick to reduce time
            catch ME
                trick = corresp_log.Var3; % trick to reduce time
            end
            trick = trick(find(~cellfun(@isempty, regexp(trick,'\d{4}-\d{2}-\d{2}'))));

            if isempty(trick)
                answer.identifier = 0;
                answer.dirpath = hitfolder;
                return
            end

            filedates = cellfun (@(x) datetime(Utl.fileparts(x),'InputFormat','yyyy-MM-dd_HH-mm-ss'),trick,'Uni',0);
            filedates =  sortrows( cell2table(filedates(~cellfun(@isempty,filedates))) ,'Var1','ascend');
            thehatchdate = datetime(hatchdate,'InputFormat','yyMMdd');
            juvsong = 0;
            adultsong = 0; %判断标识
            exception = 0;
            if  days(filedates.Var1(1) - thehatchdate) <90 %如果最早记录的song小于90dph
                juvsong = 1; % 有juv song
            end
            if  days(filedates.Var1(height(filedates))-thehatchdate) >= 90 ||isnat(thehatchdate) %如果最早记录的song大于90dph
                adultsong = 1; % 有 adult song
            end
            if ~isempty(regexpi(pure_foldername,... %如果song的记录有特殊标记
                    '|+|JK|iso|ps|female|PS|post|Surgery|old|new|mate|implant|WNS|Noise','match'))||...
                    length(regexp(bname,'\d{3}','match')) ~= 1
                exception = 1; % 一类异常，非正常song
            end


            if exception == 1
                answer.identifier = -1; %如果有异常情况，认定为什么歌都没有记录
            elseif juvsong == 0 && adultsong == 1 % 有成年歌，没有幼年歌
                answer.identifier = 1;
            elseif juvsong == 1 && adultsong == 1 % 有成年歌，同时也有幼年歌
                answer.identifier = 1.5;
            elseif juvsong == 1 && adultsong == 0 % 只有幼年歌
                answer.identifier = 0.5;
            else  % 除此以外，皆是零
                answer.identifier = 0; 
            end
            answer.dirpath = hitfolder;

        end


        function fathername = findFather(bd,birdname)

            % configs for reading the table
            birdlist = bd.birdlist;

            % get the corresponding index, set fathername
            index = find(~cellfun(@isempty, regexp(birdlist.BirdID,birdname)));
            fathername = birdlist.father(index);
            if isempty(fathername)
                fathername = string(missing);
            end
            if length(fathername)==2
                fathername = fathername{1};
            end


        end

 

        function date = getHatchdate(bd,birdname)
            %找到鸟的孵化日
            % the function to get hatch date
            birdlist = bd.birdlist;
            % get the corresponding index, set fathername
            index = find(~cellfun(@isempty, regexp(birdlist.BirdID,birdname)));
            date = birdlist.Hatchdate(index);
            disp(date);
        end




        function b = Deprecated_Bird(~) % dirpath must be a char
            diskletter = Utl.bucketletter;
            b.folders = Bird.birdsong;

            pathlog = strcat(diskletter,":\Yazaki-SugiyamaU\Bird-log_AK\Bird log2021 _ver_1.xlsx");
            pathlist = strcat(diskletter,":\Yazaki-SugiyamaU\Bird-log_AK\Bird_List_new ver_2.xlsx");

            isoid = Bird.getIso(pathlog, pathlist); % 获取 isolated birdids  
            temp = arrayfun(@(x) find(~cellfun(@isempty, regexpi({b.folders.bnames}.', x))), isoid,'UniformOutput',0);
            iso_index = vertcat(temp{:});
            morethan1bird_index = find(  cellfun(@length, regexp({b.folders.bnames}.', '[A-Za-z]+\d{3}'))~=1 );
            wierd_index = find(cellfun(@isempty, regexp({b.folders.bnames}.', '[OGBYR][a-zA-Z]*\d{3}')));
            juvonly_index = find(~[b.folders.adultsong].'); % 找到没有adult song的id

            temp_gender = cellfun(@convertStringsToChars,{b.folders.gender}.','Uni',0);
              temp_gender(cellfun(@isempty, temp_gender)) = {' '};
            female_index = find(~cellfun(@isempty, regexp(cellstr(temp_gender),'♀'))); % ♀

            temp_hatchdate = cellfun(@convertStringsToChars,{b.folders.hatchdate}.','Uni',0);
            temp_hatchdate(cellfun(@isempty, temp_hatchdate)) = {' '};

            temp_father = cellfun(@convertStringsToChars,{b.folders.father}.','Uni',0);
            temp_father(cellfun(@isempty, temp_father)) = {' '};

            bengalese_index = unique([find(~cellfun(@isempty, regexp(cellstr(temp_hatchdate),'BF')));...
                find(~cellfun(@isempty, regexp(cellstr(temp_gender),'BF')));...
                find(~cellfun(@isempty, regexp(cellstr(temp_father),'BF')))]); % 可能是Bengalese或者Beng-fostered

            all_badindex = unique([iso_index;morethan1bird_index;wierd_index;juvonly_index;female_index;bengalese_index]);

            b.adultfolders = b.folders(setdiff(1:length(b.folders),all_badindex));
% 
%             b.rand; % randomize

        end

        function b = rand(b)
            %把source文件夹的顺序打乱，但每次随机产生的顺序是一致的
            rng(1,'twister'); % 控制随机数生成器 repeatable random number
            b.folders = b.folders(randperm(numel(b.folders)));
        end


    end


    methods(Static) % 是Bird的重要底层方法，通常被其他方法调用而不会被直接使用

        function b = Deprecated_rmbad(input_dirs) % iso is a cell consists of isolated id
            input_dirs = Extract.folder("Z:\Yazaki-SugiyamaU\Bird-song").';

            numVars = 7;
            varNames = {'BirdID','Hatchdate','Gender','father','mother','parents','isolate'} ;
            varTypes = {'string','string','string','string','string','string','string'} ;
            dataStartLoc = 'A2';

            opts = spreadsheetImportOptions('Sheet',1,'NumVariables',numVars,...
                'VariableNames',varNames,...
                'VariableTypes',varTypes,...
                'DataRange', dataStartLoc);

            % preview('Z:\Yazaki-SugiyamaU\Bird-log_AK\Bird_List_new ver_2.xlsx',opts)
            birdlist = readtable('Z:\Yazaki-SugiyamaU\Bird-log_AK\Bird_List_new ver_2.xlsx',opts);



            dbstop if error
            [~,fnames,~] = cellfun(@fileparts,b.folders,'UniformOutput',false);

            for n = 1: length(fnames)

                alphabets = regexp(fnames{n},'[A-Za-z]');
                initial = fnames{n}(alphabets(1));

                number = fnames{n}(regexp(fnames{n},'\d'));

                new_fnames{n} = [initial,number];

            end

            % remove isolated birds
            for m = 1:length(iso)

                idxes = find (~cellfun(@isempty, regexpi(new_fnames, iso(m))));

                if ~isempty(idxes)
                    for bad = 1:length (idxes) % this is a bad code
                        b.folders{idxes(bad)} = NaN;
                    end
                end

            end

            % remove G123B234 like mixed marks
            not1Idx = find(  cellfun(@length, regexp(new_fnames, '[A-Za-z]+\d{3}'))~=1 );

            if ~isempty(not1Idx)

                for worst = 1:length(not1Idx)
                    b.folders{not1Idx(worst)} = NaN;
                end

            end


            % remove wierd folders, e.g. channel, test B345346 etc.

            wierdIdx = find (cellfun(@isempty, regexp(new_fnames, '[A-Za-z]+\d{3}')));

            if ~isempty(wierdIdx)

                for worse = 1:length(wierdIdx)
                    b.folders{wierdIdx(worse)} = NaN;
                end

            end

            b.folders(cellfun(@(x) any(isnan(x)),b.folders)) = [];

        end

        function all_isoids = getIso(pathlog, pathlist)
            % 获取被隔离鸟的id
            % 被隔离的birdid可以从两个部分得到，一部分是birdlog，另一个是birdlist
            diskletter = Utl.bucketletter;
            dbstop if error
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 第一部分
            pathlog =  strcat(diskletter,":\Yazaki-SugiyamaU\Bird-log_AK\Bird log2021 _ver_1.xlsx");
            [~,sheets] = xlsfinfo(pathlog);
            sheetnum = length(sheets); % number of sheet in one xls file

            isofromlog = {};
            isocollect = {};
            isocount = 0;
            for n = 1: sheetnum % birdlog按不同的日期被分成了一个个sheet，所以要对每一个单独的sheet做循环
                tempT = readtable(pathlog,'Sheet',n,'VariableNamingRule','preserve');
                vname = tempT.Properties.VariableNames; % 表格的variable names
                Racks = find (~cellfun(@isempty, regexp(vname, 'Rack'))); % 笼子的描述，即是否是isolated
                Annos =  find (~cellfun(@isempty, regexp(vname, 'Var4'))); % 找到包含Var4的那一列

                idxRack = find(~cellfun(@isempty, regexp(table2cell(tempT(:,Racks)),'iso'))); % 找到被描述为"隔离"的笼子，find rows contain 'iso'
                isoC =table2cell(tempT(idxRack, Annos)); % collection of annotations with 'iso'

                for d = 1: length(isoC) %对每一个隔离笼进行循环
                    isocount = isocount + 1;
                    candidatebids = regexp(isoC{d}, '[OGBYR]\d{3}','match');
                    isolated = candidatebids(3:end); %去掉父母鸟,这里的前提是表格的前两者总是父母鸟
                    isocount = isocount + 1;
                    isocollect{isocount} = candidatebids;
                end

            end
            isofromlog = horzcat(isocollect{:});
            isofromlog = unique(isofromlog).';

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 第二部分
            %pathlist = 'C:\Users\Zhehao\Dropbox (OIST)\My_TestData\BirdListMatlab.xlsx';
            listT = readtable(pathlist);
            idxiso = find(cellfun(@length,listT.isolate) > 1);
            isofromlist = listT.BirdID(idxiso);


           all_isoids = unique(vertcat(isofromlog, isofromlist));

        end

        function birdlist = readBirdlist(~)
            %统一的读取birdlist的方法，所有使用birdlist的函数都应该引用这个方法
            numVars = 7;
            varNames = {'BirdID','Hatchdate','Gender','father','mother','parents','isolate'} ;
            varTypes = {'string','string','string','string','string','string','string'} ;
            dataStartLoc = 'A2';

            opts = spreadsheetImportOptions('Sheet',1,'NumVariables',numVars,...
                'VariableNames',varNames,...
                'VariableTypes',varTypes,...
                'DataRange', dataStartLoc);

            bucketDriveletter = Utl.bucketletter;

            birdlist = readtable(strcat(bucketDriveletter,...
                ":\Yazaki-SugiyamaU\Bird-log_AK\Bird_List_new ver_2.xlsx"),opts);
        end

        function collectsyllables(targetfolders,file_per_folder) % generate .mat file for avgn analysis

            dbstop if error
            outdir = 'E:\automaticallyCollectedsyllables\bucket_avgn'
            mkdir(outdir);
            tic;
            if isempty(targetfolders)
                disp('Specify folders to use!!!');
                pause;
            else
                for n = 1: length(targetfolders)
                    %for n = 1: length(folders)
                    fprintf('Current Folder:%s, %n out of %n ',targetfolders{n}, n, length(targetfolders));
                    % for each folder
                    [~,inputid,~] = fileparts(targetfolders{n});
                    filenames = Extract.filename(targetfolders{1,n},'*.wav');
                    %filenames = extractAdult(filenames,inputid,pathlist);
                    %filenames = restrictCentroid(filenames); % restricit filenames by its centroid
                    filenames = filenames(randperm(numel(filenames))); % randomize filenames

                    if isempty(filenames)
                        disp('This folder contain only juvenile songs');
                        newline;
                        continue
                    end

                    if exist('file_per_folder','var')
                        minlen = min(file_per_folder, length(filenames));
                    else
                        minlen = length(filenames);
                    end

                    filenames = filenames(1:minlen);


                    %collect = cellfun(@(file) Raw(file).avgn, filenames,'UniformOutput',false);
                    % cellfun is slow!!!!!!!!!!!!!!!

                    parfor idx = 1: length(filenames)
                        collect{idx} = Raw(filenames{idx}).avgn;
                    end

                    syllables =  horzcat(collect{:});
                    parts = strsplit(targetfolders{n},'\');
                    save(sprintf('%s/%s.mat',outdir,Convert.bid(parts{end})),'syllables');



                    toc;
                    newline;
                end



            end


        end



    end

    methods(Static) % 有实用价值的方法
        
        function songstruct = birdsong(~)
            % To Extract the information about each song recording folders from bucket

            dbstop if error
            %Step-1 read birdlist
            bucketDriveletter = Utl.bucketletter;
            birdlist = Bird.readBirdlist;
            sourcedir = strcat(bucketDriveletter,":\Yazaki-SugiyamaU\Bird-song");

            folders = cellstr(Extract.folder(sourcedir).');
            songstruct =  struct('path',folders,'raw_foldernames',[],'bnames',[],'exception',0,...
                'gender',[],'hatchdate',0,'juvsong',0,'adultsong',0,'comment',[]  );
            getDateString = @(x) regexp(x,'\d{4}-\d{2}-\d{2}','match');
            datestring2Date = @(x) datetime(x,'InputFormat','yyyy-MM-dd');
            pw = PoolWaitbar(length(folders), 'Checking bucket folders');

           parfor k = 1: length(folders) % par

                try
                    songstruct(k).raw_foldernames = Utl.fileparts(folders{k});
                    songstruct(k).bnames = Convert.bid(songstruct(k).raw_foldernames);
                    corresp_index = find(~cellfun(@isempty, regexp([birdlist.BirdID].',songstruct(k).bnames)));
                    if length(corresp_index) > 1
                        corresp_index = corresp_index(1); % 有时列表里有重复的行，这种时候只取第一个
                    end
                    songstruct(k).hatchdate = birdlist.Hatchdate(corresp_index);
                    if length(songstruct(k).hatchdate) > 1
                        songstruct(k).hatchdate = unique(songstruct(k).hatchdate);
                        if length(songstruct(k).hatchdate) > 1 || ~ischar(songstruct(k).hatchdate)% 如果是cell，取第一个
                           songstruct(k).hatchdate = songstruct(k).hatchdate(1);
                        end
                    end
                    songstruct(k).gender = birdlist.Gender(corresp_index);
                    if length(songstruct(k).gender) > 1
                      songstruct(k).gender = vertcat(songstruct(k).gender{:}).'; 
                    end
                    songstruct(k).father = birdlist.father(corresp_index);
                    if ismissing(songstruct(k).father)||isempty(songstruct(k).father)
                        songstruct(k).father = '0'; % 0 means unknown
                    elseif length(songstruct(k).father)>1
                        songstruct(k).father = songstruct(k).father{1};
                    end
                    songstruct(k).mother = birdlist.mother(corresp_index);
                    songstruct(k).parents = birdlist.parents(corresp_index);
                    corresp_log = readtable(fullfile(...
                        folders{k},sprintf('%s.log',songstruct(k).raw_foldernames)));
                    trick = corresp_log.Var3([1,2,height(corresp_log)-2,height(corresp_log)-1]); % trick to reduce time
                    filedates = cellfun (@(x) datestring2Date(getDateString(char(x))),trick,'Uni',0);
                    filedates =  sortrows( cell2table(filedates(~cellfun(@isempty,filedates))) ,'Var1','ascend');
                    thehatchdate = datetime(songstruct(k).hatchdate,'InputFormat','yyMMdd');
                    songstruct(k).juvsong = 0;
                    songstruct(k).adultsong = 0;
                    if  days(filedates.Var1(1) - thehatchdate) <90
                        songstruct(k).juvsong = 1; % Juvenile song exist
                    end
                    if  days(filedates.Var1(height(filedates))-thehatchdate) >= 90 ||isnat(thehatchdate)
                        songstruct(k).adultsong = 1; % Juvenile song exist
                    end
                    if ~isempty(regexpi(songstruct(k).raw_foldernames,...
                            '|+|JK|iso|ps|female|PS|post|Surgery|old|new|mate|implant|WNS|Noise','match'))||...
                            length(regexp(songstruct(k).bnames,'\d{3}','match')) ~= 1
                        songstruct(k).exception = 1; % 一类异常，非正常song
                    end
                catch ME % No hatch date: purchased bird, adult
                    if strcmp(ME.identifier,'MATLAB:datetime:ParseErr')
                        songstruct(k).adultsong = 1;
                        songstruct(k).juvsong = 0;
                        songstruct(k).comment = 'Purchased';
                    elseif strcmp(ME.identifier,'MATLAB:table:UnrecognizedVarName')
                        songstruct(k).comment = 'Wrong_log_format';
                        songstruct(k).exception = 2;
                    elseif strcmp(ME.identifier,'MATLAB:badsubscript')
                        songstruct(k).comment = 'Folder_is_empty';
                        songstruct(k).exception = 3;
                    elseif strcmp(ME.identifier, 'MATLAB:textio:textio:FileNotFound')
                        songstruct(k).comment = 'Log_file_not_found';
                        songstruct(k).exception = 2;
                    else
                        songstruct(k).comment = ME.identifier;
                        songstruct(k).exception = 3;
                    end
                    continue
                end

                increment(pw)
            end

            % get the corresponding index, set fathername
            disp('ALl DONE')


        end

        function [fem,songstruct] = femalesong(~)
            songstruct = Bird.birdsong;
            gender = {info.gender}.';
            idx = find(cellfun(@isempty,gender)); % Find the indexes of empty cell
            gender(idx) = {'No'};

            fidx = find( ~cellfun(@isempty, regexp(cellstr(gender),'♀') ));


            fem = songstruct(fidx);

            % Replace the empty cells with '0000'
        end

        function findCandidates
            % A script to find out candidate birds for Piece experiments
            pathlog = "Z:\Yazaki-SugiyamaU\Bird-log_AK\Bird log2021 _ver_1.xlsx";
            pathlist = "Z:\Yazaki-SugiyamaU\Bird-log_AK\Bird_List_new ver_2.xlsx";
            birdlog = table2struct(readtable(pathlog, 'Sheet','today'));

            birdlog = birdlog(1:30); % 只取前30行
            hidx = find(strcmp({birdlog(:).B_H}.', 'H')); % 只取 holding cages （H）or juveniles
            juvidx = find(strcmp({birdlog(:).B_H}.', 'juv'));
            h_and_juv_idx = [hidx;juvidx];
            birdlog = birdlog(h_and_juv_idx);
            birds = horzcat(birdlog(:).Var4);

            % anno means annotation
            tokens = regexp(birds,'(?<name>[OBRGY]\d{3})♂\((?<anno>\d+/\d+)(?<name2>, ?(ZC|BC|Big|BB|BF))?\)','tokens');  % find male birds with birthday wrote
            tokens = tokens';
            candi = struct;
            for idx = 1: length(tokens)
                candi(idx).id = tokens{idx}{1};
                candi(idx).birthdate = tokens{idx}{2};
                candi(idx).reserve = tokens{idx}{3};
            end

            birdlist = table2struct(readtable(pathlist));
            s = 0; % shortlisted number
            shortlist = struct;
            for trump = 1: length(candi)
                thisbird = find(strcmp({birdlist(:).BirdID}.',candi(trump).id));
                birthdate = datetime( num2str(birdlist(thisbird).HatchDate),'InputFormat','yyMMdd');
                currentdate = datetime('today');
                dph = days(currentdate - birthdate);

                s = s + 1; % 此前限制了dph < 90,现在取消该限制
                shortlist(s).id = candi(trump).id;
                shortlist(s).age = dph;
                shortlist(s).annotation = candi(trump).reserve;
                var4 = {birdlog(:).Var4}.';
                idvar4 = find(~cellfun(@isempty,regexp(var4,candi(trump).id)));
                shortlist(s).cage =  birdlog(idvar4).Cage_;
                shortlist(s).father =  birdlist(thisbird).father;
            end

            %remove those do not have recording folders
            fdir = "Z:\Yazaki-SugiyamaU\Bird-song";
            folders = cellstr(Extract.folder(fdir).');
            %tokens = regexp(folders,'(?<color>[OBRGY][A-Za-z]+)(?<number>\d{3})','tokens');

            for b = 1: length(shortlist)
                letter = regexp(shortlist(b).id,'[OBRGY]','match');
                letter = letter{1};
                number = regexp(shortlist(b).id,'\d+','match');
                number = number{1};
                format = sprintf('%s[A-Za-z]+%s',letter, number);
                loc = find(~cellfun(@isempty, regexp(folders,format,'match')));
                adultexist = Bird.adultsongexist(shortlist(b).id);
                if adultexist == 1
                    shortlist(b).recordingState = 'Adult Yes';
                elseif adultexist == 0
                    shortlist(b).recordingState = 'No adult';
                end
            end

            for b = 1:length(shortlist)
                fprintf('Candidate%u----%s-----%u dph----Father: %s-------Cage%u--------%s---——-备注：%s\n',...
                    b,shortlist(b).id,shortlist(b).age,Bird.findFather(shortlist(b).id),shortlist(b).cage,shortlist(b).recordingState,shortlist(b).annotation);
            end

        end

        function siblings_names = findSiblings(birdname)
            % configs for reading the table
            birdlist = Bird.readBirdlist;

            % get the corresponding index, set fathername
            index = find(~cellfun(@isempty, regexp(birdlist.BirdID,birdname)));
            fathername = birdlist.father(index);

            if isempty(fathername)
                fathername = string(missing);
            end

            children_index = find(~cellfun(@isempty, regexp(birdlist.father,fathername)));
            siblings_index = setdiff(children_index, index);
            siblings_names = birdlist.BirdID(siblings_index);
            hatch_dates =birdlist.Hatchdate(siblings_index);

            for k = 1: length(siblings_names)
                disp(sprintf('sibling is:  %s, hatching date: %s, recording state: %s',siblings_names{k},hatch_dates{k},Bird.recording_state(siblings_names{k})));
                newline;
            end

        end


        function adultfilenames = getAdultSongs(input_birdname_or_dir)
            %获取adult song的文件名们
            % Step-1 : Extract hacth date from the BIRDLIST
            birdlist = Bird.readBirdlist;
            corresp_list_index = find(~cellfun(@isempty, regexp([birdlist.BirdID].',input_birdname_or_dir)));
            hatchdate = birdlist.Hatchdate(corresp_list_index);
            hatchdate = regexp(convertCharsToStrings(hatchdate),'\d+','match');
            thehatchdate = datetime(hatchdate,'InputFormat','yyMMdd');

            % Step-2 :get the earlist and the latest recording date from buckect dir storage
            fdir = "Z:\Yazaki-SugiyamaU\Bird-song";
            folders = cellstr(Extract.folder(fdir).');
            corresp_bnames = cellfun( @Convert.bid, Utl.fileparts(folders),'Uni',0);
            getDateString = @(x) regexp(x,'\d{4}-\d{2}-\d{2}','match');
            datestring2Date = @(x) datetime(x,'InputFormat','yyyy-MM-dd');

            try
                corresp_dir_index = find(~cellfun(@isempty, regexp(corresp_bnames,input_birdname_or_dir)));
                corresp_dir = convertCharsToStrings(folders{corresp_dir_index }.');

                if isempty(corresp_dir_index)
                    disp('Warning@Bird.adultsongexist: Not recorded!'); % If dir not found, return
                    return
                end

                if length(corresp_dir) > 1
                    path_len = cellfun(@length, {corresp_dir}.');
                    [~,min_ids] = min(path_len);
                    corresp_dir_index = corresp_dir_index(min_ids); % 取字符最短的，但此方法不一定对
                end
                thetargetdir = folders{corresp_dir_index };
            catch
                thetargetdir = input_birdname_or_dir;
            end

            trick = Extract.filename(thetargetdir,'*.wav');
            if isempty(trick)
                adultfilenames = [];
                return
            end
            filedates = cellfun (@(x) datestring2Date(getDateString(char(x))),trick,'Uni',0);
            filedates =  sortrows( cell2table(filedates(~cellfun(@isempty,filedates))) ,'Var1','ascend');


            if length(thehatchdate) > 1
                thehatchdate = thehatchdate(1);
            end

            % Step-3 :extract adult song names
            try
                adultids = find(days(filedates.Var1- thehatchdate) >=90);
                adultfiles = trick(adultids);
            catch % when hatching date not exist, means the bird is purcahsed
                adultfiles = trick;
                disp('Warning@Bird.adultsongexist: Purchased bird !')
            end
            adultfilenames = trick(find(~cellfun(@isempty, regexp(cellstr(trick),'WAV|wav'))));

        end

        function answer = adultsongexist(input_birdname)

            % Step-1 : Extract hacth date from the BIRDLIST
            birdlist = Bird.readBirdlist;
            corresp_list_index = find(~cellfun(@isempty, regexp([birdlist.BirdID].',input_birdname)));
            hatchdate = birdlist.Hatchdate(corresp_list_index);
            if strcmp(hatchdate,'-')
                answer = 1; % If dir not found, return
                disp('Purchased!!!');
                return
            end
            thehatchdate = datetime(hatchdate,'InputFormat','yyMMdd');

            % Step-2 :get the earlist and the latest recording date from buckect dir storage
            fdir = "Z:\Yazaki-SugiyamaU\Bird-song";
            folders = cellstr(Extract.folder(fdir).');
            corresp_bnames = cellfun( @Convert.bid, Utl.fileparts(folders),'Uni',0);
            corresp_dir_index = find(~cellfun(@isempty, regexp(input_birdname,corresp_bnames)));

            if isempty(corresp_dir_index)
                disp('Warning@Bird.adultsongexist: Not recorded!');
                answer = 0; % If dir not found, return
                return
            end

            getDateString = @(x) regexp(x,'\d{4}-\d{2}-\d{2}','match');
            datestring2Date = @(x) datetime(x,'InputFormat','yyyy-MM-dd');
            if length(corresp_dir_index) > 1
                path_len = cellfun(@length, {folders{corresp_dir_index }}.');
                [~,min_ids] = min(path_len);
                corresp_dir_index = corresp_dir_index(min_ids);
            end
            logpath =  Extract.filename(folders{corresp_dir_index },'*.log');
            corresp_log = readtable(logpath{1});
            trick = corresp_log.Var3([1,2,height(corresp_log)-2,height(corresp_log)-1]); % trick to reduce time
            filedates = cellfun (@(x) datestring2Date(getDateString(char(x))),trick,'Uni',0);
            filedates =  sortrows( cell2table(filedates(~cellfun(@isempty,filedates))) ,'Var1','ascend');

            if length(thehatchdate) > 1
                thehatchdate = thehatchdate(1);
            end
            % Step-3 :Answer the question
            if  days(filedates.Var1(1) - thehatchdate) <90
                disp('Warning@Bird.adultsongexist: Juvenile song exist !')
            end

            if  days(filedates.Var1(height(filedates))-thehatchdate) >= 90 ||isnat(thehatchdate)
                answer = 1; % Juvenile song exist
            else
                answer = 0;
            end

        end

        function answer = areTheySiblings(birdname1,birdname2)

            %判断两只鸟是否是siblings
            siblings_names = Bird.findSiblings(birdname1);
            if ismember(birdname2,siblings_names)
                answer = 1;
            else
                answer = 0;
            end
        end


        function allnames = allFathers(~)
            %找到birdlist里所有曾经产生过后代的雄鸟

            birdlist = Bird.readBirdlist;
            father = birdlist.father;
            allnames = unique(rmmissing(father));

        end

    end
end

