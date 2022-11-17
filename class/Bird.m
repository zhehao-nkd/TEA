
classdef Bird < handle
    % A class for   (1) extracting bird songs from the bucket surver, birdsong folders
    %               (2) finding candidate birds for surgery.

    properties
        folders
        selected
    end
    % from the targeted folder to Extract syllables
    methods
        
        function b = Bird(~) % dirpath must be a char
            b.folders = Extract.folder(convertStringsToChars("Z:\Yazaki-SugiyamaU\Bird-song"));
            
            if isequal(exist('isoid.mat','file'),2) % 2 means it's a file.
                display('iso exists!');
                load('isoid.mat');
            else
                display('iso mat do not exists!');
                pathlog = "Z:\Yazaki-SugiyamaU\Bird-log_AK\Bird log2021 _ver_1.xlsx"
                pathlist = "Z:\Yazaki-SugiyamaU\Bird-log_AK\Bird_List_new ver_2.xlsx"
                isoid = Bird.extriso(pathlog, pathlist);
            end
            
            b.rmbad(isoid); % remove bad folder
            b.rand; % randomize
            
        end
        
        
        function b = rand(b)
            rng(1,'twister'); % repeatable random number
            b.folders = b.folders(randperm(numel(b.folders)));
        end
        
        function b = select(b,idx)
            if exist('idx','var')
                b.selected = b.folders(idx);
            else
                b.selected = b.folders;
            end
        end
        
        function collectsyllables(b,file_per_folder) % generate .mat file for avgn analysis
            
            dbstop if error
            outdir = 'E:\automaticallyCollectedsyllables\bucket_avgn'
            mkdir(outdir);
            tic;
            if isempty(b.selected)
                disp('Specify folders to use!!!');
                pause;
            else
                for n = 1: length(b.selected)
                    %for n = 1: length(folders)
                    fprintf('Current Folder:%s, %n out of %n ',b.selected{n}, n, length(b.selected));
                    % for each folder
                    [~,inputid,~] = fileparts(b.selected{n});
                    filenames = Extract.filename(b.selected{1,n},'*.wav');
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
                    parts = strsplit(b.selected{n},'\');
                    save(sprintf('%s/%s.mat',outdir,Convert.bid(parts{end})),'syllables');
                    
                    
                    
                    toc;
                    newline;
                end
                
                
                
            end
            
            
        end
        
        
    end
    
    
    methods(Static) % 是Bird的重要底层方法，通常被其他方法调用而不会被直接使用
        
        function b = rmbad(input_dirs) % iso is a cell consists of isolated id
            input_dirs = Extract.folder("Z:\Yazaki-SugiyamaU\Bird-song").';
            
            numVars = 5;
            varNames = {'BirdID','Hatchdate','Gender','father','mother','isolate'} ;
            varTypes = {'string','string','string','string','string','string'} ;
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
        
        function isofromlog = extriso(pathlog, pathlist) % Extract id of isolated bird
            dbstop if error
            % a function for Extract id of isolated birds from birdlog
            
            pathlog = "Z:\Yazaki-SugiyamaU\Bird-log_AK\Bird log2021 _ver_1.xlsx"
      
            
            %path = 'C:\Users\Zhehao\Desktop\Bird log2018 _ver_2.xlsx'
            
            [~,sheets] = xlsfinfo(pathlog);
            sheetnum = length(sheets); % number of sheet in one xls file
            
            isofromlog = {};
            for n = 1: sheetnum
                tempT = readtable(pathlog,'Sheet',n);
                vname = tempT.Properties.VariableNames; % variable names
                
                Racks = find (~cellfun(@isempty, regexp(vname, 'Rack'))); % column# of racks
                
                Annos =  find (~cellfun(@isempty, regexp(vname, 'Var4'))); % column# of var
                
                
                idxRack = find (~cellfun(@isempty, regexp(table2cell(tempT(:,Racks))...
                    ,'iso'))); % find rows contain 'iso'
                
                
                
                isoC =table2cell(tempT(idxRack, Annos)); % collection of annotations with 'iso'
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%% for each sheet
                
                for d = 1: length(isoC)
                    
                    % For each element in isoC
                    
                    birds = regexp(isoC{d}, '[A-Z]\d\d\d') ;
                    
                    isolated = birds(3:end); % remove parents
                    
                    
                    temp = isoC{d}(sort([isolated,isolated+1,isolated+2,isolated+3]))';
                    
                    thesebirds =   cellstr(reshape (temp ,4 ,[])');
                    
                    
                    
                    isofromlog = [isofromlog;thesebirds];
                    
                end
                
            end
            
            isofromlog(cellfun(@isempty, isofromlog)) = [];
            
            isofromlog = unique(isofromlog);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % ISO FROM LIST
            %pathlist = 'C:\Users\Zhehao\Dropbox (OIST)\My_TestData\BirdListMatlab.xlsx';
            
            listT = readtable(pathlist)
            
            
            idxiso = find (cellfun(@length,listT.isolate) > 1)
            
            isofromlist = listT.BirdID(idxiso)
            
            
            ISOCELL = unique(vertcat(isofromlog, isofromlist));
            
            %{
                    No hard coding!!

                    Test = readtable('E:\0904_R612\R612_Z1_SPKC03.xlsx')
                    Test.Properties.VariableNames
                    VN = Test.Properties.VariableNames;
                    regexp(VN, '^Var\d+$')
                    ~cellfun(@isempty, regexp(VN, '^Var\d+$'))
                    idxWaveform =  ~cellfun(@isempty, regexp(VN, '^Var\d+$'));
                    unique(Test.Unit)
            %}
            
        end
        
        function birdlist = readBirdlist(~)
            numVars = 5;
            varNames = {'BirdID','Hatchdate','Gender','father','mother','isolate'} ;
            varTypes = {'string','string','string','string','string','string'} ;
            dataStartLoc = 'A2';
            
            opts = spreadsheetImportOptions('Sheet',1,'NumVariables',numVars,...
                'VariableNames',varNames,...
                'VariableTypes',varTypes,...
                'DataRange', dataStartLoc);
            
            bucketDriveletter = Utl.bucketletter;
            
            birdlist = readtable(strcat(bucketDriveletter,...
                ":\Yazaki-SugiyamaU\Bird-log_AK\Bird_List_new ver_2.xlsx"),opts);
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
            
            parfor k = 1: length(folders)
                
                try
                    songstruct(k).raw_foldernames = Utl.fileparts(folders{k});
                    songstruct(k).bnames = Convert.bid(songstruct(k).raw_foldernames);
                    corresp_index = find(~cellfun(@isempty, regexp([birdlist.BirdID].',songstruct(k).bnames)));
                    songstruct(k).hatchdate = birdlist.Hatchdate(corresp_index);
                    songstruct(k).gender = birdlist.Gender(corresp_index);
                    songstruct(k).father = birdlist.father(corresp_index);
                    if ismissing(songstruct(k).father)||isempty(songstruct(k).father)
                      songstruct(k).father = '0'; % 0 means unknown 
                    end
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
        
        function fathername = findFather(birdname)
            
            % configs for reading the table
            birdlist = Bird.readBirdlist;
            
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
        
        function answer = recording_state(birdname,mode)
            
            fdir = "Z:\Yazaki-SugiyamaU\Bird-song";
            
            folders = cellstr(Extract.folder(fdir).');
            folernames = {};
            for k = 1: length(folders)
                temp = split(folders{k},'\');
                foldernames{k} = temp{length(temp)};
                
            end
            
            color = regexp(birdname,'[A-Z]','match');
            number = regexp(birdname,'\d+','match');
            
            colorids = find(~cellfun(@isempty, regexp(foldernames,color{1})));
            numberids = find(~cellfun(@isempty, regexp(foldernames,number{1})));
            
            sharedids = intersect(colorids,numberids);
            
            if isempty(sharedids)
                answer = 'No';
            elseif length(sharedids) == 1
                answer = 'Yes';
            elseif length(sharedids) >1
                answer = 'More than one hits';
            end
            
            if exist('mode','var')&& mode == 1 % Mode1: care about whether the recorded song is adult or juvenile song 
                
            end
        end
        
        function adultfilenames = getAdultSongs(input_birdname_or_dir)
            % Step-1 : Extract hacth date from the BIRDLIST
            birdlist = Bird.readBirdlist;
            corresp_list_index = find(~cellfun(@isempty, regexp([birdlist.BirdID].',input_birdname_or_dir)));
            hatchdate = birdlist.Hatchdate(corresp_list_index);

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
        
        function date = getHatchdate(birdname)
            % the function to get hatch date
            birdlist = Bird.readBirdlist;
            % get the corresponding index, set fathername
            index = find(~cellfun(@isempty, regexp(birdlist.BirdID,birdname)));
            date = birdlist.Hatchdate(index);
            disp(date);
        end
        

        function sonnames = findSons(fathername)
               % configs for reading the table
            birdlist = Bird.readBirdlist;
            
%             % get the corresponding index, set fathername
%             index = find(~cellfun(@isempty, regexp(birdlist.BirdID,birdname)));
%             fathername = birdlist.father(index);
            
%             if isempty(fathername)
%                 fathername = string(missing);
%             end
            
            children_index = find(~cellfun(@isempty, regexp(birdlist.father,fathername)));
            male_index = find(~cellfun(@isempty, regexp(birdlist.Gender,"♂")));
            sons_index = intersect(children_index, male_index);
            sons_names = birdlist.BirdID(sons_index);
            hatch_dates =birdlist.Hatchdate(sons_index);
            
            for k = 1: length(sons_names)
                disp(sprintf('His song is:  %s, hatching date: %s, recording state: %s',sons_names{k},hatch_dates{k},Bird.recording_state(sons_names{k})));
                newline;
            end
            


        end
   
        function allnames = allFathers(~)

            birdlist = Bird.readBirdlist;
            father = birdlist.father;
            allnames = unique(rmmissing(father));

        end
    
    end
end

