% bird related function
% A class to extract song from the bucket surver, birdsong folders
classdef Bird < handle
    %All stuff about selecting birds, check their siblings, father, song
    %recording state
    
    properties
        folders
        selected
    end
    % from the targeted folder to extract syllables
    methods
        
        function b = Bird(dirpath) % dirpath must be a char
            b.folders = extract.folder(convertStringsToChars(dirpath));
            
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
                    filenames = extract.filename(b.selected{1,n},'*.wav');
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
                    save(sprintf('%s/%s.mat',outdir,convert.bid(parts{end})),'syllables');
                    
                    
                    
                    toc;
                    newline;
                end
                
                
                
            end
            
            
        end
        
        
    end
    
    
    methods(Static)
        
        
        function b = rmbad(input_dirs) % iso is a cell consists of isolated id
            input_dirs = extract.folder("Z:\Yazaki-SugiyamaU\Bird-song").';
            
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
        
        
        function isofromlog = extriso(pathlog, pathlist) % extract id of isolated bird
            dbstop if error
            % a function for extract id of isolated birds from birdlog
            
            pathlog = "Z:\Yazaki-SugiyamaU\Bird-log_AK\Bird log2021 _ver_1.xlsx"
            if ~exist(yourFolder, 'dir')
            end
            
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
        
        
    end
    
    methods(Static)
        function info = birdsong(~)
            dbstop if error
            fdir = "Z:\Yazaki-SugiyamaU\Bird-song";
            
            folders = cellstr(extract.folder(fdir).');
            
            match = regexp(folders,'(?<color>[OBRGY][A-Za-z]+)(?<number>\d{3})','match');
            
            token = regexp(folders,'(?<color>[OBRGY])[A-Za-z]+(?<number>\d{3})','tokens');
            
            not0 = find(~cellfun(@isempty,token)); % idx of tokens that are not zero
            
            token = token(not0);
            
            for k = 1: length(token)
                abbrev{k} = [token{k,1}{1}{1},token{k,1}{1}{2}];
            end
            
            
            pathlist = "Z:\Yazaki-SugiyamaU\Bird-log_AK\Bird_List_new ver_2.xlsx";
            birdlist = table2struct(readtable(pathlist));
            ID = {birdlist.BirdID}.';
            
            
            for m = 1: length(abbrev)
                thisidx = find(~cellfun(@isempty,regexp(ID,abbrev{m})));
                if ~isempty(thisidx)
                    info(m).name = abbrev{m};
                    info(m).gender = birdlist(thisidx).x__;
                end
                
            end
            
        end
        
        
        function fem = female(~)
            info = bird.birdsong;
            gender = {info.gender}.';
            idx = find(cellfun(@isempty,gender)); % Find the indexes of empty cell
            gender(idx) = {'No'};
            
            fidx = find( ~cellfun(@isempty, regexp(cellstr(gender),'♀') ));
            
            
            fem = info(fidx);
            
            % Replace the empty cells with '0000'
            
            
        end
        
        function findCandidates
            % A script to find out candidate  birds for Ephys experiments
            pathlog = "Z:\Yazaki-SugiyamaU\Bird-log_AK\Bird log2021 _ver_1.xlsx"
            
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
            %tokens = regexp(birds,'(?<name>[OBRGY]\d{3})♂\((?<anno>\d+/\d+), ZC\)','tokens')
            tokens = tokens';
            
            candi = struct;
            for idx = 1: length(tokens)
                candi(idx).id = tokens{idx}{1};
                candi(idx).birthdate = tokens{idx}{2};
                candi(idx).reserve = tokens{idx}{3};
            end
            
            pathlist = "Z:\Yazaki-SugiyamaU\Bird-log_AK\Bird_List_new ver_2.xlsx";
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
            
            folders = cellstr(extract.folder(fdir).');
            
            tokens = regexp(folders,'(?<color>[OBRGY][A-Za-z]+)(?<number>\d{3})','tokens');
            
            
            
            for b = 1: length(shortlist)
                
                %shortlist(b).id
                
                letter = regexp(shortlist(b).id,'[OBRGY]','match');
                letter = letter{1};
                number = regexp(shortlist(b).id,'\d+','match');
                number = number{1};
                
                format = sprintf('%s[A-Za-z]+%s',letter, number);
                loc = find(~cellfun(@isempty, regexp(folders,format,'match')));
                %folders{idx};
                if length(loc) == 1
                    shortlist(b).record = 'Recorded';
                elseif length(loc) == 0
                    shortlist(b).record = 'Not recorded yet';
                else
                    shortlist(b).record = '[ERROR]';
                end
                
                
                disp(sprintf('Candidate%u---------%s---------%udph--------Cage%u--------%s----——-备注：',...
                    b,shortlist(b).id,shortlist(b).age,shortlist(b).cage,shortlist(b).record,shortlist(b).annotation ));
                newline;
            end
            
            
            %%% show which bird are not recorded yet
            
            newline;
            notidx = find( strcmp({shortlist(:).record}.','Not recorded yet'));
            notyet = shortlist(notidx);
            
            for no = 1: length(notyet)
                
                
                disp(sprintf('Candidate%u------------%s---------%udph---------Cage%u----------%s----——-备注：',...
                    no,notyet(no).id,notyet(no).age,notyet(no).cage,notyet(no).record,notyet(no).annotation ));
                newline;
            end
            
            
        end
        
        function siblings_names = findSiblings(birdname)
            % configs for reading the table
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
            numVars = 5;
            varNames = {'BirdID','Hatch date','Gender','father','mother','isolate'} ;
            varTypes = {'string','string','string','string','string','string'} ;
            dataStartLoc = 'A2';
            
            opts = spreadsheetImportOptions('Sheet',1,'NumVariables',numVars,...
                'VariableNames',varNames,...
                'VariableTypes',varTypes,...
                'DataRange', dataStartLoc);
            
            % preview('Z:\Yazaki-SugiyamaU\Bird-log_AK\Bird_List_new ver_2.xlsx',opts)
            birdlist = readtable('Z:\Yazaki-SugiyamaU\Bird-log_AK\Bird_List_new ver_2.xlsx',opts);
            
            % get the corresponding index, set fathername
            index = find(~cellfun(@isempty, regexp(birdlist.BirdID,birdname)));
            fathername = birdlist.father(index);
            if isempty(fathername)
                fathername = string(missing);
            end
            
        end
        
        function answer = recording_state(birdname,mode)
            
            fdir = "Z:\Yazaki-SugiyamaU\Bird-song";
            
            folders = cellstr(extract.folder(fdir).');
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
            %             folders (sharedids)
            
            
            %             tokens = regexp(folders,'(?<color>[OBRGY][A-Za-z]+)(?<number>\d{3})','tokens');
        end
        
        
        function answer = adultsongexist(birdname)
            fdir = "Z:\Yazaki-SugiyamaU\Bird-song";
            
            folders = cellstr(extract.folder(fdir).');
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
            
            % get the corresponding index, set fathername
            index = find(~cellfun(@isempty, regexp(birdlist.BirdID,birdname)));
            date = birdlist.Hatchdate(index);
            disp(date);
        end
        
    end
end

