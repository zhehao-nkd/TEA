classdef Datafolder < handle
    
    properties
        folders
        selected
    end
    % from the targeted folder to extract syllables
    methods
        
        function d = Datafolder(dirpath) % dirpath must be a char
            d.folders = extract.folder(convertStringsToChars(dirpath));
            
            if isequal(exist('isoid.mat','file'),2) % 2 means it's a file.
                display('iso exists!');
                load('isoid.mat');
            else
               display('iso mat do not exists!');
                pathlog = 'pathlog.xls'
                pathlist = 'pathlist.xls'
                isoid = extriso(pathlog, pathlist);    
            end
            
            d.rmbad(isoid); % remove bad folder
            d.rand; % randomize
            
        end
        
        function d = rmbad(d,iso) % iso is a cell consists of isolated id
            dbstop if error
            [~,fnames,~] = cellfun(@fileparts,d.folders,'UniformOutput',false);
            
            
            
            for n = 1: length(fnames)
                
                alphabets = regexp(fnames{n},'[A-Za-z]');
                initial = fnames{n}(alphabets(1));
                
                number = fnames{n}(regexp(fnames{n},'\d'));
                
                new_fnames{n} = [initial,number];
                
            end
            
            % remove iso
            for m = 1:length(iso)
                
                idxes = find (~cellfun(@isempty, regexpi(new_fnames, iso(m))));
                
                if ~isempty(idxes)
                    for bad = 1:length (idxes) % this is a bad code
                        d.folders{idxes(bad)} = NaN;
                    end
                end
                
            end
            
            % remove G123B234 like mixed marks
            not1Idx = find(  cellfun(@length, regexp(new_fnames, '[A-Za-z]+\d{3}'))~=1 );
            
            if ~isempty(not1Idx)
                
                for worst = 1:length(not1Idx)
                    d.folders{not1Idx(worst)} = NaN;
                end
                
            end
            
            
            % remove wierd folders, e.g. channel, test B345346 etc.
            
            wierdIdx = find (cellfun(@isempty, regexp(new_fnames, '[A-Za-z]+\d{3}')));
            
            if ~isempty(wierdIdx)
                
                for worse = 1:length(wierdIdx)
                    d.folders{wierdIdx(worse)} = NaN;
                end
                
            end
            
            d.folders(cellfun(@(x) any(isnan(x)),d.folders)) = [];
            
        end
        
        function d = rand(d)
            d.folders = d.folders(randperm(numel(d.folders)));
        end
        
        function d = select(d,idx)
            if exist('idx','var')
                d.selected = d.folders(idx);
            else 
                d.selected = d.folders;
            end
        end
        
        function avgn(d,file_per_folder) % generate .mat file for avgn analysis
            
            dbstop if error
            outdir = 'bucket_avgn'
            mkdir(outdir);
            tic;
            if isempty(d.selected)
                disp('Specify folders to use!!!');
                pause;
            else
                for n = 1: length(d.selected)
                    %for n = 1: length(folders)
                    fprintf('Current Folder:%s, %n out of %n ',d.selected{n}, n, length(d.selected));
                    % for each folder
                    [~,inputid,~] = fileparts(d.selected{n});
                    filenames = extract.filename(d.selected{1,n},'*.wav');
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
                    parts = strsplit(d.selected{n},'\');
                    save(sprintf('%s/%s.mat',outdir,convert.bid(parts{end})),'syllables');
                    
                    
                    
                    toc;
                    newline;
                end
            
               
               
            end
            
            
            end
        
        
    end
    
    
    methods(Static)
        
        
    function isofromlog = extriso(pathlog, pathlist) % extract id of isolated bird
            
            % a function for extract id of isolated birds from birdlog
            
            %path = 'C:\Users\Zhehao\Desktop\Bird log2018 _ver_2.xlsx'
            
            [~,sheets] = xlsfinfo(pathlog);
            sheetnum = length(sheets); % number of sheet in one xls file
            
            isofromlog = {};
            for n = 1: sheetnum
                tempT = readtable(pathlog,'Sheet',n);
                vname = tempT.Properties.VariableNames; % variable names
                
                Racks = find (~cellfun(@isempty, regexp(vname, 'Rack'))); % column# of racks
                
                Annos =  find (~cellfun(@isempty, regexp(vname, 'Var'))); % column# of var
                
                
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
                    
                    
                    
                    isofromlog = [ isofromlog;thesebirds];
                    
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
    
end