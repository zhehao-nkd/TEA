classdef Archon
    % batch processing the raw data( .plx, .txt, .wav files sorted within
    % folders
    % 批量对raw data进行操作 ， 用以生成 Analysis object
    
    properties
        bird_folders
        zpfolders
        all_data_info
    end
    
    methods
        
        function a = Archon(dirpath)
            
            dbstop if error
            a.bird_folders = [extract.folder(dirpath)].'; %#ok<NBRAK>
            
            pure_folname = {}; % 初始化
            for k = 1: length(a.bird_folders)
                temp1 = split(a.bird_folders{k},'\');
                pure_folname{k} = cellstr(temp1{end});
            end
            
            contain_ephys_ids = find(~cellfun(@isempty, regexp([pure_folname{:}].','Ephys')));
            if ~isempty(contain_ephys_ids)
                a.bird_folders =  a.bird_folders(contain_ephys_ids);
            else
                a.bird_folders =  dirpath;
            end
            
            infocollect = {};
            for k = 1:length(a.bird_folders)
                plxfiles = extract.filename(a.bird_folders{k},'*.plx');
                txtfiles = extract.filename(a.bird_folders{k},'*.txt');
                
                if exist(sprintf('%s\\Stimuli',a.bird_folders{k}), 'dir')
                    stidirs = extract.folder(sprintf('%s\\Stimuli',a.bird_folders{k})).';
                else % 如果未把所有stimuli放进stimuli folder里面
                    stidirs = extract.folder(sprintf('%s',a.bird_folders{k})).'
                end
                
                plxname = {}; % extract .plx filenames
                for m = 1: length(plxfiles)
                    [~,plxname{m},~] = fileparts(plxfiles{m});
                end
                
                txtname = {};% extract .txt filenames
                for m = 1: length(txtfiles)
                    [~,txtname{m},~] = fileparts(txtfiles{m});
                end
                
                stiname = {}; % extract stimuli folder names
                for m = 1: length( stidirs)
                    temps = split( stidirs{m},'\');
                    stiname{m} = temps{end};
                end
                
                % Based on .plx filenames, find the corresponding .txt and
                % stimuli folder names
                datainfo = struct;
                for m = 1: length(plxname)
                    
                    txt_hit = find(~cellfun(@isempty, regexp([cellstr(txtname)].',plxname{m})));
                    
                    sti_hit = find(~cellfun(@isempty, regexp([cellstr(stiname)].',plxname{m})));
                    
                    % fill the struct
                    sadtemp = regexp(plxname{m}, '[_-]', 'split');
                    datainfo(m).animal = sadtemp{end};
                    datainfo(m).plxname = plxname{m};
                    datainfo(m).keyword = regexp(plxname{m},'Con|Deg|Detail|Frag|Repla|Mir|Rev','match');
                    datainfo(m).sectionname = regexp(plxname{m},'Z\d{2}|P\d{2}','match');
                    fuck = regexp(cellstr(datainfo(m).sectionname),'\d{2}','match');
                    datainfo(m).sectionid = str2num(fuck{1}{1});
                    datainfo(m).subfilename = regexp(plxname{m},'F\d{1}','match');
                    if ~isempty(datainfo(m).subfilename) %如果文件名里存在F标识subfile的id
                        fuck2 = regexp(datainfo(m).subfilename,'\d','match');
                        datainfo(m).subfileid = str2num(fuck{1}{1});
                    elseif ~isempty(regexp(datainfo(m).plxname,'Mrg|mrg'))% 如果不存在F标识却存在‘Mrg|mrg’标识
                        datainfo(m).subfilename = "Mrg";
                    else % 如果不存在任何标识
                        datainfo(m).subfilename = [];
                    end
                    
                    datainfo(m).plxpath = plxfiles{m};
                    if ~isempty(txt_hit)
                        datainfo(m).txtpath = txtfiles{txt_hit};
                    else
                        datainfo(m).txtpath = [];
                    end
                    
                    if ~isempty(sti_hit)
                        datainfo(m).stidirpath = stidirs{sti_hit};
                    else
                        datainfo(m).stidirpath  = [];
                    end
                    
                end
                
                infocollect{k} = datainfo;
                
            end
            
            a.all_data_info = horzcat(infocollect{:});
            
        end
        
    end
    
    
    methods(Static) % 底层方法, methods which are themselves not frequently used but are the basis of more superior methods
        function reorderFilename(targetdir, ext)
            
            files = extract.filesAllLevel(targetdir, ext);
            
            wb = waitbar(0,'Start renaming files...');
            
            for id = 1:length(files)
                
                waitbar(id/length(files),wb,sprintf('Processing %u of %u Files',id,length(files)) )
                % Get the file name
                oldname = files{id};
                
                [folder,oldfilename,ext] = fileparts(oldname);
                newfilename = strrep(oldfilename,'-','_');
                parts = split(newfilename,'_');
                
                index = find(~cellfun(@isempty, regexp(parts,'^[OBRYGX]\d{3}$'))); % ^EXP$ is perfect match
                
                % 替换index位和第一位的顺序
                if ~isempty(index)
                    temp = parts{1};
                    parts{1} = parts{index};
                    parts{index} = temp;
                end
                
                index2 = find(~cellfun(@isempty, regexp(parts,'^[PZ]\d{2}F\d{1}$'))); % ^EXP$ is perfect match
                
                
                % 替换index2位和第二位的顺序
                if ~isempty(index2)
                    temp = parts{2};
                    parts{2} = parts{index2};
                    parts{index2} = temp;
                end
                
                
                newfilename = strjoin(parts,'_');
                
                newname = fullfile(folder,strcat(newfilename,ext));
                
                if ~strcmp(oldname,newname)
                    movefile(oldname, newname,'f');
                end
            end
            
            close(wb);
        end
        % reorder filename, let bird id be the first
        function reorderDirname(targetdir)
            % reorder dirname, let bird id be the first
            dirs = extract.foldersAllLevel(targetdir);
            strlen = [];
            for m = 1: length(dirs)
                strlen(m) = length(dirs{m});
            end
            
            [newstrlen,order] = sort(strlen,'descend'); %按长度从长到短,目的是先改子文件夹名，后改父文件名
            dirs = dirs(order);
            
            wb = waitbar(0,'Start renaming folders...');
            
            for id = 1:length(dirs)-1  % 最大的folder不会改
                
                waitbar(id/length(dirs),wb,sprintf('Processing %u of %u Folders',id,length(dirs)) )
                
                % Get the file name
                oldname = dirs{id};
                FS = filesep;
                dirranks = split(oldname,FS);
                
                olddirname = dirranks{length(dirranks)};
                newdirname = strrep(olddirname,'_','-');
                parts = split(newdirname,'-');
                
                index = find(~cellfun(@isempty, regexp(parts,'^[OBRYGX]\d{3}$'))); % ^EXP$ is perfect match
                
                % 替换index位和第一位的顺序
                if ~isempty(index)
                    temp = parts{1};
                    parts{1} = parts{index};
                    parts{index} = temp;
                end
                
                index2 = find(~cellfun(@isempty, regexp(parts,'^[PZ]\d{2}F\d{1}$'))); % ^EXP$ is perfect match
                
                % 替换index2位和第二位的顺序
                if ~isempty(index2)
                    temp = parts{2};
                    parts{2} = parts{index2};
                    parts{index2} = temp;
                end
                
                
                newdirname = strjoin(parts,'-');
                
                dirranks{length(dirranks)} = newdirname;
                newname = fullfile(dirranks{:});
                
                if ~strcmp(oldname,newname)
                    movefile(oldname, newname);
                end
            end
            close(wb);
            
        end
        
        function reorderNames(targetdir)
            convert.bid_first_AllLevel(targetdir, '*.plx');
            convert.bid_first_AllLevel(targetdir, '*.txt');
            convert.bid_first_AllLevel_dirVer(targetdir);
        end
        
        function padFilename(targetdir, ext) % format data into G573_P01F2_Degs kind of format
            files = extract.filesAllLevel(targetdir, ext);
            
            wb = waitbar(0,'Start renaming files...');
            
            for id = 1:length(files)
                
                waitbar(id/length(files),wb,sprintf('Processing %u of %u Files',id,length(files)) )
                % Get the file name
                oldname = files{id};
                
                [folder,oldfilename,ext] = fileparts(oldname);
                newfilename = strrep(oldfilename,'-','_');
                parts = split(newfilename,'_');
                
                
                index = find(~cellfun(@isempty, regexp(parts,'^[OBRYGX]\d{3}$'))); % ^EXP$ is perfect match
                
                % 如果不存在BID，加上BID
                if isempty(index)
                    FS = filesep;
                    dirparts = split(targetdir,FS);
                    ephysid = find(~cellfun(@isempty, regexp(dirparts,'Ephys')));
                    if ~isempty(ephysid)
                        bid = regexp(dirparts{ephysid},'[OBRYGX]\d{3}','match');
                        bid = bid{1};
                    else
                        bid = convertStringsToChars(regexp(targetdir,'[OBRYGX]\d{3}','match'));
                    end
                    disp(parts);
                    badcodetemp = [{bid};parts];
                    parts = badcodetemp ; % parts means the parts of the filename, which is modified by this code
                    
                end
                
                
                index2 = find(~cellfun(@isempty, regexp(parts,'^[PZ]\d{2}F\d{1}$'))); % ^EXP$ is perfect match
                
                
                %如果file id格式不对，改正格式
                if isempty(index2)
                    PFid = find(~cellfun(@isempty, regexp(parts,'[PZ]\d+F\d+')));
                    
                    if ~isempty(PFid)
                        to_be_edited = parts{PFid};
                        porz = regexp(to_be_edited,'[PZ]','match'); porz = porz{1};
                        partnum = regexp(to_be_edited,'(?<=[PZ])\d+(?=F)','match'); partnum = str2num(partnum{1});
                        filenum = regexp(to_be_edited,'(?<=[PZ]\d+[F])\d+','match'); filenum = str2num(filenum{1});
                        formatSpec = '%s%02dF%1d';
                        edited = sprintf(formatSpec,porz,partnum,filenum);
                        parts{PFid} = edited;
                    end
                end
                
                
                newfilename = strjoin(parts,'_');
                
                newname = fullfile(folder,strcat(newfilename,ext));
                
                if ~strcmp(oldname,newname)
                    movefile(oldname, newname,'f');
                end
            end
            
            close(wb);
        end
    end
    
    methods(Static)
        
        % 把源数据的命名格式统一
        
        function standardizeRawname(targetdir) % 标准化文件夹里所有raw data的名称
            Archon.padFilename(targetdir, '*.txt')
            Archon.padFilename(targetdir, '*.plx')
            Archon.reorderFilename(targetdir, '*.plx');
            Archon.reorderFilename(targetdir, '*.txt');
            Archon.reorderDirname(targetdir);
            
        end
        
        % 生成analysis object
        
        function genAnalysisFromSingleZPDir(datadir)
            % generate Analysis object from a single ZP-marked raw data folder
            [neuroster,~,~,~,~] = Archon.extractAnalysisInfo(datadir);
            Archon.batch_genAnalysisSameRecordingFile(neuroster);
            
        end
        
        function  batch_genAnalysisSameRecordingFile(input_roster)
            dbstop if error
            % 从同一批原始plx数据中生成所有的analysis文件
            %采用此方法的好处是不用一次又一次地调用plexon读取方法，可以节省非常多的时间
            tic
            %[neuroster,malfunctions,cdrfroster,noSinChan_neuroster,nofeature_neuroster] = Archon.extractAnalysisInfo(folder);
            
            D = parallel.pool.DataQueue;
            h = waitbar(0, '开始生成 Analysis objects');
            num_files = length(input_roster);% Dummy call to nUpdateWaitbar to initialise
            utl.UpdateParforWaitbar(num_files, h);% Go back to simply calling nUpdateWaitbar with the data
            afterEach(D, @utl.UpdateParforWaitbar);
            
            neurondir =  unique({input_roster.neurondir}.');
            
            for mm = 1: length(neurondir) % par
                
                subroster = find(strcmp(cellstr({input_roster.neurondir}.'),neurondir{mm}));
                hit_subdir = neurondir{mm};
                plxfiles = extract.filename(hit_subdir,'*.plx');
                plxdata = {};
                for k = 1: length(plxfiles)
                    plxdata{k} = Trigger(plxfiles{k});
                end
                
                for k = 1: length(subroster) % par
                    % do stuff
                    birdname = input_roster(k).birdname;
                    ZPid = input_roster(k).neuronid;
                    channelname = input_roster(k).channelname;
                    unit = input_roster(k).unit;
                    
                    
                    %Archon.genAnalysis(birdname,ZPid,channelname,unit);
                    
                    %                    hitid = find(~cellfun(@isempty,regexp(cellstr(arch.bird_folders),birdname)));
                    %                    hitbird_folder = arch.bird_folders{hitid};
                    %
                    
                    %处理hit folder 内部的 subfolders
                    %                    subdirs = neurondir{mm};
                    %                    hitZorP = find(~cellfun(@isempty,regexp(cellstr(subdirs.'),ZPid)));
                    hit_subdir = neurondir{mm};
                    
                    %                    plxfiles = extract.filename(hit_subdir,'*.plx');
                    %                    zpids = find(~cellfun(@isempty, regexp(cellstr(plxfiles),'[ZP]\d+F\d+')));
                    %                    plxfiles = plxfiles(zpids);
                    featuredir = append(hit_subdir,'\Feature');
                    txtfiles = extract.filename(hit_subdir,'*.txt');
                    zpids = find(~cellfun(@isempty, regexp(cellstr(txtfiles),'[ZP]\d+F\d+'))); % same record
                    SPKCid = find(~cellfun(@isempty, regexp(cellstr(txtfiles),channelname))); % same channel
                    Uid = find(~cellfun(@isempty, regexp(cellstr(txtfiles),sprintf('U%d',unit)))); % same unit
                    if isempty(Uid)
                        txtfiles = txtfiles(zpids);
                    else
                        txtfiles = txtfiles(intersect(intersect(zpids,SPKCid),Uid));
                    end
                    
                    
                    stimuli_collect_dir = append(hit_subdir,'\Stimuli');
                    stimulidirs = extract.folder(stimuli_collect_dir ).';
                    
                    Archon.genAnalysisDifferentInput(plxdata,txtfiles,stimulidirs,featuredir,channelname,unit);
                    % Note we send only an "increment" for the waitbar.
                    send(D, 1);
                    
                end
                
            end
            toc
            
        end
        
        function [neuroster,malfunctions,cdrfroster,noSinChan_neuroster,nofeature_neuroster] = extractAnalysisInfo(datadir)
            
            dbstop if error
            % 分析所包含的所有Z P folder,生成Analysis对象，并且对于那些不可分析的（某些信息不够），生成记录文件以便补充这些信息
            %datadir = 'D:';
            
            allfolders = extract.foldersAllLevel(datadir).';
            index = find(~cellfun(@isempty, regexp(allfolders,'[ZP]\d+$')));
            neurondir = {allfolders{index}}.';
            neuroster = struct; % Neuron 花名册
            malfunctions = struct; % malfunctions
            count = 0;
            malcount = 0;
            
            % waitbars
            wb2 = waitbar(0,'Neuron of this bird');
            set(wb2,'doublebuffer','on');
            
            
            for kk = 1: length(neurondir)
                
                waitbar(kk/length(neurondir),wb2,sprintf('正处理 %u个神经元文件夹中的%u',length(neurondir) ,kk))
                
                
                %判断 txt 文件是否完备
                plxfiles = extract.filename(neurondir{kk},'*.plx');
                zpids = find(~cellfun(@isempty, regexp(cellstr(plxfiles),'[ZP]\d+F\d+'))); % 剔除 merge 的 plx文件
                plxfiles = plxfiles(zpids);
                
                txtfile = extract.filename(neurondir{kk},'*.txt');
                zpids = find(~cellfun(@isempty, regexp(cellstr(txtfile),'[ZP]\d+F\d+'))); %剔除额外的 plx文件
                txtfile = txtfile(zpids);
                
                stimuliparentdir = fullfile(neurondir{kk},'Stimuli');
                stimulidir = extract.folder(stimuliparentdir);
                if isempty(stimulidir)
                    malcount = malcount + 1;
                    malfunctions(malcount).birdname = regexp(neurondir{kk},'[BRGYOX]\d{3}','match');
                    malfunctions(malcount).neuronid = regexp(neurondir{kk},'[ZP]\d+','match');
                    malfunctions(malcount).reason = 'stimuli_dir的数量不对';
                    continue
                end
                zpids = find(~cellfun(@isempty, regexp(cellstr(stimulidir),'[ZP]\d+F\d+'))); %剔除额外的 plx文件
                stimulidir = stimulidir(zpids);
                
                if length(plxfiles)~= length(stimulidir)
                    malcount = malcount + 1;
                    malfunctions(malcount).birdname = regexp(neurondir{kk},'[BRGYOX]\d{3}','match');
                    malfunctions(malcount).neuronid = regexp(neurondir{kk},'[ZP]\d+','match');
                    malfunctions(malcount).reason = 'stimuli_dir的数量不对';
                    continue
                end
                
                if length(plxfiles)>length(txtfile)
                    malcount = malcount + 1;
                    malfunctions(malcount).birdname = regexp(neurondir{kk},'[BRGYOX]\d{3}','match');
                    malfunctions(malcount).neuronid = regexp(neurondir{kk},'[ZP]\d+','match');
                    malfunctions(malcount).reason = 'txt文件不足';
                    continue
                end
                
                %consistentunits = Archon.findConsistentNeurons(neurondir{kk});
                consistentunits = Archon.findConExistNeurons(neurondir{kk});
                % 判断是否存在每个文件都有的神经元
                if isempty(consistentunits)
                    malcount = malcount + 1;
                    malfunctions(malcount).birdname = regexp(neurondir{kk},'[BRGYOX]\d{3}','match');
                    malfunctions(malcount).neuronid = regexp(neurondir{kk},'[ZP]\d+','match');
                    malfunctions(malcount).reason = '不存在每个文件都有的神经元';
                    continue
                end
                
                for kkk = 1: length(consistentunits)
                    count = count + 1;
                    
                    neuroster(count).birdname = regexp(neurondir{kk},'[BRGYOX]\d{3}','match');
                    neuroster(count).neuronid = regexp(neurondir{kk},'[ZP]\d+','match');
                    neuroster(count).birdneuron = strcat(neuroster(count).birdname,neuroster(count).neuronid);
                    neuroster(count).howManyTxtFile = length(txtfile);
                    neuroster(count).neurondir = neurondir{kk};
                    
                    temp = regexp(consistentunits{kkk},'SPKC\d+','match');neuroster(count).channelname = temp{1};
                    temp = regexp(consistentunits{kkk},'(?<=_)\d+','match');neuroster(count).unit = str2num(temp{1});
                    neuroster(count).fullid = [neuroster(count).birdname,neuroster(count).neuronid,neuroster(count).channelname,neuroster(count).unit];
                    stimulidir = fullfile(neurondir{kk},'Stimuli');
                    
                    stimuli_filenames = extract.filesAllLevel(stimulidir,'*.wav');
                    
                    neuroster(count).frag = 0;  % 判断是否有frag
                    if ~isempty(find(~cellfun(@isempty,regexp(stimuli_filenames,'Frag')), 1))
                        neuroster(count).frag = 1;
                    end
                    
                    neuroster(count).deg = 0;  % 判断是否有deg
                    if ~isempty(find(~cellfun(@isempty,regexp(stimuli_filenames,'deg')), 1))
                        neuroster(count).deg = 1;
                    end
                    
                    neuroster(count).repla = 0;  % 判断是否有repla
                    if ~isempty(find(~cellfun(@isempty,regexp(stimuli_filenames,'Repla')), 1))
                        neuroster(count).repla = 1;
                    end
                    
                    neuroster(count).wns = 0;  % 判断是否有WNS
                    if ~isempty(find(~cellfun(@isempty,regexp(stimuli_filenames,'wns|WNS')), 1))
                        neuroster(count).wns = 1;
                    end
                    
                    
                    neuroster(count).allexist = 0;  % 判断是否同时存在 frag deg repla
                    if ~isempty(find(~cellfun(@isempty,regexp(stimuli_filenames,'Frag')), 1)) &&...
                            ~isempty(find(~cellfun(@isempty,regexp(stimuli_filenames,'deg')), 1)) &&...
                            ~isempty(find(~cellfun(@isempty,regexp(stimuli_filenames,'Repla')), 1))
                        neuroster(count).allexist = 1;
                    end
                    
                    neuroster(count).featureexist = 0;  % 是否存在 feature files
                    feature_dir = fullfile(neurondir{kk},'Feature');
                    if exist(feature_dir,'dir')
                        feature_files = extract.filename(feature_dir,'*.txt');
                        
                        if ~isempty(find(~cellfun(@isempty, regexp(cellstr(feature_files),'Info')))) &&...
                                ~isempty(find(~cellfun(@isempty, regexp(cellstr(feature_files),'Data'))))
                            neuroster(count).featureexist = 1;
                        end
                    end
                    
                    neuroster(count).SCDexist = 0; % 是否存在SingleChannel_Merged_Stimuli文件夹
                    subdirs = extract.folder(neurondir{kk});
                    if ~isempty(find(~cellfun(@isempty, regexp(cellstr(subdirs).','SingleChannel_Merged_Stimuli'))))
                        neuroster(count).SCDexist = 1;
                    end
                    
                end
            end
            
            if isfield(neuroster,'allexist')
                cdrfroster = neuroster([neuroster.allexist].' == 1);
                
            end
            %not_cdrfroster = neuroster([neuroster.allexist].' == 0);
            noSinChan_neuroster = neuroster([neuroster.SCDexist].' == 0);
            
            nofeature_neuroster = neuroster([neuroster.featureexist].' == 0);
            
        end
        
        function A = genAnalysisDifferentInput(plxdata,txtfiles,stimulidirs,featuredir,channelname,unit)
            
            
            Ns = {};
            
            for k = 1: length(plxdata)
                fid = regexp(plxdata{k}.inputpath,'(?<=F)\d*','match'); % id for the file rank
                
                % find corresponding txt file
                followF = regexp(cellstr(txtfiles),'(?<=F)\d*','match');
                duiying_txt_id = find(~cellfun(@isempty, regexp([followF{:}].', convertStringsToChars(fid))));
                
                if isempty(duiying_txt_id)
                    continue;
                    % 如果对应的txt不存在，说明再后续记录中，此神经元有可能因为信号太差而失去价值
                end
                followF = regexp(cellstr(stimulidirs),'(?<=F)\d*','match');
                followF = cellfun(@str2num, [followF{:}].','Uni',0);
                followF = [followF{:}].';
                duiying_stidir_id = find(followF == str2num(fid{1}));
                
                
                b = Batch(txtfiles{duiying_txt_id},plxdata{k}.inputpath,stimulidirs{duiying_stidir_id});
                
                
                this_channel = channelname;
                this_unit = unit;
                neu_list = {b.nlist.neuronname}.';
                
                channel_ids = find(~cellfun(@isempty,regexp(neu_list,this_channel)));
                unit_ids = find(~cellfun(@isempty,regexp(neu_list,sprintf('SPKC\\d{2}_%u',this_unit))));
                neuron_ids = intersect(channel_ids,unit_ids);
                
                b.select(neuron_ids);
                if isempty(b.sneu)
                    continue;
                end
                N = b.getn{1};
                % allocate sap-based feature information to each neuron
                
                
                feature_files = extract.filename(featuredir,'*.txt');
                
                datafile_id = find(~cellfun(@isempty, regexp(cellstr(feature_files),'Data')));
                infofile_id = find(~cellfun(@isempty, regexp(cellstr(feature_files),'Info')));
                % exist(feature_dir,'dir')
                if ~isempty(datafile_id) && ~isempty(infofile_id)
                    sorted_data = Neuron.extractFeaturesFromSapRawData(feature_files{datafile_id} , feature_files{infofile_id} );
                    N.setEachStimuliSapFeatures(sorted_data);
                    N.calMeanFeatures;
                end
                
                %duiying_song_folder
                Ns{k} = N;
                
            end
            
            A = Analysis(Ns);
            %A.calHarmRatio;
            
            
            
            save(A.formated_imagename,'A','-v7.3');
        end
        
        function A = genAnalysis(birdname,ZPid,channelname,unit) % generate Analysis
            arch = Archon('./');%Archon('D:/');
            if ~isempty(find(~cellfun(@isempty, regexp(cellstr(extract.filesAllLevel('./','*.mat')),...
                    sprintf('%s_%s_%s_%u.mat',birdname,ZPid,channelname,unit))))) % 如果当前folder已含有同名Analysis文件
                disp('该Neuron的分析已经存在');
                A = [];
                return
            end
            hitid = find(~cellfun(@isempty,regexp(cellstr(arch.bird_folders),birdname)));
            hitbird_folder = arch.bird_folders{hitid};
            
            
            % 处理hit folder 内部的 subfolders
            subdirs = extract.folder(hitbird_folder);
            hitZorP = find(~cellfun(@isempty,regexp(cellstr(subdirs.'),ZPid)));
            hit_subdir = subdirs {hitZorP};
            
            plxfiles = extract.filename(hit_subdir,'*.plx');
            zpids = find(~cellfun(@isempty, regexp(cellstr(plxfiles),'[ZP]\d+F\d+')));
            plxfiles = plxfiles(zpids);
            txtfiles = extract.filename(hit_subdir,'*.txt');
            zpids = find(~cellfun(@isempty, regexp(cellstr(txtfiles),'[ZP]\d+F\d+'))); % same record
            SPKCid = find(~cellfun(@isempty, regexp(cellstr(txtfiles),channelname))); % same channel
            Uid = find(~cellfun(@isempty, regexp(cellstr(txtfiles),sprintf('U%d',unit)))); % same unit
            if isempty(Uid)
                txtfiles = txtfiles(zpids);
            else
                txtfiles = txtfiles(intersect(intersect(zpids,SPKCid),Uid));
            end
            
            
            stimuli_collect_dir = append(hit_subdir,'\Stimuli');
            stimuli_dirs = extract.folder(stimuli_collect_dir ).';
            
            Ns = {};
            
            parfor k = 1: length(plxfiles)
                
                fid = regexp(plxfiles{k},'(?<=F)\d*','match'); % id for the file rank
                
                % find corresponding txt file
                followF = regexp(cellstr(txtfiles),'(?<=F)\d*','match');
                duiying_txt_id = find(~cellfun(@isempty, regexp([followF{:}].', convertStringsToChars(fid))));
                
                if isempty(duiying_txt_id)
                    continue;
                    % 如果对应的txt不存在，说明再后续记录中，此神经元有可能因为信号太差而失去价值
                end
                followF = regexp(cellstr(stimuli_dirs),'(?<=F)\d*','match');
                followF = cellfun(@str2num, [followF{:}].','Uni',0);
                followF = [followF{:}].';
                duiying_stidir_id = find(followF == str2num(fid));
                
                
                b = Batch(txtfiles{duiying_txt_id},plxfiles{k},stimuli_dirs{duiying_stidir_id});
                
                
                this_channel = channelname;
                this_unit = unit;
                neu_list = {b.nlist.neuronname}.';
                
                channel_ids = find(~cellfun(@isempty,regexp(neu_list,this_channel)));
                unit_ids = find(~cellfun(@isempty,regexp(neu_list,sprintf('SPKC\\d{2}_%u',this_unit))));
                neuron_ids = intersect(channel_ids,unit_ids);
                
                b.select(neuron_ids);
                if isempty(b.sneu)
                    continue;
                end
                N = b.getn{1};
                % allocate sap-based feature information to each neuron
                feature_dir = append(hit_subdir,'\Feature');
                
                feature_files = extract.filename(feature_dir,'*.txt');
                
                datafile_id = find(~cellfun(@isempty, regexp(cellstr(feature_files),'Data')));
                infofile_id = find(~cellfun(@isempty, regexp(cellstr(feature_files),'Info')));
                % exist(feature_dir,'dir')
                sorted_data = Neuron.extractFeaturesFromSapRawData(feature_files{datafile_id} , feature_files{infofile_id} );
                N.setEachStimuliSapFeatures(sorted_data);
                N.calMeanFeatures;
                
                %duiying_song_folder
                Ns{k} = N;
                
            end
            
            A = Analysis(Ns);
            %A.calHarmRatio;
            
            
            
            save(A.formated_imagename,'A','-v7.3'); % 保存的路径上不能存在中文
            
        end
        
        function input_roster = batch_genAnalysis(input_roster)% 生成Analysis object
            
            %clear
            tic
            D = parallel.pool.DataQueue;
            h = waitbar(0, '开始生成 Analysis objects');
            num_files = length(input_roster);
            % Dummy call to nUpdateWaitbar to initialise
            nUpdateWaitbar(num_files, h);
            % Go back to simply calling nUpdateWaitbar with the data
            afterEach(D, @nUpdateWaitbar);
            
            
            parfor k = 1: length(input_roster) % par
                
                % do stuff
                %waitbar(k/length(input_roster),bb,sprintf('此为花名册中%u个神经元中的%u',length(input_roster),k))
                birdname = input_roster(k).birdname;
                ZPid = input_roster(k).neuronid;
                channelname = input_roster(k).channelname;
                unit = input_roster(k).unit;
                try
                    Archon.genAnalysis(birdname,ZPid,channelname,unit);
                catch ME
                    input_roster(k).ME = ME;
                    disp(ME);
                end
                
                % Note we send only an "increment" for the waitbar.
                send(D, 1);
                
            end
            toc
            
            % subfunctiuon
            function p = nUpdateWaitbar(data, h)
                persistent TOTAL COUNT H
                if nargin == 2
                    % initialisation mode
                    H = h;
                    TOTAL = data;
                    COUNT = 0;
                else
                    % afterEach call, increment COUNT
                    COUNT = 1 + COUNT;
                    p = COUNT / TOTAL;
                    waitbar(p, H,sprintf('此为花名册中%u个神经元中的%u',TOTAL,COUNT));
                end
            end
            
        end
        
        function createSinChFolder(birdname,ZPid)
            arch = Archon('D:/');
            hitid = find(~cellfun(@isempty,regexp(cellstr(arch.bird_folders),birdname)));
            hitbird_folder = arch.bird_folders{hitid};
            
            
            % 处理hit folder 内部的 subfolders
            subdirs = extract.folder(hitbird_folder);
            hitZorP = find(~cellfun(@isempty,regexp(cellstr(subdirs.'),ZPid)));
            hit_subdir = subdirs {hitZorP};
            
            destineydir = convert.mergeSubfolders(fullfile(hit_subdir,'Stimuli'),'*.wav'); % 合并所有声文件
            
            convert.two2one(destineydir); % 生成单通道声文件
            
            
        end
        
        function batch_createSinChFolder(input_roster)
            
            wb = waitbar(0,'开始生成单通道声信号');
            for k = 1: length(input_roster)
                waitbar(k/length(input_roster),wb,sprintf('此为花名册中%u个神经元文件夹中的%u',length(input_roster),k))
                birdname = input_roster(k).birdname;
                ZPid = input_roster(k).neuronid;
                Archon.createSinChFolder(birdname,ZPid);
            end
            
        end
        
        function [neuroster,cdrfroster,malfunctions] = generateAnalysisAndAnalyze(datadir)% 生成Analysis object并分析
            dbstop if error
            % CDRF: Cons Degs Replas Frags
            [neuroster,malfunctions,cdrfroster,~,~] = Archon.extractAnalysisInfo(datadir);
            %[neuroster,malfunctions] = Archon.extractAnalysisInfo(datadir);
            %             cdrf_ids = find([neuroster.allexist].' == 1);
            %             cdrfroster = neuroster(cdrf_ids);
            
            for k = 1: length(neuroster)
                
                birdname = neuroster(k).birdname;
                ZPid = neuroster(k).neuronid;
                channelname =neuroster(k).channelname;
                unit = neuroster(k).unit;
                A = Archon.genAnalysis(birdname,ZPid,channelname,unit);
                A.saveSeparatedWaveform;
                for nid = 1:length(A.neurons)
                    a.neurons{nid}.rawthree;
                end
                %A.rawthree_NeuronClassVersion;
                A.sort_frags_by_response_strength_and_then_draw;
                A.drawMeanFeaturesVsRespAsLineChart;
                A.drawPairwiseFragmentsMeanFeaturesDistribution;
            end
            
        end
        
        function channel_unit_info = findSortedNeurons(path_txt) % from one file
            
            % path_txt = "D:\Ephys-O709-W\P07\O709_P07F5.txt"
            channel_unit_info = struct;
            sorted_tables = Spike.split(path_txt);
            
            for k = 1: length(sorted_tables)
                
                channel_unit_info(k).channelname = sorted_tables{k}.('channelname'){1};
                channel_unit_info(k).unit =sorted_tables{k}.('unit')(1);
                channel_unit_info(k).channel_unit = sprintf('%s_%u',channel_unit_info(k).channelname,channel_unit_info(k).unit);
            end
            
        end
        
        function shared = findConsistentNeurons(dir_txt)  % 找到在每一个文件里都有的神经元
            
            %dir_txt = "D:\Ephys-O709-W\P07"
            txtfile = extract.filename(dir_txt,'*.txt');
            
            for k = 1: length(txtfile)
                info = Archon.findSortedNeurons(txtfile{k});
                
                if length(fieldnames(info))== 0
                    shared = [];
                    return
                end
                
                consiss{k} = cellfun(@(x)convertCharsToStrings(x),{info.channel_unit}.','UniformOutput',1);
                
            end
            
            
            shared = consiss{1}; %找到所有的txt文件都拥有的一个channel+unit
            for m = 1: length(consiss)
                shared = intersect(shared,consiss{m});
            end
            
        end
        
        function conexist = findConExistNeurons(dir_txt)% 找到有Cons的神经元,这里使用F1来寻找Cons文件，后期有可能有问题
            % check which Fid is the Con session
            %dir_txt = "D:\Ephys-O709-W\P07"
            txtfile = extract.filename(dir_txt,'*.txt');
            f1txtid = find(~cellfun(@isempty,regexp(cellstr(txtfile),'F1')));
            subsetid = find(~cellfun(@isempty,regexp(cellstr(txtfile),'SPKC')));
            if ~isempty(subsetid)
                selected_txt = txtfile(intersect(f1txtid,subsetid));
                conexist = {};
                for k = 1: length(selected_txt)
                    info = Archon.findSortedNeurons(selected_txt{k});
                    conexist{k} = info.channel_unit;
                end
                
            else
                info = Archon.findSortedNeurons(txtfile{f1txtid});
                conexist = {info.channel_unit}.';
            end
            conexist = cellfun(@(x)convertCharsToStrings(x),{info.channel_unit}.','UniformOutput',1);
        end
        
        function [birddir,zpdir] = getDirpath(birdname,ZPid)
            arch = Archon('D:/');
            % birddir
            birddirs = cellstr(arch.bird_folders);
            birddir = birddirs{~cellfun(@isempty,regexp(birddirs,birdname))};
            
            % zplevel
            zpdirs = cellstr(extract.folder(birddir).');
            zpdir = zpdirs{~cellfun(@isempty,regexp(zpdirs,ZPid))};
            
        end
        
        % 处理已经生成的analysis
        
        function moveCorrespondingFiles()%basis, sourcedir, destineydir) %根据basis里的文件，移动有相同神经元编码的文件到指定文件夹
            dbstop if error
            sourcedir = "D:\CDRF_AnalysisObject";
            basis = extract.filename("C:\Users\Zhehao\Desktop\Selected_ResptoIndividual",'*.png');
            destineydir = "C:\Users\Zhehao\Desktop\Destiney";
            all_sourcefiles =cellstr( extract.filename(sourcedir,'*.*'));
            
            for k = 1: length(basis)
                
                extracted_id = regexp(basis{k},'[OBRGY]\d{3}_[ZP]\d+_SPKC\d+_\d','match');
                
                coresp_id = find(~cellfun(@isempty, regexp(all_sourcefiles,extracted_id)));
                
                for kk = 1: length(coresp_id)
                    oldpath = all_sourcefiles{coresp_id(kk)};
                    [old_dir,purename,ext] = fileparts(oldpath);
                    newpath = fullfile(destineydir,strcat(purename,ext));
                    movefile(oldpath,newpath);
                end
            end
        end
        
        
    end
    
    methods(Static) % deprecated methods, but might be useful later
        % 被弃用的一些方法
        function Deprecated_writeAnalysisObjectForMergedPlexonFiles(tablepath)
            
            %tablepath = "C:\Users\Zhehao\Desktop\AllInOne (27).xlsx";
            
            
            T = table2struct(readtable(tablepath));
            
            IDs = [T.UniqueID].';
            num_ids = rmmissing(unique(IDs(IDs~=0)));
            
            %kbad = [28,31,32,41,42,44,64,65];
            
            
            wb = waitbar(0,'Creating Neuron analysis objects');
            for k = 1: length(num_ids)
                waitbar(k/length(num_ids),wb,sprintf('%u of %u Neuron',k,length(num_ids)));
                ids_in_T = find(num_ids(k) == IDs);
                
                Ns = {};
                for i = 1:length(ids_in_T)
                    b = Batch(T(ids_in_T(i)).MergedTxtPath,T(ids_in_T(i)).MergedPlxPath,T(ids_in_T(i)).StimuliPath);
                    
                    this_channel = T(ids_in_T(i)).ChannelName;
                    this_unit = T(ids_in_T(i)).UnitName;
                    neu_list = {b.nlist.neuronname}.';
                    
                    channel_ids = find(~cellfun(@isempty,regexp(neu_list,this_channel)));
                    unit_ids = find(~cellfun(@isempty,regexp(neu_list,sprintf('_%u',this_unit))));
                    neuron_ids = intersect(channel_ids,unit_ids);
                    
                    if T(ids_in_T(i)).MergedIndex == 0 || isempty(T(ids_in_T(i)).MergedIndex ) % When the raw data is not a merged object (Or not treated as)
                        b.select(neuron_ids);
                        N = b.getn{1};
                    else
                        b.select(neuron_ids);
                        tempN = b.getn(T(ids_in_T(i)).MergedIndex);
                        N = tempN{1};
                    end
                    
                    %N.mergeIdx = T(ids_in_T(i)).mergedIdx;
                    N.signalGoodness = T(ids_in_T(i)).Goodness;
                    N.set_uniqueid(T(ids_in_T(i)).UniqueID);
                    % allocate sap-based feature information to each neuron
                    sorted_data = Neuron.extractFeaturesFromSapRawData(T(ids_in_T(i)).FeatureData, T(ids_in_T(i)).FeatureInfo);
                    N.setEachStimuliSapFeatures(sorted_data);
                    N.calMeanFeatures;
                    Ns{i} = N;
                end
                
                A = Analysis(Ns);
                A.uniqueid = num_ids(k);
                
                
                save(sprintf('%s_%u',A.birdid,A.uniqueid),'A','-v7.3');
                
                
            end
            
            close(wb);
            
            
            
            
        end
        
        function Deprecated_writeAnalysisObjectForSeparatedPlexonFiles(tablepath)
            
            T = table2struct(readtable(tablepath));
            IDs = [T.UniqueID].';
            num_ids = rmmissing(unique(IDs(IDs~=0)));
            
            %kbad = [28,31,32,41,42,44,64,65];
            
            
            wb = waitbar(0,'Creating Neuron analysis objects');
            for k = 1: length(num_ids)
                waitbar(k/length(num_ids),wb,sprintf('%u of %u Neuron',k,length(num_ids)));
                ids_in_T = find(num_ids(k) == IDs);
                
                Ns = {};
                for i = 1:length(ids_in_T)
                    b = Batch(T(ids_in_T(i)).TxtPath,T(ids_in_T(i)).PlxPath,T(ids_in_T(i)).StimuliPath);
                    
                    this_channel = T(ids_in_T(i)).ChannelName;
                    this_unit = T(ids_in_T(i)).UnitName;
                    neu_list = {b.nlist.neuronname}.';
                    
                    channel_ids = find(~cellfun(@isempty,regexp(neu_list,this_channel)));
                    unit_ids = find(~cellfun(@isempty,regexp(neu_list,sprintf('_%u',this_unit))));
                    neuron_ids = intersect(channel_ids,unit_ids);
                    
                    b.select(neuron_ids);
                    N = b.getn{1};
                    if isfield(T,'Goodness')
                        N.signalGoodness = T(ids_in_T(i)).Goodness;
                    end
                    N.set_uniqueid(T(ids_in_T(i)).UniqueID);
                    % allocate sap-based feature information to each neuron
                    if ~isnan(T(ids_in_T(i)).FeatureData)
                        sorted_data = Neuron.extractFeaturesFromSapRawData(T(ids_in_T(i)).FeatureData, T(ids_in_T(i)).FeatureInfo);
                        N.setEachStimuliSapFeatures(sorted_data);
                        N.calMeanFeatures;
                    end
                    
                    Ns{i} = N;
                end
                
                A = Analysis(Ns);
                A.uniqueid = num_ids(k);
                
                
                save(sprintf('%s_%u',A.birdid,A.uniqueid),'A','-v7.3');
                
                
            end
            
            close(wb);
        end
        
        function all_data_info = Deprecated_AnalyzeFormatedEphysData(dirpath)
            % Automatically generate plx-txt-folder table based on
            % formated Ephys folders
            % 202204211 11:44 pm paused here!
            
            dbstop if error
            bird_folders = [extract.folder(dirpath)].'; %#ok<NBRAK>
            
            pure_folname = {}; % 初始化
            for k = 1: length(bird_folders)
                temp1 = split(bird_folders{k},'\');
                pure_folname{k} = cellstr(temp1{end});
            end
            
            contain_ephys_ids = find(~cellfun(@isempty, regexp([pure_folname{:}].','Ephys')));
            if ~isempty(contain_ephys_ids)
                bird_folders =  bird_folders(contain_ephys_ids);
            else
                bird_folders =  dirpath;
            end
            
            infocollect = {};
            for k = 1:length(bird_folders)
                plxfiles = extract.filename(bird_folders{k},'*.plx');
                txtfiles = extract.filename(bird_folders{k},'*.txt');
                
                if exist(sprintf('%s\\Stimuli',bird_folders{k}), 'dir')
                    stidirs = extract.folder(sprintf('%s\\Stimuli',bird_folders{k})).';
                else % 如果未把所有stimuli放进stimuli folder里面
                    stidirs = extract.folder(sprintf('%s',bird_folders{k})).'
                end
                
                plxname = {}; % extract .plx filenames
                for m = 1: length(plxfiles)
                    [~,plxname{m},~] = fileparts(plxfiles{m});
                end
                
                txtname = {};% extract .txt filenames
                for m = 1: length(txtfiles)
                    [~,txtname{m},~] = fileparts(txtfiles{m});
                end
                
                stiname = {}; % extract stimuli folder names
                for m = 1: length( stidirs)
                    temps = split( stidirs{m},'\');
                    stiname{m} = temps{end};
                end
                
                % Based on .plx filenames, find the corresponding .txt and
                % stimuli folder names
                datainfo = struct;
                for m = 1: length(plxname)
                    
                    txt_hit = find(~cellfun(@isempty, regexp([cellstr(txtname)].',plxname{m})));
                    
                    sti_hit = find(~cellfun(@isempty, regexp([cellstr(stiname)].',plxname{m})));
                    
                    % fill the struct
                    sadtemp = regexp(plxname{m}, '[_-]', 'split');
                    datainfo(m).animal = sadtemp{end};
                    datainfo(m).plxname = plxname{m};
                    datainfo(m).keyword = regexp(plxname{m},'Con|Deg|Detail|Frag|Repla|Mir|Rev','match');
                    datainfo(m).sectionname = regexp(plxname{m},'Z\d{2}|P\d{2}','match');
                    fuck = regexp(cellstr(datainfo(m).sectionname),'\d{2}','match');
                    datainfo(m).sectionid = str2num(fuck{1}{1});
                    datainfo(m).subfilename = regexp(plxname{m},'F\d{1}','match');
                    if ~isempty(datainfo(m).subfilename) %如果文件名里存在F标识subfile的id
                        fuck2 = regexp(datainfo(m).subfilename,'\d','match');
                        datainfo(m).subfileid = str2num(fuck{1}{1});
                    elseif ~isempty(regexp(datainfo(m).plxname,'Mrg|mrg'))% 如果不存在F标识却存在‘Mrg|mrg’标识
                        datainfo(m).subfilename = "Mrg";
                    else % 如果不存在任何标识
                        datainfo(m).subfilename = [];
                    end
                    
                    datainfo(m).plxpath = plxfiles{m};
                    if ~isempty(txt_hit)
                        datainfo(m).txtpath = txtfiles{txt_hit};
                    else
                        datainfo(m).txtpath = [];
                    end
                    
                    if ~isempty(sti_hit)
                        datainfo(m).stidirpath = stidirs{sti_hit};
                    else
                        datainfo(m).stidirpath  = [];
                    end
                    
                end
                
                infocollect{k} = datainfo;
                
            end
            
            
            all_data_info = horzcat(infocollect{:});
            
            
        end
        
        function As = Deprecated_writeSingleAnalysiFile(txt,plx, folder,feature_data,feature_info)
            
            b = Batch(txt,plx,folder);
            b.select;
            Ns = b.getn;
            
            As = {};
            for k = 1: length(Ns)
                N = Ns{k};
                sorted_data = Neuron.extractFeaturesFromSapRawData(feature_data,feature_info);
                N.setEachStimuliSapFeatures(sorted_data);
                N.calMeanFeatures;
                A = Analysis(N);
                As{k} = A;
            end
            
        end
        
        function Deprecated_drawThreePlotsForAllFormattedNeurons(dir_path)
            
            dbstop if error
            
            if ~exist('dir_path','var')
                dir_path = './';
            end
            all_data_info = Sultan.AnalyzeFormatedEphysData(dir_path);
            wb = waitbar(0,'Start processing');
            
            for m = 1:length(all_data_info)
                
                waitbar(m/length(all_data_info),wb,sprintf('%u of totally %u files',m,length(all_data_info)));
                path_txt = all_data_info(m).txtpath;
                path_plx = all_data_info(m).plxpath;
                path_folder = all_data_info(m).stidirpath;
                
                if isempty(path_folder )&& ~isempty(regexp(path_plx,'Mrg|mrg'))
                    
                    sameanimal = find(all_data_info(m).animal == [all_data_info.animal].');
                    sameZfile = find(ismember([all_data_info.sectionname].',all_data_info(m).sectionname));
                    %sameZfile = find(ismember(all_data_info(m).sectionname, [all_data_info.sectionname].'));
                    sameneuronids = intersect(sameanimal,sameZfile);
                    path_folder = [all_data_info(sameneuronids).stidirpath].';
                    
                end
                
                b = Batch(path_txt,path_plx,path_folder);
                b.select;
                neuronlist = b.getn;
                for k = 1: length(neuronlist)
                    thisn = neuronlist{k};
                    thisn.rawthree;
                end
                
            end
            
            close(wb);
            
        end
        
        function Deprecated_AnalyzeNeuronsWhichIsComplete(dirpath)
            if ~exist('dirpath','var')
                dirpath = './';
            end
            final_info = Sultan.generateAFilesInputStructFromFormattedNeurons(dirpath);
            selected_info = final_info([final_info.complete] == 1);
        end
        
        function final_info = Deprecated_generateAFilesInputStructFromFormattedNeurons(dirpath)
            
            dbstop if error
            
            final_info = struct;
            count = 0;
            if ~exist('dir_path','var')
                dir_path = './';
            end
            all_data_struct = Sultan.AnalyzeFormatedEphysData(dir_path);
            
            animals = unique({all_data_struct.animal}.');
            for ani = 1: length(animals)
                
                this_animal_struct = all_data_struct(strcmp({all_data_struct.animal},  animals{ani}));
                sectionname = unique([this_animal_struct.sectionname].');
                
                for sec = 1:length(sectionname)
                    
                    this_section_struct = this_animal_struct(strcmp([this_animal_struct.sectionname],  sectionname{sec}));
                    subfilename = [this_section_struct.subfilename].';
                    
                    count = count + 1;
                    final_info(count).animal = animals{ani};
                    final_info(count).sectionname = sectionname{sec};
                    
                    
                    % 如果存在‘Mrg’标识，那么只基于Mrg文件生成Aobject
                    if ismember('Mrg',subfilename)
                        
                        [~,idx] = ismember('Mrg',subfilename);
                        final_info(count).path_txt = this_section_struct(idx).txtpath;
                        final_info(count).path_plx = this_section_struct(idx).plxpath;
                        final_info(count).path_folder = this_section_struct(idx).stidirpath;
                        final_info(count).path_pairs = [];
                    else
                        
                        for k = 1: length(this_section_struct)
                            
                            final_info(count).path_txt = [];
                            final_info(count).path_plx = [];
                            final_info(count).path_folder =[];
                            final_info(count).path_pairs{k}{1} = this_section_struct(k).txtpath;
                            final_info(count).path_pairs{k}{2} = this_section_struct(k).plxpath;
                            final_info(count).path_pairs{k}{3} = this_section_struct(k).stidirpath;
                            
                        end
                        
                    end
                    
                    % 表明这个neuron经历的stimuli都有哪些
                    final_info(count).keyword = [this_section_struct.keyword];
                    
                    disp(final_info(count).keyword)
                    
                    if ismember("Detail",final_info(count).keyword)||ismember("Frag",final_info(count).keyword)
                        final_info(count).complete = 1;
                    else
                        final_info(count).complete = 0;
                    end
                    
                    %否则，对每一个subfile生成Nobject,再基于所有Nobject生成Aobject
                end
                
                
            end
            
        end
    end
    
end

