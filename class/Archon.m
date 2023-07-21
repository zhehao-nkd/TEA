classdef Archon
    % batch processing the raw data( .pl2, .txt, .wav files sorted within
    % folders
    % 批量对raw data进行操作 ， 用以生成 Neuron object

    properties
        bird_folders
        zpfolders
        all_data_info
    end

    methods

        function a = Archon(dirpath)

            dbstop if error
            temp = Extract.foldersAllLevel(dirpath);
            a.bird_folders = {temp{find(~cellfun(@isempty,regexp(temp,'[ZP]\d+$')))}}.'; % folders ended with ZP and number

            infocollect = {};
            for k = 1:length(a.bird_folders)
                pl2files = Extract.filename(a.bird_folders{k},'*.pl2');
                txtfiles = Extract.filename(a.bird_folders{k},'*.txt');

                if exist(sprintf('%s\\Stimuli',a.bird_folders{k}), 'dir')
                    stidirs = Extract.folder(sprintf('%s\\Stimuli',a.bird_folders{k})).';
                else % 如果未把所有stimuli放进stimuli folder里面
                    stidirs = Extract.folder(sprintf('%s',a.bird_folders{k})).'
                end

                pl2name = {}; % Extract .pl2 filenames
                for m = 1: length(pl2files)
                    [~,pl2name{m},~] = fileparts(pl2files{m});
                end

                txtname = {};% Extract .txt filenames
                for m = 1: length(txtfiles)
                    [~,txtname{m},~] = fileparts(txtfiles{m});
                end

                stiname = {}; % Extract stimuli folder names
                for m = 1: length( stidirs)
                    temps = split( stidirs{m},'\');
                    stiname{m} = temps{end};
                end

                % Based on .pl2 filenames, find the corresponding .txt and
                % stimuli folder names
                datainfo = struct;
                for m = 1: length(pl2name)

                    txt_hit = find(~cellfun(@isempty, regexp([cellstr(txtname)].',pl2name{m})));

                    sti_hit = find(~cellfun(@isempty,...
                        regexp(cellfun(@Convert.fileid,[cellstr(stiname)].','Uni',0),Convert.fileid(pl2name{m}) )));

                    % fill the struct
                    sadtemp = regexp(pl2name{m}, '[_-]', 'split');
                    datainfo(m).animal = sadtemp{end};
                    datainfo(m).pl2name = pl2name{m};
                    datainfo(m).keyword = regexp(pl2name{m},'Con|Deg|Detail|Frag|Repla|Mir|Rev','match');
                    datainfo(m).sectionname = regexp(pl2name{m},'Z\d{2}|P\d{2}','match');
                    fuck = regexp(cellstr(datainfo(m).sectionname),'\d{2}','match');
                    datainfo(m).sectionid = str2num(fuck{1}{1});
                    datainfo(m).subfilename = regexp(pl2name{m},'F\d{1}','match');
                    if ~isempty(datainfo(m).subfilename) %如果文件名里存在F标识subfile的id
                        fuck2 = regexp(datainfo(m).subfilename,'\d','match');
                        datainfo(m).subfileid = str2num(fuck{1}{1});
                    elseif ~isempty(regexp(datainfo(m).pl2name,'Mrg|mrg'))% 如果不存在F标识却存在‘Mrg|mrg’标识
                        datainfo(m).subfilename = "Mrg";
                    else % 如果不存在任何标识
                        datainfo(m).subfilename = [];
                    end

                    datainfo(m).pl2path = pl2files{m};
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

            files = Extract.filesAllLevel(targetdir, ext);

            wb = waitbar(0,'Start renaming files...');

            for id = 1:length(files)

                waitbar(id/length(files),wb,sprintf('Processing %u of %u Files',id,length(files)) )
                % Get the file name
                oldname = files{id};

                [folder,oldfilename,ext] = fileparts(oldname);
                newfilename = strrep(oldfilename,'-','_');
                parts = unique(split(newfilename,'_'));

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
            dirs = Extract.foldersAllLevel(targetdir);
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
                index2 = find(~cellfun(@isempty, regexp(parts,'^[PZ]\d{2}F\d{1}$'))); % ^EXP$ is perfect match

                if isempty(index)&& ~isempty(index2)
                    temp = regexp(oldname,'[OBRYG]\d{3}','match');
                    parts{length(parts) + 1} = temp{1};
                    index = find(~cellfun(@isempty, regexp(parts,'^[OBRYGX]\d{3}$'))); % ^EXP$ is perfect match
                end
                % 替换index位和第一位的顺序
                if ~isempty(index)
                    temp = parts{1};
                    parts{1} = parts{index};
                    parts{index} = temp;
                end


                % 替换index2位和第二位的顺序
                if ~isempty(index2)
                    temp = parts{2};
                    parts{2} = parts{index2};
                    parts{index2} = temp;
                end


                newdirname = strjoin(parts,'_');

                dirranks{length(dirranks)} = newdirname;
                newname = fullfile(dirranks{:});

                if ~strcmp(oldname,newname)
                    movefile(oldname, newname);
                end
            end
            close(wb);

        end

        function reorderNames(targetdir)
            Convert.bid_first_AllLevel(targetdir, '*.pl2');
            Convert.bid_first_AllLevel(targetdir, '*.txt');
            Convert.bid_first_AllLevel_dirVer(targetdir);
        end

        function padFilename(targetdir, ext) % format data into G573_P01F2_Degs kind of format
            files = Extract.filesAllLevel(targetdir, ext);

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
                    %                     FS = filesep;
                    %                     dirparts = split(targetdir,FS);
                    %                     ephysid = find(~cellfun(@isempty, regexp(dirparts,'Piece')));
                    %                     if ~isempty(ephysid)
                    %                         bid = regexp(dirparts{ephysid},'[OBRYGX]\d{3}','match');
                    %                         bid = bid{1};
                    %                     else
                    bid = convertStringsToChars(regexp(oldname,'[OBRYGX]\d{3}','match'));
                    %                     end
                    %                     disp(parts);
                    %
                    parts = [bid;parts]; % parts means the parts of the filename, which is modified by this code

                end


                index2 = find(~cellfun(@isempty, regexp(parts,'^[PZ]\d{1,2}F\d{1}$'))); % ^EXP$ is perfect match


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

        function uncompleted_folders = verifyIfFeatureFilesAreCompleted(parentdir)
            % 在使用sap生成acoustic feature
            % files的时候，会出现随机的数据丢失，似乎多集中于以repla命名的数据里面，
            %这个function的作用就是找出有数据丢失的zpfolders
            uncompleted_folders = struct;
            count = 0;

            temp = Extract.foldersAllLevel(parentdir);
            temp = temp(find(cellfun(@isempty,regexp(cellstr(temp),'\$')))); %去除被隐藏的文件夹
            zp_folders = {temp{find(~cellfun(@isempty,regexp(temp,'[ZP]\d+$')))}}.'; % folders ended with ZP and number

            for k = 1:length(zp_folders)

                num_wav = length(Extract.filename(fullfile(zp_folders{k},'SingleChannel_Merged_Stimuli'),'*.wav'));

                txtfiles = Extract.filename(fullfile(zp_folders{k},'Feature'),'*.txt');
                infopath = txtfiles{find(~cellfun(@isempty,regexp(txtfiles,'Info')))};
                fid = fopen(infopath);
                lines =  textscan(fid,'%s','delimiter','\n');
                lines = lines{1};
                num_records = length(lines) - 1;

                if num_records < num_wav
                    count = count + 1;
                    uncompleted_folders(count).foldername = zp_folders{k};
                    uncompleted_folders(count).numfile = num_wav;
                    uncompleted_folders(count).numrecord = num_records;
                end

            end

            uncompleted_folders = uncompleted_folders.';


        end
 

        function error = reCreate(targetdir)
            % 在原位置重新生成AnalysisObject，当Analysis的生成代码发生重要变化的时候
            error = struct;
            matfiles = Extract.filesAllLevel(targetdir, '*.mat');
            unpack = @(x) x{1};
            wb = PoolWaitbar(length(matfiles),'重新生成Neurons');
            parfor k = 1:length(matfiles) % parfor
                try
                    birdname = unpack(regexp(matfiles{k},'[OGBYR]\d{3}','match'));
                    zpid = unpack(regexp(matfiles{k},'[ZP]\d{2}','match'));
                    SPKC = unpack(regexp(matfiles{k},'SPKC\d{2}','match'));
                    channel = str2num(unpack(regexp(matfiles{k},'(?<=SPKC\d{2}_)\d','match')));

                    A = Archon.genAnalysis('D:',birdname,zpid,SPKC,channel,0);
                    INF = A.getNinfo;
                    Archon.parsave_reCreate(A,INF,matfiles{k});
                    %                     save(matfiles{k},'A','INF','-v7.3'); % 保存的路径上不能存在中文
                catch ME
                    error(k).num = k;
                    error(k).ME = ME;
                    error(k).identifier = ME.identifier;
                    error(k).matfile = matfiles{k};
                end

                increment(wb);

            end

  
        end

        function parsave_reCreate(A,INF,filename)

            save(filename,'A','INF','-v7.3');

        end

        function standardizeRawname(targetdir) % 标准化文件夹里所有raw data的名称
            Archon.padFilename(targetdir, '*.txt')
            Archon.padFilename(targetdir, '*.pl2')
            Archon.reorderFilename(targetdir, '*.pl2');
            Archon.reorderFilename(targetdir, '*.txt');
            Archon.reorderDirname(targetdir);

        end

        % 生成analysis object

        function genAnalysisFromSingleZPDir(datadir)
            % generate Neuron object from a single ZP-marked raw data folder
            [neuroster,~,~] = Archon.extractAnalysisInfo(datadir);
            Archon.batch_genAnalysisSameRecordingFile(neuroster);

        end

        function [error,error2] =  Deprecated_batch_genAnalysisSameRecordingFile(input_roster, start_from)
            dbstop if error
            % 从同一批原始pl2数据中生成所有的analysis文件
            %采用此方法的好处是不用一次又一次地调用plexon读取方法，可以节省非常多的时间
            tic
            D = parallel.pool.DataQueue;
            h = waitbar(0, '开始生成 Neuron objects');
            num_files = length(input_roster);% Dummy call to nUpdateWaitbar to initialise
            Utl.UpdateParforWaitbar(num_files, h);% Go back to simply calling nUpdateWaitbar with the data
            afterEach(D, @Utl.UpdateParforWaitbar);
            neurondir =  unique({input_roster.neurondir}.');

            error = {};
            for mm = start_from: length(neurondir) % par  % 有遗漏的（因为bug，已修复，未测试）是11，21，30，20221002，ZC


                ndir_id = find(strcmp(cellstr({input_roster.neurondir}.'),neurondir{mm}));
                hit_subdir = neurondir{mm};
                pl2files = Extract.filename(hit_subdir,'*.pl2');
                pl2files = {pl2files{find(cellfun(@isempty, regexp(pl2files,'Mrg|mrg')))}}.'; % not include id with Mrg in text
                pl2data = {};
                try
                    for k = 1: length(pl2files)
                        pl2data{k} = Trigger(pl2files{k});
                    end
                catch ME
                    error2(mm).num = mm;
                    error2(mm).ME = ME;
                    continue
                end

                parfor k = 1: length(ndir_id) % par

                    try
                        % do stuff
                      

                        birdname = input_roster(ndir_id(k)).birdname;
                        ZPid = input_roster(ndir_id(k)).neuronid{1};
                        channelname = input_roster(ndir_id(k)).channelname;
                        unit = input_roster(ndir_id(k)).unit;

                        if isfile(sprintf('%s_%s_%s_%u.mat',birdname,ZPid,channelname,unit))
                            continue
                        end

                        hit_subdir = neurondir{mm};

                        featuredir = append(hit_subdir,'\Feature');
                        txtfiles = Extract.filename(hit_subdir,'*.txt');
                        zpids = find(~cellfun(@isempty, regexp(cellstr(txtfiles),'[ZP]\d+F\d+'))); % same record
                        SPKCid = find(~cellfun(@isempty, regexp(cellstr(txtfiles),channelname))); % same channel
                        Uid = find(~cellfun(@isempty, regexp(cellstr(txtfiles),sprintf('U%d',unit)))); % same unit
                        if isempty(Uid)
                            txtfiles = txtfiles(zpids);
                        else
                            txtfiles = txtfiles(intersect(intersect(zpids,SPKCid),Uid));
                        end

                        stimuli_collect_dir = append(hit_subdir,'\Stimuli');
                        stimulidirs = Extract.folder(stimuli_collect_dir ).';

                        Archon.genAnalysisDifferentInput(pl2data,txtfiles,stimulidirs,featuredir,channelname,unit);
                        % Note we send only an "increment" for the waitbar.

                    catch ME
                        send_mail_message('chengzh.nkd@gmail.com','Here comes a bug','Fix it!!!!!!')
                        error{mm,k} = struct;
                        error{mm,k}.ME =ME;
                        error{mm,k}.message =ME.message;
                        error{mm,k}.birdname =birdname ;
                        error{mm,k}.ZPid = ZPid;
                        error{mm,k}.channelname = channelname;
                        error{mm,k}.unit = unit;
                        %                         try
                       % send_mail_message('chengzh.nkd@gmail.com',ME.message,sprintf('%s_%s_%s_%u',birdname{1},ZPid{1},channelname,unit)) ;
                        %                         catch
                        %                             send_mail_message('chengzh.nkd@gmail.com','Bug','Wierd bug') ;
                        %                         end
                        error{mm,k}.unit = unit;
                        %pause
                        disp('Error!!!!!!!!!!!!')
                    end
                    
                    send(D, 1);
                end

            end


            toc

        end

        function [neuroster,malfunctions,cdrfroster] = extractAnalysisInfo(datadir)
            %[neuroster,malfunctions,cdrfroster,noSinChan_neuroster,nofeature_neuroster] = extractAnalysisInfo(datadir)
            dbstop if error
            % 分析所包含的所有Z P folder,生成Analysis对象，并且对于那些不可分析的（某些信息不够），生成记录文件以便补充这些信息

            allpl2s = Extract.filesAllLevel(datadir,'*.pl2');
            allfolders = {}; %把是否存在pl2文件当作是否是处理好的neuron原始数据的文件夹的标准
            for k = 1:length(allpl2s)
                [part1,~,~] = fileparts(allpl2s{k});
                allfolders{k} = part1;
            end
            neurondir = unique(allfolders).';
            malfunctions = struct; % malfunctions
            malcount = 0;

            % waitbars
            wb2 = waitbar(0,'Experiment of this bird');
            set(wb2,'doublebuffer','on');
            summer = {};

            unpack = @(x) x{1};
            for kk = 1: length(neurondir) % iterate 对每一个 ZP folder

                zpname = unpack(unpack(regexp(cellstr(neurondir{kk}),'[ZP]\d+','match')));

                waitbar(kk/length(neurondir),wb2,sprintf('正处理 %u个神经元文件夹中的%u',length(neurondir) ,kk))

                %判断 pl2文件是否完备
                pl2files = Extract.filename(neurondir{kk},'*.pl2');
                mergeids = find(~cellfun(@isempty, regexp(cellstr(pl2files),'merge')));
                %zpids = find(~cellfun(@isempty, regexp(cellstr(pl2files),'F\d+')));

                if ~isempty(mergeids)
                    valid_pl2files = pl2files(mergeids); % when there is a merged file
                else
                    valid_pl2files = pl2files;%(setdiff(zpids,mergeids)); % 剔除含有merge关键词的 pl2文件
                end


                stimuliparentdir = fullfile(neurondir{kk},'Stimuli');
                stimulidir = cellstr(Extract.folder(stimuliparentdir));
                stimuli_zpids = find(~cellfun(@isempty, regexp(cellstr(stimulidir),'F\d+'))); %剔除命名不合规范的stimuli文件夹
                valid_stimulidir = cellstr(stimulidir(stimuli_zpids));
                subroster = struct('fullname',{}); % Experiment 花名册
                malfunctions = struct; % malfunctions
                count = 0;


                
                for index = 1: length(valid_pl2files)  % 每一个Pl2文件

                    if  ~isempty(regexp(valid_pl2files{index},'merge'))
                        coresp_stimulidirs = valid_stimulidir;
                    else
                        hitted_stimulidir_ids = find(~cellfun(@isempty, regexp(cellstr(valid_stimulidir),...
                            unpack(regexp(valid_pl2files{index},'F\d+','match')))  ));
                        if isempty(hitted_stimulidir_ids)
                            continue
                        else
                            coresp_stimulidirs = valid_stimulidir(hitted_stimulidir_ids);
                        end
                    end

                    pl2info = PL2GetFileIndex(valid_pl2files{index});

                    cated_spikechannels = vertcat(pl2info.SpikeChannels{:}); %找到已经sorted spike channels

                    try
                        num_unit = 0;
                        names_spikechannels = {};
                        for m = 1:length(cated_spikechannels)
                            unitids = 1:1:max(find(cated_spikechannels(m).UnitCounts))-1;
                            for n = 1:length(unitids)
                                num_unit = num_unit + 1;
                                names_spikechannels{num_unit} = sprintf('%s_%s_%s_%02u',...
                                    regexp(neurondir{kk},'[BRGYOX]\d{3}','match','once'),...
                                    zpname,...
                                    unpack(unpack(regexp(cellstr({cated_spikechannels(m).Name}.'),'(?<=SPK_)SPKC\d+','match'))),...
                                    unitids(n));
                            end
                        end
                    catch
                        malcount = malcount + 1;
                        malfunctions(malcount).birdname = regexp(neurondir{kk},'[BRGYOX]\d{3}','match','once');
                        malfunctions(malcount).neuronid = regexp(neurondir{kk},'[ZP]\d+','match');
                        malfunctions(malcount).reason = 'Spike sorting 未完成';
                        continue
                    end



                    for cat = 1:length(names_spikechannels)
                        bird_channelunit = names_spikechannels{cat};
                        %                         sprintf('%s_%s',regexp(neurondir{kk},'[BRGYOX]\d{3}','match','once'),...
                        %                            names_spikechannels{cat} );

                        if  any(ismember(bird_channelunit, {subroster.fullname}.'))
                            [~,where_it_is] = ismember(bird_channelunit, {subroster.fullname}.');
                            subroster(where_it_is).pl2files{length(subroster(where_it_is).pl2files) + 1,1} = valid_pl2files(index);
                            subroster(where_it_is).stimulidirs{size(subroster(where_it_is).stimulidirs,1) + 1,1}  =  cellstr(coresp_stimulidirs);
                        else
                            disp(kk)
                            count = count + 1;

                            subroster(count).neurondir = neurondir{kk};
                            subroster(count).fullname = bird_channelunit;
                            subroster(count).pl2files = valid_pl2files(index);
                            subroster(count).stimulidirs =  cellstr(coresp_stimulidirs);
                        end

                    end

                end

                if ~isempty(subroster)
                    summer{kk} = subroster;
                end
            end

            neuroster = horzcat(summer{:});

            for k = 1:length(neuroster)
                temp = vertcat(neuroster(k).pl2files{:});
                neuroster(k).pl2files = cellstr(temp);
                try
                    temp = vertcat(neuroster(k).stimulidirs{:});
                    neuroster(k).stimulidirs = cellstr(temp);
                catch
                end
             
                neuroster(k).birdname = regexp(neuroster(k).fullname,'[BRGYOX]\d{3}','match','once');
                neuroster(k).neuronid = unpack(regexp(neuroster(k).neurondir,'[ZP]\d+','match'));
                neuroster(k).zpid =  unpack(regexp(neuroster(k).neurondir,'[ZP]\d+','match'));
                neuroster(k).birdneuron = strcat(neuroster(k).birdname,neuroster(k).neuronid);
                neuroster(k).channelname = unpack(regexp(neuroster(k).fullname,'SPKC\d+','match'));
                neuroster(k).unit = str2num(unpack(regexp(neuroster(k).fullname,'(?<=_)\d+','match')));
                [~,dirname, ~] = fileparts(neuroster(k).neurondir);
                neuroster(k).formated_name = sprintf('%s_%s_%s_%u_%s',neuroster(k).birdname,neuroster(k).zpid,...
                    neuroster(k).channelname,neuroster(k).unit,dirname);
                stimulidir = fullfile(neuroster(k).neurondir,'Stimuli');
                stimuli_filenames = Extract.filesAllLevel(stimulidir,'*.wav');

                neuroster(k).frag = 0;  % 判断是否有frag
                if ~isempty(find(~cellfun(@isempty,regexp(stimuli_filenames,'Frag')), 1))
                    neuroster(k).frag = 1;
                end

                neuroster(k).deg = 0;  % 判断是否有deg
                if ~isempty(find(~cellfun(@isempty,regexp(stimuli_filenames,'Deg|deg')), 1))
                    neuroster(k).deg = 1;
                end

                neuroster(k).repla = 0;  % 判断是否有repla
                if ~isempty(find(~cellfun(@isempty,regexp(stimuli_filenames,'Repla|repla')), 1))
                    neuroster(k).repla = 1;
                end

                neuroster(k).wns = 0;  % 判断是否有WNS
                if ~isempty(find(~cellfun(@isempty,regexp(stimuli_filenames,'wns|WNS')), 1))
                    neuroster(k).wns = 1;
                end


                neuroster(k).allexist = 0;  % 判断是否同时存在 norm frag deg repla
                if ~isempty(find(~cellfun(@isempty,regexp(stimuli_filenames,'Norm|norm')), 1)) &&...
                        ~isempty(find(~cellfun(@isempty,regexp(stimuli_filenames,'Frag|frag')), 1)) &&...
                        ~isempty(find(~cellfun(@isempty,regexp(stimuli_filenames,'deg|Deg')), 1)) &&...
                        ~isempty(find(~cellfun(@isempty,regexp(stimuli_filenames,'Repla|repla')), 1))
                    neuroster(k).allexist = 1;
                end

                neuroster(k).featureexist = 0;  % 是否存在 feature files
                feature_dir = fullfile(neurondir{kk},'Feature');
                if exist(feature_dir,'dir')
                    feature_files = Extract.filename(feature_dir,'*.txt');

                    if ~isempty(find(~cellfun(@isempty, regexp(cellstr(feature_files),'Info')))) &&...
                            ~isempty(find(~cellfun(@isempty, regexp(cellstr(feature_files),'Data'))))
                        neuroster(k).featureexist = 1;
                    end
                end

                neuroster(k).SCDexist = 0; % 是否存在SingleChannel_Merged_Stimuli文件夹
                subdirs = Extract.folder(neurondir{kk});
                if ~isempty(find(~cellfun(@isempty, regexp(cellstr(subdirs).','SingleChannel_Merged_Stimuli'))))
                    neuroster(k).SCDexist = 1;
                end
            end

            if isfield(neuroster,'allexist')
                cdrfroster = neuroster([neuroster.allexist].' == 1);
            else
                cdrfroster = [];
            end

            %             %not_cdrfroster = neuroster([neuroster.allexist].' == 0);
            %             noSinChan_neuroster = neuroster([neuroster.SCDexist].' == 0);
            %             nofeature_neuroster = neuroster([neuroster.featureexist].' == 0);

        end

        function A = genAnalysis_BasedOnRoster(roster)
            % input 是 neuroster的一行
            dbstop if error
            valid_pl2files = roster.pl2files;
            stimuli_dirs =  roster.stimulidirs;
            this_channel = roster.channelname;
            this_unit = roster.unit;
            experiments = {};
            count = 0;
            getfileid = @(x) str2num(regexp(convertCharsToStrings(x),'(?<=F)\d','match'));
            for k = 1: length(valid_pl2files) % parfor


                if ~isempty(... % 一旦目标文件是merge
                        find(~cellfun(@isempty, regexp(cellstr(valid_pl2files{k}),'merge'))))
                    b = Chorus(valid_pl2files{k},stimuli_dirs);
                else

                    fid = getfileid(valid_pl2files{k});
                    followF = cellfun(@(x) getfileid(x),stimuli_dirs);
                    duiying_stidir_id = find(followF == fid);%find(followF == str2num(fid{1}));
                    if ~isempty(duiying_stidir_id)
                        b = Chorus(valid_pl2files{k},stimuli_dirs{duiying_stidir_id});
                    else
                        continue
                    end
                end

                neu_list = {b.nlist.neuronname}.';

                channel_ids = find(~cellfun(@isempty,regexp(neu_list,this_channel)));
                unit_ids = find(~cellfun(@isempty,regexp(neu_list,sprintf('SPKC\\d{2}_%u',this_unit))));
                neuron_ids = intersect(channel_ids,unit_ids);

                b.select(neuron_ids);
                if isempty(b.sneu)
                    continue;
                end
                Exp = b.getExperiments{1};

                CALCULATE_FEATURE = 0;

                if CALCULATE_FEATURE == 1
                    Exp.setFeatures_SAT;
                    Exp.calMeanFeatures;
                end
                
                count = count + 1;
                experiments{count} = Exp;
            end
            %eleinf = load("F:\S01_GeneratedNeurons_20221209\all_eleinf.mat").eleinf;
            %A = Neuron(experiments,eleinf);
            A = Neuron(experiments);
            info = A.info;
            song = A.song;
            deg = A.deg;
            frag = A.frag;
            repla = A.repla;
            consistency = A.consistency;
            waveform = A.waveform;
            imagecache = A.imagecache;
            tags = A.tags;
            list = A.list;
% 
%             A.song = [];
%             A.deg = [];
%             A.frag = [];
%             A.repla = [];
%             A.consistency = [];
%             A.waveform = [];
%             A.imagecache = [];

            save(A.info.formated_name,'info','song','deg','frag','repla','consistency','waveform','imagecache','tags','list','-v7.3');
           % save(A.info.formated_name,'A','info');
        end

        function  genAnalysis(sourcedir,birdname,ZPid,channelname,unit,whether_to_save) % generate Neuron

            hit_subdir = sourcedir;

            pl2files = Extract.filename(hit_subdir,'*.pl2');
            mergeids = find(~cellfun(@isempty, regexp(cellstr(pl2files),'merge')));
            zpids = find(~cellfun(@isempty, regexp(cellstr(pl2files),'F\d+')));
            if ~isempty(mergeids)%判断 pl2文件是否完备
                valid_pl2files = pl2files(mergeids); % when there is a merged file
            else
                valid_pl2files = pl2files(setdiff(zpids,mergeids)); % 剔除含有merge关键词的 pl2文件
            end

            stimuli_collect_dir = append(hit_subdir,'\Stimuli');
            stimuli_dirs = Extract.folder(stimuli_collect_dir ).';

            count = 0;
           
            getfileid = @(x) str2num(regexp(convertCharsToStrings(x),'(?<=F)\d*','match'));
            for k = 1: length(valid_pl2files) % parfor

                if ~isempty(... % 一旦目标文件是merge
                        find(~cellfun(@isempty, regexp(cellstr(valid_pl2files{k}),'merge'))))
                    b = Chorus(valid_pl2files{k},stimuli_dirs);
                else
                    
                    fid = getfileid(valid_pl2files{k});
                    followF = cellfun(@(x) getfileid(x),stimuli_dirs);
                    duiying_stidir_id = find(followF == fid);%find(followF == str2num(fid{1}));
                    if ~isempty(duiying_stidir_id)
                    b = Chorus(valid_pl2files{k},stimuli_dirs{duiying_stidir_id});
                    else
                        continue
                    end
                end

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
                Exp = b.getExperiments{1};
                % allocate sap-based feature information to each neuron
                %                 try

                Exp.setFeatures_SAT;
                %                 feature_dir = append(hit_subdir,'\Feature');
                %                 feature_files = Extract.filename(feature_dir,'*.txt');
                %                 datafile_id = find(~cellfun(@isempty, regexp(cellstr(feature_files),'Data')));
                %                 infofile_id = find(~cellfun(@isempty, regexp(cellstr(feature_files),'Info')));
                %                     sorted_data = Experiment.extractFeaturesFromSapRawData(feature_files{datafile_id} , feature_files{infofile_id} );
                %                     Exp.setFeatures_SAP(sorted_data);
                Exp.calMeanFeatures;
                %                 catch
                %                 end
                count = count + 1;
                experiments{count} = Exp;

            end
            %eleinf = load("F:\S01_GeneratedNeurons_20221209\all_eleinf.mat").eleinf;
            %A = Neuron(experiments,eleinf);
            A = Neuron(experiments);
            % A.calHarmRatio;
            info = A.info;


            if ~exist('whether_to_save','var')||whether_to_save==1

                save(A.info.formated_name,'A','info','-v7.3'); % 保存的路径上不能存在中文
            end

        end


        function error = batchGenAnalysis(input_roster, from_where)% 生成Analysis object

            dbstop if error
            tic
            if ~exist('from_where','var')
                from_where = 1;
            end
            not_exist = {};
            count = 0;


            for k = 1: length(input_roster) % par
                if ~isfile(sprintf('%s.mat',input_roster(k).formated_name)) %检查当前的文件夹
                    count = count + 1;
                    not_exist{count} = input_roster(k);
                    %  Archon.genAnalysis(diskname,birdname,ZPid,channelname,unit);
                end
            end

            rest_roster = vertcat(not_exist{:});


            wb = PoolWaitbar(length(rest_roster),'开始生成 Neuron objects');

            if isempty(rest_roster)
                rest_roster = [];
               % error = [];
                return
            end

            error = struct;

           parfor k = from_where:length(rest_roster) % parfor

                tic;
                try

                    thisroster = rest_roster(k); % 否则(如果stimuli不全的话）只包括Cons（前提是F1必须是Cons）
                    Archon.genAnalysis_BasedOnRoster(thisroster);
                    %                thisroster.stimulidirs = thisroster.stimulidirs(find(~cellfun(@isempty, regexp(thisroster.stimulidirs,'F1'))));
                    %                thisroster.pl2files = thisroster.pl2files(find(~cellfun(@isempty, regexp(thisroster.pl2files,'F1'))));
                    %                A = Archon.genAnalysis_BasedOnRoster(thisroster);

                    increment(wb);

                    disp(k);
% 
                catch ME
                    disp('Error here!!!');
                    error(k).num = k;
                    error(k).ME = ME;
                    error(k).roster = rest_roster(k);
                end

                toc;

            end



        end


       
        function rest_roster = Depracated_old_batch_genAnalysis(input_roster, from_where, whichDisk)% 生成Analysis object

            dbstop if error
            tic
            if ~exist('from_where','var')
                from_where = 1;
            end
            not_exist = {};
            count = 0;


            for k = from_where: length(input_roster) % par
                birdname = input_roster(k).birdname;
                ZPid = input_roster(k).neuronid;
                channelname = input_roster(k).channelname;
                unit = input_roster(k).unit;
                diskname = regexp(input_roster(k).neurondir,'[A-Z]:','match');
                diskname = diskname{1};

                if ~isfile(sprintf('%s_%s_%s_%u.mat',birdname,ZPid,channelname,unit)) %检查当前的文件夹
                    count = count + 1;
                    not_exist{count} = input_roster(k);
                    %  Archon.genAnalysis(diskname,birdname,ZPid,channelname,unit);
                end
            end

            rest_roster = vertcat(not_exist{:});

            if exist('whichDisk','var')

                diskname = whichDisk;
            else

                diskname = regexp(input_roster(1).neurondir,'[A-Z]:','match');
                diskname = diskname{1};

            end

            wb = PoolWaitbar(length(rest_roster),'开始生成 Neuron objects');

            if length(rest_roster) == 0
                rest_roster = [];
                return
            end

            birddirs = Extract.folder(diskname);
            
           for k = from_where: length(rest_roster) % par
              % try
                   birdname = rest_roster(k).birdname;
                   ZPid = rest_roster(k).neuronid;
                   channelname = rest_roster(k).channelname;
                   unit = rest_roster(k).unit;

                   hitted_birddirs = birddirs(find(~cellfun(@isempty, regexp(birddirs,birdname))));


                   if ~isfile(sprintf('%s_%s_%s_%u.mat',birdname,ZPid,channelname,unit))
                       if neuroster(k).allexist == 1 % 如果所有stimuli都存在
                           A = Archon.genAnalysis_BasedOnRoster(neuroster(k));
                       else
                           thisroster = neuroster(k); % 前提是F1必须是Cons
                           thisroster.stimulidirs = thisroster.stimulidirs(find(~cellfun(@isempty, regexp(thisroster.stimulidirs,'F1'))));
                           thisroster.pl2files = thisroster.pl2files(find(~cellfun(@isempty, regexp(thisroster.pl2files,'F1'))));
                           A = Archon.genAnalysis_BasedOnRoster(thisroster);
                       end

                       %  Archon.genAnalysis(rest_roster(k).neurondir,birdname,ZPid,channelname,unit);
                   end
%                catch ME
%                    %                   pause
%                    rest_roster(k).ME = ME;
%                    rest_roster(k).identifier = ME.identifier;
%                    disp(ME);
%                end

                increment(wb);
            end

            toc

        end

        function createSinChFolder(zpdir)
%             arch = Archon('D:/');
%             hitbname = find(~cellfun(@isempty,regexp(cellstr(arch.bird_folders),birdname)));
%             hitzpid = find(~cellfun(@isempty,regexp(cellstr(arch.bird_folders),ZPid)));
%             hitbird_folder = arch.bird_folders{intersect(hitbname,hitzpid)};


            % 处理hit folder 内部的 subfolders
     

            destineydir = Convert.mergeSubfolders(fullfile(zpdir,'Stimuli'),'*.wav'); % 合并所有声文件

            Convert.two2one(destineydir); % 生成单通道声文件
            pause(0.1);
            status = rmdir(destineydir,'s');


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

        function channel_unit_info = findSortedNeurons(path_txt) % from  file find sorted units (unit > 0)

            % path_txt = "D:\Piece-O709-W\P07\O709_P07F5.txt"
            channel_unit_info = struct;
            sorted_tables = Spike.split(path_txt);

            for k = 1: length(sorted_tables)

                channel_unit_info(k).channelname = sorted_tables{k}.('channelname'){1};
                channel_unit_info(k).unit =sorted_tables{k}.('unit')(1);
                channel_unit_info(k).channel_unit = sprintf('%s_%u',channel_unit_info(k).channelname,channel_unit_info(k).unit);
            end

        end

        function shared = findConsistentNeurons(dir_txt)  % 找到在每一个文件里都有的神经元

            %dir_txt = "D:\Piece-O709-W\P07"
            txtfile = Extract.filename(dir_txt,'*.txt');

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


        function [birddir,zpdir] = getDirpath(birdname,ZPid)
            arch = Archon('D:/');
            % birddir
            birddirs = cellstr(arch.bird_folders);
            birddir = birddirs{~cellfun(@isempty,regexp(birddirs,birdname))};

            % zplevel
            zpdirs = cellstr(Extract.folder(birddir).');
            zpdir = zpdirs{~cellfun(@isempty,regexp(zpdirs,ZPid))};

        end

        % 处理已经生成的analysis

        function moveCorrespondingFiles()%basis, sourcedir, destineydir) %根据basis里的文件，移动有相同神经元编码的文件到指定文件夹
            dbstop if error
            sourcedir = "D:\CDRF_AnalysisObject";
            basis = Extract.filename("C:\Users\Zhehao\Desktop\Selected_ResptoIndividual",'*.png');
            destineydir = "C:\Users\Zhehao\Desktop\Destiney";
            all_sourcefiles =cellstr( Extract.filename(sourcedir,'*.*'));

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


            wb = waitbar(0,'Creating Experiment analysis objects');
            for k = 1: length(num_ids)
                waitbar(k/length(num_ids),wb,sprintf('%u of %u Experiment',k,length(num_ids)));
                ids_in_T = find(num_ids(k) == IDs);

                Ns = {};
                for i = 1:length(ids_in_T)
                    b = Chorus(T(ids_in_T(i)).MergedTxtPath,T(ids_in_T(i)).Mergedpl2Path,T(ids_in_T(i)).StimuliPath);

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
                    sorted_data = Experiment.extractFeaturesFromSapRawData(T(ids_in_T(i)).FeatureData, T(ids_in_T(i)).FeatureInfo);
                    N.setEachStimuliSapFeatures(sorted_data);
                    N.calMeanFeatures;
                    Ns{i} = N;
                end

                A = Neuron(Ns);
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


            wb = waitbar(0,'Creating Experiment analysis objects');
            for k = 1: length(num_ids)
                waitbar(k/length(num_ids),wb,sprintf('%u of %u Experiment',k,length(num_ids)));
                ids_in_T = find(num_ids(k) == IDs);

                Ns = {};
                for i = 1:length(ids_in_T)
                    b = Chorus(T(ids_in_T(i)).TxtPath,T(ids_in_T(i)).pl2Path,T(ids_in_T(i)).StimuliPath);

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
                        sorted_data = Experiment.extractFeaturesFromSapRawData(T(ids_in_T(i)).FeatureData, T(ids_in_T(i)).FeatureInfo);
                        N.setEachStimuliSapFeatures(sorted_data);
                        N.calMeanFeatures;
                    end

                    Ns{i} = N;
                end

                A = Neuron(Ns);
                A.uniqueid = num_ids(k);


                save(sprintf('%s_%u',A.birdid,A.uniqueid),'A','-v7.3');


            end

            close(wb);
        end

        function all_data_info = Deprecated_AnalyzeFormatedEphysData(dirpath)
            % Automatically generate pl2-txt-folder table based on
            % formated Piece folders
            % 202204211 11:44 pm paused here!

            dbstop if error
            bird_folders = [Extract.folder(dirpath)].'; %#ok<NBRAK>

            pure_folname = {}; % 初始化
            for k = 1: length(bird_folders)
                temp1 = split(bird_folders{k},'\');
                pure_folname{k} = cellstr(temp1{end});
            end

            contain_ephys_ids = find(~cellfun(@isempty, regexp([pure_folname{:}].','Piece')));
            if ~isempty(contain_ephys_ids)
                bird_folders =  bird_folders(contain_ephys_ids);
            else
                bird_folders =  dirpath;
            end

            infocollect = {};
            for k = 1:length(bird_folders)
                pl2files = Extract.filename(bird_folders{k},'*.pl2');
                txtfiles = Extract.filename(bird_folders{k},'*.txt');

                if exist(sprintf('%s\\Stimuli',bird_folders{k}), 'dir')
                    stidirs = Extract.folder(sprintf('%s\\Stimuli',bird_folders{k})).';
                else % 如果未把所有stimuli放进stimuli folder里面
                    stidirs = Extract.folder(sprintf('%s',bird_folders{k})).'
                end

                pl2name = {}; % Extract .pl2 filenames
                for m = 1: length(pl2files)
                    [~,pl2name{m},~] = fileparts(pl2files{m});
                end

                txtname = {};% Extract .txt filenames
                for m = 1: length(txtfiles)
                    [~,txtname{m},~] = fileparts(txtfiles{m});
                end

                stiname = {}; % Extract stimuli folder names
                for m = 1: length( stidirs)
                    temps = split( stidirs{m},'\');
                    stiname{m} = temps{end};
                end

                % Based on .pl2 filenames, find the corresponding .txt and
                % stimuli folder names
                datainfo = struct;
                for m = 1: length(pl2name)

                    txt_hit = find(~cellfun(@isempty, regexp([cellstr(txtname)].',pl2name{m})));

                    sti_hit = find(~cellfun(@isempty, regexp([cellstr(stiname)].',pl2name{m})));

                    % fill the struct
                    sadtemp = regexp(pl2name{m}, '[_-]', 'split');
                    datainfo(m).animal = sadtemp{end};
                    datainfo(m).pl2name = pl2name{m};
                    datainfo(m).keyword = regexp(pl2name{m},'Con|Deg|Detail|Frag|Repla|Mir|Rev','match');
                    datainfo(m).sectionname = regexp(pl2name{m},'Z\d{2}|P\d{2}','match');
                    just_a_temp = regexp(cellstr(datainfo(m).sectionname),'\d{2}','match');
                    datainfo(m).sectionid = str2num(just_a_temp{1}{1});
                    datainfo(m).subfilename = regexp(pl2name{m},'F\d{1}','match');
                    if ~isempty(datainfo(m).subfilename) %如果文件名里存在F标识subfile的id
                        fuck2 = regexp(datainfo(m).subfilename,'\d','match');
                        datainfo(m).subfileid = str2num(just_a_temp{1}{1});
                    elseif ~isempty(regexp(datainfo(m).pl2name,'Mrg|mrg'))% 如果不存在F标识却存在‘Mrg|mrg’标识
                        datainfo(m).subfilename = "Mrg";
                    else % 如果不存在任何标识
                        datainfo(m).subfilename = [];
                    end

                    datainfo(m).pl2path = pl2files{m};
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

        function As = Deprecated_writeSingleAnalysiFile(txt,pl2, folder,feature_data,feature_info)

            b = Chorus(txt,pl2,folder);
            b.select;
            Ns = b.getn;

            As = {};
            for k = 1: length(Ns)
                N = Ns{k};
                sorted_data = Experiment.extractFeaturesFromSapRawData(feature_data,feature_info);
                N.setEachStimuliSapFeatures(sorted_data);
                N.calMeanFeatures;
                A = Neuron(N);
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
                path_pl2 = all_data_info(m).pl2path;
                path_folder = all_data_info(m).stidirpath;

                if isempty(path_folder )&& ~isempty(regexp(path_pl2,'Mrg|mrg'))

                    sameanimal = find(all_data_info(m).animal == [all_data_info.animal].');
                    sameZfile = find(ismember([all_data_info.sectionname].',all_data_info(m).sectionname));
                    %sameZfile = find(ismember(all_data_info(m).sectionname, [all_data_info.sectionname].'));
                    sameneuronids = intersect(sameanimal,sameZfile);
                    path_folder = [all_data_info(sameneuronids).stidirpath].';

                end

                b = Chorus(path_txt,path_pl2,path_folder);
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
                        final_info(count).path_pl2 = this_section_struct(idx).pl2path;
                        final_info(count).path_folder = this_section_struct(idx).stidirpath;
                        final_info(count).path_pairs = [];
                    else

                        for k = 1: length(this_section_struct)

                            final_info(count).path_txt = [];
                            final_info(count).path_pl2 = [];
                            final_info(count).path_folder =[];
                            final_info(count).path_pairs{k}{1} = this_section_struct(k).txtpath;
                            final_info(count).path_pairs{k}{2} = this_section_struct(k).pl2path;
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

