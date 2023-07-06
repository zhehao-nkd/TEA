function [neuroster,malfunctions] = extractRAwDataInfo(datadir)
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

    waitbar(kk/length(neurondir),wb2,sprintf('Processing %u of all %u folders',kk,length(neurondir)))

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

end
