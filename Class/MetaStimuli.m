classdef MetaStimuli < handle
    % Metadata of stimuli
    % 读取stimuli的分割文件，生成eleinf
    %This class read and process the segmentation information of each stimuli to generat eth metadata of a stimuli set

    properties
        inout_dir % output dir of moveFromBucket, input dir of other functions
        bucket_dir
    end

    methods

        function bf = MetaStimuli(inout_dir)
            % bucket_dir is the path of bucket folder bird-song,collection
            % of all song folders
            if exist('inout_dir','var')
                bf.inout_dir = inout_dir;
            end
        end

    end

    methods(Static)

        function splited = split(fraginf)
            % 把集合在一起的fraginf分割成对于每个eleinf的单个fraginf
            % input is the s.sylinf, output is the sylinf for each each song
            dbstop if error
            songs = unique(cellstr({fraginf(:).songname}.'));

            parfor kk = 1: length(songs)  % I change for to par-for 1103
                idxs = find( strcmp(cellstr({fraginf(:).songname}.') ,songs{kk}) );
                if length(idxs) == 1
                    thisT = struct2table(fraginf(idxs),'AsArray',1);
                else
                    thisT =  struct2table(fraginf(idxs));
                end

                thisT = sortrows(thisT,'fragid');  % 这里可能会改   order the struct based on the syllable number
                splited{kk} = table2struct(thisT);
            end

        end

        function  moveNoiseFiles(sourcedir, destineydir)
            % 去除源文件夹里被认定为噪声的声音文件
            % to move pure noise files out of the sourcedir

            tic;

            rawfiles = Extract.filename(sourcedir,'*.wav');

            D = parallel.pool.DataQueue; % configurate waitbar in par-for
            h = waitbar(0, 'Start processing');
            num_files = length(rawfiles);
            nUpdateWaitbar(num_files, h); % Dummy call to nUpdateWaitbar to initialise
            afterEach(D, @nUpdateWaitbar); % Go back to simply calling nUpdateWaitbar with the data

            Centro_Thres = 125; % Threshold

            parfor k = 1: length(rawfiles)

                [y,fs] = audioread(rawfiles{k});% here the y is the raw y
                centroid = spectralCentroid(y,fs);
                if mean(centroid) < Centro_Thres
                    [~,fname,ext] =  fileparts(rawfiles{k});
                    new_fpath = fullfile(destineydir,strcat(fname,ext));
                    movefile(rawfiles{k},new_fpath);
                end
                send(D, 1);   % send only an "increment" for the waitbar.

            end

            toc;


            function p = nUpdateWaitbar(data, h)    % subfunctiuon for waitbar
                persistent TOTAL COUNT H
                if nargin == 2
                    H = h;TOTAL = data;COUNT = 0;   % initialisation mode
                else
                    COUNT = 1 + COUNT; p = COUNT / TOTAL;    % afterEach call, increment COUNT
                    waitbar(p, H,sprintf('Checking %u of totally %u files',COUNT,TOTAL));
                end
            end


        end

        function metadata = read2ele(file_or_dir)  % read to elements
            % 功能：读取鸟歌的分割数据到element

            if exist(file_or_dir,'dir')
                matfiles = Extract.filesAllLevel(file_or_dir,'*.mat');
            else
                matfiles = file_or_dir;
            end

            song_collect = {};

            for m = 1: length(matfiles) % parfor here

                loaded = load(matfiles{m});

                % Very dangerous bad code, temporarilary used
                %                 if isfield(loaded.segdata,'eleedge')
                %                     loaded.segdata.eleedge = loaded.segdata.eleedge -0.008;
                %                 end
                %                 if isfield(loaded.segdata,'syledge')
                %                     loaded.segdata.syledge = loaded.segdata.syledge -0.008;
                %                 end
                %                 if isfield(loaded.segdata,'motedge')
                %                     loaded.segdata.motedge = loaded.segdata.motedge -0.008;
                %                 end
                %%%%%%%%%Extremeley bad and dangerous code shown above

                if isfield(loaded.segdata,'eleedge')
                    two_eleedge = repmat(loaded.segdata.eleedge,[2,1]);
                    alledges = sort( vertcat(loaded.segdata.syledge(:),reshape(two_eleedge,[],1) ));
                else
                    alledges = sort(loaded.segdata.syledge(:));
                end

                % the following is a very bad temporily code
                try
                    fs = loaded.segdata.fs;
                catch
                    fs = 32000;
                end

                fiy = bandpass(loaded.segdata.rawy,[900 6000],fs); %% It is very important that here the y should be fiy !!!!! filtered y instead of the raw y
                I = Cal.spec(fiy,fs);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                initials = alledges(1:2:end)-0.008;
                terminals = alledges(2:2:end)+ 0.008;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Extremely extremely dangerous code
                song_eleinf = struct;

                for w = 1: length(initials) % can add a par-for
                    song_eleinf(w).initial = initials(w)*fs;
                    song_eleinf(w).terminal = terminals(w)*fs;
                    song_eleinf(w).songname = loaded.segdata.birdid;

                    try
                        song_eleinf(w).motif = MetaStimuli.findMotif(loaded.segdata.motedge,initials(w),terminals(w));
                    catch
                        song_eleinf(w).motif = 0;
                    end

                    hp_y = highpass(loaded.segdata.rawy,450,fs);
                    if isfield(loaded.segdata,'rawy')
                        song_eleinf(w).y = hp_y(max(initials(w)*fs,1):terminals(w)*fs); % originally loaded.segdata.rawy

                    end
                    if isempty(song_eleinf(w).y)
                        disp('Pause')
                    end
                    song_eleinf(w).fs = fs;

                    imgratio = size(I,2)/(length(hp_y)/fs);
                    song_eleinf(w).fragI = imresize(I(:,max(round(initials(w)*imgratio),1):round(terminals(w)*imgratio)),[257,50]);
                    song_eleinf(w).fragid = w;
                    song_eleinf(w).unifragnames = sprintf('%s-%02u',Convert.bid(song_eleinf(w).songname ),song_eleinf(w).fragid);

                    song_eleinf(w).fs = song_eleinf(w).fs;
                    song_eleinf(w).fullname = sprintf('%s-%02u',song_eleinf(w).songname,song_eleinf(w).fragid);
                end

                song_collect {m} = song_eleinf;

            end

            metadata = horzcat(song_collect{:});

        end

        function pickupMotifsFromBucket(singleFolder)

            % To Extract several song files from a target folder
            % singleFolder is the path of a single folder in bucket server
            dbstop if error

            tic
          
            adultfilenames = Bird.getAdultSongs(singleFolder); % 首先，得到 adult songs
            if isempty(adultfilenames)
                return
            end
            rawfiles = flip(sort(cellstr(adultfilenames)));
            %rawfiles = flip(sort(cellstr(Extract.filename(singleFolder,'*.wav')))); %(randperm(length(rawfiles)))
            num_tosave = 20; % 30 files to copy for each folder
            ampthres = 0.008; %声信号振幅值是否足够
            minthres = 2; %声文件最短时长
            maxthres = 15;%声文件最长时长
            shortsig = 1.3; %有声音的时间段的最短音长 % 之前使用过的threshold是1.5 seconds
            centroid_thres = 3200;
            redundancy = 0.3; % 0.5 seconds
            isi_dur_thres = 0.4; %  0.4 seconds : how long a duration will be regarded as separation of bouts
            candidates = struct([]);
            max_session_num = 10; % 最大的 running session， 为了限制运行时间

            % 为了提高运算速度，采用每次run 100 times 的方式
            num_parallel = 32;
            for session = 1:min(ceil(length(rawfiles)/num_parallel),max_session_num)

                fprintf('Running session %u',session);

                feeded_rawfiles = rawfiles(num_parallel*(session-1)+1: min(num_parallel*session,length(rawfiles)));

                multi_subsets = {};
                for n = 1: length(feeded_rawfiles) % 这一部分采用parfor
                    fprintf('Current id is %u \r',n);
                    try
                        [y,fs] = audioread(feeded_rawfiles{n});% here the y is the raw y
                    catch % if there is an error
                        continue
                    end

                    %<*>Judge by the length of the raw signal
                    if  maxthres < length(y)/fs  || length(y)/fs < minthres
                        fprintf('原始音长不足\r');
                        continue
                    end

                    fiy = highpass(abs(y),500,fs); %  to remove the noise generataed by low-frequency noise
                    ampenv = envelope(fiy,320*9,'rms'); % amplitude envelope
                    %figure; plot(ampenv); figure; Draw.spec(y,fs);% 1ms-32
                    sigpoints = find(ampenv>ampthres); % sigs means significant signal points



                    sigy = fiy(sigpoints);
                    percentile_thres = 0.3; % 30%
                    persigs = length(sigpoints)/length(y); % percentage of signals that are significant


                    if   length(sigpoints)/fs < shortsig || persigs < percentile_thres
                        fprintf('有效音长或音比不足\r');  %<*>Judge by the length of significanty signals
                        continue
                    end

                    centroid = spectralCentroid(sigy,fs); %<*>Judge by centroid to remove noise0like signal
                    mean_centroid = mean( centroid);   %<*>Judge by percentage of signaled time duration to remove calls
                    if  mean_centroid < centroid_thres
                        fprintf('音重心过低\r');
                        continue
                    end

                    %<*>Judge Inter-segment-interval: Too strict criteria

                    interval = diff(sigpoints);
                    index = find(interval> fs*isi_dur_thres); %  inter-syllable intervals 如果足够大的话就被认为是一个motif
                    inter_motif_intervals = sort([sigpoints(index);sigpoints(index+1)],'ascend');
                    edges =  [ min(sigpoints);inter_motif_intervals;max(sigpoints)];

                    subsets = struct([]);
                    internal_count = 0;
                    for kk = 1:length(edges)/2
                        sectioned_signal = y(edges(2*kk-1) :edges(2*kk));

                        if length(sectioned_signal)/fs < 0.7... %如果这一段的song太短了
                                ||mean(spectralCentroid(sectioned_signal,fs))<500
                            %fprintf('meanSpecCentroid is:%f',mean(spectralCentroid(sectioned_signal,fs)))
                            fprintf('sectioned signal not fit\r');
                            continue
                        end

                        fprintf('合格\r')
                        internal_count = internal_count + 1;

                        multi_subsets{n} = subsets;
                    end

                        % To be noticed that y and sectioned_signals are not equal
                        subsets(internal_count).y = y(max(1,edges(2*kk-1)-fs*redundancy):...
                            min(length(y), edges(2*kk)+fs*redundancy));
                        subsets(internal_count).fs = fs;
                        subsets(internal_count).initial = edges(2*kk-1);
                        subsets(internal_count).terminal = edges(2*kk);
                        subsets(internal_count).sourcefile = rawfiles{n};

%                         % 这部分是为了判断motifs，以及对song syllable做分割
%                         subset_fiy = highpass(abs(subsets(internal_count).y),500,fs); % filtered y, to remove the noise generataed by low-frequency noise
%                         subset_ampenv = envelope(subset_fiy,320*3,'rms'); % amplitude envelope
%                         %figure; plot(ampenv); figure; Draw.spec(y,fs);% 1ms-32
%                         subset_sigpoints = find(subset_ampenv>ampthres); % sigs means significant signal points
% 
% 
%                         subset_interval = diff(subset_sigpoints);
%                         subset_index = find( subset_interval> 1); % ISI is the inter-syllable intervals
%                         subset_isis = sort([subset_sigpoints(subset_index);subset_sigpoints(subset_index+1)],'ascend');
%                         subset_edges =  [ min(subset_sigpoints); subset_isis;max(subset_sigpoints)];
%                         %figure; plot(subset_fiy); hold on; xline(subset_edges)
%                         raw_onset = subset_edges(1:2:length(subset_edges));
%                         raw_offset = subset_edges(2:2:length(subset_edges));
%                         raw_off_minus_on = raw_offset - raw_onset;
%                         onset = raw_onset(find(raw_off_minus_on > 512));
%                         offset = raw_offset(find(raw_off_minus_on > 512));
%                         sim_mat = [];
%                         for t = 1:length(onset)
% 
%                             frag1_y = subsets(internal_count).y(onset(t) :offset(t));
%                             for tt = 1:length(onset)
%                                 frag2_y = subsets(internal_count).y(onset(tt) :offset(tt));
%                                 sound1 = SAT_sound(frag1_y,32000);
%                                 sound2 = SAT_sound(frag2_y,32000);
%                                 sim = SAT_similarity(sound1,sound2,0);
%                                 sim.calculate_similarity;
%                                 sim_mat(t,tt) = sim.score.accuracy;
%                             end
% 
%                         end
% 
%                      %   calculate mfcc/stft similarity matrix
%                      sat1 = SAT_sound(subset_fiy,fs);
%                      sat2 = SAT_sound(subset_fiy,fs);
%                      sim1 = SAT_similarity(sat1,sat2,0);
%                      sim1.calculate_similarity;
%                      sim1.score;
%                      figure;
%                      [coeffs,~,~,~] = mfcc(subset_fiy,fs);
%                      figure; imagesc(sim_mat)
% 
%                      sim_mfcc = [];
%                      for a = 1:size(coeffs,1)
%                          vec1 = coeffs(a,:)
%                          for aa = 1:size(coeffs,1)
%                              vec2 = coeffs(aa,:);
%                              sim_mfcc(a,aa) = 1- norm(vec1-vec2);
%                          end
%                      end
% 
%                      figure; imagesc(sim_mfcc);
%                     end


        

                toc
                candidates = horzcat(candidates,horzcat(multi_subsets{:}));

                if length(candidates) >= num_tosave
                    break
                end
                end

            end

            target_dir = '.\'; %bf.inout_dir; % "E:\WavsCollection" % destination
            [~,birdid,~] = fileparts(singleFolder);
            birdid = Convert.bid(birdid,1);
            mkdir(sprintf('%s%s',target_dir,birdid));


            for a = 1: min(num_tosave,length(candidates) )
                candidates(a).name = sprintf('%s%s\\%s-%u.wav',target_dir,birdid,birdid,a);
                audiowrite(candidates(a).name,...
                    candidates(a).y,candidates(a).fs);
            end
            candidates = rmfield(candidates,'y');
            writetable(struct2table(candidates),sprintf('%s\\%s-%s.txt',birdid,birdid,'DataInfo'),'FileType','text');
            time_used = toc;

            Running_info = struct;
            Running_info.files_processed = sprintf("%s:%u files are processed",birdid,session*num_parallel);
            Running_info.time_used = sprintf("The whole running cost %f seconds ", time_used);
            writetable(struct2table(Running_info),sprintf('%s\\%s-%s.txt',birdid,birdid,'RunInfo'),'FileType','text');

        end

        function Deprecated_pickupMotifsFromBucket2(singleFolder)

            % To Extract several song files from a target folder
            % singleFolder is the path of a single folder in bucket server
            dbstop if error

            tic
            rawfiles = flip(sort(cellstr(Extract.filename(singleFolder,'*.wav')))); %(randperm(length(rawfiles)))
            num_tosave = 10; % 30 files to copy for each folder
            ampthres = 0.008;
            minthres = 2;
            maxthres = 15;
            shortsig = 1.3; % longer than 1.5 seconds
            centroid_thres = 3200;
            redundancy = 0.3; % 0.5 seconds
            isi_dur_thres = 0.4; %  0.4 seconds : how long a duration will be regarded as separation of bouts
            %datetime('2022-06-07_10-31-03','InputFormat','yyyy-MM-dd_hh-mm-ss')
            candidates = struct([]);
            max_session_num = 10; % 最大的 running session， 为了限制运行时间

            % 为了提高运算速度，采用每次run 100 times 的方式
            num_parallel = 32;
            for session = 1:min(ceil(length(rawfiles)/num_parallel),max_session_num)

                fprintf('Running session %u',session);

                feeded_rawfiles = rawfiles(num_parallel*(session-1)+1: min(num_parallel*session,length(rawfiles)));

                multi_subsets = {};
                for n = 1: length(feeded_rawfiles) % 这一部分采用parfor
                    fprintf('Current id is %u\r',n);
                    try
                        [y,fs] = audioread(feeded_rawfiles{n});% here the y is the raw y
                    catch % if there is an error
                        continue
                    end

                    %<*>Judge by the length of the raw signal
                    if  maxthres < length(y)/fs  || length(y)/fs < minthres
                        %                         fprintf('原始音长不足\r');
                        continue
                    end

                    fiy = highpass(abs(y),500,fs); %  to remove the noise generataed by low-frequency noise
                    ampenv = envelope(fiy,320*9,'rms'); % amplitude envelope
                    %figure; plot(ampenv); figure; Draw.spec(y,fs);% 1ms-32
                    sigpoints = find(ampenv>ampthres); % sigs means significant signal points
                    sigy = fiy(sigpoints);
                    percentile_thres = 0.3; % 30%
                    persigs = length(sigpoints)/length(y); % percentage of signals that are significant


                    if   length(sigpoints)/fs < shortsig || persigs < percentile_thres
                        %                         fprintf('有效音长或音比不足\r');  %<*>Judge by the length of significanty signals
                        continue
                    end

                    centroid = spectralCentroid(sigy,fs); %<*>Judge by centroid to remove noise0like signal
                    mean_centroid = mean( centroid);   %<*>Judge by percentage of signaled time duration to remove calls
                    if  mean_centroid < centroid_thres
                        %                         fprintf('音重心过低\r');
                        continue
                    end

                    %<*>Judge Inter-segment-interval: Too strict criteria

                    interval = diff(sigpoints);
                    index = find(interval> fs*isi_dur_thres); % ISI is the inter-syllable intervals
                    isis = sort([sigpoints(index);sigpoints(index+1)],'ascend');
                    edges =  [ min(sigpoints); isis;max(sigpoints)];

                    subsets = struct;
               
                    parfor kk = 1:length(edges)/2
                        localedges = edges;
                        localy = y;


                        subsets(kk).sectioned_signal = localy(localedges(2*kk-1) :localedges(2*kk));
                        subsets(kk).duration = length(localy(localedges(2*kk-1) :localedges(2*kk)))/fs
                        %                         subsets(kk).specent = mean(spectralCentroid(localy(localedges(2*kk-1) :localedges(2*kk)),fs));
                        %
                        %                         % To be noticed that y and sectioned_signals are not equal
                        subsets(kk).y = localy(max(1,localedges(2*kk-1)-fs*redundancy):...
                            min(length(localy), localedges(2*kk)+fs*redundancy));
                        subsets(kk).fs = fs;
                        subsets(kk).initial = localedges(2*kk-1);
                        subsets(kk).terminal = localedges(2*kk);
                        subsets(kk).sourcefile = rawfiles{n};
                        try
                            subsets(kk).specent = mean(spectralCentroid(localy(localedges(2*kk-1) :localedges(2*kk)),fs)  );
                        catch
                            subsets(kk).specent =  0;
                        end

                    end

                  

                    screened_subsets = subsets(intersect(find([subsets.duration].' >=0.7), find([subsets.specent].' >=centroid_thres) ));
                    multi_subsets{n} = screened_subsets;
                end

                toc
                candidates = horzcat(candidates,horzcat(multi_subsets{:}));

                if length(candidates) >= num_tosave
                    break
                end

            end

            target_dir = '.\'; %bf.inout_dir; % "E:\WavsCollection" % destination
            [~,birdid,~] = fileparts(singleFolder);
            mkdir(sprintf('%s\\%s',target_dir,Convert.bid(birdid,1)));


            for a = 1: min(num_tosave,length(candidates) )
                candidates(a).name = sprintf('%s\\%s\\%s-%u.wav',target_dir,birdid,birdid,a);
                audiowrite(candidates(a).name,...
                    candidates(a).y,candidates(a).fs);
            end
            candidates = rmfield(candidates,'y');
            writetable(struct2table(candidates),sprintf('%s\\%s-%s.txt',birdid,birdid,'DataInfo'),'FileType','text');
            time_used = toc;

            Running_info = struct;
            Running_info.files_processed = sprintf("%s:%u files are processed",birdid,session*num_parallel);
            Running_info.time_used = sprintf("The whole running cost %f seconds ", time_used);
            writetable(struct2table(Running_info),sprintf('%s\\%s-%s.txt',birdid,birdid,'RunInfo'),'FileType','text');
            toc

        end

        function all_eleinf = getAllTwoMotifEleinf(two_motif_ele_dir)
            if ~exist('two_motif_ele_dir','var')
                two_motif_ele_dir = "G:\StimuliSource";
            end

            subdirs = Extract.folder(two_motif_ele_dir).';
            twomotif_ids = find(~cellfun(@isempty, regexp([subdirs{:}].','twoMotif')));
            twoMotif_subdirs = subdirs(twomotif_ids).';

            all_eleinf = MetaStimuli.Eleinf(twoMotif_subdirs,1,'SegData');

        end

        function which_motif = findMotif(motedges, frag_initial,frag_terminal)

            for k = 1:length(motedges)-1
                qian = motedges(k);
                hou = motedges(k+1);
                if frag_initial>qian && frag_initial<hou && frag_terminal>qian && frag_terminal <hou
                    % only when these four statement was correct, then:
                    which_motif = k;
                    break
                end
            end

        end

        function meta2fig(meta_ele)
            % 功能说明：generate figures for streamlit app classification from input_eleinf
            targetdir = 'Figs';
            mkdir(targetdir);
            for k = 1:length(meta_ele)
                figure('Position',[681 403 length(meta_ele(k).y)*187/2305 543]);
                Draw.spec(meta_ele(k).y,meta_ele(k).fs)
                saveas(gcf,fullfile(targetdir,sprintf('%s.png',meta_ele(k).fullname)))
                close(gcf)
            end

        end

        function metadata_withcatego = replenishCatego(csv_path,metadata)
            dbstop if error
            csv_info = table2struct(readtable(csv_path));
            for k = 1:length(metadata)
                if ~isempty( csv_info(find(~cellfun(@isempty, regexp({csv_info.name}.',metadata(k).fullname)))))
                    metadata(k).catego = csv_info(find(~cellfun(@isempty, regexp({csv_info.name}.',metadata(k).fullname)))).catego;
                else
                    metadata(k).catego = 0;
                end
            end

            metadata_withcatego = metadata;

        end

        function catego = categorizeByBingyang(y,fs,bingyang)

            featurespace = [];
            for k = 1:length(bingyang)
                bingyang(k).sat.meanfeatures.dur = length(bingyang(k).y); % 加上syllable duration试一试
                featurespace(k,:) = cell2mat(struct2cell(bingyang(k).sat.meanfeatures));
            end

            test = SAT_sound(y,fs);
            test.meanfeatures.dur = length(y); % 加上syllable duration
            featurespace(k+1,:) = cell2mat(struct2cell(test.meanfeatures));

            %由于每个feature绝对大小的区间不同，每个feature的比重不一样，所以要zscore一下
            z_featurespace = zscore(featurespace,0,1); % zscore比 rescale为（0，1）更靠谱
            test_coor = z_featurespace(k+1,:);

            for a = 1:length(bingyang)
                bingyang(a).coor = z_featurespace(a,:);
                bingyang(a).disfeature = norm(bingyang(a).coor - test_coor);
            end

            %计算test（被测）与每个catego的兵样的平均距离
            unique_categos = unique([bingyang.catego].');
            distances = [];
            for b = 1:length(unique_categos)
                local_BY = bingyang([bingyang.catego].'== unique_categos(b)); % 同一个catego的兵样
                distances(b) = mean([local_BY.disfeature].');
            end

            [~,catego] = min(distances); % 与哪个catego平均距离最小，就归为哪个catego

        end

    end

    methods(Static) % Deprecated but maybe still useful so don't delete!

        function all_eleinf = Deprecated_Sylinf(subdirs,id_thres,matFolder)  % 遍历 all the subdirs

            % assemble into eleinf

            % This function can be applied to both syldata or eledata
            all_collect = {};
            to_delete = [];

            for w = 1:length(subdirs)
                birdid  = split(subdirs{w},'\');
                birdid = regexp(birdid{end},'\d+','match');

                if ~isempty(birdid) % 当subidr是以birdid来命名的时候
                    birdid = str2double(birdid{1});

                    if exist('id_thres','var')
                        if birdid < id_thres
                            to_delete = [to_delete,w];
                        end
                    end
                else
                    disp('现在的folder不是以birdid命名的');
                end
            end

            subdirs(to_delete) = []; % delete those does not match the id threshold

            for r = 1:length(subdirs)  % par-for can be used here

                matfiles = Extract.filename(sprintf('%s\\%s',subdirs{r},matFolder),'*.mat'); % matFolder is the folder containing segdata.mat,
                % which can be SegData or SylData or EleData

                song_collect = {};

                for m = 1: length(matfiles)

                    loaded = load(matfiles{m});


                    alledges = sort(loaded.segdata.syledge(:));



                    % the following is a very bad temporily code
                    if ~isfield(loaded.segdata,'fs') % if there is no such field named fs
                        fs = 32000;
                    else
                        fs = loaded.segdata.fs;
                    end

                    if ~isfield(loaded.segdata,'rawy')
                        %song_eleinf = struct;
                        continue
                    end

                    fiy = bandpass(loaded.segdata.rawy,[900 6000],fs); %% It is very important that here the y should be fiy !!!!! filtered y instead of the raw y
                    I = Cal.spec(fiy,fs);

                    initials = alledges(1:2:end);
                    terminals = alledges(2:2:end);

                    song_eleinf = struct;


                    for w = 1: length(initials) % can add a par-for
                        song_eleinf(w).initial = initials(w)*fs/1000;
                        song_eleinf(w).terminal = terminals(w)*fs/1000;
                        song_eleinf(w).songname = loaded.segdata.birdid;
                        if isfield(loaded.segdata,'rawy')
                            if initials(w)== 0
                                % when the stimuli starts from the very beginning
                                song_eleinf(w).y = loaded.segdata.rawy(1:terminals(w)*fs/1000);
                            else
                                song_eleinf(w).y = loaded.segdata.rawy(initials(w)*fs/1000:terminals(w)*fs/1000);
                            end
                        end
                        song_eleinf(w).fs = fs;

                        if initials(w)== 0
                            % when the stimuli starts from the very beginning
                            song_eleinf(w).fragI = imresize(I(:,1:terminals(w)),[257,50]);
                        else
                            song_eleinf(w).fragI = imresize(I(:,initials(w):terminals(w)),[257,50]);
                        end
                    end

                    for f = 1:length(song_eleinf)
                        song_eleinf(f).fragid = f;
                    end

                    song_collect {m} = song_eleinf;

                end

                folder_eleinf = horzcat(song_collect{:});

                all_collect{r} = folder_eleinf;


            end

            all_eleinf = horzcat(all_collect{:});

            % parfor a = 1:length(all_eleinf)
            %     all_eleinf(a).uniqueid = a;
            % end


        end

        function all_eleinf = Deprecated_Fraginf(subdirs,id_thres,matFolder)  % 遍历 all the subdirs

            % If there is eleinf , then assemble into element level, while
            % if there is no eleinf but only sylinf, assemble into sylinf

            % This function can be applied to both syldata or eledata
            all_collect = {};

            to_delete = [];
            for w = 1:length(subdirs)
                birdid  = split(subdirs{w},'\');
                birdid = regexp(birdid{end},'\d+','match');
                birdid = str2double(birdid{1});
                if exist('id_thres','var')
                    if birdid < id_thres
                        to_delete = [to_delete,w];
                    end
                end
            end

            subdirs(to_delete) = []; % delete those does not match the id threshold

            for r = 1:length(subdirs)  % par-for can be used here

                matfiles = Extract.filename(sprintf('%s\\%s',subdirs{r},matFolder),'*.mat'); % matFolder is the folder containing segdata.mat,
                % which can be SegData or SylData or EleData

                song_collect = {};

                for m = 1: length(matfiles)

                    loaded = load(matfiles{m});

                    if exist('loaded.segdata.eleedge','var')
                        two_eleedge = repmat(loaded.segdata.eleedge,[2,1]);
                        alledges = sort( vertcat(loaded.segdata.syledge(:),reshape(two_eleedge,[],1) ));
                    else
                        alledges = sort( vertcat(loaded.segdata.syledge(:)));
                    end

                    % the following is a very bad temporily code
                    if ~isfield(loaded.segdata,'fs') % if there is no such field named fs
                        fs = 32000;
                    else
                        fs = loaded.segdata.fs;
                    end

                    fiy = bandpass(loaded.segdata.rawy,[900 6000],fs); %% It is very important that here the y should be fiy !!!!! filtered y instead of the raw y
                    I = Cal.spec(fiy,fs);

                    initials = alledges(1:2:end);
                    terminals = alledges(2:2:end);

                    song_eleinf = struct;


                    for w = 1: length(initials) % can add a par-for
                        song_eleinf(w).initial = initials(w)*fs/1000;
                        song_eleinf(w).terminal = terminals(w)*fs/1000;
                        song_eleinf(w).songname = loaded.segdata.birdid;
                        if isfield(loaded.segdata,'rawy')
                            song_eleinf(w).y = loaded.segdata.rawy(initials(w)*fs/1000:terminals(w)*fs/1000);
                        end
                        song_eleinf(w).fs = fs;
                        song_eleinf(w).fragI = imresize(I(:,initials(w):terminals(w)),[257,50]);
                    end

                    for f = 1:length(song_eleinf)
                        song_eleinf(f).fragid = f;
                    end

                    song_collect {m} = song_eleinf;

                end

                folder_eleinf = horzcat(song_collect{:});


                all_collect{r} = folder_eleinf;



            end

            all_eleinf = horzcat(all_collect{:});

            % parfor a = 1:length(all_eleinf)
            %     all_eleinf(a).uniqueid = a;
            % end

        end

        function HalfDeprecated_pickupRawRecordingFilesFromBucket(singleFolder)

            % To Extract several song files from a target folder
            % singleFolder is the path of a single folder in bucket server
            dbstop if error

            rawfiles = flip(sort(cellstr(Extract.filename(singleFolder,'*.wav')))); %(randperm(length(rawfiles)))
            %             rawfiles = rawfiles(randperm(length(rawfiles))); % shuffling the order of the sound files
            num_to_copy = 30; % 30 files to copy for each folder
            ids = [];
            iteration = 0; % how many turns of iteration
            LastOrNot = 0;
            keepornot = []; % judge whether a file is to be ignored ( not a good song file)

            fs = 32000; % hard-coded sampling frequency
            ampthres = 0.008;
            minthres = 2;
            maxthres = 15;
            shortsig = 1.5; % longer than 1.6 seconds
            longsig = 4.5;% shorter than 3 seconds
            centroid_thres = 3200;
            isi_dur_thres = 0.4; %  0.4 seconds : how long a duration will be regarded as separation of bouts
            num_dur = 6;  % number of bout



            %datetime('2022-06-07_10-31-03','InputFormat','yyyy-MM-dd_hh-mm-ss')
            while length(ids) < num_to_copy

                iteration = iteration + 1;
                disp('iteration update')

                if length(rawfiles) > 100* iteration
                    a = 1 + (iteration-1)* 100;
                    b = iteration*100;
                    %files = rawfiles(1 + (iteration-1)* 100 :iteration*100); % restrict the number of files
                else
                    a = 1;
                    b = length(rawfiles);
                    %files = rawfiles;
                    disp('Last 100!');
                    LastOrNot = 1;
                end

                spec = {};
                centroid = {};

                for n = a: b % here should be parfor

                    %                     try
                    [y,fs] = audioread(rawfiles{n});% here the y is the raw y
                    %                     catch
                    %                         continue
                    %                     end

                    %<*>Judge by the length of the raw signal
                    if  maxthres < length(y)/fs  || length(y)/fs < minthres
                        keepornot(n) = 0;
                        continue
                    end

                    fiy = highpass(abs(y),500,fs); %  to remove the noise generataed by low-ftrquency noise
                    ampenv = envelope(fiy,320*9,'rms'); % amplitude envelope
                    %figure; plot(ampenv); figure; Draw.spec(y,fs);% 1ms-32

                    sigs = find(ampenv>ampthres); % sigs means significant signal points
                    sigy = fiy(sigs);

                    %<*>Judge by the length of significanty signals
                    if  longsig < length(sigs)/fs  || length(sigs)/fs < shortsig
                        keepornot(n) = 0;
                        continue
                    end

                    %<*>Judge by centroid to remove noise0like signal

                    centroid{n} = spectralCentroid(sigy,fs);
                    mean_centroid = mean( centroid{n} );
                    if  mean_centroid < centroid_thres
                        keepornot(n) = 0;
                        continue
                    end

                    %<*>Judge by percentage of signaled time duration to remove calls
                    percentile_thres = 0.3; % 30%
                    persigs(n) = length(sigs)/length(y); % percentage of signals that are significant
                    if persigs(n) < percentile_thres
                        keepornot(n) = 0;
                        continue
                    end

                    %<*>Judge Inter-segment-interval: Too strict criteria,
                    % to remove songs which has too many bolts

                    interval = diff(sigs);
                    isi{n} = interval(interval~=1); % ISI is the inter-syllable intervals
                    long_isi = isi{n}(isi{n}> fs*isi_dur_thres);


                    %if length(long_isi) ~= num_dur % this criteria might be too strict
                    fprintf('长间隔的数目为: %u\n',length(long_isi));
                    if length(long_isi) > num_dur

                        keepornot(n) = 0;
                        continue
                    end

                    % calculate spec and rawspec


                    keepornot(n) = 1;

                    %diji = length(find(filejudge));
                    %fprintf('得到了第%u个合格的Song！其在文件夹中的序列是 %u\n',diji,n);
                    fprintf('得到了一个合格的Song！其在文件夹中的序列是 %u\n',n);
                    %disp('wtf')
                end


                ids = find(keepornot);
                fprintf('目前收获了_%u_首Song',length(ids));
                % Extract 5 files from each folder
                if length(ids) > num_to_copy
                    ids = ids(1: num_to_copy );
                else
                    ids = ids;
                end

                if LastOrNot == 1
                    disp('Not enough songs!');
                    break

                end

            end


            picked_files = rawfiles(ids);

            if isempty(picked_files)
                disp('This folder contains no songs that meet the requirements');
                target_dir = "E:\WavsCollection" % destination
                temptemp = split(singleFolder ,'\');
                birdid = temptemp{end};
                output_dir = sprintf('%s\\%s',target_dir,birdid);
                return
            end


            parfor w = 1: length(picked_files)
                [y,fs] = audioread(picked_files{w}); % here the y is the raw y
                fiy = highpass(abs(y),500,fs);
                ampenv = envelope(fiy,320*9,'rms'); % amplitude envelope
                sigs = find(ampenv>ampthres); % sigs means significant signal points
                %spec{w}  = imresize(Cal.spec(fiy(sigs),fs),[257 2000]);
                % rawspec{w}  = imresize(Cal.spec(fiy,fs),[257 2000]) ;
            end

            % montage ???

            target_dir = '.\'; %bf.inout_dir; % "E:\WavsCollection" % destination
            temptemp = split(singleFolder ,'\');
            birdid = temptemp{end};
            output_dir = sprintf('%s\\%s',target_dir,birdid);
            mkdir(sprintf('%s\\%s',target_dir,birdid));

            fileinfo = struct;
            for a = 1: length(picked_files)

                oldname = picked_files{a};
                temp = split(oldname,'\');
                individual = temp{end};
                newname{a} = sprintf('%s-%u.wav',birdid,a);
                destination = sprintf('%s\\%s\\%s-%u.wav',target_dir,birdid,birdid,a);
                copyfile(oldname,destination);
                fileinfo(a).newname = newname{a};
                fileinfo(a).oldname = oldname;
            end

            writetable(struct2table(fileinfo),sprintf('%s\\%s-%s.txt',birdid,birdid,'DataInfo'),'FileType','text');


        end

    end
end

