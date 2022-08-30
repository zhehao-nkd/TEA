classdef getInf < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        inout_dir % output dir of moveFromBucket, input dir of other functions
        bucket_dir
    end
    
    methods
        
        function bf = getInf(inout_dir)
            
            % bucket_dir is the path of bucket folder bird-song,collection
            % of all song folders
            if exist('inout_dir','var')
                bf.inout_dir = inout_dir;
            end
            
        end
        
        function  moveFromBucket(bf, singleFolder) 
            % To extract several song files from a target folder
            % singleFolder is the path of a single folder in bucket server
            
            rawfiles = extract.filename(singleFolder,'*.wav');
            rawfiles = rawfiles(randperm(length(rawfiles))); % shuffling the order of the sound files
            
            num_f2c = 6; % 6 files to copy for each folder
            ids = [];
            
            iteration = 0; % how many turns of iteration
            LastOrNot = 0;
            
            filejudge = []; % judge whether a file is to be ignored ( not a good song file)
            
            while length(ids) < num_f2c
                
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
                
                fs = 32000; % hard-coded sampling frequency
                ampthres = 0.008;
                shortraw = 2;
                longraw = 15;
                shortsig = 1.5; % longer than 1.6 seconds
                longsig = 4.5;% shorter than 3 seconds
                centroid_thres = 3200;
                isi_dur_thres = 0.4; %  0.4 seconds : how long a duration will be regarded as separation of bouts
                num_dur = 6;  % number of bout
                
                parfor n = a: b
                    
                    try
                        [y,fs] = audioread(rawfiles{n});% here the y is the raw y
                    catch
                        continue
                    end
                    
                    %<*>Judge by the length of the raw signal
                    if  longraw < length(y)/fs  || length(y)/fs < shortraw
                        filejudge(n) = 0;
                        continue
                    end
                    
                    fiy = highpass(abs(y),500,fs); %  to remove the noise generataed by low-ftrquency noise
                    ampenv = envelope(fiy,320*9,'rms'); % amplitude envelope
                    %figure; plot(ampenv); figure; draw.spec(y,fs);% 1ms-32
                    
                    sigs = find(ampenv>ampthres); % sigs means significant signal points
                    sigy = fiy(sigs);
                    
                    %<*>Judge by the length of significanty signals
                    if  longsig < length(sigs)/fs  || length(sigs)/fs < shortsig
                        filejudge(n) = 0;
                        continue
                    end
                    
                    %<*>Judge by centroid to remove noise0like signal
                    
                    centroid{n} = spectralCentroid(sigy,fs);
                    mean_centroid = mean( centroid{n} );
                    if  mean_centroid < centroid_thres
                        filejudge(n) = 0;
                        continue
                    end
                    
                    %<*>Judge by percentage of signaled time duration to remove calls
                    percentile_thres = 0.3; % 30%
                    persigs(n) = length(sigs)/length(y); % percentage of signals that are significant
                    if persigs(n) < percentile_thres
                        filejudge(n) = 0;
                        continue
                    end
                    
                    %<*>Judge Inter-segment-interval: Too strict criteria, to remove songs which has too many bolts
                    
                    interval = diff(sigs);
                    isi{n} = interval(interval~=1); % ISI is the inter-syllable intervals
                    long_isi = isi{n}(isi{n}> fs*isi_dur_thres);
                    %if length(long_isi) ~= num_dur % this criteria might be too strict
                    fprintf('长间隔的数目为: %u\n',length(long_isi));
                    if length(long_isi) > num_dur
                        
                        filejudge(n) = 0;
                        continue
                    end
                    
                    % calculate spec and rawspec
                    
                    
                    filejudge(n) = 1;
                    
                    %diji = length(find(filejudge));
                    %fprintf('得到了第%u个合格的Song！其在文件夹中的序列是 %u\n',diji,n);
                    fprintf('得到了一个合格的Song！其在文件夹中的序列是 %u\n',n);
                    %disp('wtf')
                end
                
                
                ids = find(filejudge);
                fprintf('目前收获了_%u_首Song',length(ids));
                % extract 5 files from each folder
                if length(ids) > num_f2c
                    ids = ids(1: num_f2c );
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
                %spec{w}  = imresize(cal.spec(fiy(sigs),fs),[257 2000]);
                % rawspec{w}  = imresize(cal.spec(fiy,fs),[257 2000]) ;
            end
            
            % draw the spec of the rest of the signals
            % figure
            % montage(spec)
            % colormap('jet')
            
            % figure
            % montage(rawspec) % raw spectrogram of selected
            % colormap('jet')
            
            target_dir = bf.inout_dir; % "E:\WavsCollection" % destination
            temptemp = split(singleFolder ,'\');
            birdid = temptemp{end};
            output_dir = sprintf('%s\\%s',target_dir,birdid);
            mkdir(sprintf('%s\\%s',target_dir,birdid));
            
            parfor a = 1: length(picked_files)
                
                oldname = picked_files{a};
                temp = split(oldname,'\');
                %     birdid = temp{end-1};
                individual = temp{end};
                %destination = sprintf('%s\\%s\\%s-%s',target_dir,birdid,birdid,individual);
                destination = sprintf('%s\\%s\\%s-%u.wav',target_dir,birdid,birdid,a);
                copyfile(oldname,destination);
                
            end
            
            
            
        end
        

    end
    
    
    methods(Static)
        
        
        function  moveNoiseFiles(sourcedir, destineydir) 
            % to move pure noise files out of the sourcedir
            tic;
            
            rawfiles = extract.filename(sourcedir,'*.wav');
            
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
        
        
        
        function all_eleinf = Fraginf(subdirs,id_thres,matFolder)  % 遍历 all the subdirs
            
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
                
                matfiles = extract.filename(sprintf('%s\\%s',subdirs{r},matFolder),'*.mat'); % matFolder is the folder containing segdata.mat,
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
                    I = cal.spec(fiy,fs);
                    
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
        
        function all_eleinf = Eleinf(subdirs,id_thres,matFolder)  % 遍历 all the subdirs
            
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
                
                matfiles = extract.filename(sprintf('%s\\%s',subdirs{r},matFolder),'*.mat'); % matFolder is the folder containing segdata.mat,
                % which can be SegData or SylData or EleData
                
                song_collect = {};
                
                for m = 1: length(matfiles)
                    
                    loaded = load(matfiles{m});
                    
                    if isfield(loaded.segdata,'eleedge')
                        two_eleedge = repmat(loaded.segdata.eleedge,[2,1]);
                        alledges = sort( vertcat(loaded.segdata.syledge(:),reshape(two_eleedge,[],1) ));
                    else
                        alledges = sort(loaded.segdata.syledge(:));
                    end
                    
                    
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
                    I = cal.spec(fiy,fs);
                    
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
        
        function all_eleinf = Sylinf(subdirs,id_thres,matFolder)  % 遍历 all the subdirs
            
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
                
                matfiles = extract.filename(sprintf('%s\\%s',subdirs{r},matFolder),'*.mat'); % matFolder is the folder containing segdata.mat,
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
                    I = cal.spec(fiy,fs);
                    
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
        
        function all_eleinf = getAllTwoMotifEleinf(two_motif_ele_dir)
            if ~exist('two_motif_ele_dir','var')
                two_motif_ele_dir = "G:\StimuliSource";
            end
            
            subdirs = extract.folder(two_motif_ele_dir).';
            twomotif_ids = find(~cellfun(@isempty, regexp([subdirs{:}].','twoMotif')));
            twoMotif_subdirs = subdirs(twomotif_ids).';
            
            all_eleinf = getInf.Eleinf(twoMotif_subdirs,1,'SegData');
            
        end
        
    end
end

