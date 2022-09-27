
function output_dir = bucketMover(goodfolder)

rawfiles = Extract.filename(goodfolder,'*.wav');
rawfiles = rawfiles(randperm(length(rawfiles))); % shuffling the order of the sound files

num_f2c = 6; % 6 files to copy for each folder
ids = [];

iteration = 0; % how many turns of iteration
LastOrNot = 0;

 filejudge = []; % judge whether a file is to be ignored ( not a good song file)

while length(ids) < num_f2c
    
    iteration = iteration + 1;
    disp('六道轮回')
    
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
        %figure; plot(ampenv); figure; Draw.spec(y,fs);% 1ms-32
        
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
    % Extract 5 files from each folder
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
    temptemp = split(goodfolder ,'\');
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

% draw the spec of the rest of the signals
% figure
% montage(spec)
% colormap('jet')

% figure
% montage(rawspec) % raw spectrogram of selected
% colormap('jet')

target_dir = "E:\WavsCollection" % destination
temptemp = split(goodfolder ,'\');
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


