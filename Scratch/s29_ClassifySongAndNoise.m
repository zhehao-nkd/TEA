
tic

b = Bird("Z:\Yazaki-SugiyamaU\Bird-song"); % bi is the Bird object
b.rand % randomize the folder paths
goodfolders = b.folders;

for f = 1: length(goodfolders)
    
    rawfiles = extract.filename(goodfolders{f},'*.wav');
    rawfiles = rawfiles(randperm(length(rawfiles))); % shuffling the order of the sound files
    if length(rawfiles) > 400
        files = rawfiles(1:400); % restrict the number of files
    end
    filejudge = [];
    spec = {};
    centroid = {};
    % entropy = {}; crest = {}; flux = {}; kurtosis = {}; rolloffPoint = {}; autcor = {};
    % $$$ problem is to restrict the files read in parfor loop to save memory
    parfor n = 1: length(files)
        [y,fs] = audioread(files{n});
        fiy = highpass(abs(y),500,fs); % filtered y
        ampenv = envelope(fiy,320*9,'rms'); % amplitude envelope
        %figure; plot(ampenv); % 1ms-32
        ampthres = 0.008;
        shorttime = 2; % longer than 1.5 seconds
        longtime = 4;
        sigs = find(ampenv>ampthres); % sigs means significant signal points
        
        %<*>Judge by the length of huge enough signal
        
        if  longtime > length(sigs)/fs  && length(sigs)/fs < shorttime
            filejudge(n) = 0;
            continue
        end
        
        %<*>Judge by centroid to remove noise0like signal
        centroid_thres = 3200;
        centroid{n} = spectralCentroid(fiy(sigs),fs);
        mean_centroid = mean( centroid{n} );
        if  mean_centroid < centroid_thres
            filejudge(n) = 0;
            continue
        end
        
        %<*>Judge by percentage of signaled time duration to remove some of the calls
        percentile_thres = 0.3; % 30%
        persigs(n) = length(sigs)/length(y); % percentage of signals that are significant
        if persigs(n) < percentile_thres
            filejudge(n) = 0;
            continue
        end
        
        %<*>Judge Inter-segment-interval: Too strict criteria, to remove songs which has too many bolts
        isi_dur_thres = 0.5; %  0.5 seconds
        interval = diff(sigs);
        isi{n} = interval(interval~=1); % ISI is the inter-syllable intervals
        this_isi = isi{n}(isi{n}> 32000*isi_dur_thres);
        if this_isi == 2
            filejudge(n) = 0;
            continue
        end
        
        
        % calculate spec and rawspec
        spec{n}  = cal.spec(fiy(sigs),fs);
        rawspec{n}  = cal.spec(fiy,fs);
        
        %     % Others
        %     crest{n} = spectralCrest(fiy(sigs),fs);
        %     flux{n} = spectralFlux(fiy(sigs),fs);
        %     entropy{n} = spectralEntropy(fiy(sigs),fs);
        %     kurtosis{n} = spectralKurtosis(fiy(sigs),fs);
        %     rolloffPoint{n} = spectralRolloffPoint(fiy(sigs),fs);
        %
        
        
        filejudge(n) = 1;
        
        
        
    end
    
    ids = find(filejudge);
    srawspec = rawspec(ids);
    sspec = spec(ids ); % s means selected mean
    %calculate mean
    scentroid = centroid(ids);
    mscentroid = cellfun(@mean,scentroid);
    
    
    spersigs = persigs(ids); % selected Percentage of significant signals
    figure; plot(spersigs)
    
    selected_isi = isi(ids);
    
    for g = 1: length(selected_isi)
        screened_selected_isi{g} = selected_isi{g}(selected_isi{g}> 32000*0.5); % 0.5s means 500ms
    end
    len_isi = cellfun(@length,screened_selected_isi);
    
    median_selected_isi = cellfun(@mean,selected_isi);
%     figure; plot(median_selected_isi)
    
    % screst = crest(ids); mscrest = cellfun(@mean,screst);
    % sflux = flux(ids); msflux = cellfun(@mean,sflux);
    % sentropy = entropy(ids); msentropy = cellfun(@mean,sentropy);
    % skurtosis = kurtosis (ids); msskurtosis = cellfun(@mean,skurtosis);
    % srolloffPoint = rolloffPoint(ids);msrolloffPoint = cellfun(@mean,srolloffPoint);
    
    for c = 1: length(scentroid)
        autcor{c} = autocorr(scentroid{c});
    end
    %ac_centroid = cellfun(@autocorr,scentroid, 'UniformOutput',0);
    for m = 1: length(sspec)
        sspec{m} = imresize(sspec{m},[257 10000]); % resizing each img to have the same size
    end
    
    % draw the spec of the rest of the signals
%     figure
%     montage(sspec)
%     colormap('jet')
    
%     figure
%     montage(srawspec) % raw spectrogram of selected
%     colormap('jet')
%     
    % figure;
    % further_srawspec = srawspec(spersigs > 0.25);
    % montage(further_srawspec);
    
    try
        figure;
        further_srawspec = srawspec(len_isi == 2);
        montage(further_srawspec);
        cateh Error
    end 
   
    
% % % % % %     % extract 10 files from each folder
% % % % % %     if length(ids) > 10
% % % % % %         restricted_ids = ids(1: 10);
% % % % % %     else
% % % % % %         restricted_ids = ids;
% % % % % %     end
% % % % % %     
% % % % % %     picked_files = files(restricted_ids);
% % % % % %     target_dir = "E:\WavsCollection" % destination
% % % % % %     
% % % % % %     for a = 1: length(picked_files)
% % % % % %         
% % % % % %         oldname = picked_files{a};
% % % % % %         temp = split(oldname,'\');
% % % % % %         birdid = temp{end-1};
% % % % % %         individual = temp{end};
% % % % % %         mkdir(sprintf('%s\\%s',target_dir,birdid));
% % % % % %         %destination = sprintf('%s\\%s\\%s-%s',target_dir,birdid,birdid,individual);
% % % % % %         destination = sprintf('%s\\%s\\%s-%u',target_dir,birdid,birdid,a);
% % % % % %         copyfile(oldname,destination);
% % % % % %         
% % % % % %     end
    
end
toc