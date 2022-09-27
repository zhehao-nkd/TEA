classdef VocClassifier
    
    methods(Static)
        
        function  extractSongs(sourceDir, destiDir,NumToExtract)
            % To Extract several song files from a target folder
            % singleFolder is the path of a single folder in bucket server
            
            rawfiles = Extract.filename(sourceDir,'*.wav');
            %rawfiles = rawfiles(randperm(length(rawfiles))); % shuffling the order of the sound files
            
            n2e = NumToExtract; % 6 files to copy for each folder
            ids = [];
            iteration = 0; % how many turns of iteration
            LastOrNot = 0;
            filejudge = []; % judge whether a file is to be ignored ( not a good song file)
            nums_per_section = 100;
            while length(ids) < n2e % evalue 100 files for each section
                
                iteration = iteration + 1;
                fprintf('# %u iterations \n',iteration);
                
                if length(rawfiles) > nums_per_section* iteration
                    first = 1 + (iteration-1)*nums_per_section;
                    last = iteration*nums_per_section;
                else
                    first = 1;
                    last = length(rawfiles);
                    disp('  Last 100 files in the folder!');
                    LastOrNot = 1;
                end

                centroid = {};
                ampthres = 0.008;
                shortraw = 2; % if the file is less than 2 seconds, ignore it
                longraw = 15; %15;
                shortsig = 1.5; % longer than 1.5 seconds
                longsig = 4.5; %4.5;% shorter than 4.5 seconds
                centroid_thres = 3700; % this thres may be too strict
                isi_dur_thres = 0.4; %  0.4 seconds : how long a duration will be regarded as separation between motifs
                max_num_dur = 6;  % maximum number of motif
                percentile_thres = 0.3; % 30%
                
                for n = first: last % par-for
                    
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
                    fprintf('  number of motifs: %u  , but not very trustable \n',length(long_isi)+1); % indicate number of motifs
                    if length(long_isi) > max_num_dur
                        filejudge(n) = 0;
                        continue
                    end
                    
                    % calculate spec and rawspec
                    filejudge(n) = 1;
                    %diji = length(find(filejudge));
                    %fprintf('得到了第%u个合格的Song！其在文件夹中的序列是 %u\n',diji,n);
                    fprintf('  You get a song file！Its index in the folder is %u\n',n);
                end
                
                
                ids = find(filejudge);
                fprintf('  Now you''ve got _%u_Songs \n',length(ids));
                % Extract n2e files from each folder
                if length(ids) > n2e
                    ids = ids(1: n2e );
                end
                
                if LastOrNot == 1
                    disp('  Not enough songs!');
                    break
                    
                end
                
            end
            
            
            picked_files = rawfiles(ids);
            
            if ~isempty(picked_files)
                temptemp = split(sourceDir,'\');
                birdid = temptemp{end};
                mkdir(sprintf('%s\\%s',destiDir,birdid));
                parfor first = 1: length(picked_files)
                    oldname = picked_files{first};
                    destination = sprintf('%s\\%s\\%s-%u.wav',destiDir,birdid,birdid,first);
                    copyfile(oldname,destination);
                end
            else
                disp('This folder contains no songs that meet the requirements');
                return
            end
            
           
        end
    end
    


end