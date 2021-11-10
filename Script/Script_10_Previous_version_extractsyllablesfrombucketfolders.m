% ver2, judge one song for each time
tic;
% for generating datasets
dbstop if error
%folders = getFolders('/Volumes/bucket/Yazaki-SugiyamaU/Bird-song');
%folders =  getFolders('/Users/Chengzh/Dropbox (OIST)/ForIdentifyingSongs');
pathlist = 'BirdListMatlab.xlsx';
pathlog = 'BirdlogMatlab.xlsx';
load('ISOCELL.mat');
%ISOCELL = extractIso(pathlog, pathlist);


folders = extract.folder('Y:\Yazaki-SugiyamaU\Bird-song');
folders = rmBadfolder(ISOCELL, folders); % remove isolated birds, and wierd things
folders = folders(randperm(numel(folders))); % random order
%load('folders.mat');

toc;
% judge bird, whether isolated or not
%[~,birdid,~] = cellfun(@fileparts, folders, 'UniformOutput', false);
%

outdir = 'matfiles';
%nr = 0; % number of records
mkdir(outdir)

disp(' About 4 hours is needed for running one trial!!!');

for n = 1141: length(folders)
    %for n = 1: length(folders)
    
    fprintf('Current Folder:%s, %u out of %u ',folders{n}, n, length(folders));
    newline;
    
    % for each folder
    [~,inputid,~] = fileparts(folders{n});
    filenames = extract.filename(folders{1,n},'*.wav');
    %filenames = extractAdult(filenames,inputid,pathlist);
    %filenames = restrictCentroid(filenames); % restricit filenames by its centroid
    filenames = filenames(randperm(numel(filenames))); % randomize filenames
    
    if isempty(filenames)
        disp('This folder contain only juvenile songs');
        newline;
        continue
    end
    

    nsongs = 0; % to count number of verified songs
    
    X = 1; % tocontrol whether to jump to next loop or not
    N = 0; % file idx
    part = 0; % save part
    
    minlen = min(50, length(filenames));
    
    parfor fileN = 1: minlen
        
        %%%%%% newly added
        [y,fs] = audioread(filenames{fileN});
        %gpuy = gpuArray(y);
        measurey = abs(y);
        measurey = highpass(measurey,400,fs);
        
        
        [yup,~] = envelope(measurey,800,'rms');
        %     figure
        %     plot(yup)
        %     figure
        [lowpks,lowlcs] = findpeaks(-yup,'MinPeakHeight',-0.005,'MinPeakDistance',500);
        %     figure
        %     findpeaks(yup,'MinPeakHeight',0.01,'MinPeakDistance',500);
        [highpks,highlcs] = findpeaks(yup,'MinPeakHeight',0.015,'MinPeakDistance',500);
        
        syllable = struct;
        syllable(1).label = 0;
        for lcsN = 2: length(lowlcs)
            
            
            
            if  ~isempty(highlcs(lowlcs(lcsN-1)<highlcs&highlcs<lowlcs(lcsN)))
                
                syllable(lcsN).label = 1;
                syllable(lcsN).initial = lowlcs(lcsN-1);
                syllable(lcsN).terminal = lowlcs(lcsN);
                syllable(lcsN).y = y(lowlcs(lcsN-1):lowlcs(lcsN));
                syllable(lcsN).path = filenames{fileN};
                
            else
                syllable(lcsN).label = 0;
            end
            
            %                 figure('Visible','off');
            %                 mySpectrogram(syllable(idx).y,fs,'Decide');
            %                 temp = getframe(gcf)
            %                 syllable(idx).I = temp.cdata;
            %                 close(gcf);
        end
        
        ind=find([syllable.label]);  % remove fake data
        if ~isempty(ind)
            syllable=syllable(ind);
        else
            syllable = [];
        end
        
        collect{fileN} = syllable;
        
        
    end
    
    
    syllables = horzcat(collect{:});
    
    try
        if ~isempty(syllables)
            eval([inputid,'_table',' = struct2table(syllables);']);
            save([outdir,'/',inputid,'_table','.mat'],[inputid,'_table']);
            eval(['clear ',inputid,'_table;']);
      
        end
    catch error
        display(error.message);
    end
    
    toc;
end


toc;

%send_mail_message('zhehao.cheng@oist.jp','FinishedSyllablesSegmentation','FinishedSyllablesSegmentation');
