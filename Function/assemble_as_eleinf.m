
function all_eleinf = assemble_as_eleinf(subdirs,id_thres,matFolder)  % 遍历 all the subdirs

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