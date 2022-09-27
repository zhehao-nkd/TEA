wav_dir = "E:\WavsCollection";
subdirs = Extract.folder(wav_dir);
fs = 32000;
tic

for r = 1:length(subdirs)
    
    mat_path = Extract.filename(sprintf('%s\\SegData',subdirs{r}),'*.mat');
    
    if isempty(mat_path)
        continue
    end
    
    for m = 1: length(mat_path)
        
        load(mat_path{m});
        
        replica_eleedge = repmat(segdata.eleedge,[2,1]);
        alledges = sort( vertcat(segdata.syledge(:),reshape(replica_eleedge,[],1) ));
        initials = alledges(1:2:end);
        terminals = alledges(2:2:end);
        fiy = bandpass(segdata.rawy,[900 6000],fs); %% It is very important that here the y should be fiy !!!!! filtered y instead of the raw y
        I = Cal.spec(fiy,fs);
       
        song_eleinf = struct;
        
        for w = 1: length(initials) % can add a par-for
            song_eleinf(w).initial = initials(w)*fs/1000;
            song_eleinf(w).terminal = terminals(w)*fs/1000;
            song_eleinf(w).songname = segdata.birdid;
            song_eleinf(w).y = segdata.rawy(initials(w)*fs/1000:terminals(w)*fs/1000);
            song_eleinf(w).fs = fs;
            song_eleinf(w).fragI = imresize(I(:,initials(w):terminals(w)),[257,50]);
            
        end
        
        parfor f = 1:length(song_eleinf)
            song_eleinf(f).fragid = f;
        end
        
        song_collect {m} = song_eleinf;
        
    end
    
    bird_eleinf = vertcat(song_collect{:});
    
    % @-----@ This section is for finding the similar fragments
    frag_img = {};
    for u = 1: length(bird_eleinf)
        frag_img{u} = bird_eleinf(u).fragI;
    end
    
    figure; montage(frag_img)
    
    heatmap_count = 0;
    heatT = table;
    img_sim = []; % image similarityu matrix
    for z = 1: length(frag_img)
        for s = 1: length(frag_img)
            img_sim(z,s) = ssim(frag_img{z},frag_img{s});
            
            % The following two lines are used to create the heat map
            heatmap_count = heatmap_count + 1;
            heatT.first(heatmap_count) = z;
            heatT.second(heatmap_count) = s;
            heatT.sim(heatmap_count) = img_sim(z,s);
            
        end
    end
   
    
    img_sim = (img_sim + img_sim.')/2; % make img_sim symmetric
    min_sim = 0.8;% minimum similarity to sort fragments together
    % figure; imagesc(img_sim); % plotting
    img_sim_tri = triu(img_sim,1); % triangular
    
    figure; heatmap(heatT,'first','second','ColorVariable','sim') 
    
    
    for c = 1:length(img_sim)
        group{c} = find(img_sim(c,:)>0.7)  %% I don't know whether this value is good or not !!!!!!!
    end
    
    %<solution-2>
    num_of_catego = 0;
    not_accessed = linspace(1,length(group),length(group));
    categos = {};
    
    for k = 1 :length(group)
        if ismember(k,not_accessed)
            
            num_of_catego = num_of_catego + 1;
            
            categos{num_of_catego} =  group{k};
            while sum(cellfun(@length,cellfun(@(x) intersect(categos{num_of_catego},x),group(not_accessed),'UniformOutput',0)))~= 0
                ids= find(cellfun(@length,cellfun(@(x) intersect(categos{num_of_catego},x),group(not_accessed),'UniformOutput',0)));
                categos{num_of_catego} = unique(horzcat(categos{num_of_catego}, horzcat( group{not_accessed(ids) })));
                not_accessed(ids) = []; % delete accessed
                disp('Loop!!!!!!!')
            end
            
        else
            continue
        end
    end
    
    
    % assign the cartego info to th eleinf
    
    for trump = 1: length(bird_eleinf)
        
        this_catego = find(cellfun(@(x) ismember(trump,x),categos));
        bird_eleinf(trump).catego = this_catego;
        
    end
   
    fprintf('Now__%u__of %u dirs are processed',r,length(subdirs));
    birdcollect{r} = bird_eleinf;
    
end

unique_eleinf = horzcat(birdcollect{:});

toc