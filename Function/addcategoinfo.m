function unique_eleinf = addcategoinfo(mergedeleinf)

songs = unique({mergedeleinf.songname}.');

% this section split the whole eleinf to separated ones for each unique
% songname
for omega = 1: length(songs)
    
    this_ids = find( strcmp({mergedeleinf.songname}.',songs{omega}));
    
    bird_ele_collect{omega} = mergedeleinf(this_ids);
end

birdcollect = {};
parfor omega = 1: length(bird_ele_collect)
 birdcollect{omega} = core( bird_ele_collect{omega});
end
    
 unique_eleinf = horzcat(birdcollect{:});   
 
 fields = {};
 for p = 1: length(birdcollect)
     fields{p} = fieldnames(birdcollect{p});
     disp(fields{p}.')
     disp(p)
 end
end



function bird_eleinf = core(bird_eleinf)
% @-----@ This section is for finding the similar fragments
frag_img = {};
for u = 1: length(bird_eleinf)
    frag_img{u} = bird_eleinf(u).fragI;
end

%figure; montage(frag_img)

heatmap_count = 0;
%heatT = table;
img_sim = []; % image similarityu matrix
for z = 1: length(frag_img)
    parfor s = 1: length(frag_img)
        img_sim(z,s) = ssim(frag_img{z},frag_img{s});
        
        % The following two lines are used to create the heat map
        %heatmap_count = heatmap_count + 1;
        %             heatT.first(heatmap_count) = z;
        %             heatT.second(heatmap_count) = s;
        %             heatT.sim(heatmap_count) = img_sim(z,s);
        
    end
end


img_sim = (img_sim + img_sim.')/2; % make img_sim symmetric
%min_sim = 0.7;% minimum similarity to sort fragments together
% figure; imagesc(img_sim); % plotting
%img_sim_tri = triu(img_sim,1); % triangular

%figure; heatmap(heatT,'first','second','ColorVariable','sim')
group = {}; %v this is very important for par-for

for c = 1:length(img_sim)
    group{c} = find(img_sim(c,:)>0.7);  %% I don't know whether this value is good or not !!!!!!!
end

if isempty(group) % in case that group is empty
    for l = 1: length(bird_eleinf)
        bird_eleinf(l).catego = l;
    end
    return
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
            disp('Loop!!!')
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


% remove field fragI
bird_eleinf = rmfield(bird_eleinf,'fragI');

%fprintf('Now__%u__of %u songs are processed',omega,length(songs));

fprintf('祇今尚有清流月，曾照高王万马过')

end