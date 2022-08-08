
data_for_origin = [];
for k = 1: length(bstruct)
    
   data_for_origin(k,:) =  [bstruct(k).songinf.neunum]
    
end

birdids = {bstruct.bname}.';

songids = {bstruct(1).songinf.songname}.';