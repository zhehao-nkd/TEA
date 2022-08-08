for x = 1: length(mergedeleinf)
    mergedeleinf(x).coor_1 = encode_z(x,1);
    mergedeleinf(x).coor_2 = encode_z(x,2);
end

eleciting_info = extract_eleciting_info(n);


% initiate wholecountsiglabel as 0
for cnm = 1: length(mergedeleinf)
    mergedeleinf(cnm).wholecountsiglabel = 0;
end

for y = 1: length(eleciting_info)
    if ~isempty( eleciting_info(y).eleciting)
        
        samesongids =  find(~cellfun(@isempty, regexp([mergedeleinf(:).songname].',eleciting_info(y).name)));
        for z = 1: length(eleciting_info(y).eleciting)
            samefrgids = find(arrayfun(@(s) ismember(s,eleciting_info(y).eleciting), [mergedeleinf(:).fragid].' ));
            ids_to_label = intersect(samesongids,samefrgids);
            for wtf = 1: length(ids_to_label)
                mergedeleinf(ids_to_label(wtf)).wholecountsiglabel = 1;
            end
        end
    end
end

