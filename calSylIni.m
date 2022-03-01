function drawAlignedNormFragTwoPlots(a)

tic

fragids = find(~cellfun(@isempty, regexp({a.list.stimuliname}.','frag|Frag|syl|Syl') ));

fraglist = a.list(fragids);

normlist = a.normlist;


for m = 1: length(fraglist)

    
    birdid = convert.bid(fraglist(m).stimuliname);

    ids_norm = find(~cellfun(@isempty, regexp({normlist.stimuliname}.',birdid) ) );
     
    if ~isempty(ids_norm)
    fraglist(m).sylIni = findIni(normlist(ids_norm).plty,fraglist(m).y);
    fraglist(m).pady = [zeros([fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext,1]);fraglist(m).plty;zeros(length(normlist(ids_norm).plty)...
        - (fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext) - length(fraglist(m).plty),1)];
    fraglist(m).padsptimes = cellfun( @(x) x + (fraglist(m).sylIni-1-fraglist(m).fs*fraglist(m).pltext)/fraglist(m).fs, fraglist(m).pltsptimes,'uni',0);
    end
end



% merge the new fraglist with the normlist
I_of_each_norm = {};
for w = 1: length(normlist)

    birdid = convert.bid(normlist(w).stimuliname);

    ids_infrag = find(~cellfun(@isempty, regexp({fraglist.stimuliname}.',birdid) ) );

    selected_fraglist = fraglist(ids_infrag);

    [~,temp_index] = sortrows([selected_fraglist.sylIni].'); 
    selected_fraglist = selected_fraglist(temp_index);


    Icollect = {};
    figure('Color','w');
    
    draw.two(normlist(w).plty,normlist(w).fs,normlist(w).pltsptimes);
    frame = getframe(gcf);
    Icollect{1} = frame.cdata;
    close(gcf)
 
    for hh = 1: length(selected_fraglist)

        figure('Color','w');
        draw.two(selected_fraglist(hh).pady,selected_fraglist(hh).fs,selected_fraglist(hh).padsptimes);
        frame = getframe(gcf);
        Icollect{1 + hh} = frame.cdata;
        close(gcf);

    end


    I_of_each_norm{w} = vertcat(Icollect{:});

   % figure; imshow(Imerged)
    %imwrite(Imerged,'TEST.png')
    %saveas(gcf,'TEST.png')

end



% padding each I based on the maximum size of local I
size1 = [];

for oo = 1: length(I_of_each_norm)
size1(oo) = size(I_of_each_norm{oo},1)
end

[max_size1,max_oo] = max(size1);

Ipad = {};
for oo = 1: length(I_of_each_norm)
    localI = I_of_each_norm{oo};
    Ibase= uint8(256*ones(size(I_of_each_norm{max_oo})));
    Ibase(1:size(localI,1),1:size(localI,2),1:size(localI,3)) = localI;
    Ipad{oo} = Ibase;
end


Iall = horzcat(Ipad{:});

imwrite(Iall,sprintf('TWO_Normfrag_%s.png',a.unique_neuronname));
toc

function SynIni = findIni(y, yfrag)
%y = B521; yfrag = B521_7;

counts = 0;
diff_info = struct;
for k = 1: length(y) - length(yfrag)
    totest = y(k:k+ length(yfrag) -1);

    if sum(totest) ~= 0

        diff = totest - yfrag;
        if diff== 0


            disp('Catch it')
            trump = 'MAGA';
            disp(k)
            SynIni = k;
            return
        else 
            counts = counts + 1;

            diff_info(counts).diff = diff;
            diff_info(counts).k = k;
        end
    end

end

[~,min_ids] = min([diff_info.diff].');
SynIni = diff_info(min_ids).k;

end


end