% scratch 10
dir = "C:\Users\Zhehao\Dropbox (OIST)\My_EphysInput\fenxi2.xls";

fea = table2struct(readtable(dir));

figure;
scatter( [fea.mean_FM].',[fea.mean_pitch].','filled');
firedidx = find([fea.label].');

onlyfire = fea(firedidx);
hold on 
scatter( [onlyfire.mean_FM].',[onlyfire.mean_pitch].','filled');

% scratch 11 
labeled = {'G548A-10','G548A-17','G548A-19','G548A-8','Y606A-10','Y606A-16','Y606A-8'};
sim = table2struct(readtable("C:\Users\Zhehao\Dropbox (OIST)\My_EphysInput\similarity.xlsx"));
%sim = sortrows(sim,'label','descend'); % sort label

for k = 1: length(sim)
    judge1 = {};
    judeg2 = {};
    for w = 1: length(labeled)
        judge1(w) =  regexp(sim(k).Sound1,labeled(w));
    end
    sum1 = [judge1{:}];
    
    if isempty(sum1)
        sim(k).label1 = 0;
    else
        sim(k).label1 = 1;
    end
    
    
    for w = 1: length(labeled)
        judge2(w) =  regexp(sim(k).Sound2,labeled(w));
    end
    sum2 = [judge2{:}];
    
    if isempty(sum2)
        sim(k).label2 = 0;
    else
        sim(k).label2 = 1;
    end
    
end
sim = struct2table(sim);
sim = sortrows(sim,'label1','descend');
%sim = sortrows(sim,{'label1','label2'},'descend');
figure;
heatmap(sim,'Sound2','Sound1','ColorVariable','Similarity')