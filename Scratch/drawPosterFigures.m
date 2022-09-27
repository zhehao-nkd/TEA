sylinf = struct;
for k = 1:length(labels)
    sylinf(k).z1 = encode_z(k,1);
    sylinf(k).z2 = encode_z(k,2);
    sylinf(k).label = labels(k);
    temp = original_birdsong_stft(k,:,:);
    sylinf(k).song = temp;
    
end

labeled = sylinf(find([sylinf.label].' == 1));
unlabeled = sylinf(find([sylinf.label].' == 0));


figure;
scatter([labeled.z2].',[labeled.z1].','filled')
hold on
scatter([unlabeled.z2].',[unlabeled.z1].','filled')
set(gca,'YTickLabel',[],'XTickLabel',[]);
%set(gca,'xtick',0.001,'ytick',0.001)

% determine position of the axes
axp = get(gca,'Position');

% determine startpoint and endpoint for the arrows 
xs=axp(1);
xe=axp(1)+axp(3)+0.04;
ys=axp(2);
ye=axp(2)+axp(4)+0.05;


% add label text

for m = 1: length(sylinf)
    dx = 0.001; dy = 0.001; % displacement so the text does not overlay the data points
    c = list(m).stimuliname;
    x = double(sylinf(m).z2);
    y = double(sylinf(m).z1);
%text(x, y, c);
hold on
end
% make the arrows
annotation('arrow', [xs xe],[ys ys]);
annotation('arrow', [xs xs],[ys ye]);
xlabel('Dim1')
ylabel('Dim2')


grid on




[y,fs] = audioread(...
"C:\Users\Zhehao\Dropbox (OIST)\My_Stimuli\SourceSongSet1@06232021\Y515A.wav");
figure; Draw.spec(y,fs);
set(gca,'FontSize',20);
% determine position of the axes
axp = get(gca,'Position');

% determine startpoint and endpoint for the arrows 
xs=axp(1);
xe=axp(1)+axp(3)+0.04;
ys=axp(2);
ye=axp(2)+axp(4)+0.05;
% make the arrows
annotation('arrow', [xs xe],[ys ys]);
annotation('arrow', [xs xs],[ys ye]);




