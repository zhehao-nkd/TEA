addpath(genpath("C:\Users\Zhehao\Dropbox (OIST)\My_Matlab\TEA"));

wavpath =  "E:\DiskDBackUp\My_NewStimuli\All_Trimmed\B346A.wav"
[y,fs] = audioread(wavpath);
[~,birdid,~] = fileparts(wavpath);

y = bandpass(y,[800 5000],fs);
img = Cal.spec(y,fs);

threshold = 0.65; % best threshold might be 0.65
high = 1;
low = 0.62;

sld1 = 1;% sliding window in ms
for k = 1+ sld: length(img)-1-sld
    sv1(k) = norm( mean(img(:,k+1:k+1+sld1),2) - mean(img(:,k-sld1:k),2) ); % variable 1 for syllable segmentation
end


ssv1 = conv(sv1,flattopwin(21),'same'); % smoothed syllable variable one


ints = polyxpoly(linspace(1,length(ssv1),length(ssv1)),ssv1,linspace(1,length(ssv1),length(ssv1)),threshold*ones(1,length(ssv1)),'unique')

% figure; plot(var1); hold on; plot(conv(var1,flattopwin(27),'same')); hold on;  smooth(var1,'rloess',30)
% hold on;

for e = 1:length(ints)
    line([ints(e),ints(e)],[0,10]);
    hold on
end

ints1 = [];
[vpks, vlocs, vw, vp] = findpeaks(-ssv1);
vpks = -vpks;
[ppks, plocs, pw, pp] = findpeaks(ssv1);

figure
findpeaks(-ssv1)
hold on;
for e = 1:length(ints)
    line([ints(e),ints(e)],[0,10]);
    hold on
end

for k = 1: length(ints)/2
    
    peaks = ppks(intersect(find(ints(2*k-1)<plocs),find(plocs<ints(2*k))) ); % peaks within this syllable duration
    
    if length(find(peaks>high)) > 1
        ints1 = [ints1, ints(2*k-1),ints(2*k)];
    end
    
end

ints2 = [];

for k = 1: length(ints1)/2-1
    valleys = vpks(intersect(find(ints1(2*k)<vlocs),find(vlocs<ints1(2*k+1))) );
    disp(length(find(valleys<low)));
    if length(find(valleys<low)) <= 0
        ints2(2*k-1) = nan;
        ints2(2*k) = nan;
    end
    
end

ints2 = rmmissing(ints) % remove nans % transfer seconds to miliseconds
ints2 = ints2*0.001;



% segment syllables into elements
syllocs = {};
for w = 1: length(ints2)/2  % fort each syllable, check whether it can be further segmented into smaller elements
    
    syl_initial = ints2(2*w-1);
    
    syl = y( ints2(2*w-1)*fs : ints2(2*w)*fs );
    sylimg = Cal.spec(syl,fs)
    
    
    
    sld = 5;
    elevar = [];
    elevar2 = [];
    elevar3 = [];
    elevar4 = []; % initialization
    
    for k = 1+ sld: size(sylimg,2)-1-sld
        
        elevar(k) = norm( mean(sylimg(:,k+1:k+1+sld),2) - mean(sylimg(:,k-sld:k),2) )/norm( mean(sylimg(:,k-sld:k),2) ); % normalized by pre
        elevar2(k) = norm( mean(sylimg(:,k+1:k+1+sld),2) - mean(sylimg(:,k-sld:k),2) )/norm( mean(sylimg(:,k+1:k+1+sld),2) ); % normalized by pro
        elevar3(k) = norm( mean(sylimg(:,k+1:k+1+sld),2) - mean(sylimg(:,k-sld:k),2) ); %v not normalized
        elevar4(k) = norm( mean(sylimg(:,k+1:k+1+sld),2) - mean(sylimg(:,k-sld:k),2) )/mean(norm( [mean(sylimg(:,k+1:k+1+sld),2),mean(sylimg(:,k-sld:k),2)] ));
        %     a = mean(img(:,k+1:k+1+sld),2)/norm(mean(img(:,k+1:k+1+sld),2));
        %     b = mean(img(:,k-sld:k),2)/norm(mean(img(:,k-sld:k),2));
        %     var(k) = norm( a-b )/norm( b );
        
    end
    
    [pks1,locs1,w1,p1] = findpeaks(elevar,"MinPeakProminence",0.15,'MinPeakHeight',0.81,'MinPeakDistance',20,'Threshold',0.05) ;
    [pks2,locs2,w2,p2] = findpeaks(elevar2,"MinPeakProminence",0.15,'MinPeakHeight',0.7,'MinPeakDistance',20,'Threshold',0.05);
    
    % plot
    
    figure
    
    imagesc(sylimg);
    yyaxis left
    
    hold on
    yyaxis right
    % plot(linspace(1+ sld,length(var)+ sld,length(var)),var,'Color','r','Linewidth',1)
    plot(linspace(1,length(elevar),length(elevar)),elevar,'Color','r','Linewidth',1)
    
    hold on
    findpeaks(elevar,"MinPeakProminence",0.15,'MinPeakHeight',0.81,'MinPeakDistance',20,'Threshold',0.05)
    
    hold on
    yyaxis right
    plot(linspace(1+ sld,length(elevar)+ sld,length(elevar)),elevar2,'Color','cyan','Linewidth',1,'LineStyle','-')
    
    hold on
    findpeaks(elevar2,"MinPeakProminence",0.15,'MinPeakHeight',0.7,'MinPeakDistance',20,'Threshold',0.05)
    hold on
    yyaxis right
    plot(linspace(1+ sld,length(elevar)+ sld,length(elevar)),elevar3,'Color','white','Linewidth',1,'LineStyle',':')
    
    % Then split the syllables into elements based on the acquired locs
    alllocs = sort([locs1,locs2],'ascend') ; % plus the time of the initial of the syllable
    
    alllocs = alllocs(alllocs>18); % 10 means 10 ms
    
    alllocs = alllocs + syl_initial*1000;
    di = diff(alllocs);
    
    ids = find(di< 15); % locs which are too close to each other
    alllocs(ids) = nan;
    alllocs = rmmissing(alllocs);
    

    
    elelocs{w} = alllocs*0.001; % colelction of the time locs for segmenting each syllable
    
    
    
   
    title(sprintf('syllable is %u',w))
    hold off
end


% check how good the segmentation is

%  for visualization of the whole song
figure
title(sprintf('num is %u',num))
imagesc(img);
yyaxis left
hold on
yyaxis right
plot(linspace(1+ sld1,length(sv1)+ sld1,length(sv1)),sv1,'Color','white','Linewidth',1,'LineStyle',':')


for k = 1: length(ints2)
    hold on
    line([ints2(k)*1000,ints2(k)*1000],[0,6],'Color','r');
end

for m = 1 : length(elelocs)
    for n = 1: length(elelocs{m})
        hold on
        line([elelocs{m}(n)*1000,elelocs{m}(n)*1000],[0,6],'Color','g');
    end
end

hold off

saveas(gcf,birdid);



% for c = 1: size(img,2) % here c represent column
%     imgnorm(:,c) = img(:,c)/norm(img(:,c));
% end
%
% figure
% title(sprintf('num is %u',num))
% imagesc(imgnorm)
%
% slid = 10 % ms
%
%
% for num = 1
% % trail 3
% sld = num;
% for k = 1+ sld: length(img)-1-sld
%
%     var(k) = norm( mean(imgnorm(:,k+1:k+1+sld),2) - mean(imgnorm(:,k-sld:k),2) )/norm( mean(imgnorm(:,k-sld:k),2) ); % normalized by pre
%
%
%     var2(k) = norm( mean(imgnorm(:,k+1:k+1+sld),2) - mean(imgnorm(:,k-sld:k),2) )/norm( mean(imgnorm(:,k+1:k+1+sld),2) ); % normalized by pro
%
%     var3(k) = norm( mean(imgnorm(:,k+1:k+1+sld),2) - mean(imgnorm(:,k-sld:k),2) ); %v not normalized
%
%
%     var4(k) = norm( mean(imgnorm(:,k+1:k+1+sld),2) - mean(imgnorm(:,k-sld:k),2) )/mean(norm( [mean(imgnorm(:,k+1:k+1+sld),2),mean(imgnorm(:,k-sld:k),2)] ));
%
% end
%
% figure
% title(sprintf('num is %u',num))
% imagesc(img);
% yyaxis left
%
% hold on
% yyaxis right
% plot(linspace(1+ sld,length(var)+ sld,length(var)),var,'Color','r','Linewidth',1)
%
% hold on
% yyaxis right
% plot(linspace(1+ sld,length(var)+ sld,length(var)),var2,'Color','cyan','Linewidth',1,'LineStyle','-')
%
% hold on
% yyaxis right
% plot(linspace(1+ sld,length(var)+ sld,length(var)),var3,'Color','white','Linewidth',1,'LineStyle',':')
%
% end



