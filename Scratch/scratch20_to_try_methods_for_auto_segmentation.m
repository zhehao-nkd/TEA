
[y,fs] = audioread(wavpath);

y = bandpass(y,[500,6000],32000);
img = Cal.spec(y,fs);

svpca = []
sld1 = 50;% sliding window size in ms
for k = 1: length(img)-sld1
    svpca(k,:) = reshape(img(:,k:k+sld1),[],1) ; % variable 1 for syllable segmentation
end
 
%[coeff,score,latent,tsquared,explained,mu] = pca(zscore(svpca,1,2),'Centered',false);

[coeff,score,latent,tsquared,explained,mu] = pca(svpca);


variable = score;



coor = [];
for y = 1: length(variable(:,1))
    coor(y,:) = variable(y,1:10); % here coor means coordinates
end

di = diff(coor);

figure 
imshow(img)
hold on

diffnorm = vecnorm(di,2,2);

% figure-1
figure;
%plot3(score(:,1),score(:,2),linspace(1,size(svpca,1),size(svpca,1)))
scatter3(variable(:,1),variable(:,2),linspace(1,size(svpca,1),size(svpca,1)))
%biplot(coeff(:,1:2),'scores',score(:,1:2));

figure
scatter3(variable(:,1),variable(:,2),variable(:,3))


% figure-2
figure;
yyaxis left
imagesc(img)
hold on 
yyaxis right
plot([zeros(sld1,1);variable(:,1)],'Color','blue') % dim 1
hold on 
plot([zeros(sld1,1);variable(:,2)],'COlor','red') % dim 2
hold off

% figure-3
figure;
yyaxis left
imagesc(img)

hold on 
yyaxis right
plot([zeros(sld1,1);diffnorm])