vae = readtable("C:\Users\Zhehao\Dropbox (OIST)\My_Python\NoSegment_VAE\otesib_similarity.xlsx")

figure

h = heatmap(vae,'Sound1','Sound2','ColorVariable','Similarity')

sim_matrix =    h.ColorData

figure; imagesc(h.ColorData)


Y = tsne(dis_matrix,'Distance','euclidean')


temp = rescale((sim_matrix + sim_matrix')/2)

temp(logical(eye(size(temp)))) = ones(length(temp),1)*max(max(temp))
temp = rescale(temp);
[x,y]=find(~(temp'==temp))

Y = mdscale(temp ,2)
figure; scatter(Y(:,1),Y(:,2))

hold on
scatter(Y(end-2:end,1),Y(end-2:end,2))