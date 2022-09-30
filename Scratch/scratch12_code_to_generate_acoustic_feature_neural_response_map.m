% scratch 12
% this code generate acoustic feature - neural response map
%addpath(genpath("C:\Users\Zhehao\Dropbox (OIST)\My_Matlab\TEA"));


for n = 1: 16
    
f = Fenxi(n);

exp_elist = f.toAcousticSpace;
all_elist = eleinf;


% construct the new all_elist with labels
for m = 1: length(all_elist)
    %nn.Three;
   thesong = find( ~cellfun(@isempty,regexp({exp_elist.stimuliname}.',all_elist(m).songname)));
   
   thesyl = find( ~cellfun(@isempty,regexp({exp_elist.stimuliname}.',sprintf('-%02d-',all_elist(m).fragid))));
   thesongsyl = intersect(thesong,thesyl);
   
   if isempty(thesongsyl)
       all_elist(m).response = -1;
       all_elist(m).maxvalue = -1;
       all_elist(m).halfsum = -1;
       all_elist(m).fullsum = -1;
   else
       all_elist(m).response = exp_elist(thesongsyl).label;
       all_elist(m).maxvalue = exp_elist(thesongsyl).maxvalue;
       all_elist(m).halfsum = exp_elist(thesongsyl).halfsum;
       all_elist(m).fullsum = exp_elist(thesongsyl).fullsum;
   end
 
    
end


% draw the scatter plot

minusone = all_elist([all_elist.response].' == -1);

zero = all_elist([all_elist.response].' == 0);

one = all_elist([all_elist.response].' == 1);

figure;

vext = 0.1 % vertical extension to esnure that the scatter is above the surface
coom1 = cell2mat({minusone.embedded_z}.'); % coordinate-minusone
scatter3(coom1(:,1),coom1(:,2),repmat([0],size(coom1(:,1)) )+vext,[],repmat([.6 .6 .6],size(coom1(:,1)) ), 'filled');

hold on
coo1 = cell2mat({one.embedded_z}.');
if ~isempty(coo1)
    scatter3(coo1(:,1),coo1(:,2),[one.maxvalue].'/max([one.maxvalue].')*10+vext,[],[one.maxvalue].'/max([one.maxvalue].')*10,'filled');
end
hold on
coo0 = cell2mat({zero.embedded_z}.');
scatter3(coo0(:,1),coo0(:,2),[zero.maxvalue].'/max([zero.maxvalue].')*10+vext,[],[zero.maxvalue].'/max([zero.maxvalue].')*10,'filled');
colorbar
xlim([-20 20])
ylim([-20 20])
view(0,-90)
pause(0.2)
saveas(gcf,sprintf('newscatter-%u.png',n))

%%%%%%%%%%%%%%%%%%%%%%%%% plot 2
tested_elist = all_elist([all_elist.maxvalue].' ~= -1);
xy = cell2mat({tested_elist.embedded_z}.');
x = double(xy(:,1));
y = double(xy(:,2));
v = [tested_elist.maxvalue].';
v = v/max(v)*10; % normalize v
% dangerous!!
%v(v==-1) = 1
[xq,yq] = meshgrid(-20:2:20, -20:2:20);
vq = griddata(x,y,v,xq,yq);
% plotting
figure
%mesh(xq,yq,vq)
surf(xq,yq,vq)
hold on
%plot3(x,y,v,'o')
%view(2)

hold on
%figure
vext = 0.1 % vertical extension to esnure that the scatter is above the surface
coom1 = cell2mat({minusone.embedded_z}.'); % coordinate-minusone
scatter3(coom1(:,1),coom1(:,2),-(repmat([0],size(coom1(:,1)) )+vext),[],repmat([.6 .6 .6],size(coom1(:,1)) ), 'filled');
hold on
coo1 = cell2mat({one.embedded_z}.');
if ~isempty(coo1)
    scatter3(coo1(:,1),coo1(:,2),-([one.maxvalue].'/max([one.maxvalue].')*10+vext),[],[one.maxvalue].'/max([one.maxvalue].')*10,'filled');
end
hold on
coo0 = cell2mat({zero.embedded_z}.');
scatter3(coo0(:,1),coo0(:,2),-([zero.maxvalue].'/max([zero.maxvalue].')*10+vext),[],[zero.maxvalue].'/max([zero.maxvalue].')*10,'filled');
colorbar
view(0,-90)
xlim([-20 20])
ylim([-20 20])


saveas(gcf,sprintf('newsurf-%u.png',n))
saveas(gcf,sprintf('fignewsurf-%u.fig',n))

end