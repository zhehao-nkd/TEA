% Scratch this code generate secondary catego-transformation map of neueral
% response to two-element sequence

% scratch 12
% this code generate acoustic feature - neural response map
%addpath(genpath("C:\Users\Zhehao\Dropbox (OIST)\My_Matlab\TEA"));

addpath(genpath("C:\Users\Zhehao\Dropbox (OIST)\My_Matlab\TEA"));

for n = 5: 12
    
f = Fenxi(n);

exp_elist = f.to2ndAcousticSpace;


if isempty(exp_elist) % if exp_elist is empty
    continue % jump to the next loop
end

% rebuild exp-NAMES  
judgeparts = strsplit(exp_elist(1).stimuliname,'-')

if length(judgeparts) < 6
    continue
end

for k = 1: length(exp_elist)
    parts = strsplit(exp_elist(k).stimuliname,'-');
    expnames(k).song = parts{3};
    expnames(k).fragid = sprintf('%02d',str2num(parts{4}));
    expnames(k).fixedsong = parts{7};
    expnames(k).fixedfragid = parts{8};
end

%对于每一个fixed song

fixed = unique({expnames.fixedsong}.');

for ff = 1: length(fixed)
    sinexpnames = expnames( ~cellfun(@isempty, regexp({expnames.fixedsong}.',fixed{ff})) );% single expnames
    all_elist = eleinf;
% construct the new all_elist with labels
for m = 1: length(all_elist)
    %nn.Three;
    
   thesong = find( ~cellfun(@isempty,regexp({sinexpnames.song}.',all_elist(m).songname)));
   
   thesyl = find( ~cellfun(@isempty,regexp({sinexpnames.fragid}.',num2str(sprintf('%02d',all_elist(m).fragid)) )));
   
   thesongsyl = intersect(thesong,thesyl);
   
   if isempty(thesongsyl)
       all_elist(m).response = -1;
       all_elist(m).maxvalue = -1;
       all_elist(m).halfsum = -1;
       all_elist(m).fullsum = -1;
       all_elist(m).beginmax = -1;
   else
       all_elist(m).response = exp_elist(thesongsyl).label;
       all_elist(m).maxvalue = exp_elist(thesongsyl).maxvalue;
       all_elist(m).halfsum = exp_elist(thesongsyl).halfsum;
       all_elist(m).fullsum = exp_elist(thesongsyl).fullsum;
       all_elist(m).testhistory = exp_elist(thesongsyl).label;
        all_elist(m).beginmax = exp_elist(thesongsyl).beginmax;
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
if ~isempty(coo0)
    scatter3(coo0(:,1),coo0(:,2),[zero.maxvalue].'/max([zero.maxvalue].')*10+vext,[],[zero.maxvalue].'/max([zero.maxvalue].')*10,'filled');
end
colorbar
xlim([-20 20])
ylim([-20 20])
view(0,-90)
pause(0.2)
saveas(gcf,sprintf('2ndscatter-%u-before-%s.png',n,fixed{ff}))

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

if ~isempty(vq)
    
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
    if ~isempty(coo0)
        scatter3(coo0(:,1),coo0(:,2),-([zero.maxvalue].'/max([zero.maxvalue].')*10+vext),[],[zero.maxvalue].'/max([zero.maxvalue].')*10,'filled');
    end
     colorbar
    view(0,-90)
    xlim([-20 20])
    ylim([-20 20])
    
    pause(0.1)
    saveas(gcf,sprintf('second-surf-%u-before-%s.png',n,fixed{ff}))
    pause(0.1)
    saveas(gcf,sprintf('fig-second-surf-%u-before-%s.fig',n,fixed{ff}))
    send_mail_message('zhehao.cheng@oist.jp',num2str(n))
    
end

 end
end

send_mail_message('zhehao.cheng@oist.jp','Done!')