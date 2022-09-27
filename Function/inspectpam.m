% h=brush;
% set(h,'Enable','on','ActionPostCallback',@brushedDataCallback);
%
% function brushedDataCallback
% disp('fuck')
% end

%sapsylinf = n.sapsylinf;
function inspectcluster(sylinf)

cluidx = [sylinf.cluidx].';
pam = {sylinf.pam}.';
pam = cell2mat(pam);
nclu = length(unique(cluidx));
uniclu = unique(cluidx);

figure('Color','w','Position',[1946 -103 1178 1086]);

scatter(pam(:,1),pam(:,2),'filled');

for id = 1: nclu
    thisclu = find(cluidx==uniclu(id));
     hold on
    scatter(pam(thisclu,1),pam(thisclu,2),'filled');
   
end

xlabel('Dim1');

ylabel('Dim2');
%title(n.neuronname);
hold off

h=brush;
set(h,'Color','cyan','Enable','on','ActionPostCallback',{@brushedDataCallback,sylinf});

    function  brushedDataCallback(~,~,sylinf)
        
        h=findobj(gca,'type','scatter');
        
        for i=1:size(h)
            idxraw=logical(get(h(i),'BrushData'));
            sidx = find(idxraw);
            % Convert to logical
            % x=get(h(i),'XData');
            %sx=x(find(idx));
            
            %y=get(h(i),'YData');
            % sy=y(find(idx));
            
            
        end
       
        
       
       
        for kk = 1: length(sidx)
            figure('Position',[30 178 352 869],'Color','w');
            Draw.spec256(sylinf(sidx(kk)).y,32000);
            sound(sylinf(sidx(kk)).y,32000);
            title(sprintf('%s-%u',sylinf(sidx(kk)).sound,sylinf(sidx(kk)).number));
           % colormap('gray')
%             bu = uibutton(gcf);
%             bu.Text = 'Plot';
        end
        
       % autoArrangeFigures()
        
    end

%set(h,'ActionPostCallback',@brushedDataCallback);
%set(h,'ActionPostCallback',@brushedDataCallback2);


end
