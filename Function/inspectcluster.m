% h=brush;
% set(h,'Enable','on','ActionPostCallback',@brushedDataCallback);
%
% function brushedDataCallback
% disp('fuck')
% end

%sapsylinf = n.sapsylinf;
function inspectcluster(sylinf)

cluidx = [sylinf.cluidx].';
dims = {sylinf.dims}.';
dims = cell2mat(dims);
nclu = length(unique(cluidx));
uniclu = unique(cluidx);

figure('Color','w','Position',[1946 -103 1178 1086]);

scatter(dims(:,1),dims(:,2),'filled');

for id = 1: nclu
    thisclu = find(cluidx==uniclu(id));
     hold on
    scatter(dims(thisclu,1),dims(thisclu,2),'filled');
   
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
            % convert to logical
            % x=get(h(i),'XData');
            %sx=x(find(idx));
            
            %y=get(h(i),'YData');
            % sy=y(find(idx));
            
            
        end
       
        
       
       
        for kk = 1: length(sidx)
            figure('Position',[30 178 352 869],'Color','w');
            draw.spec(sylinf(sidx(kk)).y,32000);
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
