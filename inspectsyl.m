% h=brush;
% set(h,'Enable','on','ActionPostCallback',@brushedDataCallback);
%
% function brushedDataCallback
% disp('fuck')
% end

%sapsylinf = n.sapsylinf;
function inspectsyl(sylinf)


idx = find ([sylinf.label] ==1);
labeled = sylinf(idx);

figure('Color','w','Position',[1946 -103 1178 1086]);
scatter([sylinf.pitch].',[sylinf.fm].','k','filled');

hold on

scatter([labeled.pitch].',[labeled.fm].','r','filled');

xlabel('Pitch');
ylabel('FM');
%title(n.neuronname);


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
            figure;
            draw.spec(sylinf(sidx(kk)).y,32000);
            title(sprintf('%s-%u',sylinf(sidx(kk)).sound,sylinf(sidx(kk)).number));
            
        end
        
        autoArrangeFigures()
        
    end

%set(h,'ActionPostCallback',@brushedDataCallback);
%set(h,'ActionPostCallback',@brushedDataCallback2);


end



