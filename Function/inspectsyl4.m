% h=brush;
% set(h,'Enable','on','ActionPostCallback',@brushedDataCallback);
%
% function brushedDataCallback
% disp('fuck')
% end

%sapsylinf = n.sapsylinf;
function inspectsyl4(sylinf)


idx = find ([sylinf.label] ==1);
labeled = sylinf(idx);

figure('Color','w','Position',[1946 -103 1178 1086]);
scatter([sylinf.pitch].',[sylinf.fm].','k','filled');

%hold on
%[C,ia,ic] = unique({sylinf(:).sound})

[~,~,ic] = unique({labeled(:).sound});
idxs = unique(ic);
for k = 1: length(idxs)
    thissong = labeled(ic == idxs(k));
    hold on
    scatter([thissong.pitch].',[thissong.fm].','filled');
end

%scatter([labeled.pitch].',[labeled.fm].','r','filled');



xlabel('Pitch');

ylabel('FM');
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
          extendPlot(sidx,sylinf);
%             bu = uibutton(gcf);
%             bu.Text = 'Plot';
        end
        
        %autoArrangeFigures()
        
    end

%set(h,'ActionPostCallback',@brushedDataCallback);
%set(h,'ActionPostCallback',@brushedDataCallback2);


end



