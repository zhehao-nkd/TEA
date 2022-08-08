


pcolormap = flip(data2);
            pcolormap = vertcat(zeros(1,size(pcolormap,2)),pcolormap);
            pcolormap = horzcat(pcolormap,zeros(size(pcolormap,1),1));
            bfig = figure; s = pcolor(pcolormap);
            s.EdgeColor = 'k';
            s.LineWidth = 1.6;
            %contour(concat_respmap,'LineColor','k');
            
              map = [ 
                1 1 1
                
                
                0 0.4470 0.7410
                
               
               0 0.4470 0.7410];
           
            map = [ 
                1 1 1
                
                
                
                0 0.4470 0.7410
                
                
               1 0 0
               
               0 0.4470 0.7410
               
               0.8500 0.3250 0.0980];
            
            colormap(bfig,map);
            %colorbar;
            set(gca,'YTick',[])
            set(gca,'XTick',[])
            xticks(1:size(concat_respmap,2));
            xticklabels(horzcat(new_cnames,spekeywords));
            yticks(1:size(concat_respmap,1));
            yticklabels(new_rnames);
            set(gca,'TickLabelInterpreter','none');
            set(gca,'TickLength',[0.001, 0.001]);
            xtickangle(45)
            savefig(bfig,'Binary_Resp_Map.fig');