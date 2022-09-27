classdef Display
    %receive input from Class Neuron
    
    properties
        info
    end
    
    methods
        function d = Display(info)
            %DISPLAY Construct an instance of this class
            %   Detailed explanation goes here
            d.info = info;
        end
        
        function showselect(d,ids,rangeratio)
            
            selected = d.info(ids);
            selected = Display.f0(selected);
            
            figure('Color','w');
            len = length(selected);
            ax = tight_subplot(2*len, 1, 0.015, 0.02, 0);
            
            for k = 1: len
                %subplot(2*len,1, 2*(k-1)+ 1);
                axes(ax(2*(k-1)+ 1));
                Draw.spec(selected(k).f0plty,selected(k).fs);
                xlabel('')
                ylabel('')
                set(gca,'TickLength',[0 .01])
                set(gca,'Yticklabel',[])
                set(gca,'Xticklabel',[])
                signallen = length(selected(k).f0plty)/selected(k).fs;% length of y signal
                xlim([rangeratio(1)*signallen,rangeratio(2)*signallen]);
                axis off
                %subplot(2*len,1, 2*k);
                axes(ax(2*k));
                Draw.raster(selected(k).f0pltsptimes,selected(k).f0plty,selected(k).fs,4.2);
                %xlim([0,length(selected(k).f0plty)/selected(k).fs]);
                xlim([rangeratio(1)*signallen,rangeratio(2)*signallen]);
                ylabel('')
                set(gca,'TickLength',[0 .01])
                if k~= len
                    set(gca,'Yticklabel',[])
                    set(gca,'Xticklabel',[])
                else
                    xlabel(selected(1).stimuliname);
                end
                axis off
                
            end
            
            
        end
        
        function showrepla(d,keyword)
             if exist('keyword','var')
                repla = d.findrepla(keyword);
            else
                repla = d.findrepla();
            end
            
            repla = Display.f0(repla);
            %deg = Display.descend(deg);
            
            figure;
            len = length(repla);
            ax = tight_subplot(2*len, 1, 0.002, 0.02, 0);
            for k = 1: len
               % subplot(2*len,1, 2*(k-1)+ 1);
                axes(ax(2*(k-1)+ 1));
                Draw.spec(repla(k).f0plty,repla(k).fs);
                xlabel('')
                ylabel('')
                
                set(gca,'TickLength',[0 .01])
                set(gca,'Yticklabel',[])
                set(gca,'Xticklabel',[])
                
               % subplot(2*len,1, 2*k);
                axes(ax(2*k));
                Draw.raster(repla(k).f0pltsptimes,repla(k).f0plty,repla(k).fs);
                ylabel('')
                set(gca,'TickLength',[0 .01])
                if k~= len
                    set(gca,'Yticklabel',[])
                    set(gca,'Xticklabel',[])
                else
                    xlabel(repla(k).stimuliname);
                end
                
            end
        end
        
        function showdeg(d,keyword)
            if exist('keyword','var')
                deg = finddeg(d,keyword);
            else 
                deg = finddeg(d);
            end
           
            deg = Display.f0(deg);
            deg = Display.descend(deg);
            
            figure; 
         
            len = length(deg);
            ax = tight_subplot(2*len, 1, 0.002, 0.02, 0);
            
            for k = 1: len
                %subplot(2*len,1, 2*(k-1)+ 1);
                axes(ax(2*(k-1)+ 1));
                Draw.spec(deg(k).f0plty,deg(k).fs);
                xlabel('')
                ylabel('')
                set(gca,'TickLength',[0 .01])
                set(gca,'Yticklabel',[])
                set(gca,'Xticklabel',[])
              
                %subplot(2*len,1, 2*k);
                axes(ax(2*k));
                Draw.raster(deg(k).f0pltsptimes,deg(k).f0plty,deg(k).fs);
                xlim([0,length(deg(k).f0plty)/deg(k).fs]);
                ylabel('')
                set(gca,'TickLength',[0 .01])
                if k~= len
                    set(gca,'Yticklabel',[])
                    set(gca,'Xticklabel',[])
                else
                    xlabel(deg(1).stimuliname);
                end
                
            end
            
           
            
            
        end
        
        function showincre(d,keyword)
            
            if exist('keyword','var')
                incre = findincre(d,keyword);
            else
                incre = findincre(d);
            end
            
            incre = Display.f0(incre);
            incre = Display.descend(incre);
            
            figure;
            len = length(incre);
            ax = tight_subplot(2*len, 1, 0.002, 0.02, 0);
            for k = 1: len
                %subplot(2*len,1, 2*(k-1)+ 1);
                axes(ax(2*(k-1)+ 1));
                Draw.spec(incre(k).f0y,incre(k).fs);
                xlabel('')
                ylabel('')
                set(gca,'TickLength',[0 .01])
                set(gca,'Yticklabel',[])
                set(gca,'Xticklabel',[])
                
                %subplot(2*len,1, 2*k);
                axes(ax(2*k));
                Draw.raster(incre(k).f0sptimes,incre(k).f0y,incre(k).fs);
                ylabel('')
                set(gca,'TickLength',[0 .01])
                if k~= len
                    set(gca,'Yticklabel',[])
                    set(gca,'Xticklabel',[])
                else
                    xlabel(incre(1).stimuliname);
                end
            end
            
            
            
        end
 
        function showcatego(d,keyword)
            
            if exist('keyword','var')
                catego = d.findcatego(keyword);
            else
                catego = d.findcatego();
            end
            
            catego = Display.f0(catego);
            %deg = Display.descend(deg);
            
            figure;
            len = length(catego);
            ax = tight_subplot(2*len, 1, 0.002, 0.02, 0);
            for k = 1: len
               % subplot(2*len,1, 2*(k-1)+ 1);
                axes(ax(2*(k-1)+ 1));
                Draw.spec(catego(k).f0plty,catego(k).fs);
                xlabel('')
                ylabel('')
                
                set(gca,'TickLength',[0 .01])
                set(gca,'Yticklabel',[])
                set(gca,'Xticklabel',[])
                
               % subplot(2*len,1, 2*k);
                axes(ax(2*k));
                Draw.raster(catego(k).f0pltsptimes,catego(k).f0plty,catego(k).fs);
                ylabel('')
                set(gca,'TickLength',[0 .01])
                if k~= len
                    set(gca,'Yticklabel',[])
                    set(gca,'Xticklabel',[])
                else
                    xlabel(catego(k).stimuliname);
                end
                
            end
        end
        
        function deg = finddeg(d,keyword)
            
            idx = find(~cellfun(@isempty,regexp(cellstr({d.info.stimuliname}.'),'deg|Deg')));
            
            deg = d.info(idx);
            
            if exist('keyword','var')
                kwidx = find(~cellfun(@isempty,regexp(cellstr({deg.stimuliname}.'),keyword)));
                deg = deg(kwidx);
            end
            
            
        end
        
        function incre = findincre(d,keyword)
            
            idx = find(~cellfun(@isempty,regexp({d.info.stimuliname}.','incre')));
            
            incre = d.info(idx);
            
            if exist('keyword','var')
                kwidx = find(~cellfun(@isempty,regexp({incre.stimuliname}.',keyword)));
                incre = incre(kwidx);
            end
            
        end
        
        function frag = findindi(d,keyword)
            
            idx = find(~cellfun(@isempty,regexp({d.info.stimuliname}.','syl')));
            
            frag = d.info(idx);
            
            if exist('keyword','var')
                kwidx = find(~cellfun(@isempty,regexp({frag.stimuliname}.',keyword)));
                frag = frag(kwidx);
            end
            
        end
        
        function catego = findcatego(d,keyword)
            idx = find(~cellfun(@isempty,regexp(cellstr({d.info.stimuliname}.'),'catego|repla')));
            catego = d.info(idx);
            if exist('keyword','var')
                
                kwidx = find(~cellfun(@isempty,regexp({catego.stimuliname}.',keyword)));
                catego = catego(kwidx);
            end
            
        end
        
        function repla = findrepla(d,keyword)
            idx = find(~cellfun(@isempty,regexp(cellstr({d.info.stimuliname}.'),'Repla|catego|repla')));
            repla = d.info(idx);
            if exist('keyword','var')
                
                kwidx = find(~cellfun(@isempty,regexp({repla.stimuliname}.',keyword)));
                repla = repla(kwidx);
            end
            
        end
        
        function trans = findtrans(d,keyword)
            idx = find(~cellfun(@isempty,regexp({d.info.stimuliname}.','trans')));
            trans = d.info(idx);
            if exist('keyword','var')
                
                kwidx = find(~cellfun(@isempty,regexp({trans.stimuliname}.',keyword)));
                trans = trans(kwidx);
            end
        end
        
        function showtrans(d,keyword)
            if exist('keyword','var')
                trans = d.findtrans(keyword);
            else
                trans = d.findtrans();
            end
            
            trans = Display.f0(trans);
            %deg = Display.descend(deg);
            
            figure;
            len = length(trans);
            ax = tight_subplot(2*len, 1, 0.002, 0.02, 0);
            for k = 1: len
                % subplot(2*len,1, 2*(k-1)+ 1);
                axes(ax(2*(k-1)+ 1));
                Draw.spec(trans(k).f0plty,trans(k).fs);
                xlabel('')
                ylabel('')
                
                set(gca,'TickLength',[0 .01])
                set(gca,'Yticklabel',[])
                set(gca,'Xticklabel',[])
                
                % subplot(2*len,1, 2*k);
                axes(ax(2*k));
                Draw.raster(trans(k).f0pltsptimes,trans(k).f0plty,trans(k).fs);
                ylabel('')
                set(gca,'TickLength',[0 .01])
                if k~= len
                    set(gca,'Yticklabel',[])
                    set(gca,'Xticklabel',[])
                else
                    xlabel(trans(k).stimuliname);
                end
                
            end
            
        end
        
        function norm = findnorm(d,keyword)
            idx = find(~cellfun(@isempty,regexp(cellstr({d.info.stimuliname}.'),'norm')));
            norm = d.info(idx);
            if exist('keyword','var')
                
                kwidx = find(~cellfun(@isempty,regexp(cellstr({norm.stimuliname}.'),keyword)));
                norm = norm(kwidx);
            end
            
        end
        
        function shownorm(d,keyword)
            
            norm = findnorm(d,keyword);  
            norm = Display.f0(norm);
            figure;
            len = length(norm);
            ax = tight_subplot(2*len, 1, 0.002, 0.02, 0);
            for k = 1: len
                %subplot(2*len,1, 2*(k-1)+ 1);
                axes(ax(2*(k-1)+ 1));
                Draw.spec(norm(k).f0plty,norm(k).fs);
                set(gca,'XTickLabel',[]);
                set(gca,'YTickLabel',[]);
                xlabel('')
                ylabel('')
                axis off
                %subplot(2*len,1, 2*k);
                axes(ax(2*k));
                Draw.raster(norm(k).f0pltsptimes,norm(k).f0plty,norm(k).fs);
                if k~= len
                    set(gca,'Yticklabel',[])
                    set(gca,'Xticklabel',[])
                else
                    xlabel(norm(k).stimuliname);
                end
                ylabel('')
                axis off
            end
        end
        
        function showfrag(d,keyword,mergedeleinf,rangeratio,ids)
            frag = alignsyl(d,keyword,mergedeleinf);
            frag = frag(ids);
            
            figure;
            len = length(frag);
            ax = tight_subplot(2*len, 1, 0.005, 0.02, 0);
            for k = 1: len
                %subplot(2*len,1, 2*(k-1)+ 1);
                axes(ax(2*(k-1)+ 1));
                Draw.spec(frag(k).f0plty,frag(k).fs);
                set(gca,'TickLength',[0 .01])
                set(gca,'XTickLabel',[]);
                set(gca,'YTickLabel',[]);
                xlabel('')
                ylabel('')
                axis off
                
                signallen = length(frag(k).f0plty)/frag(k).fs;
                if exist('rangeratio','var')
                    xlim([rangeratio(1)*signallen,rangeratio(2)*signallen]);
                end
                
                %subplot(2*len,1, 2*k);
                axes(ax(2*k));
                Draw.raster(frag(k).f0pltsptimes,frag(k).f0plty,frag(k).fs,4.2);
                
                if k~= len
                    set(gca,'TickLength',[0 .01])
                    set(gca,'Yticklabel',[])
                    set(gca,'Xticklabel',[])
                else
                    xlabel(frag(k).stimuliname);
                end
                
                if exist('rangeratio','var')
                    xlim([rangeratio(1)*signallen,rangeratio(2)*signallen]);
                end
                
                ylabel('')
                axis off
            end
            
            
        end
          
        function showsyl(d,keyword,mergedeleinf,rangeratio,ids)
            syl = alignsyl(d,keyword,mergedeleinf);
            syl = syl(ids);
            
            figure;
            len = length(syl);
            ax = tight_subplot(2*len, 1, 0.005, 0.02, 0);
            for k = 1: len
                %subplot(2*len,1, 2*(k-1)+ 1);
                axes(ax(2*(k-1)+ 1));
                Draw.spec(syl(k).f0plty,syl(k).fs);
                set(gca,'TickLength',[0 .01])
                set(gca,'XTickLabel',[]);
                set(gca,'YTickLabel',[]);
                xlabel('')
                ylabel('')
                axis off
                
                signallen = length(syl(k).f0plty)/syl(k).fs;
                if exist('rangeratio','var')
                    xlim([rangeratio(1)*signallen,rangeratio(2)*signallen]);
                end
                
                %subplot(2*len,1, 2*k);
                axes(ax(2*k));
                Draw.raster(syl(k).f0pltsptimes,syl(k).f0plty,syl(k).fs,4.2);
                
                if k~= len
                    set(gca,'TickLength',[0 .01])
                    set(gca,'Yticklabel',[])
                    set(gca,'Xticklabel',[])
                else
                    xlabel(syl(k).stimuliname);
                end
                
                if exist('rangeratio','var')
                    xlim([rangeratio(1)*signallen,rangeratio(2)*signallen]);
                end
                
                ylabel('')
                axis off
            end
            
            
        end
        
        function syls = alignsyl(d,keyword,mergedeleinf)  % the todisplay list should be a merged list between original full song and indepdenent single syllables
            
            if exist('keyword','var')
                kwidx = find(~cellfun(@isempty,regexp([d.info.stimuliname].',keyword)));
                pre = d.info(kwidx); % pre-screen
            else
                pre = d.info;
            end
            
            
            idxsyl = find(~cellfun(@isempty,regexp([pre.stimuliname].','syl')));
            
            idxnorm = find(~cellfun(@isempty,regexp([pre.stimuliname].','norm')));
            
            syls = pre(idxsyl);
            
            % restruct syls
            for a = 1: length(syls)
                parts = strsplit(syls(a).stimuliname,'-');
                purenumidx = find(cellfun(@isempty,regexp(parts,'\D')));
                syls(a).birdid = parts{2};
                syls(a).fragid = str2num(parts{purenumidx});
                
                % within mergedeleinf, find out the corresponding index
                budda = find(~cellfun(@isempty,(regexp(cellstr({mergedeleinf.songname}.'),syls(a).birdid))));
                confucius = find([mergedeleinf.fragid].'==syls(a).fragid);
                
                tao = intersect(budda,confucius);
                firstfrag = intersect(budda,find([mergedeleinf.fragid].'== 1));
                
                basal =  mergedeleinf(firstfrag).initial;
                syls(a).initial = mergedeleinf(tao).initial - basal;
                syls(a).terminal = mergedeleinf(tao).terminal - basal;
                syls(a).pregap =  mergedeleinf(tao).pregap;
                
                
            end
            
            % very bad code
            normrank = length(syls) + 1;
            syls(normrank).stimuliname  = pre(idxnorm).stimuliname;
            syls(normrank).fs  = pre(idxnorm).fs;
            syls(normrank).plxname  = pre(idxnorm).plxname;
            syls(normrank).unitname  = pre(idxnorm).unitname ;
            syls(normrank).y  = pre(idxnorm).y;
            syls(normrank).rawy  = pre(idxnorm).rawy;
            syls(normrank).plty  = pre(idxnorm).plty ;
            syls(normrank).sptimes  = pre(idxnorm).sptimes;
            syls(normrank).rawsptimes = pre(idxnorm).rawsptimes;
            syls(normrank).pltsptimes  = pre(idxnorm).pltsptimes ;
            syls(normrank).zpt  = pre(idxnorm).zpt ;
            syls(normrank).pltext  = pre(idxnorm).pltext;
            syls(normrank).birdid  = pre(idxnorm).stimuliname;
            syls(normrank).fragid  = 0;
            syls(normrank).initial  = 0;
            syls(normrank).terminal  = length(syls(normrank).y) / syls(normrank).fs;
            syls(normrank).pregap  = nan ;
            
            syls = table2struct(sortrows(struct2table(syls),'fragid'));
            
            for m = 1: length(syls)
                for k = 1: length(syls)
                    
                    syls(k).leny = length(syls(k).y);
                    syls(k).lenrawy = length(syls(k).rawy);
                    syls(k).lenplty = length(syls(k).plty);
                    
                    
                    syls(k).f0y = [zeros([int32(syls(k).initial),1]);syls(k).y;zeros([int32(syls(1).leny - syls(k).leny - syls(k).initial),1]) ];
                    syls(k).f0plty = [zeros([int32(syls(k).initial-syls(k).pltext),1]);syls(k).y;zeros([int32(syls(1).leny - syls(k).leny - syls(k).initial),1]) ];
                    %????????????????????????????
                    % f0 means front-pad-zero
                    %syls(k).f0rawy = [zeros([syls(1).leny - syls(k).lenrawy,1]);syls(k).rawy: zeros([syls(1).leny - syls(k).y - syls(k).initial*syls(k).fs,1])];
                    % syls(k).f0plty = [zeros([maxy - syls(k).leny,1]);syls(k).plty];
                    
                    diffy = syls(k).initial/syls(k).fs;
                    
                    diffrawy = diffy + syls(k).initial/syls(k).fs; %syls(1).zpt; % why use 1 insteand of k, because zpt caluclation is 100% accurate
                    diffplty = diffy - syls(k).pltext;
                    for m = 1: length(syls(k).sptimes)
                        syls(k).f0sptimes{m} = syls(k).sptimes{m} + diffy;
                        syls(k).f0rawsptimes{m} =  syls(k).rawsptimes{m} + diffrawy;
                        syls(k).f0pltsptimes{m} = syls(k).pltsptimes{m} + diffplty;  % bad code...
                    end
                end
                
            end
        end
        
        function showsimplesyl(d,mergedeleinf,keyword)
            
            idx = find(~cellfun(@isempty,regexp({d.info.stimuliname}.','syl')));
            
            syl = d.info(idx);
            
            if exist('keyword','var')
                kwidx = find(~cellfun(@isempty,regexp({syl.stimuliname}.',keyword)));
                syl = syl(kwidx);
            end
            
            
            dir = 'syllables_SVG'
            mkdir(dir);
            % get img
            for k = 1: length(syl)
                
                h = figure('Position',[522 150 368 420],'Color','white');
                
                ax = tight_subplot(2, 1, 0.002, 0.03, 0);
                %ax = tight_subplot(1, 1, 0.002, 0.03, 0);
                
                axes(ax(1));
                axis off
                Draw.spec(syl(k).plty,syl(k).fs);
                set(gca,'TickLength',[0 .01])
                set(gca,'XTickLabel',[]);
                set(gca,'YTickLabel',[]);
                xlabel('')
                ylabel('')
                axis off
                
                saveas(gca,sprintf('%s/%s.svg',dir,syl(k).stimuliname));
                
                % tight_subplot(2*len,1, 2*k);
                axes(ax(2));
                
                Draw.rasterBeta(syl(k).pltsptimes,syl(k).plty,syl(k).fs,4.5);  %%% tick label !!!!!!!! need adding
                
                ylabel('')
                set(gca,'TickLength',[0 .01])
                set(gca,'Yticklabel',[])
                set(gca,'Xticklabel',[])
                
                xlabel(syl(k).stimuliname);
                axis off
                
                
                frame = getframe(gcf);
                img{k} = frame.cdata;
                close(h);
            end
            
            % create figure
            
            column = 20;
            row = ceil(length(syl)/column);
            
            rest = column*row - length(img);
            white = uint8(255*ones(size(img{1})));
            
            if rest > 0
                for k = 1:rest
                    img = [img,white];
                end
            end
            
            reshapedI = reshape(img, column,[])';
            IMG = cell2mat(reshapedI);
            
            figure;
            imshow(IMG);
            nowtime = datestr(datetime('now'),'yyyy-mmm-dd-HH-MM-SS');
            channelname = syl(1).channelname;
            plxname = syl(1).plxname;
            unitname = syl(1).unitname;
            imwrite(IMG,sprintf('%s-%s-%u-%s.png',plxname,channelname,unitname,nowtime))
            
            %             figure;
            %             bx = tight_subplot(column,row,0,0,0);
            %
            %             for m = 1:length(img)
            %                 axes(bx(m));
            %                 imshow(img{m});
            %             end
        end
        
        function target = findtarget(d,keyword,mergedeleinf)
            syl = alignsyl(d,keyword,mergedeleinf);
            deg = finddeg(d,keyword);
            incre = findincre(d,keyword);
            
            
        end
        
        function showtarget(d,keyword)
        end
        
    end
    
    methods(Static)
        
        function info = f0(info)
            
            leny = cellfun(@length,{info.y}.');
            lenrawy = cellfun(@length,{info.rawy}.');
            lenplty = cellfun(@length,{info.plty}.');
            
            maxy = max(leny);
            maxrawy = max(lenrawy);
            maxplty = max(lenplty);
            
            for k = 1: length(info)
                
                info(k).leny = leny(k);
                info(k).lenrawy = lenrawy(k);
                info(k).lenplty = lenplty(k);
                
                info(k).f0y = [zeros([maxy - info(k).leny,1]);info(k).y]; % f0 means front-pad-zero
                info(k).f0rawy = [zeros([maxy - info(k).lenrawy,1]);info(k).rawy];
                info(k).f0plty = [zeros([maxy - info(k).leny,1]);info(k).plty];  
                
                diffy = (maxy - info(k).leny)/info(k).fs;
                diffrawy = (maxrawy - info(k).lenrawy)/info(k).fs;
                diffplty = (maxplty - info(k).lenplty)/info(k).fs;
                
                for m = 1: length(info(k).sptimes)
                    info(k).f0sptimes{m} = info(k).sptimes{m} + diffy;
                    info(k).f0rawsptimes{m} = info(k).rawsptimes{m} + diffrawy;
                    info(k).f0pltsptimes{m} = info(k).pltsptimes{m} + diffplty;  % bad code...
                end
            end
         
        end  % f0 means front pad zero
        
        function info = descend(info)
            if length(info) == 1
                t = struct2table(info,'AsArray',1);
            else
                t = struct2table(info);
            end
            t = sortrows(t,'leny','descend');
            info = table2struct(t);
        end
        
    end
end

