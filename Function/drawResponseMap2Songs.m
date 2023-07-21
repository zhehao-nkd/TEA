

        function [fig,persong] = drawResponseMap2Songs(neuronfiles)
            % 为了计算Neuron对每一个song的反应的大小
            
            dbstop if error

%             conkeywords = {'B346','B512','B521','B554','B606','G429','G506','G518','G573',...
%                 'G578','O331','O507','O509','O540','Y515','Y606','Y616','G699'}; % 如果查看其他song,那就在此处添加!!!
           num_files = length(neuronfiles);
            wb = PoolWaitbar(num_files,'loading...');


            responseinfo = struct;

            summer = {};
            %    findit = @(x) find();

            perneuron = struct;
            for k = 1:length(neuronfiles) % should be par-for
                %                 try
                loaded = load(neuronfiles(k).filepath,'list','info');
                list18 = loaded.list;%song.normlist(~cellfun(@isempty, regexp(cellstr({loaded.song.normlist.stimuliname}.'),strjoin(conkeywords,'|'))));
                list18  = Song.judgeConResp(list18);
                perneuron(k).nsbs = neuronfiles(k).nsbs;
                for kk = 1:length(list18)
                    list18(kk).neuronname = loaded.info.formated_name;
                    list18(kk).nsbs = neuronfiles(k).nsbs;
          
                end
                summer{k} = list18;
                perneuron(k).numpositive = length(find([list18.label].' == 1));
                perneuron(k).neuronname = loaded.info.formated_name;
              

                %perneuron(k).respsong18 = neuronfiles(k).RepSong18;
                perneuron(k).resp2any = sum([list18.label].');
                increment(wb);
                %
            end

            [~,index] = sortrows([perneuron.numpositive].'); perneuron = perneuron(index(end:-1:1)); clear index

            columnnames = cellstr({perneuron.neuronname}.');

            respinfo = horzcat(summer{:});
            for k = 1:length(respinfo)
                if respinfo(k).nsbs == 1 && respinfo(k).label == 1
                    respinfo(k).colorstate = 1;
                elseif respinfo(k).nsbs == 2 && respinfo(k).label == 1
                    respinfo(k).colorstate = 2;
                else
                    respinfo(k).colorstate = 0;
                end
            end

            respinfoT = struct2table(respinfo);


            % calculate the stimuli information
            unique_songname = unique(cellstr({respinfo.stimuliname}.'));

            persong = struct;
            for k = 1:length(unique_songname)

                hitted = find(~cellfun(@isempty, regexp(cellstr({respinfo.stimuliname}.'),unique_songname{k})));
                sublist = respinfo(hitted);
                persong(k).songname = unique_songname{k};
                persong(k).numpositive =  length(find([sublist.label].' == 1));

            end


            perneuron = table2struct(sortrows(struct2table(perneuron),{'nsbs','numpositive'},{'ascend','descend'}));

            selected_perneuron = perneuron(find([perneuron.resp2any].' ~= 0));
            columnnames = cellstr({selected_perneuron.neuronname}.');

            persong = table2struct(sortrows(struct2table(persong),{'numpositive'},{'descend'}));
            
            rownames = cellstr({persong.songname}.');



            fig = figure;
            title(neuronfiles(1).neuronname);
            h = heatmap(respinfoT,'neuronname','stimuliname','ColorVariable','colorstate','CellLabelColor','none');

            % Define the custom colormap
            customColormap = [0.99 0.99 0.99; 0.1 0.1 1; 1 0.1 0.1]; % White, Blue, Red

            % Assign the custom colormap to the heatmap object
            h.Colormap = customColormap;

            h.YDisplayData = rownames;
            h.XDisplayData = columnnames;

            % title('Diff')
            %h = heatmap(T_sum_perstimulus,'fragname','neuronname','ColorVariable','diff','CellLabelColor','none')
            h.YDisplayData = rownames;
            h.XDisplayData = columnnames;
            textObjs = findobj(h,'Type','Text');
            set(textObjs,'Interpreter','none');
          
           



        end

