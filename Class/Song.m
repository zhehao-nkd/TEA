classdef Song < handle
    %@ 类说明 Member class of Neuron 专门process Songs
    
    properties
        normlist
        formated_name
        list18
        info
        num_resp_to_18
        name_resp_to_18
    end
    
    methods
        
        function sg = Song(input_list,formated_name)

            sg.formated_name = formated_name;

            % 首先把list根据Fid分裂开来，然后再对每一个分裂的list算num of
            % norms的个数，取放了最多的stimuli的分裂list作为normlist
            %如果有重复，那么就接受有多个normlist？ 或者产生一个替补normlist，因为之所以会有重复原因就在于替补
            %此外还需要建立一个名为sibling的class
            separated ={};
            unique_fids = unique({input_list.Fid}.');
            for k = 1:length(unique_fids) %有多少file，对应多少experiment
                temp = input_list(find(~cellfun(@isempty, regexp(cellstr({input_list.Fid}.'),unique_fids{k}))));
                separated{k}= temp(find(~cellfun(@isempty, regexp(cellstr({temp.stimuliname}.'),'norm'))));
            end

     
            [~,maxloc] = max(cellfun(@length,separated));
            sg.normlist = separated{maxloc};

            sg.normlist = Song.judgeConResp(sg.normlist);


%             sg.normlist = input_list(find(~cellfun(@isempty, regexp(cellstr({input_list.stimuliname}.'),'norm'))));
            sg.labelSongs;
            sg.getInfo;
            
        end

        function sg = getInfo(sg)

            respond_18_ids = find([sg.normlist.label].' ==1);
            sg.num_resp_to_18 = length(respond_18_ids);
            if ~(sg.num_resp_to_18 == 0)
                sg.name_resp_to_18 = cellfun(@Convert.bid,cellstr({sg.normlist(respond_18_ids).stimuliname}.'),'Uni',0);
            else
                sg.name_resp_to_18 =[];
            end

        end

        function sg = labelSongs(sg)
            %由于不同neuron的stimuliname和tutbos不一样，这个方法是用来标记出每个stimuli是cons，het还是其他

            bd = Bird; %使用Bird类判断father song
            bosname = Convert.bid(sg.formated_name);
            fathername =  bd.findFather(bosname);


            names18 = {'B346','B512','B521','B554','B606','G429','G506','G518','G548','G573',...
                'G578','O331','O507','O509','O540','Y515','Y606','Y616','G699'};

            for k = 1:length(sg.normlist)

                 if ~isempty(regexp(sg.normlist(k).stimuliname,strjoin(names18,'|'))) 
                     sg.normlist(k).type = 1;
                     sg.normlist(k).name = regexp(convertCharsToStrings(sg.normlist(k).stimuliname),strjoin(names18,'|'),'match');
                 elseif ~isempty(regexp(sg.normlist(k).stimuliname,'WNS|Fcall|Het|Mcall'))
                     sg.normlist(k).type = 3;
                     sg.normlist(k).name = regexp(convertCharsToStrings(sg.normlist(k).stimuliname),'WNS|Fcall|Het|Mcall','match');
                 elseif ~isempty(regexp(sg.normlist(k).stimuliname,bosname))
                     sg.normlist(k).type = 2;
                     sg.normlist(k).name = "BOS";
                 elseif ~isempty(regexp(sg.normlist(k).stimuliname,fathername))
                      sg.normlist(k).type = 2;
                     sg.normlist(k).name = "TUT";         
                 end

            end



%             findit = @(x) find(~cellfun(@isempty, regexp(cellstr({sg.normlist.stimuliname}.'),x)));

          

        end
        
        function drawThree(sg,outdir)

            I = {};
            for idx = 1: length(sg.normlist)
                figure('Color','w','Position', [216 590 1274 447]);
                Draw.three(sg.normlist(idx).plty,sg.normlist(idx).fs,sg.normlist(idx).pltsptimes);
                xlabel(sprintf('%s-%s',sg.normlist(idx).stimuliname,sg.normlist(idx).Fid));
                frame = getframe(gcf);
                I{idx} = frame.cdata;
                saveas(gcf,sprintf('%s\\Song_%s_%s.fig',outdir,sg.formated_name,sg.normlist(idx).stimuliname));
                saveas(gcf,sprintf('%s\\Song_%s_%s.eps',outdir,sg.formated_name,sg.normlist(idx).stimuliname),'epsc');
                saveas(gcf,sprintf('%s\\Song_%s_%s.svg',outdir,sg.formated_name,sg.normlist(idx).stimuliname),'svg');
                close(gcf);
            end

%             figure('Color','w','Position', [1933 672 673 497]);
%             neu.drawFirstWaveform;     % draw waveform
%             frame = getframe(gcf);
%             I{length(I)+ 1} = frame.cdata;
%             close(gcf);

            % draw blank white
            lieshu = 3;
            hangshu = ceil(length(I)/lieshu);
            rest = lieshu*hangshu - length(I);
            white = uint8(255*ones(size(I{1})));

            if rest > 0
                for k = 1:rest
                    I = [I,white];
                end
            end

            reshapedI = reshape(I, lieshu,[])';
            clear I
            img = cell2mat(reshapedI);
            imwrite(img,sprintf('%s\\临安府NormThree_%s.png',outdir,sg.formated_name));


        end

        function img = saveDrawSortedRespToCons(neu) %做三分析图，把使反应的neuron标记为黄色，不反应的不标记
            %非常难以找到这个function

            neu.judgeConResp;
            %             finalfig = getframe(gcf).cdata;
            %             close(gcf)
            ids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'norm')));
            normlist = neu.list(ids);
            if isempty(normlist)
                return
            end
            %sorted_fraglist = table2struct(sortrows( struct2table(fraglist) ,'maxvalue','descend'));
            sorted_normlist = table2struct(sortrows( struct2table(normlist) ,'maxsdf','descend'));


            I = {}; % collection of frag-response-three images
            for k = 1: length(sorted_normlist)

                if sorted_normlist(k).label == 0
                    h =  figure('Position',[1606 287 1343 660],'Color','w');
                elseif sorted_normlist(k).label == 1
                    h =  figure('Position',[1606 287 1343 660],'Color','y');
                end

                %h.WindowState = 'maximized';
                Draw.two(sorted_normlist(k).plty,sorted_normlist(k).fs,sorted_normlist(k).pltsptimes);
                xlabel(sprintf('%s-maxsdf: %f',sorted_normlist(k).stimuliname,sorted_normlist(k).maxsdf));
                temp = getframe(gcf);
                I{k} = temp.cdata;


                close(h)
            end

            %             I{k+1} = finalfig; % performance figure

            %             figure;
            %             n.draw_waveform;     % draw waveform
            %             frame = getframe(gcf);
            %             I{length(I)+ 1} = frame.cdata;
            %             close(gcf);

            % draw blank white
            lieshu = 10;
            hangshu = ceil(length(I)/lieshu);
            rest = lieshu*hangshu - length(I);
            white = uint8(255*ones(size(I{1})));

            if rest > 0
                for k = 1:rest
                    I = [I,white];
                end
            end

            reshapedI = reshape(I, lieshu,[])';
            clear I;
            img = cell2mat(reshapedI);
            imwrite(img,sprintf('秦RespToSongs_%s.png',neu.info.formated_name));

        end

    
    end
    
    methods(Static)



        function outputlist = judgeConResp(inputlist,mode)

            if exist("mode",'var')  % 原来的备选方法


                % ids = find(~cellfun(@isempty, regexp(cellstr({inputlist.stimuliname}.'),'norm|deg|repla'))); % find all norms
                % ’syl'可以兼容旧的stimuli命名规则

                for k = 1: length(inputlist)
                    %  thisi = ids(n);

                    presdf = Cal.sdf(inputlist(k).prejudgerespsptimes,zeros(length(inputlist(k).judgerespy),1),inputlist(k).fs,0.001,0.02);
                    sdf = Cal.sdf(inputlist(k).judgerespsptimes,inputlist(k).judgerespy,inputlist(k).fs,0.001,0.02); % 0.001,0.004
                    %                 sdf = Cal.sdf(neu.list(thisi).judgerespsptimes,neu.list(thisi).judgerespy,neu.list(thisi).fs,0.001,0.02);
                    [maxpresdf,~] = max(presdf);
                    [maxsdf,maxidx] = max(sdf);

                    pre_frs = Cal.eachTrialFiringRate(inputlist(k).prejudgerespsptimes,length(inputlist(k).judgerespy)/inputlist(k).fs);
                    sti_frs = Cal.eachTrialFiringRate(inputlist(k).judgerespsptimes,length(inputlist(k).judgerespy)/inputlist(k).fs);

                    [h,p] = ttest(sti_frs,pre_frs,'Tail','Right','Alpha',0.05);
                    inputlist(k).pvalue = p;
                    inputlist(k).label = 0; % 初始化
                    inputlist(k).sti_frs = mean(sti_frs);
                    inputlist(k).pre_frs = mean(pre_frs);
                    if h == 1
                        inputlist(k).label = 1;
                        %                     if (num_of_not_empty_trials/length(neu.list(thisi).pltsptimes)<0.5)||(num_of_not_empty_trials<5) % if not 60% of the trails are not empty
                        %                         neu.list(thisi).label = 0;
                        %                     end
                    elseif maxsdf > 17 && maxsdf > maxpresdf % neu rescue
                        inputlist(k).label = 1;

                    end
                end

                return


            end
            %             ids = find(~cellfun(@isempty, regexp(cellstr({neu.list.stimuliname}.'),'norm|deg'))); % find all norms
            %             % ’syl'可以兼容旧的stimuli命名规则

            for k = 1: length(inputlist)
                % thisi = ids(n);
              %  presptimes = Extract.sptimes

                presdf = Cal.sdf(inputlist(k).presptimes,zeros(length(inputlist(k).y),1),inputlist(k).fs,0.001,0.02);
                sdf = Cal.sdf(inputlist(k).sptimes,inputlist(k).y,inputlist(k).fs,0.001,0.02); % 0.001,0.004
                %figure; plot(sdf);
                % figure; Draw.three(neu.list(thisi).judgerespy,neu.list(thisi).fs,neu.list(thisi).judgerespsptimes);
                % figure; Draw.three(neu.list(thisi).plty,neu.list(thisi).fs,neu.list(thisi).pltsptimes);
                [maxpresdf,~] = max(presdf);
                [maxsdf,maxidx] = max(sdf);
                percentage_max = maxidx/length(sdf);
                time_max = length(inputlist(k).y)/inputlist(k).fs*percentage_max;
                % check whether the surroding are has spikes in most of the
                % trials
                extracted_sptimes = Extract.sptimes(inputlist(k).sptimes,time_max - 0.15, time_max + 0.15); % 前后 100ms
                num_of_not_empty_trials = length(find(~cellfun(@isempty, extracted_sptimes)));

                % ttest2 to test whether sti_frs are significantly higher
                % than pre_frs
                pre_frs = Cal.eachTrialFiringRate(inputlist(k).presptimes,length(inputlist(k).y)/inputlist(k).fs);
                sti_frs = Cal.eachTrialFiringRate(inputlist(k).sptimes,length(inputlist(k).y)/inputlist(k).fs);

                [p,h] = signrank(sti_frs,pre_frs,'Tail','Right','Alpha',0.05);
                inputlist(k).pvalue = p;
                % mean_maxsdf = maxsdf/length(neu.list(thisi).judgerespsptimes);
                inputlist(k).pre_fr = mean(pre_frs);
                inputlist(k).sti_fr = mean(sti_frs);
                inputlist(k).rs = inputlist(k).sti_fr - inputlist(k).pre_fr;

                inputlist(k).maxsdf = maxsdf;
                inputlist(k).label = h; % 初始化

            


                if maxsdf > 17 && maxsdf > maxpresdf %如果是 time-locked response
                    inputlist(k).label = 1;
                end

                if (num_of_not_empty_trials/length(inputlist(k).pltsptimes)<0.5)||(num_of_not_empty_trials<5) % if not 60% of the trails are not empty
                    inputlist(k).label = 0;
                end



                %                 figure;  % 测试用代码
                %                 Draw.three(neu.list(thisi).plty,neu.list(thisi).fs,neu.list(thisi).pltsptimes);
                %                 title(sprintf('Label is %u',neu.list(thisi).label) );
                %                 close(gcf)

            end

            zscored_rs = [inputlist.rs].';

            for k = 1:length(inputlist)
                inputlist(k).zscored_rs = zscored_rs(k);
            end

            outputlist = inputlist;

        end


        
        function namelist = get_all_elements_involved(rawdir_path)
            % 功能说明： 从dir里包含的所用音频文件里提取参与到replacement的音节的名称，称为namelist
            % from all ZP folder, know all the syllable transitions presented as repla stimuli
            filenames = cellfun(@Utl.fileparts,Extract.filesAllLevel(rawdir_path, '*.wav'),'Uni',0);
            repla_filenames = filenames(find(~cellfun(@isempty, regexp(filenames,'Repla|repla'))));
            %上面那行可以用来获取ZPdir里面的Individual Fragments的filename，以期后续对其分类
            
            namelist = struct;
            locY = @(x,y) x{y}; % get element in location y
            for k = 1:length(repla_filenames)
                splited = split(repla_filenames{k},{'-before-','-gapis-'});
                namelist(k).bname1 = Convert.bid(splited{1});
                namelist(k).fid1 = str2num(locY(split(splited{1},'-'),4));
                namelist(k).bname2 = Convert.bid(splited{2});
                namelist(k).fid2 = str2num(locY(split(splited{2},'-'),3));
                namelist(k).concat1 = sprintf('%s-%02u',namelist(k).bname1,namelist(k).fid1);
                namelist(k).concat2 = sprintf('%s-%02u',namelist(k).bname2,namelist(k).fid2);
                namelist(k).fullrepla = sprintf('%s-%02u-%s-%02u',namelist(k).bname1,namelist(k).fid1,namelist(k).bname2,namelist(k).fid2);
            end
            
            [~,uniqueids,~] = unique({namelist.fullrepla}.'); % 去除重复的rows
            namelist = namelist(uniqueids);
            
        end
        
        
        function metadata = namelist2meta(namelist,refer_list)
            % @功能说明： convert namelist to standard metadata, based on refer_list
            unifragnames = unique(vertcat({namelist.concat1}.',{namelist.concat2}.')); % unique fragnames in ZPdirs of repla stimuli
            %reference_eleinf = Utl.load('con_and_allspe_eleinf.mat');
            
            bnames = cellfun(@Convert.bid,{refer_list.fullname}.','UniformOutput',false);
            func1 = @(x) regexp(x,'(?<=-)\d+','match'); % reference_eleinf(1).fullname
            func2 = @(x) x{1};
            fragids = cellfun(@(x) func2(func1(x)),{refer_list.fullname}.','UniformOutput',false) ;
            bname_fragid = cellfun(@(x,y) horzcat(x,'-',y), bnames,fragids,'Uni',0); % reference_eleinf
            
            metadata = struct;
            for k = 1:length(unifragnames)
                ids = find(strcmp(unifragnames{k},bname_fragid));
                metadata(k).unifragnames = unifragnames{k};
                if ~isempty(ids)
                    metadata(k).y = refer_list(ids).y;
                    metadata(k).fs= refer_list(ids).fs;
                    metadata(k).fullname = refer_list(ids).fullname;
                else
                    metadata(k).y = [];
                    metadata(k).fs= [];
                    metadata(k).fullname = [];
                end
            end
            
            metadata = metadata(~cellfun(@isempty,{metadata.fullname}.'));
            
        end
        
        
        function [transition_statics,relation_map] = get_trans_statics(sequential_song_meta)
            % @功能说明： 生成transition statics
            % remove those with catego num as 0
            idsnot0 = find([sequential_song_meta.catego].' ~= 0);
            sequential_song_meta = sequential_song_meta(idsnot0);
            disp('Remove those with catego num as 0');
            splited = MetaStimuli.split(sequential_song_meta); % 把大的struct分解成对单个song的
            collect_transfering = {};
            for k = 1: length(splited) % 对每个song计算transfering的计数
                transfering = {};
                parfor kk = 1: length(splited{k})-1
                    transfering{kk} = [splited{k}(kk).catego,splited{k}(kk+1).catego];
                end
                collect_transfering{k} = transfering;
            end
            all_transfering = horzcat(collect_transfering{:}); % Concat这些song里面所有的transfering
            
            % 集合所有存在的transition事件
            transitioncollect = struct;
            for k = 1:length(all_transfering)
                transitioncollect(k).first = all_transfering{k}(1);
                transitioncollect(k).second = all_transfering{k}(2);
            end
            
            % 对每对可能的catego transition统计对应的频数
            all_categos = unique([sequential_song_meta .catego].');
            relation_map = zeros(length(all_categos),length(all_categos)); % initialization
            transition_statics = struct;
            
            count = 0;
            for k = 1:length(all_categos)
                for kk = 1:length(all_categos)
                    count = count + 1;
                    relation_map(k,kk) = length(intersect(...
                        find([transitioncollect.first] == all_categos(k)),...
                        find([transitioncollect.second] == all_categos(kk))));
                    
                    transition_statics(count).firstsyllable = k;
                    transition_statics(count).secondsyllable = kk;
                    transition_statics(count).firstsecond = [k,kk];
                    transition_statics(count).counts = relation_map(k,kk);
                end
            end
            
             all_counts = sum(sum(relation_map));
            
            % sort the transition_statics and assign the rank
            transition_statics = table2struct(sortrows(struct2table(transition_statics),'counts','descend'));
            for k = 1:length(transition_statics)
                transition_statics(k).rank = k; % rank越小， transition越常见
                transition_statics(k).chance =  transition_statics(k).counts/all_counts;        
            end
            
        end
        
        
        function drawFragSpec(all_eleinf,mode)
            % A function to Draw spectrograms of already categorized fragments
            
            fs = 32000;
            img = {};
            
            switch mode
                case 'by categos' % by number of categories
                    categos = unique([all_eleinf.catego].'); % 计算categos的数量
                    concated_y = {};
                    for k = 1: length(categos)
                        
                        ids = find([all_eleinf.catego].' == categos(k));
                        gap = zeros(0.3*fs,1);
                        categoys = cellfun(@(x) vertcat(x,gap),{all_eleinf(ids).y}.','Uni',0);
                        concated_y{k} = vertcat(categoys{:});
                        
                    end
                    
                    max_y_length = max(cellfun(@length,concated_y));
                    
                    
                    for k = 1: length(categos)
                        
                        
                        concated_y{k} = [concated_y{k};zeros(max_y_length - length(concated_y{k}),1)];
                        
                        figure('Position',[8 310 6000 716]);
                        Draw.spec(concated_y{k},fs);
                        img{k} = getframe(gcf).cdata;
                        close(gcf);
                    end
                    
                    
                case 'by rows' % by number of rows
                    
                    [~,index] = sortrows([all_eleinf.catego].'); all_eleinf = all_eleinf(index); clear index % sort all_eleinf
                    column_num = 10; % I should use fixed syllables for replacement
                    row_num = round(length(all_eleinf)/column_num);
                    
                    for k = 1: row_num
                        
                        ids = column_num*(k-1) + 1: min(column_num*k,length(all_eleinf));
                        gap = zeros(0.3*fs,1);
                        categoys = cellfun(@(x) vertcat(x,gap),{all_eleinf(ids).y}.','Uni',0);
                        concated_y = vertcat(categoys{:});
                        
                        figure('Position',[14 829 1902 280]);
                        Draw.spec(concated_y,fs);
                        img{k} = getframe(gcf).cdata;
                        close(gcf);
                        
                    end
            end
            
            whole_img = vertcat(img{:});
            imwrite(whole_img,'New_Frags_separated_by_frags.png')
        end

        
        % from csv sort figures


        function csv2Sortfig(csv_path,figdir)
            %csvpath = "E:\My_Analysis_by_date\0925_CheckNeuralResponseToReplas\FigOTE.csv";
            
            otestruct = table2struct(readtable(csv_path))
            %figdir = "E:\My_Analysis_by_date\0925_CheckNeuralResponseToReplas\Figs_NewCon_OTE_For_TransStatastics\Figs_OTE";
            [~,figdirname,~]= fileparts(figdir)
            which_categos = unique([otestruct.catego].');
            outdir = sprintf('Sorted_%s_figs',figdirname);
            
            for k = 1:length(which_categos)
                subdirs{k} = sprintf('%s\\%u',outdir,which_categos(k));
                mkdir(subdirs{k});
            end
            
            if min(which_categos) == 0
                method = 0;
            else
                methods = 1;
            end
              
            
            for k = 1:length(otestruct)
                if method == 0
                    targetdir = subdirs{otestruct(k).catego + 1};
                    
                else
                    targetdir = subdirs{otestruct(k).catego};
                end
                source = fullfile(figdir,otestruct(k).name);
                destiney = fullfile(targetdir,otestruct(k).name);
                copyfile(source,destiney);
            end
            
        end

        % from sorted fig regenrate csv
        function figdir2csv(figdir)
            %figdir = "E:\My_Analysis_by_date\0925_CheckNeuralResponseToReplas\Sorted_ote_repla_figs";
            
            subdirs = Extract.folder(figdir);
            
            figstruct = struct;
            count = 0;
            for k = 1:length(subdirs)
                
                [~,categoname,~] = fileparts(subdirs{k});
                categonum = int16(str2double(categoname));
                files = Extract.filename(subdirs{k},'*.png');
                for kk = 1:length(files)
                    count = count + 1;
                    [~,name,ext] = fileparts(files{kk});
                    figstruct(count).name = strcat(name, ext);
                    figstruct(count).catego = categonum;
                end
                
            end
            
            figtable = struct2table(figstruct);
            [~,figdirname,~] = fileparts(figdir);
            
            writetable(figtable,sprintf('%s.csv',figdirname));
        end

        
    end






end

