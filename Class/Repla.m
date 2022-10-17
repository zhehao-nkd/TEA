classdef Repla < handle
    %@ 类说明 Member class of Neuron
    
    properties
        replalist
        list
        formated_name
    end
    
    methods
        
        function rp = Repla(list)
            % 从 Neuron 文件中提取neuron对replas的反应，已经具体的replas的名称
            rp.list = list;
            rp.judgeReplaResp;
            replalist = rp.list(find(~cellfun(@isempty, regexp(cellstr({list.stimuliname}.'),'Repla'))));
            
            Replst_fnames = [replalist.stimuliname].';
            locY = @(x,y) x{y}; % get element in location y
            for k = 1:length(replalist)
                splited = split(Replst_fnames{k},{'-before-','-gapis-'});
                replalist(k).bname1 = Convert.bid(splited{1});
                replalist(k).fid1 = str2num(locY(split(splited{1},'-'),4));
                replalist(k).bname2 = Convert.bid(splited{2});
                replalist(k).fid2 = str2num(locY(split(splited{2},'-'),3));
                replalist(k).concat1 = sprintf('%s-%02u',replalist(k).bname1,replalist(k).fid1);
                replalist(k).concat2 = sprintf('%s-%02u',replalist(k).bname2,replalist(k).fid2);
                replalist(k).fullrepla = sprintf('%s-%02u-%s-%02u',replalist(k).bname1,replalist(k).fid1,replalist(k).bname2,replalist(k).fid2);
            end
            
            rp.replalist = replalist;
            
            
        end
        
        function  rp = judgeReplaResp(rp)
            dbstop if error
            
            % evaluate the responsiveness of song replacements
            replaids = find( ~cellfun(@isempty, regexp([rp.list.stimuliname].','Repla') ));
            
            if isempty(replaids)
                return
            end
            % Find out that the repla stimuliname correspond to how many natrual songs
            bnames = unique(cellfun(@(x)Convert.bid(x,2),[rp.list(replaids).stimuliname].','Uni',0));
            
            findCorrespConID = @(x) intersect(find( ~cellfun(@isempty, regexp([rp.list.stimuliname].','norm') )),...
                find( ~cellfun(@isempty, regexp([rp.list.stimuliname].',x) )) ); % function handle to get corresponding norm ids
            
            findCorrespReplaID = @(x) intersect(find( ~cellfun(@isempty, regexp([rp.list.stimuliname].','Repla') )),...
                find( ~cellfun(@isempty, regexp([rp.list.stimuliname].',x) )) ); % function handle to get corresponding norm ids
            
            for k = 1:length(bnames)
                nid = findCorrespConID(bnames{k}); % 有可能得到多个nid
                % 如果有多个nid
                if length(nid) > 1
                    nid = nid(1);
                    disp('Dangerous Code: @Neuron.judgeReplaResp');
                end
                rids = findCorrespReplaID(bnames{k});
                
                for kk = 1:length(rids)
                    [Ini_y,Ini_replay] = Neuron.findConergentPointBetwenNormAndRepla(...
                        rp.list(nid).y,...
                        rp.list(rids(kk)).y );

                    num_of_zeros = find(rp.list(nid).y(Ini_y:end),1)-1; % 在分歧点后多长一段是0
                    Ini_y = Ini_y + num_of_zeros;
                    Ini_replay = Ini_replay + num_of_zeros;
                    try
                        rp.list(rids(kk)).targety = rp.list(rids(kk)).y(Ini_replay:Ini_replay + 0.2*32000); % 截取200ms
                        rp.list(rids(kk)).replaceparty = rp.list(rids(kk)).y(1:Ini_replay ); % 被替换的那一部分的y值
                    catch % if the second index overceed
                        rp.list(rids(kk)).targety = rp.list(rids(kk)).y(Ini_replay:end); % 截取到 end
                        rp.list(rids(kk)).replaceparty = rp.list(rids(kk)).y(1:Ini_replay ); % 被替换的那一部分的y值
                        disp(' MEssage@Neuron.judgeReplaResp :Overceed');
                    end
                    rp.list(rids(kk)).targetsptimes = Extract.sptimes_resetSP(rp.list(rids(kk)).sptimes,Ini_replay/32000,(Ini_replay + 0.2*32000)/32000);
                    rp.list(rids(kk)).replacepartsptimes = Extract.sptimes_resetSP(... % replacepart, sptimes
                        rp.list(rids(kk)).sptimes,1,Ini_replay/32000);
                    % corresponding pre (targety) data
                    rp.list(rids(kk)).pretargety = zeros(length(rp.list(rids(kk)).targety),1);%rp.list(nid).prey(end - 0.1*32000:end); % 截取100ms
                    rp.list(rids(kk)).pretargetsptimes = Extract.sptimes_resetSP(...
                        rp.list(rids(kk)).presptimes,length(rp.list(rids(kk)).pretargety)/32000 -0.2*32000,length(rp.list(rids(kk)).pretargety)/32000 );%length(rp.list(rids(kk)).prey)/32000 -0.1*32000
                    %calculate and judege whether the neuron respond to the target area or not
                    [rp.list(rids(kk)).replaresp,rp.list(rids(kk)).replaPvalue] = Neuron.UseTtestToJudegeRespOrNot(...
                        rp.list(rids(kk)).targety,rp.list(rids(kk)).targetsptimes,...
                        rp.list(rids(kk)).pretargety ,rp.list(rids(kk)).pretargetsptimes ,32000);
                    rp.list(rids(kk)).targetfr = length(vertcat(rp.list(rids(kk)).targetsptimes{:}))/200; % per milisecond
                end
                
                
            end
            
            
        end
        
        function toCheck(rp)
            % 功能说明：把replalist的神经元反应结果作图，进而探明到底哪些replacement可以使神经元反应
           dbstop if error
            replalist = rp.replalist; % let us make replalist as an internal variable
            if isempty(replalist)
                return
            end
            replalist = table2struct(sortrows(struct2table(replalist),'targetfr','descend')); % rank by firing rate
            
           
            replafig = {};
            prefig = {};
            replacepartfig = {};
            for k = 1:length(replalist)
                
                figure('Position',[2189 255 560 1000])
                Draw.two(replalist(k).targety, 32000,replalist(k).targetsptimes); % something wrong with the getReplallist
                 xlabel(sprintf(...
                    'mFR (target dur): %f, Result: % u ',replalist(k).targetfr,replalist(k).replaresp));%???? big problem here
                replacepartfig {k} = getframe(gcf).cdata;
                replafig{k} = getframe(gcf).cdata;
                
                close(gcf)
                
                figure('Position',[2189 255 560 1000])  %basline level prestimuli
                Draw.two(replalist(k).pretargety, 32000,replalist(k).pretargetsptimes);
                prefig{k} = getframe(gcf).cdata;
                close(gcf)
                
                
%                 figure('Position',[2189 255 200 840]) % the
%                 Draw.two(replalist(k).replaceparty, 32000,replalist(k).replacepartsptimes);
%                 replacepartfig {k} = getframe(gcf).cdata;
%                 close(gcf)
                
                figure('Position',[2189 255 400 1000])
               
                Draw.two(replalist(k).replaceparty, 32000,replalist(k).replacepartsptimes);
                replacepartfig {k} = getframe(gcf).cdata;
                close(gcf)
            end
            
            sum_replafig = vertcat(replafig{:});
            sum_prefig = vertcat(prefig{:});
            sum_replacepartfig = vertcat(replacepartfig{:});
            sum_all = horzcat(sum_prefig,sum_replacepartfig,sum_replafig);
            
            %imwrite(sum_two,sprintf('%s_Sum_repla_fig.png',A.formated_name));
            
            imwrite(sum_all,sprintf('%s_Sum_preNreplaceNtarget.png',rp.formated_name));
            
        end
        
        
        function replalist = does_neuron_respond_to_specific_transition(rp,transition_statics,syllable_theircatego)
            % 根据Catego Trans Probility，解答最终的问题，即是否不常见的combination不会使neuron反应，或者正相反
            replalist = rp.replalist;
            for k = 1:length(replalist)
                
                ids1 = find(~cellfun(@isempty, regexp({syllable_theircatego.unifragnames}.', replalist(k).concat1)));
                
                if ~isempty(ids1)
                    replalist(k).catego1 =  syllable_theircatego(ids1).catego;
                else
                    replalist(k).catego1 = 0;
                end
                
                ids2 = find(~cellfun(@isempty, regexp({syllable_theircatego.unifragnames}.', replalist(k).concat2)));
                if ~isempty(ids2)
                    replalist(k).catego2 =  syllable_theircatego(ids2).catego;
                else
                    replalist(k).catego2 = 0;
                end
                
                replalist(k).catego1catego2 = [replalist(k).catego1,replalist(k).catego2];
                
                transID = find(cellfun(@(x)isequal(replalist(k).catego1catego2,x),{transition_statics.firstsecond}.','Uni',1));
                if ~isempty(transID)
                    replalist(k).transrank = transition_statics(transID).rank;
                    replalist(k).count = transition_statics(transID).counts;
                     replalist(k).chance = transition_statics(transID).chance;
                else
                    replalist(k).transrank = 0;
                    replalist(k).count = 0;
                    replalist(k).chance = 0;
                end
                
            end
            
        end
        
    end
    
    methods(Static)
        
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

