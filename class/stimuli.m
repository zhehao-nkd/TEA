%%% a class to generate Stimuli

% segment songs to get a syllable table
% generate normalized songs
% generate reversed songs
% generate single syllables
% generate mirrored songs
% generate context replaced songs
% select syllable from different clusters
% add trigger channel

% the name of this class will be Stimuli??? ort other names???

classdef Stimuli < handle

    properties
%         raw % infment information struct
%         fs
%         splited % splited is a struct for each song, by modifying the splited, like delete a syllable or, replace the y to the normalize y, or even replace the syllable to other syllables, ne wStimuli will be gnerated
%         %%% That is to say, function normalize should be renewed to
%         %%% become a general function. Used by all other fucntions to
%         %%% assemble the syllable/elemnt within songs into a single song
%         responded  % it is the collection of songs that can activate a neurons
%         selected11 % selected songs that can activate a neuron
%         selected01    %not responded  ??????
%         selectedsplited %% splited data, but not full, instead , it is selected
        rawfraglist
        preprocessed% preprocessed info, which is firstly highpass-filtered and then normalized
       % separated_rawfraglist % separated processed
        separated
        fs
%         songnames  % song names within info
%         songnum  % number of songs within info
%         outdir
%         num_of_near
%         num_of_far
%         num_of_sim
%         repla_pregap % pre_gap for replacement
%         initial_pregap  %pre_gap of the first elements % seconds
% 

    end

    methods  % 底层算法，不为产生wav文件

        function st = Stimuli(info)
    
            st.fs = unique([info.fs].');  % a single value
            %st.preprocessed = Stimuli.getPrepro(info);
            summer = {};
            unique_songnames = unique(cellstr({info.songname}.'));
            for k = 1: length(unique_songnames)  % I change for to par-for 1103
                hitted = find( strcmp(cellstr({info.songname}.') ,unique_songnames{k}) );
                temp =  info(hitted);
                [~,index] = sortrows([temp.fragid].'); temp = temp(index); clear index
                %st.separated_rawfraglist{k} = temp;
                st.separated{k} = Substimuli(temp);
                summer{k} = st.separated{k}.preprocessed;
            end
            st.preprocessed = vertcat(summer{:});

%             s.splited = Stimuli.split(info);

%             s.initial_pregap = 0.05;

            % judge whether those infos are from a single song or multiple
            % song

%             s.songnames = unique(cellstr( {info.songname}.'));
%             s.songnum = length(s.songnames);
%             s.outdir = './';
%             s.set_near_far(10,10);


        end

        function st = setoutdir(st,outdir) % output directory of generated Stimuli
            st.outdir = outdir;
        end


    end

    methods  % 直接产生wav文件的方法

        function writeSongs(st,outdir) %Stimuli set 1

            if ~exist('outdir','var')
                outdir = 'Songs'; % default value
                mkdir(outdir);
            end


            for k = 1:length(st.separated)
                singlesong = st.separated{k}.getSong
                audiowrite(sprintf('%s\\norm-%s.wav',outdir,singlesong.songname),singlesong.y,32000);
            end
        end

        function writeDegsAndSamesongSyllables(st,outdir) % Stimuli set 2


            if ~exist('outdir','var')
                outdir  = 'Degs&SamesongFrags'; % default value
                mkdir(outdir);
            end


            for k = 1:length(st.separated)

                subdir = sprintf('%s\\%s',outdir,st.separated{k}.songname);
                mkdir(subdir);
                this_deg = st.separated{k}.getDegs;%(length(st.separated{k}.preprocessed));
                for kk = 1:length(this_deg)
                    audiowrite(sprintf('%s\\Deg-%s-%02u.wav',subdir,this_deg(kk).songname,this_deg(kk).rank),this_deg(kk).y,32000);
                end
            end

            for k = 1:length(st.separated)
                subdir = sprintf('%s\\%s',outdir,st.separated{k}.songname);
                % mkdir(subdir);
                this_samesong = st.separated{k}.getSamesongSyllables;
                for kk = 1:length(this_samesong)
                    audiowrite(sprintf('%s\\Frag-%s-samesong.wav',subdir,this_samesong(kk).unifragnames),this_samesong(kk).y,32000);
                end
            end

        end

        function writeReplacements(st,categolist) %set 3 Testing of the replacement
            dbstop if error

            wb = PoolWaitbar(length(st.separated),'For each CON song');
            outdir  = 'ReplasZeroGap';
            mkdir(outdir);
            for k = 1:length(st.separated)
                this_separated = st.separated{k};

                parfor kk = 1:length(this_separated.preprocessed) % par

                    if kk~= 1 % 这里是为了把原来的combination生成出来，且gap为0
                        local_categolist = categolist;
                        local_categolist(length(categolist) +1).initial = this_separated.preprocessed(kk-1).initial;
                        local_categolist(length(categolist) +1).terminal = this_separated.preprocessed(kk-1).terminal;
                        local_categolist(length(categolist) +1).songname = this_separated.preprocessed(kk-1).songname;
                        local_categolist(length(categolist) +1).motif = this_separated.preprocessed(kk-1).motif;
                        local_categolist(length(categolist) +1).y = this_separated.preprocessed(kk-1).y;
                        local_categolist(length(categolist) +1).fs = this_separated.preprocessed(kk-1).fs;
                        local_categolist(length(categolist) +1).fragI = this_separated.preprocessed(kk-1).fragI;
                        local_categolist(length(categolist) +1).fragid = this_separated.preprocessed(kk-1).fragid;
                        local_categolist(length(categolist) +1).unifragnames = this_separated.preprocessed(kk-1).unifragnames;
                        local_categolist(length(categolist) +1).fullname = this_separated.preprocessed(kk-1).fullname;
                        local_categolist(length(categolist) +1).catego = this_separated.preprocessed(kk-1).catego;
                        %local_categolist(length(categolist) +1).distance = this_separated.preprocessed(kk-1).distance;
                        local_categolist(length(categolist) +1).sat = this_separated.preprocessed(kk-1).sat;
                    else
                        local_categolist = categolist;
                    end
                    replalist = this_separated.getReplas(local_categolist,kk,0);

                    subdir = sprintf('%s\\Target-%s',outdir,replalist(1).targetname);
                    mkdir(subdir);
                    for kkk = 1:length(replalist)
                        %或许最好改成ms
                        audiowrite(sprintf('%s\\Repla-%s-Type%02u-before-%s-Type%02u-gapis-%.0f.wav',...
                            subdir,replalist(kkk).contextname,replalist(kkk).catego,replalist(kkk).targetname,...
                            this_separated.preprocessed(kk).catego,replalist(kkk).gapduration*1000),...
                            replalist(kkk).y,replalist(kkk).fs); % 毫秒
                    end

                end

                increment(wb);
            end

        end

        function writeChangeGapReplacements(st,categolist,gap_durations) % Stimuli set 4 确认replacement没什么影响

            outdir  = 'ReplasChangeGap';
            mkdir(outdir);
           % gap_durations = [nan,0.05,0.1,0.2,0.4,0.8,1.6];% 需要设置一个梯度的gap durations
            wb = PoolWaitbar(length(gap_durations),'对于每个gap duration');

            for gp = 1:length(gap_durations)

                for k = 1:length(st.separated)

                    this_separated = st.separated{k};

                    for kk = 1:length(this_separated.preprocessed) % par

                        if kk~= 1 % 这里是为了把原来的combination生成出来，且gap为0
                            local_categolist = categolist;
                            local_categolist(length(categolist) +1).initial = this_separated.preprocessed(kk-1).initial;
                            local_categolist(length(categolist) +1).terminal = this_separated.preprocessed(kk-1).terminal;
                            local_categolist(length(categolist) +1).songname = this_separated.preprocessed(kk-1).songname;
                            local_categolist(length(categolist) +1).motif = this_separated.preprocessed(kk-1).motif;
                            local_categolist(length(categolist) +1).y = this_separated.preprocessed(kk-1).y;
                            local_categolist(length(categolist) +1).fs = this_separated.preprocessed(kk-1).fs;
                            local_categolist(length(categolist) +1).fragI = this_separated.preprocessed(kk-1).fragI;
                            local_categolist(length(categolist) +1).fragid = this_separated.preprocessed(kk-1).fragid;
                            local_categolist(length(categolist) +1).unifragnames = this_separated.preprocessed(kk-1).unifragnames;
                            local_categolist(length(categolist) +1).fullname = this_separated.preprocessed(kk-1).fullname;
                            local_categolist(length(categolist) +1).catego = this_separated.preprocessed(kk-1).catego;
                            %local_categolist(length(categolist) +1).distance = this_separated.preprocessed(kk-1).distance;
                            local_categolist(length(categolist) +1).sat = this_separated.preprocessed(kk-1).sat;
                        else
                            local_categolist = categolist;
                        end

                        replalist = this_separated.getReplas(local_categolist,kk,gap_durations(gp));


                        subdir = sprintf('%s\\Target-%s',outdir,replalist(1).targetname);
                        mkdir(subdir);
                        for kkk = 1:length(replalist)
                            %或许最好改成ms
                            audiowrite(sprintf('%s\\Repla-%s-Type%02u-before-%s-Type%02u-gapis-%.0f.wav',...
                            subdir,replalist(kkk).contextname,replalist(kkk).catego,replalist(kkk).targetname,...
                            this_separated.preprocessed(kk).catego,replalist(kkk).gapduration*1000),...
                            replalist(kkk).y,replalist(kkk).fs); % 毫秒
                        end
                    end
                end

                increment(wb);


            end
       
         

        end




    end

    methods % 作图方法
        function draw_frag_scatter_from_folder(s,dirpath)


        end

        function draw_spec_in_space(s,imagename)
            % 把频谱图按照其坐标展现在二维空间中, 速度非常快，new version

            mergedeleinf = s.prepro;

            isize = 12000;
            kuanchangbi = 0.6;

            coor_1 = [mergedeleinf.coor_1].'; coor_2 = [mergedeleinf.coor_2].';
            coor_1 = rescale(coor_1, isize*0.03*kuanchangbi, isize*0.97*kuanchangbi);
            coor_2 = rescale(coor_2, isize*0.03, isize*0.97);

            %             parentimg = ones([isize,isize]);
            %             figure('Position',[202 155 1400 1000]);
            %             h = image(parentimg);
            %             hold on
            %             scale = 30; % scale = 50;
            %             for w = 1:length(mergedeleinf)%1: 100
            %                 img = rescale(Cal.spec(mergedeleinf(w).y,mergedeleinf(w).fs),0,1);
            %                 %                 [alpha,beta] = size(img);
            %                 %                 img = imresize(flip(img,1),[alpha,beta*2]);
            %                 image([coor_1(w),coor_1(w) + size(img,2)], [coor_2(w),coor_2(w) + size(img,1)], img, 'CDataMapping', 'scaled');
            %
            %             end

            parentimg = zeros([isize*kuanchangbi,isize]);
            %             scale = 30; % scale = 50;
            for w = 1:length(mergedeleinf)%1: 100
                img = rescale(Cal.spec(mergedeleinf(w).y,mergedeleinf(w).fs),0,1);
                [alpha,beta] = size(img);
                img = imresize(img,[alpha*0.6,beta*1.8]);
                parentimg(round(coor_1(w)):round(coor_1(w))+ size(img,1)-1,round(coor_2(w)):round(coor_2(w)) + size(img,2)-1) = img;
                % image([coor_1(w),coor_1(w) + size(img,2)], [coor_2(w),coor_2(w) + size(img,1)], img, 'CDataMapping', 'scaled');

            end

            parentimg = im2uint8(parentimg);

            cmap = hot(128); % Or whatever one you want.
            parentimg = ind2rgb(parentimg, cmap);
            parentimg = im2uint8(parentimg);

            % configure tiff
            t = Tiff(imagename,'w8');
            setTag(t,'ImageWidth',isize);
            setTag(t,'ImageLength',isize*(kuanchangbi+0.005));
            setTag(t,'Photometric',Tiff.Photometric.RGB);
            setTag(t,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
            setTag(t,'BitsPerSample',8);
            setTag(t,'SamplesPerPixel',3);

            % write data
            write(t,parentimg);
            close(t);

        end % plot spectrogram in scatter space

        function draw_songtrace_in_space(s, inputnames,filename)
            mergedeleinf = s.prepro;

            isize = 12000;
            kuanchangbi = 0.6;

            coor_1 = [mergedeleinf.coor_1].'; coor_2 = [mergedeleinf.coor_2].';
            coor_1 = rescale(coor_1, isize*0.03*kuanchangbi, isize*0.97*kuanchangbi);
            coor_2 = rescale(coor_2, isize*0.03, isize*0.97);

            for k = 1:length(mergedeleinf)
                mergedeleinf(k).coor_1 = coor_1(k);
                mergedeleinf(k).coor_2 = coor_2(k);
            end

            parentimg = ones([isize*kuanchangbi,isize]);
            %             figure('Position', get(0, 'Screensize'));
            figure('Position', [15,15,1910,1190]);
            h = image(parentimg);
            hold on

            axis off
            ax = gca;
            outerpos = ax.OuterPosition;
            ti = ax.TightInset;
            left = outerpos(1) + ti(1);
            bottom = outerpos(2) + ti(2);
            ax_width = outerpos(3) - ti(1) - ti(3);
            ax_height = outerpos(4) - ti(2) - ti(4);
            ax.Position = [left bottom ax_width ax_height];



            %             scale = 30; % scale = 50;
            for w = 1:length(mergedeleinf)%1: 100
                img = rescale(Cal.spec(mergedeleinf(w).y,mergedeleinf(w).fs),0,1);
                [alpha,beta] = size(img);
                img = imresize(img,[alpha*0.6,beta*1.8]);
                image([coor_2(w)-size(img,1)/2,coor_2(w) + size(img,1)/2], [coor_1(w)-size(img,2)/2,coor_1(w) + size(img,2)/2], img, 'CDataMapping', 'scaled');

            end

            linecmap = colormap('lines');
            for k = 1:length(inputnames)
                correspid = find(~cellfun(@isempty, regexp(cellstr({mergedeleinf.songname}.'),inputnames(k))));
                thissong_eleinf = mergedeleinf(correspid);
                [~,index] = sortrows([thissong_eleinf.fragid].');
                thissong_eleinf = thissong_eleinf(index); clear index % sort by fragid
                plot([thissong_eleinf.coor_2].',[thissong_eleinf.coor_1].','Color',linecmap(k,:),'LineWidth',2);

            end
            colormap('hot');
            saveas(gcf,filename);
            close(gcf)
            %             parentimg = zeros([isize*kuanchangbi,isize]);
            %             %             scale = 30; % scale = 50;
            %             for w = 1:length(mergedeleinf)%1: 100
            %                 img = rescale(Cal.spec(mergedeleinf(w).y,mergedeleinf(w).fs),0,1);
            %                 [alpha,beta] = size(img);
            %                 img = imresize(img,[alpha*0.9,beta*1.8]);
            %                 parentimg(round(coor_1(w)):round(coor_1(w))+ size(img,1)-1,round(coor_2(w)):round(coor_2(w)) + size(img,2)-1) = img;
            %                 % image([coor_1(w),coor_1(w) + size(img,2)], [coor_2(w),coor_2(w) + size(img,1)], img, 'CDataMapping', 'scaled');
            %
            %             end

            %             parentimg = im2uint8(parentimg);
            %
            %             cmap = jet(128); % Or whatever one you want.
            %             parentimg = ind2rgb(parentimg, cmap);
            %             parentimg = im2uint8(parentimg);
            %
            %             % configure tiff
            %             t = Tiff(imagename,'w8');
            %             setTag(t,'ImageWidth',isize);
            %             setTag(t,'ImageLength',isize*(kuanchangbi+0.005));
            %             setTag(t,'Photometric',Tiff.Photometric.RGB);
            %             setTag(t,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
            %             setTag(t,'BitsPerSample',8);
            %             setTag(t,'SamplesPerPixel',3);
            %
            %             % write data
            %             write(t,parentimg);
            %             close(t);

        end

        function draw_scatter(s,ids_green,ids_blue,ids_red)

            mergedeleinf = s.prepro;
            coor_1 = [mergedeleinf.coor_1].'; coor_2 = [mergedeleinf.coor_2].';
            coordinate = horzcat(coor_1,coor_2);

            % set initial state as black
            for h = 1:length(mergedeleinf)
                mergedeleinf(h).imagehandle = 0; % 2 means braod
            end

            if exist('ids_green','var')
                for h = ids_green(:).'
                    mergedeleinf(h).imagehandle = -1;
                end
            end

            if exist('ids_blue','var')
                for h = ids_blue(:).'
                    mergedeleinf(h).imagehandle = 2; % 2 means braod
                end
            end

            if exist('ids_red','var')
                for h = ids_red(:).'
                    mergedeleinf(h).imagehandle = 1;
                end
            end

            figure;
            hold on

            for f = 1: length(mergedeleinf)

                if mergedeleinf(f).imagehandle == 0
                    scatter(mergedeleinf(f).coor_1,mergedeleinf(f).coor_2,[],'k','filled'); % not used
                elseif mergedeleinf(f).imagehandle == 2
                    scatter(mergedeleinf(f).coor_1,mergedeleinf(f).coor_2,[],'b','filled'); % broad-blue
                elseif mergedeleinf(f).imagehandle == 1
                    scatter(mergedeleinf(f).coor_1,mergedeleinf(f).coor_2,[],'r','filled'); % near-red
                elseif mergedeleinf(f).imagehandle == -1
                    scatter(mergedeleinf(f).coor_1,mergedeleinf(f).coor_2,[],'g','filled'); % target_ele-green
                end
                %drawnow
            end
        end

        function draw_conscatter(s,ids_green,ids_blue,ids_red) % These ids are global
            mergedeleinf = s.prepro;

            con_ids = find(~cellfun(@isempty,regexp([mergedeleinf(:).songname].','CON|SPE')));

            con_eleinf = mergedeleinf(con_ids);


            coor_1 = [con_eleinf.coor_1].'; coor_2 = [con_eleinf.coor_2].';
            coordinate = horzcat(coor_1,coor_2);

            % set initial state as black
            for h = 1:length(con_eleinf)
                con_eleinf(h).imagehandle = 0; % 2 means broad
            end

            if exist('ids_green','var')
                for h = ids_green(:).'
                    target_id = find(con_ids == h);
                    con_eleinf( target_id).imagehandle = -1;
                end
            end

            if exist('ids_blue','var')
                for h = ids_blue(:).'
                    target_id = find(con_ids == h);
                    con_eleinf(target_id).imagehandle = 2; % 2 means braod
                end
            end

            if exist('ids_red','var')
                for h = ids_red(:).'
                    target_id = find(con_ids == h);
                    con_eleinf(target_id).imagehandle = 1;
                end
            end

            figure;
            hold on

            for f = 1: length(con_eleinf)

                if con_eleinf(f).imagehandle == 0
                    scatter(con_eleinf(f).coor_1,con_eleinf(f).coor_2,[],'k','filled'); % not used
                elseif con_eleinf(f).imagehandle == 2
                    scatter(con_eleinf(f).coor_1,con_eleinf(f).coor_2,[],'b','filled'); % broad-blue
                elseif con_eleinf(f).imagehandle == 1
                    scatter(con_eleinf(f).coor_1,con_eleinf(f).coor_2,[],'r','filled'); % near-red
                elseif con_eleinf(f).imagehandle == -1
                    scatter(con_eleinf(f).coor_1,con_eleinf(f).coor_2,[],'g','filled'); % target_ele-green
                end
                %drawnow
            end




        end

        function draw_conspace(s,global_target)

            % get random
            global_eleinf = s.prepro;
            senatus_ids = find(~cellfun(@isempty,regexp([global_eleinf(:).songname].','CON|SPE'))); % ids of CON/SPE songs in all_eleinf
            senatus_eleinf = global_eleinf(senatus_ids);

            %randomize the senatus
            rng(1);
            random_ids = senatus_ids(randperm(length(senatus_ids))); % randomize mergedeleinf


            %get ids
            far_ids = random_ids(1:s.num_of_far);
            nearby_ids = Stimuli.findnearby(senatus_eleinf,s.num_of_near,global_target);
            global_far_ids = [global_eleinf(far_ids).uniqueid].';
            global_nearby_ids = [global_eleinf(nearby_ids).uniqueid].';
            %draw
            %             s.draw_conscatter(far_ids);
            %             s.draw_scatter();
            s.draw_conscatter(global_far_ids,global_nearby_ids,global_target);

        end

        function draw_allspace(s,global_target)

            % get random
            global_eleinf = s.prepro;
            rng(1);
            random_ids = randperm (numel(global_eleinf));

            %get ids
            far_ids = random_ids(1:s.num_of_far);
            nearby_ids = Stimuli.findnearby(global_eleinf,s.num_of_near,global_target);

            %draw
            s.draw_scatter(far_ids,nearby_ids,global_target);

        end


    end

    methods(Static)

        function writeSyllablesFromDifferentCategories(categolist) % Stimuli set 5 有可能4和5能放在一起
            outdir = 'DifferentCategoFrags';
            mkdir(outdir);
            for k = 1:length(categolist)
                audiowrite(sprintf('%s\\Frag-%s-Type%02u.wav',outdir,categolist(k).unifragnames,categolist(k).catego),categolist(k).y,32000);
            end

        end

        function preprocessed = getPrepro(info)

            for k = 1: length(info)
                info(k).dur = length(info(k).y)/info(k).fs;

                info(k).uniqueid = k;

                if k==1 || ~strcmp(info(k-1).songname,info(k).songname)
                    info(k).pregap = inf;   %%%%%% This is a super stupid solution !!!!!!!!!!!!!!!!!!!!!!!!!
                else
                    info(k).pregap = (info(k).initial - info(k-1).terminal)/info(k).fs;
                end
            end

            temp = Substimuli.highpass(info,450);
            temp = Substimuli.normalize(temp, 0.05);

            if length(temp) == 1
                preprocessed = temp;
            else
                preprocessed = table2struct(sortrows(struct2table(temp),'uniqueid','ascend')); % re-order
            end
            % already highpass-filtered and normalized for each song

        end


        function fromConRespCopyResponsiveFrags(conresp_siginfo)
            % 这个方法是为了快速产生下一阶段的测试需要的stimuli
            dbstop if error
            % 找到所有在EachSongFrags里的文件的path
            outdir = '.\RespElicitingFrags';
            mkdir(outdir);

            sharedDepot = "G:\SharedDepot";
            specialDepot = "G:\SpecialDepot";

            temp1 = Extract.foldersAllLevel(sharedDepot).';
            temp2 = Extract.foldersAllLevel(specialDepot).';
            temp = vertcat(temp1,temp2);

            eachfragids = find(~cellfun(@isempty,regexp(cellstr(temp),'EachSongFrags')));
            fragsfolders = temp(eachfragids);

            for k = 1: length(conresp_siginfo)

                dirids = find(~cellfun(@isempty, regexp(fragsfolders,conresp_siginfo(k).songid)));
                contain_frag_dir = fragsfolders(dirids);

                files_in_dir = Extract.filename(contain_frag_dir,'*.wav');
                regexpression =  sprintf('%s\\S+%02u',conresp_siginfo(k).songid,conresp_siginfo(k).fragid);

                fileid = find(~cellfun(@isempty, regexp(files_in_dir,regexpression)));
                target_file_path = files_in_dir{fileid};

                [~,nameonly,ext] = fileparts(target_file_path);

                purename = strcat(nameonly,ext);
                destiney_file_path = fullfile(outdir,purename);
                copyfile(target_file_path,destiney_file_path);
            end

        end

        function  wavTransformer(audiopath)
            % 这个方法生成改变了某些feature的stimuli
            outdir = './transformed'

            [y,fs] = audioread(audiopath);
            [dir,name,ext] = fileparts(audiopath);

            mkdir(outdir);

            % strength_change
            for scale = [0.5^4,0.125,0.25,0.5,1,2,4]
                newy = DQ.transform(y,fs,'strength_change',scale);
                audiowrite(sprintf('%s/trans-strength-%g-%s.wav',outdir,scale,name),newy,fs);
            end

            % frequency change

            for scale = [0.4,0.6,0.8,1.2,1.4,1.6]
                newy = DQ.transform(y,fs,'frequency_change',scale);

                audiowrite(sprintf('%s/trans-freq-%g-%s.wav',outdir,scale,name),newy,fs);
            end

            % duration change


            for scale = [0.2,0.4,0.6,0.8,1.2,1.4,1.6,1.8,2]
                newy = DQ.transform(y,fs,'duration_change',scale);

                audiowrite(sprintf('%s/trans-dur-%g-%s.wav',outdir,scale,name),newy,fs);
            end



        end


        function finalelen = timeSpent(dir_path)
            % 估算播放一个文件夹中的stimuli所需要的时间
            % calculate the time it takes for Stimuli presentation

            files = Extract.filename(dir_path,'*.wav');

            len = [];
            for k = 1: length(files)

                [y,fs] = audioread(files{k});

                len(k) = size(y,1)/fs;
            end

            finallen = sum(len*10)/60; % mins

            fprintf('It takes %f mins\n',finallen);

        end

        function copyStimuli(varargin)
            % copy  all the Stimuli of specific birdname as inputs
            % 后续应该手动删掉不需要的Stimuli,这样比在copy阶段筛选更快

            outdir = './Copied';   mkdir(outdir);
            Stimuli_DEPOT = "E:\Stimuli_Synthesis_Source\StimuliDepot";

            files = Extract.filesAllLevel(Stimuli_DEPOT,'*.wav');
            for k = 1:length(varargin)
                selected_filenames = {files{find(~cellfun(@isempty, regexp(files,varargin{k})))}}.';
                for kk = 1:length(selected_filenames)
                    [~,name,ext] = fileparts(selected_filenames{kk});
                    copyfile(selected_filenames{kk},fullfile(outdir,strcat(name,ext)) );
                end
            end

        end


    end

    methods % Frozen or Deprecated

        function s = Deprecated_setpregap(pregap)
            s.repla_pregap = pregap;
        end

        function s = Deprecated_set_near_far(s,num_of_near,num_of_far) % near is how many closest element to be used, same idea to far
            s.num_of_near = num_of_near;
            s.num_of_far = num_of_far;
        end

        function replaced = Deprecated_replace(s,radish,pit)    % 一个萝卜一个坑
            % I have a variable of syl info which indicate whether this
            % sylinf can activate a neuron within a asong and another
            % variable indicate that it activate multiple songs

            pitsongname = s.prepro(pit).songname;
            pitid = s.prepro(pit).fragid;

            radsongname = s.prepro(radish).songname;
            radid = s.prepro(radish).fragid;

            songnames = {};
            for k = 1: length(s.splited)
                songnames{k} = unique(cellstr({s.splited{k}.songname}.'));
                songnames{k} = songnames{k}{1};
            end


            oldidx = find(~cellfun(@isempty,regexp(songnames,pitsongname)));
            old = s.splited{oldidx};



            %replacement
            former = old;
            pit_splited = find(pitid==[former.fragid].');
            former(pit_splited :height(old)) = [];   % delete later part


            % add later parts of syllables from the radish into the new

            radsongidx = find(~cellfun(@isempty,regexp(songnames,radsongname)));
            radsong = s.splited{radsongidx}; % which splitted struct
            radish_splited = find(radid == [radsong.fragid].');
            later = radsong(radish_splited : height(radsong));


            % change the later part
            later(1).pregap = old(pit_splited).pregap;
            timediff = radsong(radish_splited).initial - old(pit_splited).initial;

            for mm = 1: length(later)
                later(mm).initial =  later(mm).initial + timediff;
                later(mm).terminal =  later(mm).terminal + timediff;

            end


            replaced = vertcat(former, later);



            for obama = 1: length(replaced)
                if obama == pit_splited   % if it is the replacement location
                    replaced(obama).repla = 1;
                else
                    replaced(obama).repla = 0;
                end
            end


        end


        function replaced = Deprecated_New_replace(s,radish,pit)    % 一个萝卜一个坑
            % I have a variable of syl info which indicate whether this
            % sylinf can activate a neuron within a asong and another
            % variable indicate that it activate multiple songs

            pitsongname = s.prepro(pit).songname;
            pitid = s.prepro(pit).fragid;

            radsongname = s.prepro(radish).songname;
            radid = s.prepro(radish).fragid;

            songnames = {};
            for k = 1: length(s.splited)
                songnames{k} = unique(cellstr({s.splited{k}.songname}.'));
                songnames{k} = songnames{k}{1};
            end


            oldidx = find(~cellfun(@isempty,regexp(songnames,pitsongname)));
            old = s.splited{oldidx};



            %replacement
            former = old;
            pit_splited = find(pitid==[former.fragid].');
            former(pit_splited :height(old)) = [];   % delete later part


            % add later parts of syllables from the radish into the new

            radsongidx = find(~cellfun(@isempty,regexp(songnames,radsongname)));
            radsong = s.splited{radsongidx}; % which splitted struct
            radish_splited = find(radid == [radsong.fragid].');
            later = radsong(radish_splited : height(radsong));


            % change the later part
            later(1).pregap = old(pit_splited).pregap;
            timediff = radsong(radish_splited).initial - old(pit_splited).initial;

            for mm = 1: length(later)
                later(mm).initial =  later(mm).initial + timediff;
                later(mm).terminal =  later(mm).terminal + timediff;

            end


            replaced = vertcat(former, later);



            for obama = 1: length(replaced)
                if obama == pit_splited   % if it is the replacement location
                    replaced(obama).repla = 1;
                else
                    replaced(obama).repla = 0;
                end
            end


        end



        function incre = Deprecared_incremental(s,singlesyl,start)
            incre = struct;
            for k = 1: start
                rank =k -1;
                smallsylinf = singlesyl(start-rank:end);
                incre(k).y = s.assemble(smallsylinf);
                incre(k).rank = rank;
            end

        end

        function deg = Deprecated_degressive2(s,singlesyl, initial, terminal)  % initial-terminal
            deg = struct;
            diff = terminal - initial + 1;
            for k = 1: diff
                rank =k - 1;
                smallsylinf = singlesyl(initial+rank:terminal);
                deg(k).y = s.assemble(smallsylinf);
                deg(k).rank = rank;
            end

            % tester
            %             figure
            %             for mm = 1: length(ans)
            %                 subplot(9,1,mm);
            %                Draw.spec( ans(mm).y, 32000);
            %             end

        end


        function pregap = Deprecated_getpregap(s,num)

            pregap = s.raw(num).pregap;
        end

        function Deprecated_find_close_10(s)
            eleinf = s.prepro;
            coor_1 = [eleinf.coor_1].'; coor_2 = [eleinf.coor_2].';
            coordinate = horzcat(coor_1,coor_2);

            distance = [];

            for indi = 1: length(coordinate)

                for a = 1: length(coordinate)
                    distance(a) = norm(coordinate(a,:) - coordinate(indi,:));
                end

                [B,I] = mink(distance,16);


                %cloest = coordinate(I,:);

                %uni_id_cloest = [eleinf(I).fragid].';  % the unique id of element which is most cloest 10 element to current element
                eleinf(indi).closest10 = I(2:15);
                %try to verift this part
                %
                %             figure;
                %             scatter(coordinate(:,1),coordinate(:,2));
                %             hold on
                %             scatter(cloest(:,1),cloest(:,2));

            end
            s.prepro = eleinf;
        end

        function [nearall,nearsena,farall,farsena] = Deprecated_get_near_far_ids(s,global_id)

            global_eleinf = s.prepro;
            senatus_ids = find(~cellfun(@isempty,regexp([global_eleinf(:).songname].','CON|SPE'))); % ids of CON/SPE songs in all_eleinf
            senatus_eleinf = global_eleinf(senatus_ids);

            nearall = Stimuli.findnearby(global_eleinf,s.num_of_near,global_id);

            dijige = find([senatus_eleinf.uniqueid].'==global_id);
            temp = Stimuli.findnearby(senatus_eleinf,s.num_of_near,dijige);
            nearsena = [senatus_eleinf(temp).uniqueid].';

            rng(1);
            temp = randperm (numel(global_eleinf));
            farall = temp(1:s.num_of_far);
            rng(1);
            temp = randperm (numel(senatus_eleinf));
            temp = temp(1:s.num_of_far);
            farsena = [senatus_eleinf(temp).uniqueid].';


        end


        function new = Deprecated_normalize2(s,old,tarrms)
            % this normalize all element/syllable to have the same rms
            %old = s.splited;
            new = {};

            if ~exist('tarrms','var')
                tarrms = 0.05; % default rms
            end

            for k = 1: length(old)

                new{k} = TempStimcore(old{k}).normalize2(tarrms);


            end

        end

        function Deprecated_set_num_of_sim(s,num_of_sim)
            s.num_of_sim = num_of_sim; % set the number of similar frags(just the same elements in different motifs
        end

        function all_near = Deprecated_from_senatus_near_find_all_near(s,senatus_near,global_id) % as explained by the name
            % find global index of nearby elements
            global_eleinf = s.prepro;

            %%%%%% Some thing extremely wrong!!!!!!!!!
            disp('Wrong !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            senatus_ids = find(~cellfun(@isempty,regexp([global_eleinf(:).songname].','CON|SPE'))); % ids of CON/SPE songs in all_eleinf
            senatus_eleinf = global_eleinf(senatus_ids);

            % find the corresponding ids of global_id in senatus_eleinf
            dijige = find([senatus_eleinf.uniqueid].'==global_id);
            nearbyids = Stimuli.findnearby(senatus_eleinf,senatus_near,dijige);


            allnearby_rank = Stimuli.findnearby(global_eleinf,length(global_eleinf) -1,global_id);
            unique_id_all_nearby = [global_eleinf(allnearby_rank).uniqueid].';
            %allnearby_rank = Stimuli.findnearby(global_eleinf,length(global_eleinf) -1,global_id);
            all_near = find(unique_id_all_nearby == senatus_eleinf(nearbyids(end)).uniqueid);
            % all_near = find(allnearby_rank == senatus_eleinf(357).uniqueid);

            % To test whether this function works well
            %             figure;
            %             for w = 1: length(nearbyids)
            %                 edge(w) = find(unique_id_all_nearby == senatus_eleinf(nearbyids(w)).uniqueid);
            %             end
            %             plot(edge);

        end



        function Deprecated_writeDeg_target(s,target)
            if s.songnum > 1
                disp('Hey! This dataset contains multiple song, do you really want to generate degressive songs based on it?');

            end


            if ~exist('initial','var')
                initial =  1;
            end

            if ~exist('terminal','var')
                terminal =  length(s.prepro);
            end
            deg = Stimuli.degressive(s.prepro, initial, terminal);

            nowtime = datestr(datetime('now'),'yyyy-mmm-dd-HH-MM-SS');
            degname = s.songnames{1};

            if isempty(s.outdir)

                dir = sprintf('%s/degressive_%s',PM.TEAOutputFolder,nowtime);
                mkdir(dir);

                for nn = 1: length(deg)
                    audiowrite(sprintf('%s/deg-%s-%u.wav',dir,degname,deg(nn).rank),deg(nn).y,s.fs);
                end
            else
                dir = sprintf('%s/Deg-%s',s.outdir,degname);
                mkdir(dir);

                for nn = 2: length(deg)
                    audiowrite(sprintf('%s/deg-%s-%02u.wav',dir,degname,deg(nn).rank),deg(nn).y,s.fs);
                end

            end

        end

        function Deprecated_writedegressive(s, initial, terminal)
            % initial 是deg起始处，terminal是deg结束处

            if s.songnum > 1
                disp('Hey! This dataset contains multiple song, do you really want to generate degressive songs based on it?');
            end


            if ~exist('initial','var')
                initial =  1;
            end

            if ~exist('terminal','var')
                terminal =  length(s.prepro);
            end
            deg = Stimuli.degressive(s.prepro, initial, terminal);

            nowtime = datestr(datetime('now'),'yyyy-mmm-dd-HH-MM-SS');
            degname = s.songnames{1};

            if isempty(s.outdir)

                dir = sprintf('%s/degressive_%s',PM.TEAOutputFolder,nowtime);
                mkdir(dir);

                for nn = 1: length(deg)
                    audiowrite(sprintf('%s/deg-%s-%u.wav',dir,degname,deg(nn).rank),deg(nn).y,s.fs);
                end
            else
                dir = sprintf('%s/Deg-%s',s.outdir,degname);
                mkdir(dir);

                for nn = 2: length(deg)
                    audiowrite(sprintf('%s/deg-%s-%02u.wav',dir,degname,deg(nn).rank),deg(nn).y,s.fs);
                end

            end
        end

        function Deprecated_New_writereplace(s,radish,pit)
            % 需要一个名为terminal的parameter
            nowtime = datestr(datetime('now'),'yyyy-mmm-dd-HH-MM-SS');
            dir = sprintf('%s/replaced_%s',PM.TEAOutputFolder,nowtime);

            mkdir(dir);
            replaced = s.replace(radish,pit);
            %[~, ~,beta] = unique(cellstr({replaced.songname}.'),'stable');
            start = find([replaced(:).repla].');
            %start = beta(2);
            sylname = sprintf('%s%u-follow-%s',replaced(start).songname,replaced(start).fragid,replaced(start-1).songname);
            incre = s.incremental(replaced,start);

            for k = 1: length(incre)
                audiowrite(sprintf('%s/incre%u-%s.wav',dir,incre(k).rank,sylname),incre(k).y,s.fs);
            end


        end

        function Deprecated_writereplace(s,radish,pit)
            % 需要一个名为terminal的parameter
            nowtime = datestr(datetime('now'),'yyyy-mmm-dd-HH-MM-SS');
            dir = sprintf('%s/replaced_%s',PM.TEAOutputFolder,nowtime);

            mkdir(dir);
            replaced = s.replace(radish,pit);
            %[~, ~,beta] = unique(cellstr({replaced.songname}.'),'stable');
            start = find([replaced(:).repla].');
            %start = beta(2);
            sylname = sprintf('%s%u-follow-%s',replaced(start).songname,replaced(start).fragid,replaced(start-1).songname);
            incre = s.incremental(replaced,start);

            for k = 1: length(incre)
                audiowrite(sprintf('%s/incre%u-%s.wav',dir,incre(k).rank,sylname),incre(k).y,s.fs);
            end


        end

        function Deprecated_writeintofrags(s,number) % write into segments

            if isempty(s.outdir) % if s.outdir is empty
                newinf = s.prepro;

                nowtime = datestr(datetime('now'),'yyyy-mmm-dd-HH-MM-SS');
                dir = sprintf('%s/frags_%s',PM.TEAOutputFolder,nowtime);
                mkdir(dir);

                for xx = 1: length(newinf)

                    audiowrite(sprintf('%s/syl-%s-%u.wav',dir,newinf(xx).songname,newinf(xx).fragid), newinf(xx).y,s.fs);
                end

            else


                newinf = s.prepro;

                songname = unique([newinf.songname].');

                if length(songname) > 1
                    songname = strcat(songname{:});
                end

                targetfragid = number;
                dir = sprintf('%s/Detail-%s-%u',s.outdir,songname,targetfragid);

                mkdir(dir);

                for xx = 1: length(newinf)

                    audiowrite(sprintf('%s/Single-%s-%u.wav',dir,newinf(xx).songname,newinf(xx).fragid), newinf(xx).y,s.fs);
                end



            end

        end

        function Deprecated_writenorm(s,~,~)

            if isempty(s.outdir)
                nowtime = datestr(datetime('now'),'yyyy-mmm-dd-HH-MM-SS');
                dir = sprintf('%s/norm_%s',PM.TEAOutputFolder,nowtime);
            else
                dir = sprintf('%s/ConsMirsRevs',s.outdir);
            end

            mkdir(dir);

            source = Stimuli.split(s.prepro);

            wb = waitbar(0,'In Progress');
            for m = 1: length(source)
                waitbar(m/length(source),wb,'In Progress');
                summer = s.assemble(source{m});
                audiowrite(sprintf('%s/norm-%s.wav',dir,string(unique(cellstr({source{m}.songname}.')))),summer,s.fs);

                disp(string(unique(cellstr({source{m}.songname}.'))));
            end
            close(wb);


        end

        function Deprecated_writereverses(s,~,~)


            dir = sprintf('%s/ConsMirsRevs',s.outdir);
            mkdir(dir);

            fraginf = s.prepro;
            for k = 1: length(fraginf)
                fraginf(k).y = flip(fraginf(k).y); % first flip
            end



            source = Stimuli.split(fraginf);

            wb = waitbar(0,'In Progress');
            for m = 1: length(source)
                waitbar(m/length(source),wb,'In Progress');
                %                 Comment.initial = [source{m}.initial].';
                %                 Comment.terminal = [source{m}.terminal].';
                %                 commentstr = jsonencode(Comment);
                summer = flip(s.assemble(source{m})); % second flip
                audiowrite(sprintf('%s/reverse-%s.wav',dir,string(unique(cellstr({source{m}.songname}.')))),summer,s.fs);

                disp(string(unique(cellstr({source{m}.songname}.'))));
            end
            close(wb);


        end

        function Deprecated_writemirrors(s,~,~)


            dir = sprintf('%s/ConsMirsRevs',s.outdir);
            mkdir(dir);

            fraginf = s.prepro;

            source = Stimuli.split(fraginf);

            wb = waitbar(0,'In Progress');
            for m = 1: length(source)
                waitbar(m/length(source),wb,'In Progress');
                %                 Comment.initial = [source{m}.initial].';
                %                 Comment.terminal = [source{m}.terminal].';
                %                 commentstr = jsonencode(Comment);
                summer = flip(s.assemble(source{m})); % second flip
                audiowrite(sprintf('%s/mirror-%s.wav',dir,string(unique(cellstr({source{m}.songname}.')))),summer,s.fs);

                disp(string(unique(cellstr({source{m}.songname}.'))));
            end
            close(wb);


        end

        function Deprecated_writecategories(s,number,gap) % write different categories
            dbstop if error

            mergedeleinf = s.prepro;

            [samesongeles,targetfragid] = Stimuli.samesong(mergedeleinf,number);

            targetsongname = mergedeleinf(number).songname;

            samesongeles(targetfragid).pregap = gap;
            samesongeles(1:targetfragid-1) = [];
            if isfield(samesongeles,'similar')
                samesongeles = rmfield(samesongeles,'similar');
            end


            fraginf = s.prepro;
            labels = unique([fraginf(:).label].');

            nowtime = datestr(datetime('now'),'yyyy-mmm-dd-HH-MM-SS');
            dir = sprintf('%s/catego_precede_%s_%u_%s',PM.TEAOutputFolder,targetsongname,targetfragid,nowtime);
            mkdir(dir);

            if ~isempty(find(labels == -1))  % remove not element rows
                index = find(labels == -1)
                labels(index) = [];
            end

            for kk = 1: length(labels)
                idxs = find([fraginf(:).label].' == labels(kk));
                if length(idxs) == 1
                    thisT = struct2table(fraginf(idxs),'AsArray',1);
                else
                    thisT =  struct2table(fraginf(idxs));
                end

                thisT = sortrows(thisT,'fragid');  % 这里可能会改   order the struct based on the syllable number
                splited{kk} = table2struct(thisT);
            end

            for m = 1: length(splited)
                number = 5;  %%% edit here ! how many replacements in each catego will take place
                thiscatego = splited{m};
                catego = thiscatego(1).label;
                thiscatego = rmfield(thiscatego,'label');
                len = min(length(thiscatego),number);


                for jj = 1: len

                    summedele = vertcat(thiscatego(jj), samesongeles);

                    summedy = Stimuli.assemble(summedele);
                    % figure; Draw.spec(summedy,32000)
                    audiowrite(sprintf('%s/catego-%u-%s-%u-before-%s-%u-%s.wav',dir,catego,thiscatego(jj).songname,thiscatego(jj).fragid,targetsongname,targetfragid,num2str(gap)),summedy,s.fs);

                    % audiowrite(sprintf('%s/catego-%u-%s-%s.wav',dir,catego,string(unique(cellstr({thiscatego(jj).songname}.'))),num2str(gap)),summedy,s.fs);
                end

            end
        end

        function Deprecated_writetransform(s,number)  % number is thr uniqueid



            [same,localidx] = Stimuli.samesong(s.prepro,number);


            raw = same(localidx:end);





            nowtime = datestr(datetime('now'),'yyyy-mmm-dd-HH-MM-SS');
            dir = sprintf('%s/transform_%s-%u_%s',PM.TEAOutputFolder,raw(1).songname,raw(1).fragid,nowtime);
            mkdir(dir);



            % strength_change
            for scale = [0.25,0.5,1,2,4]
                temp = raw;
                temp(1).y = DQ.transform(temp(1).y,s.fs,'strength_change',scale);
                summedy = Stimuli.assemble(temp);
                audiowrite(sprintf('%s/trans-%s-%u-strength-%g.wav',dir,temp(1).songname,temp(1).fragid,scale),summedy,s.fs);
            end

            % frequency change

            for scale = [0.6,0.8,1.2,1.4]
                temp = raw;
                temp(1).y = DQ.transform(temp(1).y,s.fs,'frequency_change',scale);
                summedy = Stimuli.assemble(temp);
                audiowrite(sprintf('%s/trans-%s-%u-freq-%g.wav',dir,temp(1).songname,temp(1).fragid,scale),summedy,s.fs);
            end

            % duration change


            for scale = [0.6,0.8,1.2,1.4]
                temp = raw;
                temp(1).y = DQ.transform(temp(1).y,s.fs,'duration_change',scale);
                summedy = Stimuli.assemble(temp);
                audiowrite(sprintf('%s/trans-%s-%u-dur-%g.wav',dir,temp(1).songname,temp(1).fragid,scale),summedy,s.fs);
            end



        end

        % to use the following functions, the eleinf must contain
        % coordinate information, which might be acquired from
        % Experiment-object

        function Deprecated_writeReplace_fixedEleinfToReplace(s,number,replace_eleinf)

            dbstop if error

            % Extract the samesong_eleinf after the target_index ele
            mergedeleinf = s.prepro;
            [local_eleinf,local_id] = Stimuli.samesong(mergedeleinf,number);
            targetsongname = mergedeleinf(number).songname;
            if ~isempty(s.repla_pregap)
                local_eleinf(local_id).pregap = s.repla_pregap;
            end
            local_eleinf(1:local_id-1) = [];

            % assemble the eleinf to continuous y values
            collect = [];
            for k = 1: length(local_eleinf)
                if local_eleinf(k).pregap == inf
                    local_eleinf(k).pregap = 0.05; % Dangerous code!!
                end
                collect = [collect;zeros([round(local_eleinf(k).pregap*local_eleinf(k).fs),1]);local_eleinf(k).y];
            end

            %dir/path issues
            dir = sprintf('%s/Detail-%s-%02u',s.outdir,targetsongname,local_id);
            mkdir(dir);

            % write nearby-replace songs
            rng(1);
            random_ids = randperm (numel(replace_eleinf));
            %far_ids = random_ids(1:s.num_of_far);

            for h = 1: length(replace_eleinf)

                to_be_summed = replace_eleinf(h).y;
                summedy = [to_be_summed;collect];
                % figure; Draw.spec(summedy,32000)
                audiowrite(sprintf('%s/Repla-%s-%02u-before-%s-%02u-gapis-%s.wav',dir,replace_eleinf(h).songname,replace_eleinf(h).fragid,targetsongname,local_id,num2str(mergedeleinf(number).pregap)),summedy,s.fs);

            end

            %             % draw to verify
            %             s.draw_scatter(far_ids,number);
            %             title('Repla_far_from_all');
            %             saveas(gcf,sprintf('%s/Repla-far-from-all.png',dir));
            %             close(gcf)



        end

        function Deprecated_writeFrag_samesong(s,global_id) % given the index of a simuli, write every single element from the same song

            % s.prepro include filtering of eleinf
            global_eleinf = s.prepro;

            % Extract the inf for song contain number_indexed element
            this_songname = global_eleinf(global_id).songname;
            this_song_ids =  find(~cellfun(@isempty, regexp(this_songname,[global_eleinf.songname].'))) ;
            local_eleinf = global_eleinf(this_song_ids);

            % define dir,path,name
            songname = unique([local_eleinf.songname].');
            local_fradid = global_eleinf(global_id).fragid;
            dir = sprintf('%s/Detail-%s-%02u',s.outdir,songname,local_fradid);
            mkdir(dir);

            %write wavs
            for xx = 1: length(local_eleinf)
                audiowrite(sprintf('%s/Frag-%s-%02u.wav',dir,local_eleinf(xx).songname,local_eleinf(xx).fragid), local_eleinf(xx).y,s.fs);
            end


        end

        function Deprecated_writeFrag_near_from_all(s,global_id)% function to generat indivcidual element song with the new method

            % near---same catego eles
            % all---all eles in the global_eleinf is used

            % find global index of nearby elements
            global_eleinf = s.prepro;
            nearbyids = Stimuli.findnearby(global_eleinf,s.num_of_near,global_id)


            % define path and name
            songname = global_eleinf(global_id).songname;
            local_fragid = global_eleinf(global_id).fragid;
            dir = sprintf('%s/Detail-%s-%02u',s.outdir,songname,local_fragid);
            mkdir(dir);

            % construct nearbyinf and write wavs
            nearbyinf = global_eleinf(nearbyids);

            for xx = 1: length(nearbyinf)
                audiowrite(sprintf('%s/Frag-%s-%02u.wav',dir,nearbyinf(xx).songname,nearbyinf(xx).fragid), nearbyinf(xx).y,s.fs);
            end

            %             %draw to verify
            %             s.draw_scatter(nearbyids,nearbyids);
            %             title('Frag_near_from_all');
            %             saveas(gcf,sprintf('%s/Frag-near-from-all.png',dir));
            %             close(gcf)

        end

        function Deprecated_writeFrag_near_from_senatus(s,global_id)% function to generat indivcidual element song with the new method

            % near---same catego eles
            % all---all eles in the global_eleinf is used

            % find global index of nearby elements
            global_eleinf = s.prepro;
            senatus_ids = find(~cellfun(@isempty,regexp([global_eleinf(:).songname].','CON|SPE'))); % ids of CON/SPE songs in all_eleinf
            senatus_eleinf = global_eleinf(senatus_ids);

            % find the corresponding ids of global_id in senatus_eleinf
            dijige = find([senatus_eleinf.uniqueid].'==global_id);
            nearbyids = Stimuli.findnearby(senatus_eleinf,s.num_of_near,dijige)


            % define path and name
            songname = global_eleinf(global_id).songname;
            local_fragid = global_eleinf(global_id).fragid;
            dir = sprintf('%s/Detail-%s-%02u',s.outdir,songname,local_fragid);
            mkdir(dir);

            % construct nearbyinf and write wavs
            nearbyinf = senatus_eleinf(nearbyids);

            for xx = 1: length(nearbyinf)
                audiowrite(sprintf('%s/Frag-%s-%02u.wav',dir,nearbyinf(xx).songname,nearbyinf(xx).fragid), nearbyinf(xx).y,s.fs);
            end

            %             %draw to verify
            %             s.draw_scatter([senatus_eleinf(nearbyids).uniqueid].',global_id);
            %             title('Frag_near_from_senatus');
            %             saveas(gcf,sprintf('%s/Frag-near-from-senatus.png',dir));
            %             close(gcf)

        end

        function Deprecated_writeFrag_far_from_all(s)% function to generat indivcidual element song with the new method

            % construct and randmize eleinf
            global_eleinf = s.prepro;
            rng(1);
            random_ids = randperm (numel(global_eleinf)); % randomize mergedeleinf

            % path/name
            dir = sprintf('%s/FarFrag',s.outdir);
            mkdir(dir);

            % write wavs
            far_ids = random_ids(1:s.num_of_far);
            for xx = far_ids(:).'
                audiowrite(sprintf('%s/Frag-%s-%02u.wav',dir,global_eleinf(xx).songname,global_eleinf(xx).fragid),global_eleinf(xx).y,s.fs);
            end

            %             %draw to verify
            %             s.draw_scatter(far_ids);
            %             title('Frag_far_from_all');
            %             saveas(gcf,sprintf('%s/Frag-far-from-all.png',dir));
            %             close(gcf)
        end

        function Deprecated_writeFrag_far_from_senatus(s)% function to generat indivcidual element song with the new method

            % Extract senatus_eleinf
            global_eleinf = s.prepro;
            senatus_ids = find(~cellfun(@isempty,regexp([global_eleinf(:).songname].','CON'))); % ids of CON/SPE songs in all_eleinf
            senatus_eleinf = global_eleinf(senatus_ids);

            %randomize the senatus
            rng(1);
            random_ids = randperm (numel(senatus_eleinf)); % randomize mergedeleinf

            % path/name
            dir = sprintf('%s/FarFrag',s.outdir);
            mkdir(dir);

            % write wavs
            far_ids = random_ids(1:s.num_of_far);
            disp('Far_IDS_Of_Frag');
            disp(far_ids);
            disp([senatus_eleinf(far_ids).songname]);
            for xx = far_ids(:).'
                audiowrite(sprintf('%s/Frag-%s-%02u.wav',dir,senatus_eleinf(xx).songname,senatus_eleinf(xx).fragid),senatus_eleinf(xx).y,s.fs);
            end

            %             %draw to verify
            %             s.draw_scatter([senatus_eleinf(far_ids).uniqueid].');
            %             title('Frag_far_from_senatus');
            %             saveas(gcf,sprintf('%s/Frag-far-from-senatus.png',dir));
            %             close(gcf)

        end

        function Deprecated_writeRepla_near_from_all(s,global_id) % function to generate catego_replace song with the new method

            dbstop if error

            % Extract the samesong_eleinf after the target_index ele
            mergedeleinf = s.prepro;
            [local_eleinf,local_id] = Stimuli.samesong(mergedeleinf,global_id);
            targetsongname = mergedeleinf(global_id).songname;
            if ~isempty(s.repla_pregap)
                local_eleinf(local_id).pregap = s.repla_pregap;
            end

            local_eleinf(1:local_id-1) = [];

            % assemble the eleinf to continuous y values
            collect = [];
            for k = 1: length(local_eleinf)
                collect = [collect;zeros([round(local_eleinf(k).pregap*local_eleinf(k).fs),1]);local_eleinf(k).y];
            end

            %dir/path issues
            dir = sprintf('%s/Detail-%s-%02u',s.outdir,targetsongname,local_id);
            mkdir(dir);

            % write nearby-replace songs
            nearby_ids = Stimuli.findnearby(mergedeleinf,s.num_of_near,global_id-1);
            for h = 1: length(nearby_ids)

                rp_id = nearby_ids(h);  % mergedeleinf id for replacement
                to_be_summed = mergedeleinf(rp_id).y;
                summedy = [to_be_summed;collect];
                % figure; Draw.spec(summedy,32000)
                audiowrite(sprintf('%s/Repla-%s-%02u-before-%s-%02u-gapis-%s.wav',dir,mergedeleinf(rp_id).songname,mergedeleinf(rp_id).fragid,targetsongname,local_id,num2str(mergedeleinf(global_id).pregap)),summedy,s.fs);

            end

            %             %draw to verify
            %             s.draw_scatter(nearby_ids,global_id);
            %             title('Repla_near_from_all');
            %             saveas(gcf,sprintf('%s/Repla-near-from-all.png',dir));
            %             close(gcf)

        end

        function Deprecated_writeRepla_near_from_senatus(s,number)
            dbstop if error

            % Extract the samesong_eleinf after the target_index ele
            mergedeleinf = s.prepro;
            [local_eleinf,local_id] = Stimuli.samesong(mergedeleinf,number);
            targetsongname = mergedeleinf(number).songname;
            if ~isempty(s.repla_pregap)
                local_eleinf(local_id).pregap = s.repla_pregap;
            end
            local_eleinf(1:local_id-1) = [];

            % Extract senatus
            senatus_ids = find(~cellfun(@isempty,regexp([mergedeleinf(:).songname].','CON|SPE'))); % ids of CON/SPE songs in all_eleinf
            senatus_eleinf = mergedeleinf(senatus_ids);

            % assemble the eleinf to continuous y values
            collect = [];
            for k = 1: length(local_eleinf)
                collect = [collect;zeros([round(local_eleinf(k).pregap*local_eleinf(k).fs),1]);local_eleinf(k).y];
            end

            %dir/path issues
            dir = sprintf('%s/Detail-%s-%02u',s.outdir,targetsongname,local_id);
            mkdir(dir);

            % write nearby-replace songs
            dijige = find([senatus_eleinf.uniqueid].'== number);
            nearby_ids = Stimuli.findnearby(senatus_eleinf,s.num_of_near,dijige-1); % !!! might be wrong
            for h = 1: length(nearby_ids)

                rp_id = nearby_ids(h);  % mergedeleinf id for replacement
                to_be_summed = senatus_eleinf(rp_id).y;
                summedy = [to_be_summed;collect];
                % figure; Draw.spec(summedy,32000)
                audiowrite(sprintf('%s/Repla-%s-%02u-before-%s-%02u-gapis-%s.wav',dir,mergedeleinf(rp_id).songname,mergedeleinf(rp_id).fragid,targetsongname,local_id,num2str(mergedeleinf(number).pregap)),summedy,s.fs);

            end

            %             % draw to verify
            %             s.draw_scatter([senatus_eleinf(nearbyids).uniqueid].',number);
            %             title('Repla_near_from_senatus');
            %             saveas(gcf,sprintf('%s/Repla-near-from-senatus.png',dir));
            %             close(gcf)

        end

        function Deprecated_writeRepla_far_from_all(s,number) % function to generate catego_replace song with the new method
            dbstop if error

            % Extract the samesong_eleinf after the target_index ele
            mergedeleinf = s.prepro;
            [local_eleinf,local_id] = Stimuli.samesong(mergedeleinf,number);
            targetsongname = mergedeleinf(number).songname;
            if ~isempty(s.repla_pregap)
                local_eleinf(local_id).pregap = s.repla_pregap;
            end
            local_eleinf(1:local_id-1) = [];

            % assemble the eleinf to continuous y values
            collect = [];
            for k = 1: length(local_eleinf)
                collect = [collect;zeros([round(local_eleinf(k).pregap*local_eleinf(k).fs),1]);local_eleinf(k).y];
            end

            %dir/path issues
            dir = sprintf('%s/Detail-%s-%02u',s.outdir,targetsongname,local_id);
            mkdir(dir);

            % write nearby-replace songs
            rng(1);
            random_ids = randperm (numel( mergedeleinf));
            far_ids = random_ids(1:s.num_of_far);

            for h = 1: length(far_ids)
                rp_id = far_ids(h);  % mergedeleinf id for replacement
                to_be_summed = mergedeleinf(rp_id).y;
                summedy = [to_be_summed;collect];
                % figure; Draw.spec(summedy,32000)
                audiowrite(sprintf('%s/Repla-%s-%02u-before-%s-%02u-gapis-%s.wav',dir,mergedeleinf(rp_id).songname,mergedeleinf(rp_id).fragid,targetsongname,local_id,num2str(mergedeleinf(number).pregap)),summedy,s.fs);

            end

            %             % draw to verify
            %             s.draw_scatter(far_ids,number);
            %             title('Repla_far_from_all');
            %             saveas(gcf,sprintf('%s/Repla-far-from-all.png',dir));
            %             close(gcf)

        end

        function Deprecated_writeRepla_far_from_senatus(s,number)
            dbstop if error

            % Extract the samesong_eleinf after the target_index ele
            mergedeleinf = s.prepro;
            [local_eleinf,local_id] = Stimuli.samesong(mergedeleinf,number);
            targetsongname = mergedeleinf(number).songname;
            if ~isempty(s.repla_pregap)
                local_eleinf(local_id).pregap = s.repla_pregap;
            end
            local_eleinf(1:local_id-1) = [];

            % Extract senatus
            senatus_ids = find(~cellfun(@isempty,regexp([mergedeleinf(:).songname].','CON'))); % ids of CON/SPE songs in all_eleinf
            senatus_eleinf = mergedeleinf(senatus_ids);

            % if local_eleinf(1) Then set local_eleinf(1).pregap = 0.05
            if local_eleinf(1).pregap == inf
                local_eleinf(1).pregap = s.initial_pregap;
            end
            % assemble the eleinf to continuous y values
            collect = [];
            for k = 1: length(local_eleinf)
                collect = [collect;zeros([round(local_eleinf(k).pregap*local_eleinf(k).fs),1]);local_eleinf(k).y];
            end

            %dir/path issues
            dir = sprintf('%s/Detail-%s-%02u',s.outdir,targetsongname,local_id);
            mkdir(dir);



            % write nearby-replace songs
            rng(1);
            random_ids = randperm (numel(senatus_eleinf));
            far_ids = random_ids(1:s.num_of_far);
            disp('Far_IDS_Of_Repla');
            disp([senatus_eleinf(far_ids).songname]);

            for h = 1: length(far_ids)
                rp_id = far_ids(h);  % mergedeleinf id for replacement
                to_be_summed = senatus_eleinf(rp_id).y;
                summedy = [to_be_summed;collect];
                % figure; Draw.spec(summedy,32000)
                audiowrite(sprintf('%s/Repla-%s-%02u-before-%s-%02u-gapis-%s.wav',dir,senatus_eleinf(rp_id).songname,senatus_eleinf(rp_id).fragid,targetsongname,local_id,num2str(senatus_eleinf(local_id).pregap)),summedy,s.fs);
            end

            %             %draw to verify
            %             s.draw_scatter([senatus_eleinf(far_ids).uniqueid].',number);
            %             title('Repla_far_from_senatus');
            %             saveas(gcf,sprintf('%s/Repla-far-from-senatus.png',dir));
            %             close(gcf)

        end




        function Deprecated_with_sampling_writeRepla_near_from_all(s,global_id,edge_converted_from_senatus,mode)
            dbstop if error

            % Extract the samesong_eleinf after(include) the target_index ele
            mergedeleinf = s.prepro;
            [local_eleinf,local_id] = Stimuli.samesong(mergedeleinf,global_id);
            targetsongname = mergedeleinf(global_id).songname;
            if ~isempty(s.repla_pregap)
                local_eleinf(local_id).pregap = s.repla_pregap;
            end

            local_eleinf(1:local_id-1) = []; % 把local_id之前的都删掉

            % assemble the eleinf to continuous y values
            collect = [];
            for k = 1: length(local_eleinf)
                collect = [collect;zeros([round(local_eleinf(k).pregap*local_eleinf(k).fs),1]);local_eleinf(k).y];
            end

            %dir/path issues
            dir = sprintf('%s/Detail-%s-%02u',s.outdir,targetsongname,local_id);
            mkdir(dir);

            % sampling, but this might not be good
            nearby_ids = Stimuli.findnearby(mergedeleinf,edge_converted_from_senatus,global_id-1);

            if exist('mode','var')
                if strcmp(mode,'random')

                    if length(nearbyids) >= s.num_of_near
                        rng(1);
                        sampled_nearby_ids = nearby_ids(sort(randsample(length(nearby_ids),s.num_of_near),'ascend'));
                    else
                        rng(1);
                        sampled_nearby_ids = nearby_ids(sort(randsample(length(nearby_ids),length(nearby_ids)),'ascend'));
                    end

                elseif strcmp(mode,'descend')

                    sampled_nearby_ids = nearby_ids(1:s.num_of_near); % which means the most nearest element will be used

                end
            else % default mode is random
                if length(nearbyids) >= s.num_of_near
                    rng(1);
                    sampled_nearby_ids = nearby_ids(sort(randsample(length(nearby_ids),s.num_of_near),'ascend'));
                else
                    rng(1);
                    sampled_nearby_ids = nearby_ids(sort(randsample(length(nearby_ids),length(nearby_ids)),'ascend'));
                end

            end

            %sampled_nearby_ids = sort(randsample(length(nearby_ids),s.num_of_near),'ascend');


            % write nearby-replace songs
            for h = 1: length(sampled_nearby_ids )

                rp_id = sampled_nearby_ids (h);  % mergedeleinf id for replacement
                to_be_summed = mergedeleinf(rp_id).y;
                summedy = [to_be_summed;collect];
                % figure; Draw.spec(summedy,32000)
                audiowrite(sprintf('%s/Repla-%s-%02u-before-%s-%02u-gapis-%s.wav',dir,mergedeleinf(rp_id).songname,mergedeleinf(rp_id).fragid,targetsongname,local_id,num2str(mergedeleinf(global_id).pregap)),summedy,s.fs);

            end

            %             %draw to verify
            %             s.draw_scatter(nearby_ids,global_id);
            %             title('Repla_near_from_all');
            %             saveas(gcf,sprintf('%s/Repla-near-from-all.png',dir));
            %             close(gcf)
        end

        function Deprecated_with_sampling_writeFrag_near_from_all(s,global_id,edge_converted_from_senatus,mode)
            % near---same catego eles
            % all---all eles in the global_eleinf is used

            % find global index of nearby elements
            global_eleinf = s.prepro;
            nearbyids = Stimuli.findnearby(global_eleinf,edge_converted_from_senatus,global_id);

            % sampling, but this might not be good
            if exist('mode','var')
                if strcmp(mode,'random')

                    if length(nearbyids) >= s.num_of_near
                        rng(1);
                        sampled_nearby_ids = nearbyids(sort(randsample(length(nearbyids),s.num_of_near),'ascend'));
                    else
                        rng(1);
                        sampled_nearby_ids = nearbyids(sort(randsample(length(nearbyids),length(nearbyids)),'ascend'));
                    end

                elseif strcmp(mode,'descend')

                    if length(nearbyids) >= s.num_of_near
                        sampled_nearby_ids = nearbyids(1:s.num_of_near);
                    else
                        sampled_nearby_ids = nearbyids(1:length(nearbyids));
                    end

                end
            else % default mode is random
                if length(nearbyids) >= s.num_of_near
                    rng(1);
                    sampled_nearby_ids = nearbyids(sort(randsample(length(nearbyids),s.num_of_near),'ascend'));
                else
                    rng(1);
                    sampled_nearby_ids = nearbyids(sort(randsample(length(nearbyids),length(nearbyids)),'ascend'));
                end
            end



            % define path and name
            songname = global_eleinf(global_id).songname;
            local_fragid = global_eleinf(global_id).fragid;
            dir = sprintf('%s/Detail-%s-%02u',s.outdir,songname,local_fragid);
            mkdir(dir);

            % construct nearbyinf and write wavs
            nearbyinf = global_eleinf(sampled_nearby_ids);

            for xx = 1: length(nearbyinf)
                audiowrite(sprintf('%s/Frag-%s-%02u.wav',dir,nearbyinf(xx).songname,nearbyinf(xx).fragid), nearbyinf(xx).y,s.fs);
            end

            %             %draw to verify
            %             s.draw_scatter(nearbyids,nearbyids);
            %             title('Frag_near_from_all');
            %             saveas(gcf,sprintf('%s/Frag-near-from-all.png',dir));
            %             close(gcf)
        end

        function Deprecated_writeEachSongFrag(s) % in each folder

            source = Stimuli.split(s.prepro);

            wb = waitbar(0,'In Progress');
            for m = 1: length(source)

                %wait bar
                local_eleinf = source{m};
                waitbar(m/length(source),wb,sprintf('In Progress...%s',local_eleinf(1).songname));

                %write
                dir = sprintf('%s/EachSongFrags-%s',s.outdir,local_eleinf(1).songname);
                mkdir(dir);
                for k = 1: length(local_eleinf)
                    audiowrite(sprintf('%s/Frag-%s-%02u.wav',dir,local_eleinf(k).songname,local_eleinf(k).fragid), local_eleinf(k).y,s.fs);
                end
                disp(string(unique(cellstr({source{m}.songname}.'))));

            end
            close(wb);

        end

        function Deprecated_writeFlippedEachSongFrag(s)
            source = Stimuli.split(s.prepro);

            wb = waitbar(0,'In Progress');
            for m = 1: length(source)

                %wait bar
                local_eleinf = source{m};
                waitbar(m/length(source),wb,sprintf('In Progress...%s',local_eleinf(1).songname));

                %write
                dir = sprintf('%s/EachSongFrags-%s',s.outdir,local_eleinf(1).songname);
                mkdir(dir);
                for k = 1: length(local_eleinf)
                    % Garf means fliped fragments
                    audiowrite(sprintf('%s/Garf-%s-%02u.wav',dir,local_eleinf(k).songname,local_eleinf(k).fragid), flip(local_eleinf(k).y),s.fs);
                end
                disp(string(unique(cellstr({source{m}.songname}.'))));

            end
            close(wb);

        end

        function Deprecated_writeEachSongFragInOneFolder(s) % in each folder

            source = Stimuli.split(s.prepro);

            wb = waitbar(0,'In Progress');
            for m = 1: length(source)

                %wait bar
                local_eleinf = source{m};
                waitbar(m/length(source),wb,sprintf('In Progress...%s',local_eleinf(1).songname));

                %write
                dir = sprintf('%s/SongFragsFromAllSongs',s.outdir);
                mkdir(dir);
                for k = 1: length(local_eleinf)
                    audiowrite(sprintf('%s/Frag-%s-%02u.wav',dir,local_eleinf(k).songname,local_eleinf(k).fragid), local_eleinf(k).y,s.fs);
                end
                disp(string(unique(cellstr({source{m}.songname}.'))));

            end
            close(wb);

        end


        function Deprecated_writeSimilarFragsButDifferentMotif(s)
            mergedeleinf = s.prepro;

            conids = find(~cellfun(@isempty, regexp([mergedeleinf.songname].','CON')));

            for k = 1: length(conids)
                target_norm_name = mergedeleinf(conids(k)).normed_name;
                target_catego = mergedeleinf(conids(k)).catego;

                samesong = find(strcmp ([mergedeleinf.normed_name].',target_norm_name));
                samecatego =  find([mergedeleinf.catego].'==target_catego);
                samefrag = setdiff(intersect(samesong,samecatego),conids(k));

                subfoldername = sprintf('%s/SameFragDiffMotif-%s-%02u',s.outdir, target_norm_name , target_catego);

                mkdir(subfoldername);
                for w = 1:length(samefrag)
                    audiowrite(sprintf('%s//%s-%02u.wav',subfoldername,mergedeleinf(samefrag(w)).songname,w),mergedeleinf(samefrag(w)).y,mergedeleinf(samefrag(w)).fs);
                end

                % I guess that there might be something wrong with the
                % code, as usually I will only get two samefragdiffmoitif
                % songs !!!!!!!! By Zhehao Cheng @ 11.16

            end

        end



        function s = Frozen_set_near(s,num_of_near) % number of near_Stimuli to generate
            s.num_of_near = num_of_near;
        end

        function s = Frozen_set_far(s,num_of_far) % number of far_Stimuli to generate
            s.num_of_far = num_of_far;
        end

        function Deprecated_draw_complex_scatter(s,number)
            % 把频谱图按照其坐标展现在二维空间中, old version, suitable for small datasets

            mergedeleinf = s.prepro;

            coor_1 = [mergedeleinf.coor_1].'; coor_2 = [mergedeleinf.coor_2].';
            coordinate = horzcat(coor_1,coor_2);


            for indi = 1: length(coordinate)

                distance = [];
                for a = 1: length(coordinate)
                    distance(a) = norm(coordinate(a,:) - coordinate(indi,:));
                end
                [B,I] = mink(distance,16);
                mergedeleinf(indi).closest10 = I(2:16);

            end



            for h0 = 1: length(mergedeleinf)
                mergedeleinf(h0).imagehandle = 0; % 2 means braod
            end

            whole_ids = randperm (numel(mergedeleinf));
            broad_ids = whole_ids(1:45);
            for h = broad_ids
                mergedeleinf(h).imagehandle = 2; % 2 means braod
            end

            close_ids = mergedeleinf(number).closest10;
            for hh = close_ids
                mergedeleinf(hh).imagehandle = 1; % 2 means braod
            end

            mergedeleinf(number).imagehandle = -1;

            cmap1 = jet;
            cmap2 = summer;
            cmapm1 = winter;


            figure
            set(gca,'Color','w')
            hold on
            scale = 30; % scale = 50;
            axis([min([mergedeleinf.coor_1]*scale) max([mergedeleinf.coor_1]*scale) min([mergedeleinf.coor_2]*scale) max([mergedeleinf.coor_2]*scale)])

            for w = 1:length(mergedeleinf)%1: 100

                img = Cal.spec(mergedeleinf(w).y,mergedeleinf(w).fs);
                [alpha,beta] = size(img);
                img = imresize(flip(img,1),[alpha,beta*2]);

                try
                    if length(mergedeleinf(w).y)  > 50

                        if mergedeleinf(w).imagehandle == 0

                            h = subimage(mergedeleinf(w).coor_1*scale,mergedeleinf(w).coor_2*scale,img*1);
                            set(h, 'AlphaData', 0.7);
                        elseif mergedeleinf(w).imagehandle == 1

                            h = subimage(mergedeleinf(w).coor_1*scale,mergedeleinf(w).coor_2*scale,img*10,cmap1(1:12:end,:));
                            set(h, 'AlphaData', 0.7);

                        elseif mergedeleinf(w).imagehandle == 2
                            h = subimage(mergedeleinf(w).coor_1*scale,mergedeleinf(w).coor_2*scale,img*10,cmap2(1:12:end,:));
                            set(h, 'AlphaData', 0.7);

                        elseif mergedeleinf(w).imagehandle == -1
                            h = subimage(mergedeleinf(w).coor_1*scale,mergedeleinf(w).coor_2*scale,img*10,cmapm1(1:12:end,:));
                            set(h, 'AlphaData', 0.7);
                        end
                        %colormap(ax,'jet')
                        %scatter(all_eleinf(w).xfake,all_eleinf(w).yfake);
                        %drawnow
                    end
                catch error
                end

            end
        end % plot spectrogram in scatter space


        function Deprecated_draw_Stimuli_space_unique_ele_only(s,number)

            mergedeleinf = s.prepro;

            % here I need to do some screening to keep only unique
            % elements!!!!!!!!!!!!!!!!!!!!!!!

            % split element based on their catego information

            % how to do that??? by concatenate the songname with the
            % catgeo
            concat = {};
            for biden = 1: length(mergedeleinf)
                concat{biden} = sprintf('%s-%u',mergedeleinf(biden).songname, mergedeleinf(biden).catego);
            end

            [~,ia,~] = unique(concat);

            mergedeleinf = mergedeleinf(ia); % this is to screen the mergedeleinf so only those unique will maintain

            coor_1 = [mergedeleinf.coor_1].'; coor_2 = [mergedeleinf.coor_2].';
            coordinate = horzcat(coor_1,coor_2);


            for indi = 1: length(coordinate)

                distance = [];

                for a = 1: length(coordinate)
                    distance(a) = norm(coordinate(a,:) - coordinate(indi,:));
                end

                [B,I] = mink(distance,16);

                mergedeleinf(indi).closest10 = I(2:16);

            end



            for h0 = 1: length(mergedeleinf)
                mergedeleinf(h0).imagehandle = 0; % 2 means braod
            end

            whole_ids = randperm (numel(mergedeleinf));
            broad_ids = whole_ids(1:45);
            for h = broad_ids
                mergedeleinf(h).imagehandle = 2; % 2 means braod
            end

            close_ids = mergedeleinf(number).closest10;
            for hh = close_ids
                mergedeleinf(hh).imagehandle = 1; % 2 means braod
            end

            mergedeleinf(number).imagehandle = -1;

            figure;
            hold on

            for f = 1: length(mergedeleinf)

                if mergedeleinf(f).imagehandle == 0
                    scatter(mergedeleinf(f).coor_1,mergedeleinf(f).coor_2,[],'k','filled');
                elseif mergedeleinf(f).imagehandle == 2
                    scatter(mergedeleinf(f).coor_1,mergedeleinf(f).coor_2,[],'b','filled');
                elseif mergedeleinf(f).imagehandle == 1
                    scatter(mergedeleinf(f).coor_1,mergedeleinf(f).coor_2,[],'r','filled');
                elseif mergedeleinf(f).imagehandle == -1
                    scatter(mergedeleinf(f).coor_1,mergedeleinf(f).coor_2,[],'g','filled');
                end
                %drawnow
            end


            % draw spectrum in space
            cmap1 = jet;
            cmap2 = summer;
            cmapm1 = winter;


            figure
            set(gca,'Color','w')
            hold on
            scale = 30;
            axis([min([mergedeleinf.coor_1]*scale) max([mergedeleinf.coor_1]*scale) min([mergedeleinf.coor_2]*scale) max([mergedeleinf.coor_2]*scale)])

            for w = 1: 100

                img = Cal.spec(mergedeleinf(w).y,mergedeleinf(w).fs);
                [alpha,beta] = size(img);
                img = imresize(flip(img,1),[alpha,beta*2]);

                try
                    if length(mergedeleinf(w).y)  > 50

                        if mergedeleinf(w).imagehandle == 0

                            h = subimage(mergedeleinf(w).coor_1*scale,mergedeleinf(w).coor_2*scale,img*1);
                            set(h, 'AlphaData', 0.7);
                        elseif mergedeleinf(w).imagehandle == 1

                            h = subimage(mergedeleinf(w).coor_1*scale,mergedeleinf(w).coor_2*scale,img*10,cmap1(1:12:end,:));
                            set(h, 'AlphaData', 0.7);

                        elseif mergedeleinf(w).imagehandle == 2
                            h = subimage(mergedeleinf(w).coor_1*scale,mergedeleinf(w).coor_2*scale,img*10,cmap2(1:12:end,:));
                            set(h, 'AlphaData', 0.7);

                        elseif mergedeleinf(w).imagehandle == -1
                            h = subimage(mergedeleinf(w).coor_1*scale,mergedeleinf(w).coor_2*scale,img*10,cmapm1(1:12:end,:));
                            set(h, 'AlphaData', 0.7);
                        end
                        %colormap(ax,'jet')
                        %scatter(all_eleinf(w).xfake,all_eleinf(w).yfake);
                        %drawnow
                    end
                catch error
                end

            end



        end

        function Deprecated_write_syl_unique_only_con_first(s,number)

            % It is necessary to find a way to restrict the same-catego
            % elements to mainly come from those CONs presented previously
            % rather than those not presented.

            % Besides, when trying to find the element clusters, it is
            % still important to include as much elements as possible to
            % get the clear structure/border of element space.
            %
            %s.outdir % Necessary

            mergedeleinf = unique_it(s.prepro);

            closeids = Stimuli.findclose(mergedeleinf,15,number);

            songname = mergedeleinf(number).songname;
            target_fragid = mergedeleinf(number).fragid;
            newinf = mergedeleinf(closeids);
            dir = sprintf('%s/Second-%s-%u',s.outdir,songname,target_fragid);
            mkdir(dir);

            for xx = 1: length(newinf)
                audiowrite(sprintf('%s/OteSyl-%s-%u.wav',dir,newinf(xx).songname,newinf(xx).fragid), newinf(xx).y,s.fs);
            end

        end

        function Deprecated_write_repla_unique_only_con_first(s,number)

            % It is necessary to find a way to restrict the same-catego
            % elements to mainly come from those CONs presented previously
            % rather than those not presented.

            % Besides, when trying to find the element clusters, it is
            % still important to include as much elements as possible to
            % get the clear structure/border of element space.
            %
            %s.outdir % Necessary

            dbstop if error
            gap = 0.010;
            s.find_close_10;
            mergedeleinf = unique_it(s.prepro);
            [samesongeles,targetfragid] = Stimuli.samesong(mergedeleinf,number);
            targetsongname = mergedeleinf(number).songname;
            samesongeles(targetfragid).pregap = gap;
            samesongeles(1:targetfragid) = [];
            how_many_close = 15; %%% edit here ! how many replacements in each catego will take place

            collect = [];

            for k = 1: length(samesongeles)
                collect = [collect;zeros([round(samesongeles(k).pregap*samesongeles(k).fs),1]);samesongeles(k).y];
            end


            dir = sprintf('%s/Second-%s-%u',s.outdir,targetsongname,targetfragid);
            mkdir(dir);

            close_points = mergedeleinf(number).closest10;

            for h = 1: how_many_close
                rp_id = Stimuli.findclose(mergedeleinf,how_many_close,number) % mergedeleinf id for replacement
                to_be_summed = mergedeleinf(rp_id).y;

                summedy = [to_be_summed;collect];
                % figure; Draw.spec(summedy,32000)
                audiowrite(sprintf('%s/catego-%s-%u-before-%s-%u-%s.wav',dir,mergedeleinf(rp_id).songname,mergedeleinf(rp_id).fragid,targetsongname,targetfragid,num2str(gap)),summedy,s.fs);

                % audiowrite(sprintf('%s/catego-%u-%s-%s.wav',dir,catego,string(unique(cellstr({thiscatego(jj).songname}.'))),num2str(gap)),summedy,s.fs);
            end


        end


    end

    methods(Static) %Frozen or Deprecated




        function sylinf = Deprecated_syl(dir) % this function Convert songs from a folder into a single syl struct

            names = Extract.filename(dir,'*.wav');
            sylinf = struct;
            count = 0;
            for idx = 1: length(names)
                [y,fs] = audioread(names{idx});

                [~,soundname,~] = fileparts(names{idx});

                info = segment(y,fs).seg3;

                for k = 1: length(info)
                    count = count + 1;
                    sylinf(count).songname = soundname;
                    sylinf(count).y = y(info(k).initial:info(k).terminal);
                    sylinf(count).initial = info(k).initial;
                    sylinf(count).terminal = info(k).terminal;
                    if k ~= 1
                        sylinf(count).pregap = (info(k).initial - info(k-1).terminal)/fs;
                    else
                        sylinf(count).pregap = info; % pre-gap duration
                    end
                    sylinf(count).dur = length(sylinf(count).y)/fs;

                end

            end
        end

        function Deprecated_writeSingleWNS(length_in_seconds)
            data = randn(length_in_seconds*32000,1);
            raw_rms = rms(data);
            normalized_data = data*0.05/raw_rms;
            audiowrite('norm-WNS-rms0.05.wav',normalized_data,32000);
            disp(rms(normalized_data));


        end


        function sylinf = Deprecated_addpregap(sylinf)
            for k = 1: length(sylinf)
                if k ~= 1
                    if strcmp(sylinf(k-1).songname,sylinf(k).songname)
                        sylinf(k).pregap = (sylinf(k).initial - sylinf(k-1).terminal)/sylinf(k).fs;
                    else
                        sylinf(k).pregap = inf;
                    end
                else
                    sylinf(k).pregap = inf;   %%%%%% This is a super stupid solution !!!!!!!!!!!!!!!!!!!!!!!!!
                end
            end
        end

        function newinfo = Deprecated_highpass(fraginf,hpf) % hpf is the high pass frequency
            % this normalize all element/syllable to have the same rms
            if ~exist('hpf','var')
                hpf = 450;
            end
            % old = s.splited;
            splited = Stimuli.split(fraginf);
            newsplited = {};

            parfor k = 1: length(splited) % highpass is very slow
                newsplited{k} = TempStimcore.highpass(splited{k},hpf);
            end

            newinfo = Stimuli.merge(newsplited);


        end

        function newinfo = Deprecated_normalize(fraginf,trms)  % this normalize works for every song

            if ~exist('trms','var')
                trms = 0.05;
            end
            %old = s.splited;
            splited = Stimuli.split(fraginf);
            newsplited = {};

            for k = 1: length(splited )
                newsplited{k} = TempStimcore.normalize(splited {k},trms);
            end

            newinfo = Stimuli.merge(newsplited);

        end


        function splited = Deprecated_split(fraginf) % input is the s.sylinf, output is the sylinf for each each song
            dbstop if error
            songs = unique(cellstr({fraginf(:).songname}.'));

            parfor kk = 1: length(songs)  % I change for to par-for 1103
                idxs = find( strcmp(cellstr({fraginf(:).songname}.') ,songs{kk}) );
                if length(idxs) == 1
                    thisT = struct2table(fraginf(idxs),'AsArray',1);
                else
                    thisT =  struct2table(fraginf(idxs));
                end

                thisT = sortrows(thisT,'fragid');  % 这里可能会改   order the struct based on the syllable number
                splited{kk} = table2struct(thisT);

            end

        end

        function newinf = Deprecated_merge(splited)  % merge the splited and processed info back
            newinf = vertcat(splited{:})';
        end

        function summed = Deprecated_assemble(fraginf) % a function to assemble constitute together


            ys = {fraginf.y}.';
            gaps = [fraginf(:).pregap].';
            gaps(1) = 0; % replace the first pregap to 0,which is a info
            fs = unique([fraginf.fs].');

            summed = [];

            for haris = 1: length(ys)
                if gaps(haris) >= 0
                    summed = [summed;zeros(int64(gaps(haris)*fs),1);ys{haris}];
                else
                    konoy =  ys{haris};
                    summed = [summed;konoy(abs(int64(gaps(haris)*fs)) + 1:end)];  %%%%% compensate when gap is minus
                end

            end


        end

        function deg = Deprecated_degressive(fraginf, initial, terminal) % initial-end
            deg = struct;
            diff = terminal - initial + 1;
            for k = 1: diff
                rank =k - 1;
                shortensylinf = fraginf(initial+rank:end);
                deg(k).y = Stimuli.assemble(shortensylinf);
                deg(k).rank = rank;
            end

        end

        function [same,localidx] = Deprecated_samesong(fraginf,idx) % get all the syllables from the same song

            splited = Stimuli.split(fraginf);

            songnames = {};
            for k = 1: length(splited)
                songnames{k} = splited{k}(1).songname;
            end

            splitidx = find(strcmp(cellstr(songnames),fraginf(idx).songname));

            same = splited{splitidx};

            localidx = find([same.fragid].'== fraginf(idx).fragid);


        end

        function closeids = Deprecated_findnearby(mergedeleinf,how_many_close,number)
            eleinf = mergedeleinf;
            coor_1 = [eleinf.coor_1].'; coor_2 = [eleinf.coor_2].';
            coordinate = horzcat(coor_1,coor_2);

            distance = [];

            for a = 1: length(coordinate)
                distance(a) = norm(coordinate(a,:) - coordinate(number,:));
            end

            [~,I] = sort(distance,'ascend');
            %[B,I] = mink(distance,how_many_close + 1); %Another option
            closeids = I(2:how_many_close + 1); % close id 是从相似度由高到低排列的

            for trump = 1: length(eleinf)
                eleinf(trump).distance = distance(trump);
            end





            sorted_eleinf = eleinf(I);
            %try to verift this part

            %                         figure;
            %                         scatter(coordinate(:,1),coordinate(:,2),[],'k','filled');,
            %                         hold on
            %                         c = distance(closeids);
            %                         for k = 1:length(closeids)
            %                         scatter(eleinf(closeids(k)).coor_1,eleinf(closeids(k)).coor_2,[],c(k),'filled');
            %                         end
            %                         scatter(eleinf(number).coor_1,eleinf(number).coor_2,[],'r','filled','h');

        end

        function newinf = Deprecated_unique_it(mergedeleinf)
            % 每个category取一个syllable时用此方法unique一下
            % based on categories to unique the mergedeleinf
            concat = {};

            for k = 1: length(mergedeleinf)
                concat{k} = sprintf('%s-%u',mergedeleinf(k).songname,mergedeleinf(k).catego);
            end

            [~,ia,~] = unique(concat);

            newinf = mergedeleinf(ia);

        end


        function yFaded = Frozen_fadeOut(y,fs)
            % Create a fading envelope signal.
            numSamples = length(y);
            fadenum = 500;
            % Create a fading envelope
            fadeSignal = [ones(1, numSamples-fadenum), linspace(1, 0, fadenum)]';
            yFaded = y .* fadeSignal;
            % Plot fading signal.

            if exist('fs','var')
                figure
                subplot(3, 1, 1);
                Draw.spec(y,fs);
                grid on;

                subplot(3, 1, 2);
                plot(fadeSignal, '-', 'LineWidth', 2);
                %title('Fading Envelope Signal', 'FontSize', fontSize);
                grid on;
                % Apply the fading by multiplying

                % Plot fading signal.
                subplot(3, 1, 3);
                Draw.spec(yFaded,fs);
                %title('Faded Waveform', 'FontSize', fontSize);
                grid on;
                % audiowrite('2.wav' ,y,Fs);

            end

        end

        function yFaded = Frozen_fadeIn(y,fs)
            % Create a fading envelope signal.
            numSamples = length(y);
            fadenum = 500;
            % Create a fading envelope
            fadeSignal = [ linspace(1, 0, fadenum),ones(1, numSamples-fadenum)]';
            yFaded = y .* fadeSignal;
            % Plot fading signal.

            if exist('fs','var')
                figure
                subplot(3, 1, 1);
                Draw.spec(y,fs);
                grid on;

                subplot(3, 1, 2);
                plot(fadeSignal, '-', 'LineWidth', 2);
                %title('Fading Envelope Signal', 'FontSize', fontSize);
                grid on;
                % Apply the fading by multiplying

                % Plot fading signal.
                subplot(3, 1, 3);
                Draw.spec(yFaded,fs);
                %title('Faded Waveform', 'FontSize', fontSize);
                grid on;
                % audiowrite('2.wav' ,y,Fs);

            end

        end

        function sylinf = Frozen_polish(sylinf)
            %似乎仅仅是补充了fs和duration
            for k = 1: length(sylinf)
                sylinf(k).fs = 32000;
                sylinf(k).dur = length(sylinf(k).y)/sylinf(k).fs;
            end
        end

        function sylinf = Frozen_addfs(sylinf)
            for k = 1: length(sylinf)
                sylinf(k).fs = 32000;
            end
        end

        function sylinf = Frozen_adduniqueid(sylinf)
            for k = 1: length(sylinf)
                sylinf(k).uniqueid = k;
            end
        end

        function sylinf = Frozen_adddur(sylinf)
            if ~isfield(sylinf,'dur')
                parfor k = 1: length(sylinf)
                    sylinf(k).dur = length(sylinf(k).y)/sylinf(k).fs;
                end
            end
        end


        function fraginf = Frozen_fadeout_all(fraginf)  % this normalize works for every song

            for k = 1: length(fraginf)
                fraginf(k).y = Stimuli.fadeOut(fraginf(k).y);
            end

        end

        function fraginf = Frozen_fadein_all(fraginf)  % this normalize works for every song

            for k = 1: length(fraginf)
                fraginf(k).y = Stimuli.fadeIn(fraginf(k).y);
            end

        end


        function new_conresp_siginfo = Frozen_expand_siginfo(conresp_siginfo,before)

            % 这个方法原先似乎是给没有fragid信息的eleinf补上fragid的，现在貌似是没有用了
            % 2022.12.20 暂时冻结
            % 扩展之后更不容易漏掉 response-eliciting elements
            count = 0;
            extra = struct;
            for k = 1:length(conresp_siginfo)
                for a = 1:before
                    if conresp_siginfo(k).fragid -a > 0
                        count = count + 1;
                        extra(count).songname = conresp_siginfo(k).songname;
                        extra(count).songid = conresp_siginfo(k).songid;
                        extra(count).fragid = conresp_siginfo(k).fragid - a;
                    end
                end
            end

            new_conresp_siginfo = horzcat(conresp_siginfo, extra);

        end


    end

end
