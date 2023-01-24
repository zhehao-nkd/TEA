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

classdef StimCore < handle

    properties
        raw % infment information struct
        fs
        %splited % splited is a struct for each song, by modifying the splited, like delete a syllable or, replace the y to the normalize y, or even replace the syllable to other syllables, ne wStimuli will be gnerated
        %%% That is to say, function normalize should be renewed to
        %%% become a general function. Used by all other fucntions to
        %%% assemble the syllable/elemnt within songs into a single song
        responded  % it is the collection of songs that can activate a neurons
        selected11 % selected songs that can activate a neuron
        selected01    %not responded  ??????
        selectedsplited %% splited data, but not full, instead , it is selected

        prepro % preprocessed info, which is firstly highpass-filtered and then normalized
        songnames  % song names within info
        songnum  % number of songs within info
        outdir
        num_of_near
        num_of_far
        num_of_sim
        repla_pregap % pre_gap for replacement
        initial_pregap  %pre_gap of the first elements % seconds


    end

    methods  % 底层算法，不为产生wav文件

        function s = StimCore(info)
            %initialize the eleinf's uniqueid        
            s.fs = unique([info.fs].');  % a single value
            s.prepro = StimCore.getprepro(info);
            s.initial_pregap = 0.05;

            % judge whether those infos are from a single song or multiple
            % song

            s.songnames = unique(cellstr( {info.songname}.'));
            s.songnum = length(s.songnames);



        end

        function new = normalize2(s,old,tarrms)
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

        function set_num_of_sim(s,num_of_sim)
            s.num_of_sim = num_of_sim; % set the number of similar frags(just the same elements in different motifs
        end


    end

    methods  % 直接产生wav文件的方法

        function summer = getnormy(s)
            %得到 purified song data
            summer = s.assemble(s.prepro);
        end

        function motify = getmotify(s, motifindex)
            %得到 purified motif data
            thismotif = s.prepro([s.prepro.motif].' == motifindex);
            motify = s.assemble(thismotif);
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
    
    methods % Frozen or Deprecated

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

    methods(Static)



        function prepro = getprepro(info)

            for k = 1: length(info)
                info(k).dur = length(info(k).y)/info(k).fs;

                info(k).uniqueid = k;

                if k==1 || ~strcmp(info(k-1).songname,info(k).songname)
                    info(k).pregap = inf;   %%%%%% This is a super stupid solution !!!!!!!!!!!!!!!!!!!!!!!!!
                else
                    info(k).pregap = (info(k).initial - info(k-1).terminal)/info(k).fs;
                end
            end

            temp = StimCore.highpass(info,450);
            temp = StimCore.normalize(temp, 0.05);

            prepro = table2struct(sortrows(struct2table(temp),'fragid','ascend')); % re-order
            % already highpass-filtered and normalized for each song
        
        end  


        function new = highpass(fraginf,hpf) % 重要hpf is the high pass frequency
            % this normalize all element/syllable to have the same rms
            if ~exist('hpf','var')
                hpf = 450; % 似乎之前也有使用过420，需要知道哪个更适宜
            end
            % old = s.splited;
            %                 newinfo = TempStimcore.highpass(fraginf,hpf);
            %                      new = highpass(fraginf,hpf) % hpf is the high pass frequency
            % this normalize all element/syllable to have the same rms

           old = fraginf;
           new = old;

           % Before, 08.30.2022, the highpass filter was applied for each
           % fragment individually, which is not correct because it will
           % cause sharp moving noise

           % after 08,30,2022 the high-pass filter was applied for the
           % whole song

           for k = 1:length(old)
               temp = {old(1:k-1).y}.';
               old(k).firstdatapoint = length(vertcat(temp{:}))+1;
               temp = {old(1:k).y}.';
               old(k).lastdatapoint = length(vertcat(temp{:}));
           end

           sumy = {old(:).y}.';
           sumy = vertcat(sumy{:});
           hp_sumy = highpass(sumy, hpf,old(1).fs);

           for k = 1: length(old)
               new(k).y = hp_sumy(old(k).firstdatapoint:old(k).lastdatapoint);
           end


        end

        function new = normalize(fraginf,trms)  % 重要this normalize works for every song

            if ~exist('trms','var')
                trms = 0.05;
            end

            %             new= TempStimcore.normalize(fraginf,trms);
            % normalize to make the song have the same rms
            if ~exist('tarrms','var')
                tarrms = 0.05; % default high-pass filter
            end

            old = fraginf;
            new = old;
            ys = {old(:).y}.';
            sumy = vertcat(ys{:});

            sumrms = rms(sumy);
            ratio = tarrms/sumrms;
            for k = 1: length(old)
                new(k).y = ratio*old(k).y;
            end

        end

        function summed = assemble(fraginf) % 重要a function to assemble constitute together
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

        function deg = degressive(fraginf, initial, terminal) % initial-end
            deg = struct;
            diff = terminal - initial + 1;
            for k = 1: diff
                rank =k - 1;
                shortensylinf = fraginf(initial+rank:end);
                deg(k).y = Stimuli.assemble(shortensylinf);
                deg(k).rank = rank;
            end

        end

        function [same,localidx] = samesong(fraginf,idx) % get all the syllables from the same song

            splited = Stimuli.split(fraginf);

            songnames = {};
            for k = 1: length(splited)
                songnames{k} = splited{k}(1).songname;
            end

            splitidx = find(strcmp(cellstr(songnames),fraginf(idx).songname));

            same = splited{splitidx};

            localidx = find([same.fragid].'== fraginf(idx).fragid);


        end

        function closeids = findnearby(mergedeleinf,how_many_close,number)
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

        function newinf = unique_it(mergedeleinf)
            % 每个category取一个syllable时用此方法unique一下
            % based on categories to unique the mergedeleinf
            concat = {};

            for k = 1: length(mergedeleinf)
                concat{k} = sprintf('%s-%u',mergedeleinf(k).songname,mergedeleinf(k).catego);
            end

            [~,ia,~] = unique(concat);

            newinf = mergedeleinf(ia);

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
   
    end




end
