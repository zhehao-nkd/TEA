classdef Simatrix < handle
    % similarity matrix of feature-feature comparision, subclass of class Neuron

    properties
        fraglist
        features
        matrix
        yesids
        compareids
        deviations % deviation between responsive frags' feature difference and random selected frags'
        formated_name
        pvalues

    end

    methods
        function s = Simatrix(fraglist)


            s.features = struct;
            if isempty(fraglist)
                return
            end

            % all to all
            if length(find(cellfun(@isempty,{fraglist.meanfeatures}.'))) <3 % 如果小于三，姑且放过
                fraglist = fraglist(find(~cellfun(@isempty,{fraglist.meanfeatures}.'))); % 大于三就不能姑息，且让它报错
            end

            for k = 1:length(fraglist)
                s.features.stimuliname{k} = fraglist(k).stimuliname;
                thosefeatures = fieldnames(fraglist(k).meanfeatures);
                for kk = 1:length(thosefeatures)
                    eval(['s.features.',thosefeatures{kk},'(k)','= fraglist(k).meanfeatures.',thosefeatures{kk},';']);
                end
                s.features.dur = length(fraglist(k).y); % 计算input的长度，并把它作为一个feature
            end


            % calculate distance matrixes

            featurenames = setdiff(fieldnames(s.features),'stimuliname');

            for kk = 1:length(featurenames)
                eval(['s.matrix.',featurenames{kk},'= rescale(squareform( pdist(s.features.',featurenames{kk},''')));']);  % 验证过，与写一个nested循环得到的结果相同
            end


            %claculate global matrixes : accuracy and similarity
            % 目前这个太费时间了，暂时不要使用，等meeting过了有时间了以后再使用
%             for v = 1:length(fraglist)
%                 for vv = 1:length(fraglist)
%                     sim = SAT_similarity(SAT_sound(fraglist(v).y,fraglist(v).fs),SAT_sound(fraglist(vv).y,fraglist(vv).fs),0);
%                     sim.calculate_similarity;
%                     s.matrix.accuarcy(v,vv)= sim.score.accuracy;
%                     s.matrix.similarity(v,vv) = sim.score.similarity; % or accuracy
%                 end
%             end
            featurenames = fieldnames(s.matrix);   
            global_matrix = zeros(size(s.matrix.(featurenames{1}))); % 初始化
            for k = 1:length(featurenames)
                global_matrix = global_matrix + s.matrix.(featurenames{k});
            end

            s.matrix.global = rescale(global_matrix/length(featurenames));
            


            % yes ids and compare ids

            s.yesids= find([fraglist.label].' == 1);

            s.compareids = cell(120,1);
            for k = 1:120
                s.compareids{k} = randsample(1:length(fraglist),length(s.yesids));
            end


            % 下面是关于best stimuli 的初始化
            if isfield(fraglist,'maxsdf')
                s.fraglist = table2struct( sortrows(struct2table(fraglist),'maxsdf','descend') ); % 这一步确定谁是best stimuli，十分重要
                %
                beststimuli = SAT_sound(s.fraglist(1).y,s.fraglist(1).fs);

                try
                    [mfcc_beststimuli,~,~,~] = mfcc(s.fraglist(1).y,s.fraglist(1).fs);
                catch % 如果fraglist y 太短的话，在后面补零试一试，补零之后再去掉补零部分的mfcc,但因为只有一帧所以没办法补
                    [mfcc_beststimuli,~,~,~] = mfcc([s.fraglist(1).y;zeros(s.fraglist(1).fs*0.03 -length(s.fraglist(1).y),1)],fraglist(1).fs);
                end


                % For calculating similarities towards the best stimuli
                for k = 1:length(s.fraglist)
                    sim = SAT_similarity(beststimuli,SAT_sound(s.fraglist(k).y,s.fraglist(k).fs),0);
                    sim.calculate_similarity;
                    s.fraglist(k).accuracy = sim.score.accuracy;
                    s.fraglist(k).similarity = sim.score.similarity; % or accuracy
                    %calculate distance in different features

                    feature_names = setdiff(fieldnames(s.features),'stimuliname');

                    for kk = 1:length(feature_names)

                        eval(['s.fraglist(k).',feature_names{kk},' = norm(fraglist(1).meanfeatures.',feature_names{kk},' - fraglist(k).meanfeatures.',feature_names{kk},');']);

                    end


                    try
                        [mfcc_specific_stimuli,~,~,~] = mfcc(s.fraglist(k).y,s.fraglist(k).fs);
                    catch % 如果fraglist y 太短的话，在后面补零试一试，补零之后再去掉补零部分的mfcc,但因为只有一帧所以没办法补
                        [mfcc_specific_stimuli,~,~,~] = mfcc([s.fraglist(k).y;zeros(s.fraglist(k).fs*0.03 -length(s.fraglist(k).y),1)],fraglist(k).fs);
                    end

                    mfcc_matrix = [];
                    for a = 1:size(mfcc_beststimuli,1)
                        vec1 = mfcc_beststimuli(a,:);
                        for aa = 1:size(mfcc_specific_stimuli,1)
                            vec2 = mfcc_specific_stimuli(aa,:);
                            mfcc_matrix(a,aa) =  norm(vec1-vec2);
                        end
                    end
                    mfcc_matrix = 1 - rescale(mfcc_matrix); % transfer dissimilarity matrix to similarity matrix
                    fraglist(k).mfcc_meansim = mean(mean(mfcc_matrix));
                    fraglist(k).perc_highsim = length(find(mfcc_matrix > 0.78))/(size(mfcc_matrix,1)*size(mfcc_matrix,2));

                    %                 figure; imagesc(mfcc_matrix);
                end

            end


        end

        function drawBoxPlot(s)

           
            if length(s.yesids) == 0  %如果根本没有syllables which trigger significant response
                return
            end
            s.deviations = struct; % deviation between responsive frags' feature difference and random selected frags'

            yestemp = sort(reshape(s.matrix.amplitude(s.yesids,s.yesids),[],length(s.yesids).^2));
            collect_temp = {};
            for k = 1:length(s.compareids)
                comparetemp = sort(reshape(s.matrix.amplitude(s.compareids{k},s.compareids{k}),[],length(s.compareids{k}).^2));
                collect_temp{k}= yestemp - comparetemp;
            end
            s.deviations.amplitude = horzcat(collect_temp{:})/std(horzcat(collect_temp{:}),'omitnan');
            [~,s.pvalues(1)] = ttest(rmoutliers(s.deviations.amplitude),0,'Tail','left');

            yestemp = sort(reshape(s.matrix.mfa(s.yesids,s.yesids),[],length(s.yesids).^2));
            collect_temp = {};
            for k = 1:length(s.compareids)
                comparetemp = sort(reshape(s.matrix.mfa(s.compareids{k},s.compareids{k}),[],length(s.compareids{k}).^2));
                collect_temp{k}= yestemp - comparetemp;
            end
            s.deviations.mfa = horzcat(collect_temp{:})/std(horzcat(collect_temp{:}),'omitnan');
            [~,s.pvalues(2)] = ttest(rmoutliers(s.deviations.mfa),0,'Tail','left')

            yestemp = sort(reshape(s.matrix.pitch(s.yesids,s.yesids),[],length(s.yesids).^2));
            collect_temp = {};
            for k = 1:length(s.compareids)
                comparetemp = sort(reshape(s.matrix.pitch(s.compareids{k},s.compareids{k}),[],length(s.compareids{k}).^2));
                collect_temp{k}= yestemp - comparetemp;
            end
            s.deviations.pitch = horzcat(collect_temp{:})/std(horzcat(collect_temp{:}),'omitnan');
            [~,s.pvalues(3)] = ttest(rmoutliers(s.deviations.pitch),0,'Tail','left')

            yestemp = sort(reshape(s.matrix.mf(s.yesids,s.yesids),[],length(s.yesids).^2));
            collect_temp = {};
            for k = 1:length(s.compareids)
                comparetemp = sort(reshape(s.matrix.mf(s.compareids{k},s.compareids{k}),[],length(s.compareids{k}).^2));
                collect_temp{k}= yestemp - comparetemp;
            end
            s.deviations.mf = horzcat(collect_temp{:})/std(horzcat(collect_temp{:}),'omitnan');
            [~,s.pvalues(4)] = ttest(rmoutliers(s.deviations.mf),0,'Tail','left')

            yestemp = sort(reshape(s.matrix.FM(s.yesids,s.yesids),[],length(s.yesids).^2));
            collect_temp = {};
            for k = 1:length(s.compareids)
                comparetemp = sort(reshape(s.matrix.FM(s.compareids{k},s.compareids{k}),[],length(s.compareids{k}).^2));
                collect_temp{k}= yestemp - comparetemp;
            end
            s.deviations.FM = horzcat(collect_temp{:})/std(horzcat(collect_temp{:}),'omitnan');
            [~,s.pvalues(5)] = ttest(rmoutliers(s.deviations.FM),0,'Tail','left')

            yestemp = sort(reshape(s.matrix.AM(s.yesids,s.yesids),[],length(s.yesids).^2));
            collect_temp = {};
            for k = 1:length(s.compareids)
                comparetemp = sort(reshape(s.matrix.AM(s.compareids{k},s.compareids{k}),[],length(s.compareids{k}).^2));
                collect_temp{k}= yestemp - comparetemp;
            end
            s.deviations.AM = horzcat(collect_temp{:})/std(horzcat(collect_temp{:}),'omitnan');
            [~,s.pvalues(6)] = ttest(rmoutliers(s.deviations.AM),0,'Tail','left')

            yestemp = sort(reshape(s.matrix.goodness(s.yesids,s.yesids),[],length(s.yesids).^2));
            collect_temp = {};
            for k = 1:length(s.compareids)
                comparetemp = sort(reshape(s.matrix.goodness(s.compareids{k},s.compareids{k}),[],length(s.compareids{k}).^2));
                collect_temp{k}= yestemp - comparetemp;
            end
            s.deviations.goodness = horzcat(collect_temp{:})/std(horzcat(collect_temp{:}),'omitnan');
            [~,s.pvalues(7)] = ttest(rmoutliers(s.deviations.goodness),0,'Tail','left')

            yestemp = sort(reshape(s.matrix.entropy(s.yesids,s.yesids),[],length(s.yesids).^2));
            collect_temp = {};
            for k = 1:length(s.compareids)
                comparetemp = sort(reshape(s.matrix.entropy(s.compareids{k},s.compareids{k}),[],length(s.compareids{k}).^2));
                collect_temp{k}= yestemp - comparetemp;
            end
            s.deviations.entropy = horzcat(collect_temp{:})/std(horzcat(collect_temp{:}),'omitnan');
            [~,s.pvalues(8)] = ttest(rmoutliers(s.deviations.entropy),0,'Tail','left')

            yestemp = sort(reshape(s.matrix.pf(s.yesids,s.yesids),[],length(s.yesids).^2));
            collect_temp = {};
            for k = 1:length(s.compareids)
                comparetemp = sort(reshape(s.matrix.pf(s.compareids{k},s.compareids{k}),[],length(s.compareids{k}).^2));
                collect_temp{k}= yestemp - comparetemp;
            end
            s.deviations.pf = horzcat(collect_temp{:})/std(horzcat(collect_temp{:}),'omitnan');
            [~,s.pvalues(9)] = ttest(rmoutliers(s.deviations.pf),0,'Tail','left')

            yestemp = sort(reshape(s.matrix.ct(s.yesids,s.yesids),[],length(s.yesids).^2));
            collect_temp = {};
            for k = 1:length(s.compareids)
                comparetemp = sort(reshape(s.matrix.ct(s.compareids{k},s.compareids{k}),[],length(s.compareids{k}).^2));
                collect_temp{k}= yestemp - comparetemp;
            end
            s.deviations.ct = horzcat(collect_temp{:})/std(horzcat(collect_temp{:}),'omitnan');
            [~,s.pvalues(10)] = ttest(rmoutliers(s.deviations.ct),0,'Tail','left')

            yestemp = sort(reshape(s.matrix.cf(s.yesids,s.yesids),[],length(s.yesids).^2));
            collect_temp = {};
            for k = 1:length(s.compareids)
                comparetemp = sort(reshape(s.matrix.cf(s.compareids{k},s.compareids{k}),[],length(s.compareids{k}).^2));
                collect_temp{k}= yestemp - comparetemp;
            end
            s.deviations.cf = horzcat(collect_temp{:})/std(horzcat(collect_temp{:}),'omitnan');
            [~,s.pvalues(11)] = ttest(rmoutliers(s.deviations.cf),0,'Tail','left')

            allDatacollect = {};
            allCats = [];
            allFields = string(fieldnames(s.deviations));
            ttest_allFields = allFields;
            for k = 1:length(s.pvalues)
                if s.pvalues(k)>0.5
                    ttest_allFields(k) = sprintf("%s(n.s.)",allFields(k));
                elseif s.pvalues(k)>0.01 && s.pvalues(k)<=0.5
                    ttest_allFields(k) = sprintf("%s(*)",allFields(k));
                elseif s.pvalues(k)>0.001 && s.pvalues(k)<=0.01
                    ttest_allFields(k) = sprintf("%s(**)",allFields(k));
                elseif s.pvalues(k)>0.0001 &&s.pvalues(k)<=0.001
                    ttest_allFields(k) = sprintf("%s(***)",allFields(k));
                elseif s.pvalues(k)<=0.0001
                    ttest_allFields(k) = sprintf("%s(****)",allFields(k));
                end
            end
            for iField = 1:length(allFields)
                allDatacollect{iField} = s.deviations.(allFields(iField));
                allCats = [allCats,repelem(ttest_allFields(iField),length(s.deviations.(allFields(iField))));]
            end
            allData = vertcat(allDatacollect{:});
            figure('Position',[1906 -46 946 1217],'Color','w');
            hold on
            boxplot(allData,allCats,'Symbol','');
            ylim([-2.5,2.5])
            for k = 1:length(allData(:,1))
                plot(k,mean(rmoutliers(allData(1,:))), 'dk')
            end
            set(gca,'FontSize',20)
            set(findobj(gca,'type','line'),'linew',2)
            yline(0,':k','LineWidth',1.5)
            saveas(gcf,sprintf('BoxPlot_%s.png',s.formated_name));
            close(gcf)


        end


        function Relation_DistBestStimuli_NeuralResponse(s)

            if isempty(s.fraglist)
                return
            end

            ranking_by = {'accuracy','similarity','amplitude','mfa','pitch','mf','AM','goodness','entropy','pf','ct','cf'};
            figure('Color','w');
            hold on

            notordered_labels = [s.fraglist.label].';
            num_of_random = 100;
            for ww = 1:num_of_random
                notreordered = notordered_labels(randperm(length(s.fraglist)));
                pers = getPercentageOfResponse(notreordered);
                hold on
                plot(pers,'Color',[0.42,0.42,0.42]);
            end

            cmap = colormap('lines');
            for k = 1:length(ranking_by)
                local_fraglist = table2struct( sortrows(struct2table(s.fraglist),ranking_by{k}) ); % by what kind of similarity to the best stimuli
                ordered_labels = [local_fraglist.label].';
                pers_ordered = getPercentageOfResponse(ordered_labels);
                plot (pers_ordered,'Color',cmap(k,:));
            end

            positivecontrol = sort(ordered_labels, 'descend' ) ;
            pers_pcontrol = getPercentageOfResponse(positivecontrol);
            plot (pers_pcontrol,'Color','k') ;
            hold off
            title(sprintf('%s.png',s.formated_name));
            saveas(gcf,sprintf('PercenExplained_%s.png',s.formated_name));
            close(gcf);
            function pers = getPercentageOfResponse(labels)
                for dd= 1: length(labels)
                    pers(dd) = length(find(find(labels==1) <=dd))/length(find(labels==1));
                end
            end

        end
        

        function tempDrawGlobalMatrix(s)



            
            [~,neworder] = sortrows(struct2table(s.fraglist),{'label','maxsdf'});
        

            figure('Color','w','Position',[680 393 873 705]);  % 暂时性的，后期要删掉的
            hold on
            newmatrix = Simatrix.reorderMatrix(temp_globalmatrix,neworder);
            newmatrix = 1 - rescale(newmatrix,0,1);
            imagesc(newmatrix);
            xlabel(sprintf('Global similarity'),'Fontsize',23);

            bound = length(find([s.fraglist.label].' == 1));

            plot([1-0.5 bound+0.5], [1-0.5 1-0.5], 'r', 'LineWidth', 2);
            plot( [1-0.5 1-0.5], [1-0.5 bound+0.5],'r', 'LineWidth', 2);
            plot(  [bound+0.5 bound+0.5],[1-0.5 bound+0.5],'r', 'LineWidth', 2);
            plot(  [1-0.5 bound+0.5],[bound+0.5 bound+0.5],'r', 'LineWidth', 2);

            set(gca,'XTick',[],'YTick',[]);

           
            frame = getframe(gcf);
            I = frame.cdata;
            close(gcf);
            imwrite(I,sprintf('郑GlobalSimilarityMatrix_%s.png',s.formated_name));


        end

        function drawMatrix(s)% draw those similarity matrixes
             % 注意matrix是row distance matrix，作图需要的是similarity matrix

            [~,neworder] = sortrows(struct2table(s.fraglist),{'label','maxsdf'});

             names = {'amplitude','pitch','AM','FM','mf','pf','entropy','goodness','mfa','ct','cf'};
            
            I = {}; % inatialize the Image collection
            f = waitbar(0,'Start');
            
            for na = 1: length(names)
                waitbar(na/length(names),f,replace(sprintf('Generating %s similarity matrix ...',names{na}),'_','-'));
                figure('Color','w','Position',[680 531 715 567]);
                hold on
                newmatrix = Simatrix.reorderMatrix(s.matrix.(names{na}),neworder);
                newmatrix = 1 - rescale(newmatrix,0,1);
                imagesc(newmatrix);
                xlabel(sprintf('Similarity in %s',names{na}),'Fontsize',20); % ordered by maxsdf
                bound = length(find([s.fraglist.label].' == 1));

                plot([1-0.5 bound+0.5], [1-0.5 1-0.5], 'r', 'LineWidth', 2);
                plot( [1-0.5 1-0.5], [1-0.5 bound+0.5],'r', 'LineWidth', 2);
                plot(  [bound+0.5 bound+0.5],[1-0.5 bound+0.5],'r', 'LineWidth', 2);
                plot(  [1-0.5 bound+0.5],[bound+0.5 bound+0.5],'r', 'LineWidth', 2);

                set(gca,'XTick',[],'YTick',[]);
                frame = getframe(gcf);
                I{na} = frame.cdata;
                close(gcf);
            end
            close(f);


            lieshu = 4;
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
            IMG = cell2mat(reshapedI);
            imwrite(IMG,sprintf('鲁SimilarityMatrix_%s.png',s.formated_name));

        end
    
   
    
    end

    methods(Static)

        function newmatrix = reorderMatrix(oldmatrix, neworder)
            % 把一个矩阵，比如相似度矩阵重新排序

            newmatrix = zeros(size(oldmatrix));

            for k = 1:length(oldmatrix)

                for kk = 1:length(oldmatrix)
                    newmatrix(neworder(k),neworder(kk)) = oldmatrix(k,kk);
                end

            end



        end

    end
end