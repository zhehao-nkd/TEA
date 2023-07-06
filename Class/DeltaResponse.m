classdef DeltaResponse
    %R(A+B) - RA - RB

    properties
        perstimulus
        neuronname
        perneuron
        percatego
    end

    methods
        function dr = DeltaResponse(frag,repla,consistency,newcategolist,neuronname)

            allfraglist = frag.allfraglist;
            hitted_repla = repla;
            allreplalist = hitted_repla.replalist{1};
            hitted_consistency = consistency;
            temptargets = hitted_repla.knowTarget;
            targets = regexp(convertCharsToStrings(temptargets),'[OGBYR]\d+-\d+','match');
            temp_responseB =  hitted_consistency.sublist(intersect(find(~cellfun(@isempty, regexp(cellstr({ hitted_consistency.sublist.stimuliname}.'),targets))),...
                find(~cellfun(@isempty, regexp(cellstr({ hitted_consistency.sublist.stimuliname}.'),'frag|Frag'))) ) );
            responseB = temp_responseB(end);

            %deltaresponse = DeltaResponse(responseB,allfraglist,allreplalist, newcategolist,neuronname);
            perstimulus = struct;
            for k = 1:length(allreplalist)
                parts = strsplit(allreplalist(k).stimuliname,'before');
                allreplalist(k).part1 = parts{1};
            end

            for kk = 1:length(newcategolist)
                hid_frag = intersect( find(~cellfun(@isempty, regexp(cellstr({allfraglist.stimuliname}.'),newcategolist(kk).fullname))),...
                    find(~cellfun(@isempty, regexp(cellstr({allfraglist.stimuliname}.'),'Type'))) );
                hid_repla = find(~cellfun(@isempty, regexp(cellstr({allreplalist.part1}.'),newcategolist(kk).fullname)));

                responseAB = allreplalist(hid_repla);
                responseA = allfraglist(hid_frag);

                [istart,istop,dist] = findsignal(responseAB.rawy, responseB.y);


                dur = istart/responseAB.fs - responseAB.zpt; % 第二个syllable的onset和第一个syllable的onset之间的距离
                % dur = responseAB.onsetB/responseAB.fs - responseAB.zpt

                tocompare_sylB_sptimes = Extract.sptimes_resetSP(responseB.rawsptimes, responseB.zpt - dur, responseB.zpt - dur + 0.5);
                tocompare_sylB_presptimes = Extract.sptimes_resetSP(responseB.rawsptimes, responseB.zpt - dur-0.5, responseB.zpt - dur);

                frag_zpt = responseA.zpt;
                tocompare_sylA_sptimes = Extract.sptimes_resetSP(responseA.rawsptimes, frag_zpt, frag_zpt + 0.5); % 取500ms的间距
                tocompare_sylA_presptimes = Extract.sptimes_resetSP(responseA.rawsptimes, frag_zpt -0.5, frag_zpt ); % 取500ms的间距

                sylA300_sptimes = Extract.sptimes_resetSP(responseA.rawsptimes, frag_zpt, frag_zpt + 0.3); % 取500ms的间距
                sylA300_presptimes = Extract.sptimes_resetSP(responseA.rawsptimes, frag_zpt -0.3, frag_zpt ); % 取500ms的间距


                repla_zpt = responseAB.zpt;
                sptimes_AB = Extract.sptimes_resetSP(responseAB.rawsptimes, repla_zpt , repla_zpt + 0.5); % 取500ms的间距
                presptimes_AB = Extract.sptimes_resetSP(responseAB.rawsptimes, repla_zpt -0.5, repla_zpt ); % 取500ms的间距

                perstimulus(kk).neuronname = neuronname;
                perstimulus(kk).fragname = newcategolist(kk).fullname;
                perstimulus(kk).RS_B =  length(vertcat(tocompare_sylB_sptimes{:})) - length(vertcat(tocompare_sylB_presptimes{:}));
                perstimulus(kk).yB = responseB.rawy(int64((responseB.zpt - dur)*responseB.fs): int64((responseB.zpt - dur + 0.5)*responseB.fs));
                perstimulus(kk).sptimesB = tocompare_sylB_sptimes;
                perstimulus(kk).eachtrail_RS_B = cellfun(@length,perstimulus(kk).sptimesB);

                perstimulus(kk).RS_A = length(vertcat(tocompare_sylA_sptimes{:})) - length(vertcat(tocompare_sylA_presptimes{:}));
                perstimulus(kk).RS_A_300 = length(vertcat(sylA300_sptimes{:})) - length(vertcat(sylA300_presptimes{:}));
                perstimulus(kk).response_A_300 = length(vertcat(sylA300_sptimes{:}))/(0.3*10);
                perstimulus(kk).response_A = length(vertcat(tocompare_sylA_sptimes{:}));
                perstimulus(kk).yA = responseA.rawy(int64(responseA.zpt*responseA.fs): int64((responseA.zpt + 0.5)*responseA.fs));
                perstimulus(kk).sptimesA = tocompare_sylA_sptimes;
                perstimulus(kk).eachtrail_RS_A = cellfun(@length,perstimulus(kk).sptimesA)

                perstimulus(kk).fid_A = responseA.Fid;
                perstimulus(kk).fid_B = responseB.Fid;
                perstimulus(kk).fid_A = responseAB.Fid;


                perstimulus(kk).linearsum_RS_A_RS_B = perstimulus(kk).eachtrail_RS_A  + perstimulus(kk).eachtrail_RS_B;


                perstimulus(kk).RS_AB = length(vertcat(sptimes_AB{:})) - length(vertcat(presptimes_AB{:}));
                perstimulus(kk).yAB = responseAB.rawy(int64(responseAB.zpt*responseAB.fs): int64((responseAB.zpt + 0.5)*responseAB.fs));
                perstimulus(kk).sptimesAB = sptimes_AB;
                perstimulus(kk).eachtrail_RS_AB = cellfun(@length,perstimulus(kk).sptimesAB);

                [h1,p1] = ttest(perstimulus(kk).linearsum_RS_A_RS_B,perstimulus(kk).eachtrail_RS_AB);
                perstimulus(kk).hvalue = nan; % 初始化
                if h1 ==1
                    [h2,p2] = ttest(perstimulus(kk).linearsum_RS_A_RS_B,perstimulus(kk).eachtrail_RS_AB,'Tail','right');
                    if h2 == 0
                        perstimulus(kk).hvalue = 1;
                    else
                        perstimulus(kk).hvalue = -1;
                    end
                else
                    perstimulus(kk).hvalue = 0;
                end

                perstimulus(kk).ratio = perstimulus(kk).RS_AB/ ( perstimulus(kk).RS_A + perstimulus(kk).RS_B  );
                perstimulus(kk).log_ratio = log(perstimulus(kk).RS_AB/ ( perstimulus(kk).RS_A + perstimulus(kk).RS_B  ));
                perstimulus(kk).diff = perstimulus(kk).RS_AB - ( perstimulus(kk).RS_A + perstimulus(kk).RS_B  );
                perstimulus(kk).stimuliname = newcategolist(kk).fullname;
                perstimulus(kk).catego = newcategolist(kk).catego;
                perstimulus(kk).dRS_RS = perstimulus(kk).diff/( perstimulus(kk).RS_A + perstimulus(kk).RS_B);
                perstimulus(kk).diff_sum_ratio = perstimulus(kk).diff/( abs(perstimulus(kk).RS_A + perstimulus(kk).RS_B) +abs(perstimulus(kk).RS_AB));

                if abs(perstimulus(kk).diff_sum_ratio) <=0.5
                    perstimulus(kk).bi_contexteffect = 0.0;
                elseif perstimulus(kk).diff_sum_ratio > 0.5
                    perstimulus(kk).bi_contexteffect = 1.0;
                elseif perstimulus(kk).diff_sum_ratio <-0.5
                    perstimulus(kk).bi_contexteffect = -1.0;
                else
                    perstimulus(kk).bi_contexteffect = nan;
                end
                
                perstimulus(kk).log_dRS_RS = sign(perstimulus(kk).dRS_RS)*log(1+abs(perstimulus(kk).dRS_RS)/10^1);
%                 try
                    perstimulus(kk).eachtrail_RS_diff = perstimulus(kk).eachtrail_RS_AB -(perstimulus(kk).eachtrail_RS_A + perstimulus(kk).eachtrail_RS_B);
%                 catch
% 
%                     deltainfo(kk).eachtrail_RS_diff = [];
% 
%                 end
            end


            % 计算那些需要zscored的数据
            normalized_diff = ([perstimulus.diff].'- mean([perstimulus.diff].'))/std([perstimulus.diff].');
            normalized_dRS_RS = ([perstimulus.dRS_RS].'- mean([perstimulus.dRS_RS].'))/std([perstimulus.dRS_RS].');
            normalized_RS_A = ([perstimulus.RS_A].' - mean([perstimulus.RS_A].'))/std([perstimulus.RS_A].');
            normalized_response_A = ([perstimulus.response_A].' - mean([perstimulus.response_A].'))/std([perstimulus.response_A].');

            normalized_diff_sum_ratio = ([perstimulus.diff_sum_ratio].' - mean([perstimulus.diff_sum_ratio].','omitnan'))/std([perstimulus.diff_sum_ratio].','omitnan');
            zscored_diff_sum_ratio = zscore([perstimulus.diff_sum_ratio].');
            zscored_RS_A = zscore([perstimulus.RS_A].');
            zscored_response_A = zscore([perstimulus.response_A].');

            for k = 1:length(perstimulus)
                perstimulus(k).normalized_diff = normalized_diff(k);
                perstimulus(k).normalized_dRS_RS = normalized_dRS_RS(k);
                perstimulus(k).normalized_RS_A = normalized_RS_A(k);
                perstimulus(k).zscored_RS_A = zscored_RS_A(k);
                perstimulus(k).normalized_response_A = normalized_response_A(k);
                perstimulus(k).zscored_response_A = zscored_response_A(k);
                perstimulus(k).normalized_diff_sum_ratio = normalized_diff_sum_ratio(k);
                perstimulus(k).zscored_diff_sum_ratio = zscored_diff_sum_ratio(k);
            end


            %计算Neurons关于每个类别的syllables的反应的d prime value，但是不一定对
            percatego = struct;

            categos = 1:8;

            for k = 1:length(categos)
                hitted = find(categos(k) == [perstimulus.catego].');
                percatego(k).response_A = [perstimulus(hitted).response_A].';
                percatego(k).RS_A_300 = [perstimulus(hitted).RS_A_300].';
                percatego(k).response_A_300 = [perstimulus(hitted).response_A_300].';
                percatego(k).mean_response_A = mean([perstimulus(hitted).response_A].');
                percatego(k).normalized_diff_sum_ratio = mean([perstimulus(k).normalized_diff_sum_ratio].');
                percatego(k).diff_sum_ratio = [perstimulus(hitted).diff_sum_ratio].';
                percatego(k).bi_contexteffect = [perstimulus(hitted).bi_contexteffect].';
                percatego(k).catego = categos(k);
                percatego(k).neuronname = neuronname;
            end

            for k = 1:length(percatego) % 计算d prime values

                count = 0;
                for kk = setdiff(1:length(percatego),k)
                    count = count + 1;
                    mean1 = mean([percatego(k).response_A].');
                    mean2 = mean([percatego(kk).response_A].');
                    variance1 = var([percatego(k).response_A].');
                    variance2 = var([percatego(kk).response_A].');
                    percatego(k).dprime(count) = 2*(mean1-mean2)/sqrt(variance1 + variance2);


                    %计算 diff_sum_ratio 的 dprime values
                    mean1 = mean([percatego(k).diff_sum_ratio].');
                    mean2 = mean([percatego(kk).diff_sum_ratio].');
                    variance1 = var([percatego(k).diff_sum_ratio].');
                    variance2 = var([percatego(kk).diff_sum_ratio].');
                    percatego(k).contextdprime(count) = 2*(mean1-mean2)/sqrt(variance1 + variance2);


                    % 计算 binary 的 d prime values
                    mean1 = mean([percatego(k).bi_contexteffect].');
                    mean2 = mean([percatego(kk).bi_contexteffect].');
                    variance1 = var([percatego(k).bi_contexteffect].');
                    variance2 = var([percatego(kk).bi_contexteffect].');
                    percatego(k).bidprime(count) = 2*(mean1-mean2)/sqrt(variance1 + variance2);

                end

                percatego(k).mindprime = min(percatego(k).dprime);
                percatego(k).mincontextdprime = min(percatego(k).contextdprime);
                percatego(k).maxcontextdprime = max(percatego(k).contextdprime);

                percatego(k).minbidprime = min(percatego(k).bidprime);
                percatego(k).maxbidprime = max(percatego(k).bidprime);

            end



            % 计算per-neuron 的 global info
            % 计算ANOVA
            perneuron.selectivecatego = [percatego(find([percatego.mindprime].' > 0.5)).catego].';
            perneuron.selective_positive_contextcatego = [percatego(find([percatego.mincontextdprime].' > 0.5)).catego].';
            perneuron.selective_negative_contextcatego = [percatego(find([percatego.maxcontextdprime].' < -0.5)).catego].';
            perneuron.selective_positive_bicatego = [percatego(find([percatego.minbidprime].' > 0.5)).catego].';
            perneuron.selective_negative_bicatego = [percatego(find([percatego.maxbidprime].' < -0.5)).catego].';

            if isempty(perneuron.selectivecatego )
                perneuron.selectivecatego = 0;
            end

            if isempty(perneuron.selective_positive_contextcatego )
                perneuron.selective_positive_contextcatego = 0;
            end

            if isempty(perneuron.selective_negative_contextcatego )
                perneuron.selective_negative_contextcatego = 0;
            end

             if isempty(perneuron.selective_positive_bicatego )
                perneuron.selective_positive_bicatego = 0;
             end

             if isempty(perneuron.selective_negative_bicatego )
                 perneuron.selective_negative_bicatego = 0;
             end


            perneuron.selectivemindprime =[percatego(find([percatego.mindprime].' > 0.5)).mindprime].';
            perneuron.selective_positive_mindprime = [percatego(find([percatego.mincontextdprime].' > 0.5)).catego].';
            perneuron.selective_negative_mindprime = [percatego(find([percatego.mincontextdprime].' < -0.5)).catego].';

          %  globalinfo.selectivepositivemindprime =[percategoinfo(find([percategoinfo.mindprime].' > 0.5)).mindprime].';

            for k = 1:length(perstimulus) % 方差齐性检测
                [perstimulus(k).bartlett_p,perstimulus(k).bartlett_stats] = vartestn(perstimulus(k).eachtrail_RS_diff.','Display','off');
            end
            

            [perneuron.anovap,tbl,stats] = anova1(vertcat(perstimulus.eachtrail_RS_diff).',[],'off');
            perneuron.anovaf = tbl{2,5};
            perneuron.anovatbl = tbl;

            post_hoc = multcompare(stats,"Display","off");

            perneuron.anova_stats = stats;
            perneuron.post_hoc = post_hoc;
            perneuron.post_hoc_report = array2table(post_hoc,"VariableNames", ...
                ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);

            unique_categos = unique([perstimulus.catego].');
            categoinfo = struct;
            for k = 1:length(unique_categos)
                categoinfo(k).catego = unique_categos(k);
                subinfo = perstimulus(find([perstimulus.catego].' == unique_categos(k)));
                categoinfo(k).mean_drs_rs = mean([subinfo.dRS_RS].');
                categoinfo(k).mean_A = mean([subinfo.RS_A].');
                categoinfo(k).mean_diff = mean([subinfo.diff].');
            end

            % global information
            [perneuron.bestresponseA, temp] = max([categoinfo.mean_A].');
            perneuron.bestcategoA = categoinfo(temp).catego;

            perneuron.neuronname = neuronname;
            [perneuron.bestresponse_drsrs, index] = max([categoinfo.mean_drs_rs].');
            [perneuron.worstresponse__drsrs, index2] = min([categoinfo.mean_drs_rs].');
            perneuron.bestcatego = categoinfo(index).catego;
            perneuron.worstcatego = categoinfo(index2).catego;


            perneuron.overalleffect = sum([perstimulus.bi_contexteffect].');




            parts = strsplit(allreplalist(1).stimuliname,'-before-');
            part2 = parts{2};
            perneuron.secondsyl_catego = str2num(regexp(convertCharsToStrings(part2),'(?<=Type)\d+','match'));
       

            [perneuron.bestresponse_diff, index3] = max([categoinfo.mean_diff].');
            perneuron.bestcatego_by_diff = categoinfo(index3).catego;

            dr.perneuron = perneuron;
            dr.neuronname = neuronname;

            dr.perstimulus = perstimulus;
            dr.percatego = percatego;

        

        
        end

        function  DrawSorted(dr, outdir)

            if ~exist('outdir','var')

                outdir = '.';

            end

            goodlist = dr.perstimulus;


            [~,index] = sortrows([goodlist.catego].'); goodlist = goodlist(index); clear index

            I = {};

            for k = 1: length(goodlist)

                figure('Color','w','Position',[2172 -413 696 1515]);
                Draw.six_forPoster(goodlist(k).yA,goodlist(k).sptimesA,...
                    goodlist(k).yB,goodlist(k).sptimesB,...
                    goodlist(k).yAB,goodlist(k).sptimesAB,...
                    32000);
                frame = getframe(gcf);
                temp = frame.cdata;
                saveas(gcf,sprintf('%s\\Delta_%s.fig',outdir,goodlist(k).stimuliname));
                saveas(gcf,sprintf('%s\\Delta_%s.eps',outdir,goodlist(k).stimuliname),'epsc');
                saveas(gcf,sprintf('%s\\Delta_%s.svg',outdir,goodlist(k).stimuliname),'svg');
                close(gcf);
                % imwrite(temp,'sss.tiff')

                I{k}= Convert.colorEdge(temp,'k',7);


            end


            % draw blank white

            lieshu = 10;

            hangshu = ceil(length(I)/lieshu);
            rest = lieshu*hangshu - length(I);
            white = uint8(255*ones(size(I{1})));

            if rest > 0
                for k = 1:rest
                    I = [I,white];
                    %                     ax = gcf;
                    %                     ax.Position(3) = 560;
                    %                     ax.Position(4) = 420;
                end
            end

            reshapedI = reshape(I, lieshu,[])';
            %clear I
            IMG = cell2mat(reshapedI);
            imwrite(IMG,sprintf('%s\\差值DeltaResponse_%s.png',outdir,dr.neuronname));

        end

    
    end
end