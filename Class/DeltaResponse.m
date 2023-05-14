classdef DeltaResponse
    %R(A+B) - RA - RB

    properties
        deltainfo
        neuronname
        globalinfo
        percategoinfo
    end

    methods
        function dr = DeltaResponse(responseB,allfraglist,allreplalist,newcategolist, neuronname)


            deltainfo = struct;

            for kk = 1:length(newcategolist)
                hid_frag = intersect( find(~cellfun(@isempty, regexp(cellstr({allfraglist.stimuliname}.'),newcategolist(kk).fullname))),...
                    find(~cellfun(@isempty, regexp(cellstr({allfraglist.stimuliname}.'),'Type'))) );
                hid_repla = find(~cellfun(@isempty, regexp(cellstr({allreplalist.stimuliname}.'),newcategolist(kk).fullname)));

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


                repla_zpt = responseAB.zpt;
                sptimes_AB = Extract.sptimes_resetSP(responseAB.rawsptimes, repla_zpt , repla_zpt + 0.5); % 取500ms的间距
                presptimes_AB = Extract.sptimes_resetSP(responseAB.rawsptimes, repla_zpt -0.5, repla_zpt ); % 取500ms的间距

                deltainfo(kk).neuronname = neuronname;
                deltainfo(kk).fragname = newcategolist(kk).fullname;
                deltainfo(kk).RS_B =  length(vertcat(tocompare_sylB_sptimes{:})) - length(vertcat(tocompare_sylB_presptimes{:}));
                deltainfo(kk).yB = responseB.rawy(int64((responseB.zpt - dur)*responseB.fs): int64((responseB.zpt - dur + 0.5)*responseB.fs));
                deltainfo(kk).sptimesB = tocompare_sylB_sptimes;
                deltainfo(kk).eachtrail_RS_B = cellfun(@length,deltainfo(kk).sptimesB);

                deltainfo(kk).RS_A = length(vertcat(tocompare_sylA_sptimes{:})) - length(vertcat(tocompare_sylA_presptimes{:}));
                deltainfo(kk).response_A = length(vertcat(tocompare_sylA_sptimes{:}));
                deltainfo(kk).yA = responseA.rawy(int64(responseA.zpt*responseA.fs): int64((responseA.zpt + 0.5)*responseA.fs));
                deltainfo(kk).sptimesA = tocompare_sylA_sptimes;
                deltainfo(kk).eachtrail_RS_A = cellfun(@length,deltainfo(kk).sptimesA);


                deltainfo(kk).linearsum_RS_A_RS_B = deltainfo(kk).eachtrail_RS_A  + deltainfo(kk).eachtrail_RS_B;


                deltainfo(kk).RS_AB = length(vertcat(sptimes_AB{:})) - length(vertcat(presptimes_AB{:}));
                deltainfo(kk).yAB = responseAB.rawy(int64(responseAB.zpt*responseAB.fs): int64((responseAB.zpt + 0.5)*responseAB.fs));
                deltainfo(kk).sptimesAB = sptimes_AB;
                deltainfo(kk).eachtrail_RS_AB = cellfun(@length,deltainfo(kk).sptimesAB);

                [h1,p1] = ttest(deltainfo(kk).linearsum_RS_A_RS_B,deltainfo(kk).eachtrail_RS_AB);
                if h1 ==1
                    [h2,p2] = ttest(deltainfo(kk).linearsum_RS_A_RS_B,deltainfo(kk).eachtrail_RS_AB,'Tail','right');
                    if h2 == 0
                        deltainfo(kk).hvalue = 1;
                    else
                        deltainfo(kk).havlue = -1;
                    end
                else
                    deltainfo(kk).hvalue = 0;
                end



                deltainfo(kk).ratio = deltainfo(kk).RS_AB/ ( deltainfo(kk).RS_A + deltainfo(kk).RS_B  );
                deltainfo(kk).log_ratio = log(deltainfo(kk).RS_AB/ ( deltainfo(kk).RS_A + deltainfo(kk).RS_B  ));
                deltainfo(kk).diff = deltainfo(kk).RS_AB - ( deltainfo(kk).RS_A + deltainfo(kk).RS_B  );
                deltainfo(kk).stimuliname = newcategolist(kk).fullname;
                deltainfo(kk).catego = newcategolist(kk).catego;
                deltainfo(kk).dRS_RS = deltainfo(kk).diff/( deltainfo(kk).RS_A + deltainfo(kk).RS_B);
                deltainfo(kk).diff_sum_ratio = deltainfo(kk).diff/( abs(deltainfo(kk).RS_A + deltainfo(kk).RS_B) +abs(deltainfo(kk).RS_AB));
                deltainfo(kk).log_dRS_RS = sign(deltainfo(kk).dRS_RS)*log(1+abs(deltainfo(kk).dRS_RS)/10^1);
%                 try
                    deltainfo(kk).eachtrail_RS_diff = deltainfo(kk).eachtrail_RS_AB -(deltainfo(kk).eachtrail_RS_A + deltainfo(kk).eachtrail_RS_B);
%                 catch
% 
%                     deltainfo(kk).eachtrail_RS_diff = [];
% 
%                 end
            end


            % 计算那些需要zscored的数据
            normalized_diff = ([deltainfo.diff].'- mean([deltainfo.diff].'))/std([deltainfo.diff].');
            normalized_dRS_RS = ([deltainfo.dRS_RS].'- mean([deltainfo.dRS_RS].'))/std([deltainfo.dRS_RS].');
            normalized_RS_A = ([deltainfo.RS_A].' - mean([deltainfo.RS_A].'))/std([deltainfo.RS_A].');
            normalized_response_A = ([deltainfo.response_A].' - mean([deltainfo.response_A].'))/std([deltainfo.response_A].');
            zscored_RS_A = zscore([deltainfo.RS_A].');
            zscored_response_A = zscore([deltainfo.response_A].');

            for k = 1:length(deltainfo)
                deltainfo(k).normalized_diff = normalized_diff(k);
                deltainfo(k).normalized_dRS_RS = normalized_dRS_RS(k);
                deltainfo(k).normalized_RS_A = normalized_RS_A(k);
                deltainfo(k).zscored_RS_A = zscored_RS_A(k);
                deltainfo(k).normalized_response_A = normalized_response_A(k);
                deltainfo(k).zscored_response_A = zscored_response_A(k);
            end


            %计算Neurons关于每个类别的syllables的反应的d prime value，但是不一定对
            percategoinfo = struct;

            categos = 1:8;

            for k = 1:length(categos)
                hitted = find(categos(k) == [deltainfo.catego].');
                percategoinfo(k).response_A = [deltainfo(hitted).response_A].';
                percategoinfo(k).mean_response_A = mean([deltainfo(hitted).response_A].');
                percategoinfo(k).catego = categos(k);
                percategoinfo(k).neuronname = neuronname;
            end

            for k = 1:length(percategoinfo) % 计算d prime values

                count = 0;
                for kk = setdiff(1:length(percategoinfo),k)
                    count = count + 1;
                    mean1 = mean([percategoinfo(k).response_A].');
                    mean2 = mean([percategoinfo(kk).response_A].');
                    variance1 = var([percategoinfo(k).response_A].');
                    variance2 = var([percategoinfo(kk).response_A].');
                    percategoinfo(k).dprime(count) = 2*(mean1-mean2)/sqrt(variance1 + variance2);
                end

                percategoinfo(k).mindprime = min(percategoinfo(k).dprime);

            end



            % 计算per-neuron 的 global info
            % 计算ANOVA
            globalinfo.selectivecatego = [percategoinfo(find([percategoinfo.mindprime].' > 0.5)).catego].';

            if isempty(globalinfo.selectivecatego )
                globalinfo.selectivecatego = 0;
            end
            globalinfo.selectivemindprime =[percategoinfo(find([percategoinfo.mindprime].' > 0.5)).mindprime].';

            for k = 1:length(deltainfo) % 方差齐性检测
                [deltainfo(k).bartlett_p,deltainfo(k).bartlett_stats] = vartestn(deltainfo(k).eachtrail_RS_diff.','Display','off');
            end
            

            [globalinfo.anovap,tbl,stats] = anova1(vertcat(deltainfo.eachtrail_RS_diff).',[],'off');
            globalinfo.anovaf = tbl{2,5};
            globalinfo.anovatbl = tbl;

            post_hoc = multcompare(stats,"Display","off");

            globalinfo.anova_stats = stats;
            globalinfo.post_hoc = post_hoc;
            globalinfo.post_hoc_report = array2table(post_hoc,"VariableNames", ...
                ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);

            unique_categos = unique([deltainfo.catego].');
            categoinfo = struct;
            for k = 1:length(unique_categos)
                categoinfo(k).catego = unique_categos(k);
                subinfo = deltainfo(find([deltainfo.catego].' == unique_categos(k)));
                categoinfo(k).mean_drs_rs = mean([subinfo.dRS_RS].');
                categoinfo(k).mean_A = mean([subinfo.RS_A].');
                categoinfo(k).mean_diff = mean([subinfo.diff].');
            end

            % global information
            [globalinfo.bestresponseA, temp] = max([categoinfo.mean_A].');
            globalinfo.bestcategoA = categoinfo(temp).catego;

            globalinfo.neuronname = neuronname;
            [globalinfo.bestresponse_drsrs, index] = max([categoinfo.mean_drs_rs].');
            [globalinfo.worstresponse__drsrs, index2] = min([categoinfo.mean_drs_rs].');
            globalinfo.bestcatego = categoinfo(index).catego;
            globalinfo.worstcatego = categoinfo(index2).catego;


            parts = strsplit(allreplalist(1).stimuliname,'-before-');
            part2 = parts{2};
            globalinfo.secondsyl_catego = str2num(regexp(convertCharsToStrings(part2),'(?<=Type)\d+','match'));
       

            [globalinfo.bestresponse_diff, index3] = max([categoinfo.mean_diff].');
            globalinfo.bestcatego_by_diff = categoinfo(index3).catego;

            dr.globalinfo = globalinfo;
            dr.neuronname = neuronname;

            dr.deltainfo = deltainfo;
            dr.percategoinfo = percategoinfo;

        

        
        end

        function  DrawSorted(dr)


            goodlist = dr.deltainfo;


            [~,index] = sortrows([goodlist.catego].'); goodlist = goodlist(index); clear index

            I = {};

            for k = 1: length(goodlist)

                subI = {};
                figure('Color','w','Position',PM.size1);
                Draw.two(goodlist(k).yA,32000,goodlist(k).sptimesA);
                frame = getframe(gcf);
                subI{1} = frame.cdata;
                close(gcf);

                figure('Color','w','Position',PM.size1);
                Draw.two(goodlist(k).yB,32000,goodlist(k).sptimesB);
               % xlabel(goodlist(k).fragname);
                frame = getframe(gcf);
                subI{2} = frame.cdata;
                close(gcf);

                figure('Color','w','Position',PM.size1);
                Draw.two(goodlist(k).yAB,32000,goodlist(k).sptimesAB);
                xlabel(goodlist(k).fragname);
                frame = getframe(gcf);

                subI{3} = frame.cdata;
                close(gcf);

                temp = vertcat(subI{:});

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
            imwrite(IMG,sprintf('肃州DeltaResponse_%s.png',dr.neuronname));

        end

    
    end
end