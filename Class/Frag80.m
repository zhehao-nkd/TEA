classdef Frag80
    %R(A+B) - RA - RB

    properties
        perstimulus
        neuronname
        perneuron
        percatego
    end

    methods
        function dr = Frag80(frag,consistency,newcategolist,neuronname)
            dbstop if error

            allfraglist = frag.allfraglist;
            hitted_consistency = consistency;
            %deltaresponse = DeltaResponse(responseB,allfraglist,allreplalist, newcategolist,neuronname);
            perstimulus = struct;
            for kk = 1:length(newcategolist)
                hid_frag = intersect( find(~cellfun(@isempty, regexp(cellstr({allfraglist.stimuliname}.'),newcategolist(kk).fullname))),...
                    find(~cellfun(@isempty, regexp(cellstr({allfraglist.stimuliname}.'),'Type'))) );
              
                responseA = allfraglist(hid_frag);



                frag_zpt = responseA.zpt;
                sylA500_sptimes = Extract.sptimes_resetSP(responseA.rawsptimes, frag_zpt, frag_zpt + 0.5); % 取500ms的间距
                sylA500_presptimes = Extract.sptimes_resetSP(responseA.rawsptimes, frag_zpt -0.5, frag_zpt ); % 取500ms的间距

                sylA200_sptimes = Extract.sptimes_resetSP(responseA.rawsptimes, frag_zpt, frag_zpt + 0.2); % 取500ms的间距
                sylA200_presptimes = Extract.sptimes_resetSP(responseA.rawsptimes, frag_zpt -0.2, frag_zpt ); % 取500ms的间距

                sylA300_sptimes = Extract.sptimes_resetSP(responseA.rawsptimes, frag_zpt, frag_zpt + 0.3); % 取500ms的间距
                sylA300_presptimes = Extract.sptimes_resetSP(responseA.rawsptimes, frag_zpt -0.3, frag_zpt ); % 取500ms的间距

                sylA250_sptimes = Extract.sptimes_resetSP(responseA.rawsptimes, frag_zpt, frag_zpt + 0.25); % 取500ms的间距
                sylA250_presptimes = Extract.sptimes_resetSP(responseA.rawsptimes, frag_zpt -0.25, frag_zpt ); % 取500ms的间距


               
                perstimulus(kk).neuronname = neuronname;
                perstimulus(kk).fragname = newcategolist(kk).fullname;
               

               
                perstimulus(kk).RS_A_200 = length(vertcat(sylA200_sptimes{:})) - length(vertcat(sylA200_presptimes{:}));
                perstimulus(kk).RS_A_250 = length(vertcat(sylA250_sptimes{:})) - length(vertcat(sylA250_presptimes{:}));
                perstimulus(kk).RS_A_300 = length(vertcat(sylA300_sptimes{:})) - length(vertcat(sylA300_presptimes{:}));
                perstimulus(kk).RS_A_500 = length(vertcat(sylA500_sptimes{:})) - length(vertcat(sylA500_presptimes{:}));
                perstimulus(kk).response_A_300 = length(vertcat(sylA300_sptimes{:}))/(0.2*10);
                perstimulus(kk).response_A = length(vertcat(sylA500_sptimes{:}));
                perstimulus(kk).yA = responseA.rawy(int64(responseA.zpt*responseA.fs): int64((responseA.zpt + 0.5)*responseA.fs));
                perstimulus(kk).sptimesA = sylA500_sptimes;
                perstimulus(kk).eachtrail_RS_A = cellfun(@length,perstimulus(kk).sptimesA)

                perstimulus(kk).fid_A = responseA.Fid;

                perstimulus(kk).stimuliname = newcategolist(kk).fullname;
                perstimulus(kk).catego = newcategolist(kk).catego;

     
%                 end
            end


            % 计算那些需要zscored的数据
            normalized_RS_A_200 = ([perstimulus.RS_A_200].' - mean([perstimulus.RS_A_200].'))/std([perstimulus.RS_A_200].');
            normalized_RS_A_250 = ([perstimulus.RS_A_250].' - mean([perstimulus.RS_A_250].'))/std([perstimulus.RS_A_250].');
            normalized_RS_A_300 = ([perstimulus.RS_A_300].' - mean([perstimulus.RS_A_300].'))/std([perstimulus.RS_A_300].');
            normalized_RS_A_500 = ([perstimulus.RS_A_500].' - mean([perstimulus.RS_A_500].'))/std([perstimulus.RS_A_500].');

            normalized_response_A = ([perstimulus.response_A].' - mean([perstimulus.response_A].'))/std([perstimulus.response_A].');

         
            zscored_RS_A_200 = zscore([perstimulus.RS_A_200].');
            zscored_RS_A_250 = zscore([perstimulus.RS_A_250].');
            zscored_RS_A_300 = zscore([perstimulus.RS_A_300].');
            zscored_RS_A_500 = zscore([perstimulus.RS_A_500].');
            zscored_response_A = zscore([perstimulus.response_A].');

            for k = 1:length(perstimulus)
                perstimulus(k).normalized_RS_A_200 = normalized_RS_A_200(k);
                perstimulus(k).normalized_RS_A_250 = normalized_RS_A_250(k);
                perstimulus(k).normalized_RS_A_300 = normalized_RS_A_300(k);
                perstimulus(k).normalized_RS_A_500 = normalized_RS_A_500(k);
                perstimulus(k).zscored_RS_A_200 = zscored_RS_A_200(k);
                perstimulus(k).zscored_RS_A_250 = zscored_RS_A_250(k);
                perstimulus(k).zscored_RS_A_300 = zscored_RS_A_300(k);
                perstimulus(k).zscored_RS_A_500 = zscored_RS_A_500(k);
                perstimulus(k).normalized_response_A = normalized_response_A(k);
                perstimulus(k).zscored_response_A = zscored_response_A(k);

            end


            %计算Neurons关于每个类别的syllables的反应的d prime value，但是不一定对
            percatego = struct;

            categos = 1:8;

            for k = 1:length(categos)
                hitted = find(categos(k) == [perstimulus.catego].');
                percatego(k).response_A = [perstimulus(hitted).response_A].';
                percatego(k).RS_A_200 = [perstimulus(hitted).RS_A_200].';
                percatego(k).response_A_300 = [perstimulus(hitted).response_A_300].';
                percatego(k).mean_response_A = mean([perstimulus(hitted).response_A].');
                percatego(k).catego = categos(k);
                percatego(k).neuronname = neuronname;
            end

            for k = 1:length(percatego) % 计算d prime values

                count = 0;
                for kk = setdiff(1:length(percatego),k)

                    %计算 responseA 的 dprime values
                    count = count + 1;
                    mean1 = mean([percatego(k).RS_A_200].');
                    mean2 = mean([percatego(kk).RS_A_200].');
                    variance1 = var([percatego(k).RS_A_200].');
                    variance2 = var([percatego(kk).RS_A_200].');
                    percatego(k).dprime(count) = 2*(mean1-mean2)/sqrt(variance1 + variance2);


   

                end

                percatego(k).mindprime = min(percatego(k).dprime);
 

            end



            % 计算per-neuron 的 global info
            % 计算ANOVA
            perneuron.selectivecatego = [percatego(find([percatego.mindprime].' > 0.5)).catego].';

            if isempty(perneuron.selectivecatego )
                perneuron.selectivecatego = 0;
            end

    


            perneuron.selectivemindprime =[percatego(find([percatego.mindprime].' > 0.5)).mindprime].';
    
          %  globalinfo.selectivepositivemindprime =[percategoinfo(find([percategoinfo.mindprime].' > 0.5)).mindprime].';

         

            unique_categos = unique([perstimulus.catego].');
            categoinfo = struct;
            for k = 1:length(unique_categos)
                categoinfo(k).catego = unique_categos(k);
                subinfo = perstimulus(find([perstimulus.catego].' == unique_categos(k)));
                categoinfo(k).mean_A = mean([subinfo.RS_A_200].');
            end

            % global information
            [perneuron.bestresponseA, temp] = max([categoinfo.mean_A].');
            perneuron.bestcategoA = categoinfo(temp).catego;

            perneuron.neuronname = neuronname;
     

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