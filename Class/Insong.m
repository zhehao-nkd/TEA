classdef Insong
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        insonglist
    end

    methods
        function obj = Insong(normlist)
            % 根据normlist构造出neuron对insong syllables的反应
            
        end


        function insonglist = getInsongFragRespList(neu)
            % get In-songFrag Response List
            cons_neuron = neu.experiments{neu.song_id};
            insonglist = cons_neuron.toList_Insong;

        end

        
        function saveDrawInsong_lineChart_2D_Cumulative(neu)
            img1 = neu.Insong_drawPitchHarmoLineChart;
            img2 = neu.Insong_draw2DPitchVsHarmo;
            img3 = neu.Insong_drawCumulativePitchDiff;
            img = horzcat(img1,img2,img3);
            imwrite(img,sprintf('汉Insong_lineChart_2D_Cumulative_%s.png',neu.formated_name));
        end

        function saveDrawMeanFeaturesInSongVsRespAsLineChart(neu) % draw the distribution of mean features

            fraglist = neu.calResponseToWithinSongFragsFromEleinf;
            if isempty(fraglist)
                return;
            end
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).maxvalue;
            end

            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).amplitude = fraglist(u).mean_amplitude;
                toshow(u).pitch = fraglist(u).mean_pitch;
                toshow(u).AM = fraglist(u).mean_AM;
                toshow(u).mean_frequency = fraglist(u).mean_mean_frequency;
                toshow(u).FM = fraglist(u).mean_FM;
                toshow(u).peak_frequency = fraglist(u).mean_peak_frequency;
                toshow(u).entropy = fraglist(u).mean_entropy;
                toshow(u).resp = fraglist(u).responseMeasure;
            end

            toshow  = table2struct(sortrows(struct2table(toshow), 'resp','ascend'));


            num_toshow = [toshow.amplitude; toshow.pitch; toshow.AM; toshow.mean_frequency; toshow.FM; ...
                toshow.peak_frequency; toshow.entropy].';

            znum_toshow = zscore(num_toshow,0,1);
            %            figure
            %            plot(znum_toshow.');

            figure('Position',[1997 233 1388 658],'Color','w');
            hold on
            c = 1- rescale([toshow.resp].',0.1,1);
            for r = 1: length(znum_toshow)
                plot(znum_toshow(r,:),'Color',repmat(c(r),3,1));
            end

            colormap(flip(repmat(unique(c),1,3)))
            colorbar
            xlim([0,8]);
            xticks([0 1 2 3 4 5 6 7 8])
            xticklabels({'','amplitude','pitch','AM','mean_frequency','FM','peak-frequency','entropy',''});

            ylabel('Zscored Feature(averaged)');
            title(sprintf('Totally %u song elements',length(fraglist)));
            saveas(gcf,sprintf('V1-WithinSongsLineChartMeanFeaturesVsResp-%s.png',neu.formated_name));
            close(gcf);
        end

        

    end
end