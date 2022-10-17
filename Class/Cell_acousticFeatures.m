classdef Cell_acousticFeatures < handle
    %Another axillary class of class Neuron
    
    methods
        function [dists1,dists0,featurename] = calCumulativeFeatureDiff(a,featurename)
            %fraglist = a.judgeFragResponse;
            a.judgeFragResp_FR;
            fragids = find(~cellfun(@isempty, regexp(cellstr({a.list.stimuliname}.'),'Frag|frag|syl|ele')));
            fraglist = a.list(fragids);
            labeled_frags = fraglist([fraglist.label].' == 1);
            fe1 = []; for k = 1:length(labeled_frags); eval(sprintf('fe1(k) = labeled_frags(k).meanfeatures.%s',featurename)); end
            if isempty(fe1)
                dists1 = [];
                dists0 = [];
                featurename = featurename;
                return
            end
            pairs1 = (nchoosek(fe1,2));
            dists1 = abs(pairs1(:,1)-pairs1(:,2));
            
            
            unlabeled_frags = fraglist([fraglist.label].' == 0);
            fe0 = [];  error = struct;
            for k = 1:length(unlabeled_frags)
                try
                    eval(sprintf('fe0(k) = unlabeled_frags(k).meanfeatures.%s',featurename));
                catch ME
                    error(k).ME = ME;
                    fe0(k) = nan;
                    disp('Feature Info Missing!!')
                end
            end
            fe0 = rmmissing(fe0);
            pairs0 = (nchoosek(fe0,2));
            dists0 = abs(pairs0(:,1)-pairs0(:,2));
            % fig = figure('Position',a.fig_size1,'Color','w');  hold on;
            
            %             if isempty(dists1)
            %                 frame = getframe(gcf);
            %                 img = frame.cdata;
            %                 close(gcf);
            %                 return
            %             end
            %             cdfplot(dists1);
            %             cdfplot(dists0);
            %             legend('Response-eliciting elements','Not eliciting elements','Location','best')
            %             xlabel(sprintf('Difference between Mean %s (Hz)',featurename),'interpreter', 'none');
            %             ylabel('Cumulative %');
            %             hold off
            %             %             [f0,x0]= ecdf(dists0)
            %             %             [f1,x1]= ecdf(dists1)
            %             [h,p] = kstest2(dists1,dists0);
            %             %set(fig,'defaultTextInterpreter','none')
            %             title(sprintf('%s P-value : %.8f',a.formated_name,p),'interpreter', 'none');
            %             %saveas(gcf,sprintf('CDF-%s.png',a.formated_name));
            %             frame = getframe(gcf);
            %             img = frame.cdata;
            %             close(gcf);
            
        end
        function img = exportDrawPitchHarmoLineChart(a)
            % draw the distribution of mean features
            a.calHarmRatio;
            fraglist = a.judgeFragResponse;
            if isempty(fraglist)
                return;
            end
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).maxvalue;
            end
            
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).pitch = fraglist(u).meanfeatures.pitch;
                toshow(u).harmratio = fraglist(u).meanfeatures.harmratio;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            
            toshow  = table2struct(sortrows(struct2table(toshow), 'resp','ascend'));
            
            
            num_toshow = [ toshow.pitch; toshow.harmratio].';
            
            znum_toshow = zscore(num_toshow,0,1);
            %
            %            figure
            %            plot(znum_toshow.');
            
            figure('Position',a.fig_size1,'Color','w');
            hold on
            c = 1- rescale([toshow.resp].',0.1,1);
            for r = 1: size(znum_toshow,1)
                plot(znum_toshow(r,:),'Color',repmat(c(r),3,1));
                %drawnow
                % pause(0.5)
            end
            colormap(flip(repmat(unique(c),1,3)))
            colorbar
            xlim([0,3]);
            xticks([0 1 2 3])
            xticklabels({'','Pitch','Harmonic ratio',''});
            title(sprintf('%s---%u song elements',a.formated_name,length(fraglist)),'interpreter','none');
            ylabel('Zscored Feature(averaged)');
            
            saveas(gcf,sprintf('PitchHarmoLineChart-%s.png',a.formated_name));
            
            frame = getframe(gcf);
            img = frame.cdata;
            close(gcf);
        end
        function img = exportDrawInsong_drawPitchHarmoLineChart(a)
            % draw the distribution of mean features
            fraglist = a.getInsongFragRespList;
            if isempty(fraglist)
                return;
            end
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).label;
            end
            
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).pitch = fraglist(u).meanfeatures.pitch;
                toshow(u).harmratio = fraglist(u).meanfeatures.harmratio;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            
            toshow  = table2struct(sortrows(struct2table(toshow), 'resp','ascend'));
            
            [~,index] = sortrows([toshow.resp].'); toshow = toshow(index); clear index
            num_toshow = [ toshow.pitch; toshow.harmratio].';
            
            znum_toshow = zscore(num_toshow,0,1);
            %
            %            figure
            %            plot(znum_toshow.');
            
            figure('Position',a.fig_size1,'Color','w');
            hold on
            c = 1- rescale([toshow.resp].',0.1,1);
            for r = 1: size(znum_toshow,1)
                plot(znum_toshow(r,:),'Color',repmat(c(r),3,1));
                %drawnow
                % pause(0.5)
            end
            colormap(flip(repmat(unique(c),1,3)))
            colorbar
            xlim([0,3]);
            xticks([0 1 2 3])
            xticklabels({'','Pitch','Harmonic ratio',''});
            title(sprintf('%s---%u song elements',a.formated_name,length(fraglist)),'interpreter','none');
            ylabel('Zscored Feature(averaged)');
            
            saveas(gcf,sprintf('Insong-PitchHarmoLineChart-%s.png',a.formated_name));
            
            frame = getframe(gcf);
            img = frame.cdata;
            close(gcf);
        end
        function img = exportDrawdraw2DPitchVsHarmo(a)
            
            fraglist = a.judgeFragResponse;
            if isempty(fraglist)||isempty(fraglist(1).meanfeatures)
                return
            end
            
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).rs;
            end
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).pitch = fraglist(u).meanfeatures.pitch;
                toshow(u).harmratio = fraglist(u).meanfeatures.harmratio;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            [~,index] = sortrows([toshow.resp].'); toshow = toshow(index); clear index
            
            name = {'pitch','harmratio'};
            ncomb = nchoosek(name,2);
            
            I = {};
            for idx = 1:size(ncomb,1) % 此处可用parfor
                
                %subplot(hangshu,lieshu,idx);
                figure('Position',a.fig_size1,'Color','w');
                hold on
                cmap = 1 - rescale([toshow.resp].');
                for  md = 1: length(toshow)
                    scatter(toshow(md).(ncomb{idx,1})', toshow(md).(ncomb{idx,2}),[],repmat(cmap(md,:),1,3),'filled');
                end
                xlabel(replace(ncomb{idx,1},'_','-'));
                ylabel(replace(ncomb{idx,2},'_','-'));
                title(a.formated_name,'interpreter', 'none')
                frame = getframe(gcf);
                I{idx} = frame.cdata;
                close(gcf);
                
            end
            
            reshapedI = reshape(I, 1,[])';
            clear I
            img = cell2mat(reshapedI);
            imwrite(img,sprintf('2DPitchHarmo_%s.png',a.formated_name));
            
        end
        
        function img = exportDrawInsong_draw2DPitchVsHarmo(a)
            
            fraglist = a.judgeFragResponse;
            if isempty(fraglist)||isempty(fraglist(1).meanfeatures)
                return
            end
            
            for k = 1: length(fraglist)
                fraglist(k).responseMeasure = fraglist(k).label; % 二分法
            end
            toshow = struct;
            for u = 1: length(fraglist)
                toshow(u).pitch = fraglist(u).meanfeatures.pitch;
                toshow(u).harmratio = fraglist(u).meanfeatures.harmratio;
                toshow(u).resp = fraglist(u).responseMeasure;
            end
            
            [~,index] = sortrows([toshow.resp].'); toshow = toshow(index); clear index
            name = {'pitch','harmratio'};
            ncomb = nchoosek(name,2);
            
            I = {};
            for idx = 1:size(ncomb,1) % 此处可用parfor
                
                %subplot(hangshu,lieshu,idx);
                figure('Position',a.fig_size1,'Color','w');
                hold on
                cmap = 0.9 - rescale([toshow.resp].')*0.9;
                for  md = 1: length(toshow)
                    scatter(toshow(md).(ncomb{idx,1})', toshow(md).(ncomb{idx,2}),[],repmat(cmap(md,:),1,3),'filled');
                end
                xlabel(replace(ncomb{idx,1},'_','-'));
                ylabel(replace(ncomb{idx,2},'_','-'));
                title(a.formated_name,'interpreter', 'none')
                frame = getframe(gcf);
                I{idx} = frame.cdata;
                close(gcf);
                
            end
            
            reshapedI = reshape(I, 1,[])';
            clear I
            img = cell2mat(reshapedI);
            imwrite(img,sprintf('Insong_2DPitchHarmo_%s.png',a.formated_name));
            
        end
        
        function img = exportDrawCumulativePitchDiff(a)
            fraglist = a.judgeFragResponse;
            labeled_frags = fraglist([fraglist.label].' == 1);
            pitch1 = []; for k = 1:length(labeled_frags); pitch1(k) = labeled_frags(k).meanfeatures.pitch; end
            pairs1 = (nchoosek(pitch1,2));
            dists1 = abs(pairs1(:,1)-pairs1(:,2));
            
            
            unlabeled_frags = fraglist([fraglist.label].' == 0);
            pitch0 = []; for k = 1:length(unlabeled_frags); pitch0(k) = unlabeled_frags(k).meanfeatures.pitch; end
            pairs0 = (nchoosek(pitch0,2));
            dists0 = abs(pairs0(:,1)-pairs0(:,2));
            fig = figure('Position',a.fig_size1,'Color','w');  hold on;
            
            if isempty(dists1)
                frame = getframe(gcf);
                img = frame.cdata;
                close(gcf);
                return
            end
            cdfplot(dists1);
            cdfplot(dists0);
            legend('Response-eliciting elements','Not eliciting elements','Location','best')
            xlabel('Difference between Mean Pitch (Hz)');
            ylabel('Cumulative %');
            hold off
            %             [f0,x0]= ecdf(dists0)
            %             [f1,x1]= ecdf(dists1)
            [h,p] = kstest2(dists1,dists0);
            %set(fig,'defaultTextInterpreter','none')
            title(sprintf('%s P-value : %.8f',a.formated_name,p),'interpreter', 'none');
            saveas(gcf,sprintf('CDF-%s.png',a.formated_name));
            frame = getframe(gcf);
            img = frame.cdata;
            close(gcf);
            
        end
        
        function img = exportdrawCumulativeFeatureDiff(a,featurename)
            fraglist = a.judgeFragResponse;
            labeled_frags = fraglist([fraglist.label].' == 1);
            fe1 = []; for k = 1:length(labeled_frags); eval(sprintf('fe1(k) = labeled_frags(k).meanfeatures.%s',featurename)); end
            pairs1 = (nchoosek(fe1,2));
            dists1 = abs(pairs1(:,1)-pairs1(:,2));
            
            
            unlabeled_frags = fraglist([fraglist.label].' == 0);
            fe0 = []; for k = 1:length(unlabeled_frags); eval(sprintf('fe0(k) = unlabeled_frags(k).meanfeatures.%s',featurename)); end
            pairs0 = (nchoosek(fe0,2));
            dists0 = abs(pairs0(:,1)-pairs0(:,2));
            fig = figure('Position',a.fig_size1,'Color','w');  hold on;
            
            if isempty(dists1)
                frame = getframe(gcf);
                img = frame.cdata;
                close(gcf);
                return
            end
            cdfplot(dists1);
            cdfplot(dists0);
            legend('Response-eliciting elements','Not eliciting elements','Location','best')
            xlabel(sprintf('Difference between Mean %s (Hz)',featurename),'interpreter', 'none');
            ylabel('Cumulative %');
            hold off
            %             [f0,x0]= ecdf(dists0)
            %             [f1,x1]= ecdf(dists1)
            [h,p] = kstest2(dists1,dists0);
            %set(fig,'defaultTextInterpreter','none')
            title(sprintf('%s P-value : %.8f',a.formated_name,p),'interpreter', 'none');
            %saveas(gcf,sprintf('CDF-%s.png',a.formated_name));
            frame = getframe(gcf);
            img = frame.cdata;
            close(gcf);
            
        end
        
        function img = exportdrawInsongCumulativePitchDiff(a)
            fraglist = a.getInsongFragRespList;
            labeled_frags = fraglist([fraglist.label].' == 1);
            pitch1 = []; for k = 1:length(labeled_frags); pitch1(k) = labeled_frags(k).meanfeatures.pitch; end
            if isempty(pitch1)||length(pitch1) == 1
                figure('Position',[1999 84 872 891],'Color','w');  hold on;
                title('malfunction')
                frame = getframe(gcf);
                img = frame.cdata;
                close(gcf);
                return
            end
            pairs1 = (nchoosek(pitch1,2));
            dists1 = abs(pairs1(:,1)-pairs1(:,2));
            
            
            unlabeled_frags = fraglist([fraglist.label].' == 0);
            pitch0 = []; for k = 1:length(unlabeled_frags); pitch0(k) = unlabeled_frags(k).meanfeatures.pitch; end
            pairs0 = (nchoosek(pitch0,2));
            dists0 = abs(pairs0(:,1)-pairs0(:,2));
            fig = figure('Position',a.fig_size1,'Color','w');  hold on;
            
            if isempty(dists1)
                frame = getframe(gcf);
                img = frame.cdata;
                close(gcf);
                return
            end
            cdfplot(dists1);
            cdfplot(dists0);
            legend('Response-eliciting elements','Not eliciting elements','Location','best')
            xlabel('Pairwise Mean Pitch Difference (Hz)');
            ylabel('Cumulative %');
            hold off
            %             [f0,x0]= ecdf(dists0)
            %             [f1,x1]= ecdf(dists1)
            [h,p] = kstest2(dists1,dists0);
            %set(fig,'defaultTextInterpreter','none')
            title(sprintf('%s P-value : %.8f',a.formated_name,p),'interpreter', 'none');
            saveas(gcf,sprintf('Insong-CDF-%s.png',a.formated_name));
            frame = getframe(gcf);
            img = frame.cdata;
            close(gcf);
            
        end
           
    end
end

