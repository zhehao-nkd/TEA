classdef Ephys < handle
    %EPHYS2 Summary of this class goes here
    %   Detailed explanation goes here
    properties(Constant)
        pltext = 0.15 % previosuly it was 0.1
        
    end
    properties
        spike
        trigger
        sound
        fs
        rawy  % rawy is the y of the original wav file
        rawfeatures
        y % y is rawy remove pre and post silence
        sptimes
        features
        meanfeatures
        presptimes
        plty  % plty is y plus pre / post pltext, respectivley
        pltsptimes
        pltfeatures
        zpt
        npt
        rawsptimes
        latency   
        frate % firing rate
        fnum % firing number
        label % responsive(1) or not(0)
        sigloc % location (time) of the significant neural response
        eliciting_element_counting_whole %  index of the response-elciting 
        %element by taking neural response to all stimuli into consideration
        % this Paramter are calculated outside the Ephys class
        trigger_onset % onset of trigger stimuli
    end
    
    methods
        function e = inilatency(e)
            e.latency = 0.1;
        end
        function e = Ephys(spike, trigger, sound)
            
            
            e.spike = spike;
            e.trigger = trigger;
            e.sound = sound;
            e.rawy = e.sound.rawy; % rawy
            e.y = e.sound.corey; % corey is y
            % disp('MAGA!');
            e.fs = sound.fs;
            e.zpt = e.sound.initial(1)/e.fs;
            e.npt = e.sound.terminal(end)/e.fs;
            e.allocate;
            e.updateplt;
            if ~isempty(regexp(e.sound.name,'norm'))  %%% This code has not been completed yet!!!
                [~,e.sigloc] = e.getSigRespnse;
                % calculate this value only when stimuli is not a syllable
                % THis is a very bad solution. As the real problem here is
                % that when I calculate the sdf for syllables, I onlly
                % include the time duration the same as the stimuli, didn;t
                % take the latency into consideration
            end
            
            
        end
        
        function e = updateplt(e)
            e.allocate;
            e.plty = [zeros(e.pltext*e.fs,1);e.y;zeros(e.pltext*e.fs,1)];
        end
        function e = allocate(e) % here the response is the collection of response
            
            n = e.sound.trigger;
            index = find([e.trigger.info.name].'== n);
            initials = e.trigger.info(index).time; % initials of all the repeats
            e.trigger_onset = initials;
            ylength = length(e.y)/e.sound.fs;
            
            for k = 1:length(initials)
                e.sptimes{k} = e.spike.time( initials(k)<e.spike.time & e.spike.time<initials(k) + ylength)...
                    - initials(k);
                e.presptimes{k} = e.spike.time( initials(k)- ylength<e.spike.time & e.spike.time<initials(k))- double((initials(k)-ylength));
                e.pltsptimes{k} = e.spike.time( initials(k)- e.pltext<e.spike.time& e.spike.time<initials(k) + ylength + e.pltext)...
                    - (initials(k)- e.pltext);
                e.rawsptimes{k} = e.spike.time( initials(k)- e.zpt<e.spike.time& e.spike.time<initials(k) + ylength + length(e.rawy)/e.fs- e.npt)...
                    - (initials(k)- e.zpt);  % 这里可能有一个巨大的bug
               
            end
           
            disp('Here might be a huge bug in Ephys.allocate!!!!!!');
        end
        
        function e = raster(e)% draw raster
            color = 'k';
            title = ' ';
            %figure;
            draw.raster(e.pltsptimes,e.plty,e.fs,color);
            %draw.raster(e.sound.corey, e.fs, e.sptimes);
        end
        
         function e = rasterpure(e)% draw raster pure means signal only
            color = 'k';
            title = ' ';
            %figure;
            draw.raster(e.sptimes,e.y,e.fs,color);
            %draw.raster(e.sound.corey, e.fs, e.sptimes);
        end
        
        
        function e = spec(e) % draw spectrogram
            
            figure;
            draw.spec(e.plty,e.fs);
            %draw.raster(e.sound.corey, e.fs, e.sptimes);
        end
        
        function e = sdf(e)    % draw sdf
          %  figure;
            draw.sdf(e.plty,e.fs,e.pltsptimes);
        end
        
        function e = threeWithFreqRelatedFeatures(e) % three plots for single syllable/element
            e.updateplt;
            
            figure('Visible','off','color','w');
            %figure;
            %subplot(3,1,1);
            t = tiledlayout(3,1);
            ax1 = axes(t);
            draw.spec(e.plty,e.fs);
            ylims = get(gca,'YLim');
            
            if ~isempty(e.pltfeatures)
                ax2 = axes(t);
                plot(e.pltfeatures.pitch,'-r')
                ax2.Color = 'none';
                ylim(ylims);
                
                
                ax3 = axes(t);
                plot(e.pltfeatures.mean_frequency,'-w')
                ax3.Color = 'none';
                ylim(ylims);
                
                ax4 = axes(t);
                plot(e.pltfeatures.peak_frequency,'-g')
                ax4.Color = 'none';
                ax1.Box = 'off';
                ax2.Box = 'off';
                ax3.Box = 'off';
                ax4.Box = 'off';
                ylim(ylims);
                
            end
            
            nexttile
            draw.raster(e.pltsptimes,e.plty,e.fs);
            
            nexttile
            draw.sdf(e.plty,e.fs,e.pltsptimes);
            xlabel(e.sound.name);
            
            
        end
        
        function e = threesingle(e) % three plots for single syllable/element
             e.updateplt;
            figure('Visible','off','color','w')
            draw.three(e.plty,e.fs,e.pltsptimes);
            
            xlabel(e.sound.name);
        end
        
        function e = slendertwo(e) % slender version of two plots
             e.updateplt;
            figure('Visible','off','color','w','Position',[2104 35 468 1013])
            draw.two(e.plty,e.fs,e.pltsptimes);
            xlabel(e.sound.name);
        end
        
        function e = three(e)% draw three plots
            e.updateplt;
            figure('color','w')
            draw.three(e.plty,e.fs,e.pltsptimes);
            xlabel(e.sound.name);
            subplot(3,1,1);
            locs = e.getIndexOfSigSyllable;
            for idx = 1: length(locs)
                initial = (e.sound.initial(locs(idx))-e.sound.initial(1))/e.fs + e.pltext;
                terminal = (e.sound.terminal(locs(idx))-e.sound.initial(1))/e.fs+ e.pltext;
                line([initial,terminal],[10,10],'color','r');
            end
        end
        
        function e = three_with_sig_resp(e)% draw three plots with marks on the significant response
            e.updateplt;
            figure('color','w')
            draw.three(e.plty,e.fs,e.pltsptimes);
            xlabel(e.sound.name);
            subplot(3,1,2);
            [pks,locs] = e.getSigRespnse;
            for idx = 1: length(locs) 
              xline(locs(idx) + e.pltext,'--hr');
            end
            subplot(3,1,3);
            hold on
            
            for idx = 1: length(locs) 
              scatter(locs(idx)+ e.pltext,pks(idx),[],'r','*')
            end
            hold off
            
        end
        
        function e = four(e)
            e.updateplt;
            figure('color','w')
            draw.four(e.plty,e.fs,e.pltsptimes);
            xlabel(e.sound.name);
            subplot(4,1,2);
            locs = e.sig;
            for idx = 1: length(locs)
                initial = (e.sound.initial(locs(idx))-e.sound.initial(1))/e.fs + e.pltext;
                terminal = (e.sound.terminal(locs(idx))-e.sound.initial(1))/e.fs+ e.pltext;
                line([initial,terminal],[10,10],'color','r');
            end
        end
    
        function e = rawthree(e)% draw three plots
            figure('color','w')
            draw.three(e.rawy,e.fs,e.rawsptimes);
            xlabel(e.sound.name);
            subplot(3,1,1);
            [~,locs] = e.getSigRespnse;
            locs = e.getIndexOfSigSyllable(0) % arbitritary set latency as 0 !!!!!!
            %locs = e.sig;
            for idx = 1: length(locs)
                initial = e.sound.initial(locs(idx))/e.fs;
                terminal = e.sound.terminal(locs(idx))/e.fs;
                line([initial,terminal],[10,10],'color','r');
            end
        end
        
        function [sdf_values,times_from_onset] = getSigRespnse(e) 
            % This function find out the response peak events above threshold,
            % regarding them as significant response event of the signal
            resolution = .001;
            gausswidth = .02;
            FORCEDHEIGHT= 15;
            minpeakprominence = 5;
            minpeakdistance = 30;
            range = 0.05;
            presdf = cal.sdf(e.presptimes,e.y,e.fs,resolution,gausswidth);
            thres = cal.thres(presdf,range);
            thres = max(thres, FORCEDHEIGHT);
            sdf = cal.sdf(e.sptimes,e.y,e.fs,resolution,gausswidth);
            % for debugging
            %             figure;
            %             draw.raster(e.sptimes,e.y,e.fs);
            %             figure
            %             findpeaks(sdf,"MinPeakHeight",thres,...
            %                 'MinPeakProminence',minpeakprominence,'MinPeakDistance',minpeakdistance);
            [sdf_values,times_from_onset,~,~]  = findpeaks(sdf,"MinPeakHeight",thres,... %Error using findpeaks Expected MinPeakDistance to be a scalar with value < 20.??????
                'MinPeakProminence',minpeakprominence,'MinPeakDistance',min(minpeakdistance,length(sdf) -2));  % very bad code

            
            % This section is used for screening out those response events
            % for which spikes do not happen for each trail
            width_thres = 50;
            for j = 1: length(times_from_onset)
                low_t = times_from_onset(j)-width_thres;
                high_t = times_from_onset(j)+ width_thres; 
                num_sig_trials = length(find(~cellfun(@isempty, cellfun(@(x)find(low_t<x<high_t),e.sptimes,'UniformOutput',0)))); % in how many trails there are spikes in this duration
                if num_sig_trials < length(e.sptimes)/2 % if num_sig_trails do not contain half of the trials
                    times_from_onset(j) = nan;
                end
            end
            ids = find(~isnan(times_from_onset)); % find ids which is not nan
            times_from_onset = times_from_onset(ids);
            sdf_values = sdf_values(ids);
%             w = w(ids);
%             p = p(ids);
            
            times_from_onset = times_from_onset*0.001;% convert miliseconds to seconds 
            
        end
       
        
        function sylidx = getIndexOfSigSyllable(e,latency) % index of significant syllable
            % This is based on neural respons to a simgle stimuli, which is
            % usually not accurate, so this function is not very useful
            % --10.13,2021,Zhehao
        
            if ~exist('latency','var')
                latency = 0; % default latency is 0
            end
            
            dbstop if error
            initial = ( e.sound.initial -e.sound.initial(1) )/e.fs;  % this is for converting initial to seconds relative to the first syllable
            [~,locs] = e.getSigRespnse;
            if ~isempty(locs)
                for idx = 1: length(locs)
                    sylidx(idx) = length(find(locs(idx)-latency>initial)); % return the idx of the significant syllables
                end
            else
                sylidx = [];
                return
            end
            sylidx = horzcat(sylidx(:));
            
            sylidx = unique(sylidx,'stable'); % remove repeated element
            
        end
        
        
        function respinfo = JudgeRespToWIthinSongFrags_fromStimuliToResponse(e,latency) % index of significant syllable
            % based on latency information to set hysteretic window for
            % delineate corresponding sptimes
        
            if ~exist('latency','var')
                latency = 0; % default latency is 0
            end
            
            respinfo = e.sound.fragment;
            
            dbstop if error
            
            for k = 1: length(e.sound.initial)
                this_initial = e.sound.initial(k) /e.fs;
                this_terminal = e.sound.terminal(k) /e.fs;
                len_rawy = length(e.sound.rawy)/e.fs;
                per_onset = this_initial/len_rawy; % percentage of onset
                per_offset = this_terminal/len_rawy; % percentage of offset
                
                featurenames = fieldnames(e.rawfeatures);
                featurenames = setdiff(featurenames,{'file_index','file_name'});
                for gg = 1: length(featurenames)
                    temp = e.rawfeatures.(featurenames{gg});
                    
                    if round(length(temp)*per_onset) == 0
                        respinfo(k).(featurenames{gg}) = temp(1:round(length(temp)*per_offset)); % if it is zero, set the initial to 1
                    else
                        respinfo(k).(featurenames{gg}) = temp(round(length(temp)*per_onset):round(length(temp)*per_offset));
                    end
                    respinfo(k).(sprintf('mean_%s',featurenames{gg})) = mean(respinfo(k).(featurenames{gg}));
                end

                respinfo(k).songname = e.sound.name;
                respinfo(k).fragnum = k;
                respinfo(k).fs = e.fs;
                respinfo(k).lagging_sptimes = convert.sptimesOnset2Zero(extract.sptimes(e.rawsptimes,this_initial + latency,this_terminal + latency), this_initial + latency);
                disp(respinfo(k).lagging_sptimes);
                respinfo(k).ypluslatency = [respinfo(k).y;zeros(fix(latency*e.fs),1)];
                
                
                temp_of_sum = cal.psth_frag(respinfo(k).y,respinfo(k).fs,respinfo(k).lagging_sptimes);
                halfsum = sum(temp_of_sum(end/2:end));
                fullsum = sum(temp_of_sum);
                maxvalue = max(cal.psth_frag(respinfo(k).y,respinfo(k).fs,respinfo(k).lagging_sptimes));
                respinfo(k).maxvalue = maxvalue;
                respinfo(k).halfsum = halfsum;
                respinfo(k).fullsum = fullsum;
                if maxvalue > 6 % here the threshold is very important % originally set as 8
                    respinfo(k).label = 1;
                else
                    respinfo(k).label = 0;
                end  
            end
          
        end
        
        
        
        function sylidx = sapsig(e) % index of significant syllable
            dbstop if error
            initial = ( e.sound.sapinitial -e.sound.sapinitial(1) )/e.fs;  % this is for cvonverting initial to seconds relative to the first syllable
            [~,locs] = e.getSigRespnse;
            if ~isempty(locs)
                for idx = 1: length(locs)
                    sylidx(idx) = length(find(~(locs(idx)<initial))); % return the idx of the significant syllables
                end
            else
                sylidx = [];
            end
            sylidx = horzcat(sylidx(:));
            
            sylidx = unique(sylidx,'stable'); % remove repeated element
            
        end
        
        function siginf = siginf(e) % full info of sig syllable
            dbstop if error
            if isempty(e.sig)
                siginf = [];
                return
            end
            idx = find([e.sylinf.label] == 1);
            temp = e.sylinf;
            siginf = temp(idx); % derive from sylinf
            
        end
        
        function num = sylnum(e)
            num = length(e.sound.initial);
        end
        
        function syl = sylinf(e) % full info of all syllables
            dbstop if error
            number = e.sig;
            syl = struct;
            for n = 1:length(e.sound.initial)
                
                
                syl(n).sound = convertStringsToChars(e.sound.name);
                syl(n).number = n;
                syl(n).plx = convertStringsToChars(e.trigger.plxname);
                syl(n).channel = e.spike.channel;
                syl(n).unit = e.spike.unit;
                syl(n).y = e.sound.fragment(n).y;
                %syl(n).image = cal.img(syl(n).y,e.fs); % store the image matrix
                temp = extract.feature(syl(n).y,e.fs);
                syl(n).goodness = temp.goodness;
                syl(n).meanf = temp.mean_frequency;
                syl(n).fm = temp.fm;
                syl(n).amplitude = temp.amplitude;
                syl(n).entropy = temp.entropy;
                syl(n).pitch = temp.pitch;
                syl(n).rawpitch = temp.rawpitch;
                syl(n).am = temp.AM;
                syl(n).initial = e.sound.fragment(n).initial;
                syl(n).terminal = e.sound.fragment(n).terminal;
                if n ~= 1
                    syl(n).pregap = (syl(n).initial - syl(n-1).terminal)/e.fs;
                else
                    syl(n).pregap = inf; % pre-gap duration
                end
                syl(n).dur = length(syl(n).y)/e.fs;
                
                if ismember(n,number)
                    syl(n).label = 1; % significant
                else
                    syl(n).label = 0; % not significant
                end
                
                
            end
            
        end
        
        function syl = sapsylinf(e)
            dbstop if error
            number = e.sapsig;
            syl = struct;
            thissap = e.sound.getsap;
            for n = 1:length(e.sound.sapinitial)
                
                
                syl(n).sound = convertStringsToChars(e.sound.name);
                syl(n).number = n;
                syl(n).plx = convertStringsToChars(e.trigger.plxname);
                syl(n).channel = e.spike.channel;
                syl(n).unit = e.spike.unit;
                syl(n).y = e.sound.sapfragment(n).y;
               
                %syl(n).image = cal.img(syl(n).y,e.fs); % store the image matrix
                %temp = extract.feature(syl(n).y,e.fs);
                syl(n).goodness = thissap(n).mean_goodness;
                syl(n).meanf = thissap(n).mean_freq;
                syl(n).fm = thissap(n).mean_FM;
                syl(n).amplitude = thissap(n).mean_amplitude;
                syl(n).entropy = thissap(n).mean_entropy;
                syl(n).pitch = thissap(n).mean_pitch;
               % syl(n).rawpitch = temp.rawpitch;
                syl(n).am = thissap(n).mean_AM;
                syl(n).initial = e.sound.sapfragment(n).initial;
                syl(n).terminal = e.sound.sapfragment(n).terminal;
                syl(n).dur = (syl(n).terminal - syl(n).initial)/e.fs;
                syl(n).yplt = [e.sound.sapfragment(n).y;zeros(e.latency*e.fs,1)];
                e.inilatency;
                syl(n).sptimesplt = extract.cutspt(e.rawsptimes,syl(n).initial/e.fs, syl(n).dur + e.latency);
                
                %syl(n).sptimes = e.sptimes(syl(n).initial<e.sptimes< syl(n).terminal); % !~!!!!!!!!!!!
                if n ~= 1
                    syl(n).pregap = (syl(n).initial - syl(n-1).terminal)/e.fs;
                else
                    syl(n).pregap = inf; % pre-gap duration
                end
                syl(n).dur = length(syl(n).y)/e.fs;
                
                if ismember(n,number)
                    syl(n).label = 1; % significant
                else
                    syl(n).label = 0; % not significant
                end
                
                
            end
        end
        
        function syl = sylinfq(e) % quick but less
            dbstop if error
            number = e.sig;
            syl = struct;
            for n = 1:length(e.sound.initial)
                
                
                syl(n).sound = convertStringsToChars(e.sound.name);
                syl(n).number = n;
                syl(n).plx = convertStringsToChars(e.trigger.plxname);
                syl(n).channel = e.spike.channel;
                syl(n).unit = e.spike.unit;
                syl(n).y = e.sound.fragment(n).y;
                syl(n).fs = e.fs;
                syl(n).initial = e.sound.fragment(n).initial;
                syl(n).terminal = e.sound.fragment(n).terminal;
                if n ~= 1
                     syl(n).pregap = (syl(n).initial - syl(n-1).terminal)/e.fs;
                else
                    syl(n).pregap = inf; % pre-gap duration
                end
                syl(n).dur = length(syl(n).y)/e.fs;
                
                if ismember(n,number)
                    syl(n).label = 1; % significant
                else
                    syl(n).label = 0; % not significant
                end
                
                
            end
            
        end    
        
        function syl = avgn(e) % full info of all syllables
            dbstop if error
            number = e.sig;
            syl = struct;
            for n = 1:length(e.sound.initial)
                syl(n).birdid = convert.bid(convertStringsToChars(e.sound.name));
                syl(n).filename = convertStringsToChars(e.sound.name);
                syl(n).number = n;
                syl(n).plx = convertStringsToChars(e.trigger.plxname);
                syl(n).channel = e.spike.channel;
                syl(n).unit = e.spike.unit;
                syl(n).y = e.sound.fragment(n).y;
                syl(n).hpy = highpass(syl(n).y,450,e.fs); % high passed y, threshold is 400
                syl(n).image = cal.img(syl(n).hpy,e.fs); % store the image matrix

                if ismember(n,number)
                    syl(n).label = 1; % significant
                else
                    syl(n).label = 0; % not significant
                end
            end
           
            
        end
        
        function response = resp(e) % mimic the old response function
            
            response.plxname = e.trigger.plxname;
            response.channel_name = e.spike.channel;
            response.unit_name = e.spike.unit;
            response.stimuli_name = convertStringsToChars(e.sound.name);
            response.wav_Y = e.sound.y;
            response.initial = e.sound.initial;
            response.terminal = [e.sound.fragment.terminal].';
            response.sigsylidx = e.sig;
            
        end
        
        function pre = preinf(e, duration)% info of pre-peak duration
            
            if ~exist('duration','var')
                dur = 0.200 ;% second , this is default
            else
                dur = duration;
            end
            
            locs =  e.peak;
            if ~isempty(locs)
                
                for n = 1: length(locs)
                    pre(n).y = Sound.intercept(e.y,e.fs,locs(n)-dur,locs(n));
                    pre(n).fs = e.fs;
                    temp = extract.feature(pre(n).y,e.fs);
                    pre(n).goodness = temp.goodness;
                    pre(n).meanf = temp.mean_frequency;
                    pre(n).fm = temp.fm;
                    pre(n).amplitude = temp.amplitude;
                    pre(n).entropy = temp.entropy;
                    pre(n).pitch = temp.pitch;
                    pre(n).am = temp.AM;
                    pre(n).plx = e.trigger.plxname;
                    pre(n).channel = e.spike.channel{1};
                    pre(n).unit = e.spike.unit;
                    pre(n).sound = e.sound.name;
                end
            else
                pre = [];
            end
            
        end
        
        function fea = featurecurve(e)
            temp = SAT_hijack(e.y,e.fs);
            fea.rawpitch = temp.features.pitch;
            newpitch = temp.features.pitch;
            newpitch(newpitch > 2000) = nan;
            newpitch(newpitch < 400) = nan;
            fea.pitch = newpitch;
            fea.amplitude = temp.features.amplitude;
            fea.FM = temp.features.FM;
            fea.AM = temp.features.AM;
            fea.goodness = temp.features.goodness;
            fea.entropy = temp.features.entropy;
            fea.mfreq = temp.features.mean_frequency;
            fea.y = e.y;
        end
        
        function drawfeaturecurve(e)
            fea = featurecurve(e);
            figure('color','w');
            name ={'pitch', 'amplitude', 'FM', 'AM', 'goodness', 'entropy', 'mfreq'};
            subplot(length(name)+2,1,1);
            draw.spec(e.y,e.fs);
            ylabel('Hz');
            subplot(length(name)+2,1,2);
            draw.raster(e.sptimes,e.y,e.fs,'k');
            ylabel('trails');
            for idx = 1: length(name)
                subplot(length(name)+2,1,idx+ 2);
                eval(['plot(fea.',name{idx},');']);
                xlim([0,length(fea.mfreq)]);
                ylabel(name{idx});
                hold on
            end
        end
         
    end
end








