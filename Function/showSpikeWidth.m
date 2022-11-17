
function spike_info = showSpikeWidth(path_pl2)


spike_info = struct;
count = 0;


pl2file = readpl2FileC(convertStringsToChars(path_pl2),'all')

fs_pl2 = 40000;

for k = 1:length(pl2file.SpikeChannels)
    
    units = setdiff(unique(pl2file.SpikeChannels(k).Units),[0]);
    
     for u = 1: length(units)
         
         ids = find(pl2file.SpikeChannels(k).Units==units(u));
         
         spikewfs = pl2file.SpikeChannels(k).Waves(:,ids);
         
         
        [~,peak] =  max(spikewfs);
         [~,valley] =  min(spikewfs);
         
         widths_in_datapoints = peak - valley;
         
         
         
         width_in_time =  widths_in_datapoints/fs_pl2*1000; % ms
         
         mean_width = mean(width_in_time);
         
         disp(sprintf('Channel: %s, Unit: %u, Spike width: %f ms, NumOfSpikes: %u',...
             pl2file.SpikeChannels(k).Name, units(u),...
             mean_width,length(ids)));
         
         count = count + 1;
         spike_info(count).channel = pl2file.SpikeChannels(k).Name;
         spike_info(count).unit = units(u);
         spike_info(count).meanspikewidth = mean_width;
         spike_info(count).numofspikes = length(ids);
        
     end
  
end

[~,pl2name,~] = fileparts(path_pl2);
save(sprintf('Spike_Info_%s.mat',strrep(pl2name,'-','_')),'spike_info');

end