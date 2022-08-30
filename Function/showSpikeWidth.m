
function spike_info = showSpikeWidth(path_plx)


spike_info = struct;
count = 0;


plxfile = readPLXFileC(convertStringsToChars(path_plx),'all')

fs_plx = 40000;

for k = 1:length(plxfile.SpikeChannels)
    
    units = setdiff(unique(plxfile.SpikeChannels(k).Units),[0]);
    
     for u = 1: length(units)
         
         ids = find(plxfile.SpikeChannels(k).Units==units(u));
         
         spikewfs = plxfile.SpikeChannels(k).Waves(:,ids);
         
         
        [~,peak] =  max(spikewfs);
         [~,valley] =  min(spikewfs);
         
         widths_in_datapoints = peak - valley;
         
         
         
         width_in_time =  widths_in_datapoints/fs_plx*1000; % ms
         
         mean_width = mean(width_in_time);
         
         disp(sprintf('Channel: %s, Unit: %u, Spike width: %f ms, NumOfSpikes: %u',...
             plxfile.SpikeChannels(k).Name, units(u),...
             mean_width,length(ids)));
         
         count = count + 1;
         spike_info(count).channel = plxfile.SpikeChannels(k).Name;
         spike_info(count).unit = units(u);
         spike_info(count).meanspikewidth = mean_width;
         spike_info(count).numofspikes = length(ids);
        
     end
  
end

[~,plxname,~] = fileparts(path_plx);
save(sprintf('Spike_Info_%s.mat',strrep(plxname,'-','_')),'spike_info');

end