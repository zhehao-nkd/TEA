% a class

classdef Trigger
    
    properties
        equipment
        info
        plxname
        raw
    end
    
    methods
        
        function t = Trigger(path_plx)
            
            [~,t.plxname,~] = fileparts(path_plx);
            path_plx = convertStringsToChars(path_plx);
            plexonfile = readPLXFileC(path_plx,'events'); % read .plx files as a structure
            
            PULSE_LOCATION = 0;
            %
            % Find the trigger
            
            % Find the trigger plexonfile.EventChannels(1).Timestamps
            for n = 1:length(plexonfile.EventChannels)
                
                if  strcmp(plexonfile.EventChannels(n).Name, 'DIG01')
                    PULSE_LOCATION = n;
                    t.equipment = 'ZEUS';
                    t.raw = plexonfile.EventChannels(PULSE_LOCATION);
                    t.info = t.digextract; % digital
                    clear plexonfile
                    return
                end
                
            end
            
            plexonfile = readPLXFileC(path_plx,'continuous');
            
            for m = 1:length(plexonfile.ContinuousChannels)
                
                if  strcmp(plexonfile.ContinuousChannels(m).Name, 'AI01')
                    PULSE_LOCATION = m;
                    t.equipment = 'PLEXON';
                    raw = plexonfile.ContinuousChannels(PULSE_LOCATION);
                    t.info = t.extract; % not digital
                    clear plexonfile
                    return
                end
            end
            
            
            
        end
        
        
        function time = getTime(t, pulse_number)
            if exist('pulse_number','var')
                time = t.info(pulse_number).time;
            else
                time = {t.info.time}.';
            end
        end
     
        
        function properties = getsingle(t, which)
            properties.equipment = t.equipment;
            properties.info = t.info(which);
            properties.plxname = t.plxname;
        end
        
        function info = digextract(t)
            
            
            fs = 30000;
            pulse_initials = t.raw.Timestamps;
            
            
            
            % line([pulse_initials,pulse_initials],[0,5]); %for debug
            
            m = 1;
            groups ={};
            coder.varsize(groups);
            
            % As different stimuli have different number of square pulses, this works for grouping pulse clusters
            
            WINDOW = fs; % The length of searching window, 30000samples means 1 second for Zeus
            
            
            for n = 1: length(pulse_initials)
                
                if n == length(pulse_initials)
                    groups{m}= [groups{m},pulse_initials(n)];
                    
                    
                    
                else
                    
                    if pulse_initials(n+1)-pulse_initials(n)<WINDOW
                        if isempty(groups)
                            groups{m} = pulse_initials(n);
                        else
                            groups{m}=[ groups{m},pulse_initials(n)];
                        end
                        
                    else
                        if isempty(groups)
                            groups{m}= pulse_initials(n);
                            m = m+1;
                            groups{m} =[];
                        else
                            groups{m}= [groups{m},pulse_initials(n)];
                            m = m+1;
                            groups{m} =[];
                        end
                    end
                    
                    
                    
                end
                
            end
            
            
            
            % For collecting .time and .number
            detection.samplepoint = [];
            detection.name= [];
            coder.varsize(detection.name,detection.samplepoint);
            
            for o = 1:length(groups)
                
                detection.samplepoint(o) = groups{o}(1);
                detection.name(o) = length(groups{o});
            end
            
            detection.time = detection.samplepoint/fs;
            
            
            name = unique(detection.name);
            
            %sort triggers with same number of pulses
            
            for n = 1:length(name) % for write name into the struct 'info'
                info(n).name = name(n);
            end
            
            for n = 1:length(name)
                info(n).repeats = length(find(detection.name==name(n)));
                info(n).time = detection.time(find(detection.name==name(n)));
            end
            
            
        end
        
        function info = extract(t)
            
            fs = 20000;
            %dbstop if error
            pulse_channel = t.raw.Values;
            %detecting triggers and recording their time
            
            % debug pulse_channel = trigger.Values;
            
            MIN_THRESHOLD = 5000; % May vary in different situations
            MAX_THRESHOLD = 8500; % May vary in different situations
            
            
            pulse_channel(pulse_channel < MIN_THRESHOLD | pulse_channel > MAX_THRESHOLD) = 0; %Convert singals of pulse channel to 0 and 1
            %figure; plot(pulse_channel);
            
            pulse_channel(pulse_channel ~= 0) = 1; % based on the threshold
            %figure; plot(pulse_channel);
            
            
            
            k = find(pulse_channel);   % Find all the 1
            
            if k(1) == 1 % incase that k(1) = 1, when there are noises interrrput the pulse channel
                tempinitial = find(pulse_channel == 0);
                newinitial = tempinitial(1);
                pulse_channel = pulse_channel(newinitial:end);
                k = find(pulse_channel);
            end
            
            %    difference =  k(2:end)- k(1:end-1) ;
            %    idxDif = find(difference ~= 1);
            
            
            pulse_initials = k(pulse_channel(k - 1)==0);      % Find the initial of all the pulses
            pulse_ends = k(pulse_channel(k + 1)==0);          % Find the ends of all the pulses
            
            
            
            % Remove noises
            for n = 1:length(pulse_initials)
                
                if pulse_ends(n) - pulse_initials(n) < 10e-3 / 2 * fs % 20000 is the fs of AI channel
                    pulse_initials(n) = NaN;
                    pulse_ends(n) = NaN;
                end
                
            end
            
            
            pulse_initials (isnan(pulse_initials)) = [];
            pulse_ends (isnan(pulse_ends)) = [];
            
            
            % line([pulse_initials,pulse_initials],[0,5]); %for debug
            
            m = 1;
            groups ={};
            coder.varsize(groups);
            
            % As different stimuli have different number of square pulses, this works for grouping pulse clusters
            
            WINDOW = fs; % The length of searching window, 20000samples means 1 second
            
            
            for n = 1: length(pulse_initials)
                
                if n == length(pulse_initials)
                    groups{m}= [groups{m},pulse_initials(n)];
                    
                    
                    
                else
                    
                    if pulse_initials(n+1)-pulse_initials(n)<WINDOW
                        if isempty(groups)
                            groups{m} = pulse_initials(n);
                        else
                            groups{m}=[ groups{m},pulse_initials(n)];
                        end
                        
                    else
                        if isempty(groups)
                            groups{m}= pulse_initials(n);
                            m = m+1;
                            groups{m} =[];
                        else
                            groups{m}= [groups{m},pulse_initials(n)];
                            m = m+1;
                            groups{m} =[];
                        end
                    end
                    
                    
                    
                end
                
            end
            
            
            
            % For collecting .time and .number
            detection.samplepoint = [];
            detection.name= [];
            coder.varsize(detection.name,detection.samplepoint);
            
            for o = 1:length(groups)
                
                detection.samplepoint(o) = groups{o}(1);
                detection.name(o) = length(groups{o});
            end
            
            detection.time = detection.samplepoint/fs;
            
            
            name = unique(detection.name);
            
            %sort triggers with same number of pulses
            
            for n = 1:length(name) % for write name into the struct 'info'
                info(n).name = name(n);
            end
            
            for n = 1:length(name)
                info(n).repeats = length(find(detection.name==name(n)));
                info(n).time = detection.time(find(detection.name==name(n)));
            end
            
            
        end
        
    end
    
    
    
end




