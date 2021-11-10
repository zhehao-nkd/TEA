%  This is a script for processing normalizede stimuli series.
% It works by adding silent paddings together with trigger pulses to
% another channel

% By Zhehao Cheng, 0820,2020


%　This is compitible for generating binary-code triggers


addpath(genpath("C:\Users\Zhehao\Dropbox (OIST)\My_Matlab\TEA"));
dbstop if error
 
pre_length = 1; post_length = 3; % 224
GAP_DURATION = 1*8e-3; PULSE_DURATION = 5*8e-3; AMPLITUDE = 1;
from_which = 1;
%initial = 148;

    

dirpath = uigetdir();
files = extract.filename(dirpath, '*.wav');
files = flip(files,1);
outdir = 'AAAAAA'
mkdir (outdir);

for n = 1:length(files)
       
   [y,fs] = audioread(files{n}); % y means original y
    
    binary_code = de2bi(n);
   
   zero_pulse = [zeros(GAP_DURATION*fs,1);AMPLITUDE*ones(GAP_DURATION*fs,1)];
   one_pulse = [zeros(GAP_DURATION*fs,1);zeros(GAP_DURATION*fs,1);AMPLITUDE*ones(GAP_DURATION*fs,1)];
   
   % write data channel

   
   yD = [zeros(pre_length*fs,1);y;zeros(post_length*fs,1)];
   
   
   % write trigger channel
    pulses = AMPLITUDE*ones(GAP_DURATION*fs,1); % Here the initial state of variable-pulses is a pulse
    coder.varsize(pulses); % yT means y Trigger channel
    
    for m = 1:length(binary_code)
        if binary_code(m) == 0
            pulses = [pulses;zero_pulse];
        elseif binary_code(m) == 1
            pulses = [pulses;one_pulse];
        end
        
    end
    
    yT = [zeros(pre_length*fs,1);pulses;zeros(length(yD) - length(pulses) - pre_length*fs,1)]; % padding zeros
    
    yOut = [yT,yD];
    
    [~,name,~] = fileparts(files{n});
    
    audiowrite(sprintf('%s/%s-%uPulses.wav',outdir,name,from_which),yOut,fs);
    
    from_which = from_which + 1;
end

