% To test spectrogram parameters


[y,fs] = audioread('example 1.wav');



figure;
draw.spec(y,fs);
entropy = spectralEntropy(y,fs );

figure;
plot(entropy)


fname = "E:\FeatureOfStimuli\B659_2.json"'; 
fid = fopen(fname); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
val = jsondecode(str);

figure;
plot([val.amplitude].')