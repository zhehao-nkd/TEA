%addpath(genpath("C:\Users\Zhehao\Dropbox (OIST)\My_Matlab\TEA"));

function ConDirs = fromWavFindParentFolder(dir_of_wavs)
wavs = extract.filename(dir_of_wavs,'*.wav');

wavfile = {};

for n = 1: length(wavs)
    
    [~,wavfile{n},~] = fileparts(wavs{n});
    
    
end

subdir = extract.folder("Z:\Yazaki-SugiyamaU\Bird-song");

dirend = {};
for k = 1: length(subdir)
    
    temp = split(subdir{k},'\');
    dirend{k} = convert.bid(temp{end});
    
end

ids = {};
for m = 1: length(wavfile)
    ids{m} = find(~cellfun(@isempty,regexp( wavfile{m},dirend)));
    
end

ids = vertcat(ids{:});

ConDirs = subdir(ids);
end