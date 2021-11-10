
% -------------- FT options --------------
% labeled_data = load('C:\Code\birdsong_new\list.mat');
% dataDir = "C:\code\birdsong_new\";
function preprocessVrae(eleinf)

datafilename = "eleConOnly.mat";
%labeled_data = load(".\" + datafilename);
dataSaveDir = ".\";

labeled_data.eleinf = eleinf;

wlen = 300;
hop = 75;
nfft = 300;

fs = 32000;

anal_win = blackmanharris(wlen, 'periodic');
synth_win = hamming(wlen, 'periodic');


try
    list = labeled_data.list;
catch
    list = labeled_data.eleinf;

end

data = zeros(length(list), 1, 151);

% figure
songnames = [];
%figure
for n_s = 1 :  length(list)
    
    syllable = list(n_s);
    y = syllable.y / max(abs(syllable.y));
    
    [STFT, F, T] = DQ.stft(y, anal_win, hop, nfft, fs);
    
    S_softplused = bdsoftplus(abs(STFT));
    
    labeled_data.list(n_s).S_softplused = S_softplused';
    labeled_data.list(n_s).STFT = STFT;
    labeled_data.list(n_s).T = T;
    labeled_data.list(n_s).F = F;
    
    data(n_s, 1:size(S_softplused, 2), :) = S_softplused';
    
    %     contourf(T, F, S_softplused, 30, 'LineStyle', 'None')
    %     title(list(n_s).songname)
    %     drawnow
    %     pause(2)
    
    try
        labels(n_s) = list(n_s).label;
    catch
        try
            if isempty(songnames) ||  sum(ismember(songnames, list(n_s).songname)) < 1
                songnames = [songnames, list(n_s).songname];
                labels(n_s) =find(ismember(songnames, list(n_s).songname)) ;
            else
                labels(n_s) =find(ismember(songnames, list(n_s).songname)) ;
            end
        catch
            labels(n_s) = 0;
        end
    end
    
    try
        fr = 0;
        spiketimes = [];
        for k = 1:10
            fr = fr + size(list(n_s).sptimes{k}, 1);
            spiketimes = [spiketimes; list(n_s).sptimes{k}];
        end
        if isempty(spiketimes)
            spiketimes_std(n_s) = -1;
        else
            spiketimes_std(n_s) = std(spiketimes);
        end
        frs(n_s) = fr;
    catch
        frs = 0;
        spiketimes_std = 0;
    end
    
    labeled_data.list(n_s).label = labels(n_s);
    
end


save(dataSaveDir + "labeled_" + datafilename, "labeled_data", "data", "labels", "frs", "spiketimes_std")

end

function y = bdsoftplus(x)
    
   y = (1 - 2*(x < 0)).* log(1 + (1 - 2*(x < 0)) .* x);

end


function x = inversebdsoftplus(y)
    
   x = (1 - 2*(y < 0)).* exp(y) - (1 - 2*(y < 0));

end