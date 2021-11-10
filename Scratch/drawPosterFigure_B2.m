% drawPosterFigure_B2


info = f6.list;


wlen = 300;
hop = 75;
nfft = 300;

fs = 32000;


for k = 1: length(info)
    [info(k).stft,~,~] = DQ.stft(info(k).y,wlen, hop,nfft,fs);
    [S,F,T] = spectrogram(y,hamming(512),512-round(fs/1e3),512,fs);
    info(k).specdata = log(1+abs(S));
end

sima = []

for a = 1:length(info)
    for b = 1: length(info)
        alpha = info(a).specdata;
        beta = info(b).specdata;
        sima(a,b) = norm(alpha(:)-beta(:));
        
        %sima(a,b) = ssim(info(a).stft,info(b).stft);
    end
end





