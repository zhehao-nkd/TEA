for idx = 1: length(ssi)
    figure;
    thisy = ssi(idx).y;
    Draw.spec2( thisy ,32000)
    hold on
    len = length( thisy )
    plot(pitch( thisy ,32000,'Method','SRH','WindowLength',round(len/20),'OverlapLength',round(len/60)))
end

autoArrangeFigures;