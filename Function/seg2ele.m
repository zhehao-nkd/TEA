function eleinf = seg2ele(segdata)

all_collect = {};
song_collect = {};


two_eleedge = repmat(segdata.eleedge,[2,1]);
alledges = sort( vertcat(segdata.syledge(:),reshape(two_eleedge,[],1) ));

% the following is a very bad temporily code
if ~isfield(segdata,'fs') % if there is no such field named fs
    fs = 32000;
else
    fs = segdata.fs;
end

initials = alledges(1:2:end);
terminals = alledges(2:2:end);

song_eleinf = struct;

for w = 1: length(initials) % can add a par-for
    song_eleinf(w).initial = initials(w)*fs/1000;
    song_eleinf(w).terminal = terminals(w)*fs/1000;
    song_eleinf(w).songname = segdata.birdid;
    
    if isfield(segdata,'y')
        song_eleinf(w).y = segdata.y(initials(w)*fs/1000:terminals(w)*fs/1000);  % very bad code, but meibanfa
    elseif isfield(segdata,'rawy')
        song_eleinf(w).y = segdata.rawy(initials(w)*fs/1000:terminals(w)*fs/1000);
    end
    
    song_eleinf(w).fs = fs;
end

parfor f = 1:length(song_eleinf)
    song_eleinf(f).fragid = f;
end

eleinf = song_eleinf;

end