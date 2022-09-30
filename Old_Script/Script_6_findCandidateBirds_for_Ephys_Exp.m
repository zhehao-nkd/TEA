% A script to find out candidate  birds for Piece experiments
pathlog = "Z:\Yazaki-SugiyamaU\Bird-log_AK\Bird log2021 _ver_1.xlsx"
birdlog = table2struct(readtable(pathlog, 'Sheet','today'));

birdlog = birdlog(1:30); % 只取前30行
hidx = find(strcmp({birdlog(:).B_H}.', 'H')); % 只取 holding cages （H）
birdlog = birdlog(hidx);

birds = horzcat(birdlog(:).Var4);

tokens = regexp(birds,'(?<name>[OBRGY]\d{3})♂\((?<anno>\d+/\d+)\)','tokens');  % find male birds with birthday wrote
tokens = tokens';

candi = {};
for idx = 1: length(tokens)
    candi{idx} = tokens{idx}{1};
    %candi(idx).anno = tokens{idx}{2};
end
    
pathlist = "Z:\Yazaki-SugiyamaU\Bird-log_AK\Bird_List_new ver_2.xlsx";
birdlist = table2struct(readtable(pathlist));


s = 0; % shortlisted number
shortlist = struct;
for trump = 1: length(candi)
    
    thisbird = find(strcmp({birdlist(:).BirdID}.',candi{trump}));
    
    
    birthdate = datetime( num2str(birdlist(thisbird).HatchDate),'InputFormat','yyMMdd');
    currentdate = datetime('today');
    dph = days(currentdate - birthdate);
    
    if dph > 90
        s = s + 1;
        shortlist(s).id = candi{trump};
        shortlist(s).age = dph;
        
        var4 = {birdlog(:).Var4}.';
        idvar4 = find(~cellfun(@isempty,regexp(var4,candi{trump})));
        shortlist(s).cage =  birdlog(idvar4).Cage_;
        shortlist(s).father =  birdlist(thisbird).father;
    end
    
end


%remove those do not have recording folders
fdir = "Z:\Yazaki-SugiyamaU\Bird-song";

folders = cellstr(Extract.folder(fdir).');

tokens = regexp(folders,'(?<color>[OBRGY][A-Za-z]+)(?<number>\d{3})','tokens');

 

for b = 1: length(shortlist)
    
    %shortlist(b).id
    
    letter = regexp(shortlist(b).id,'[OBRGY]','match');
    letter = letter{1};
    number = regexp(shortlist(b).id,'\d+','match');
    number = number{1};
    
    format = sprintf('%s[A-Za-z]+%s',letter, number);
    loc = find(~cellfun(@isempty, regexp(folders,format,'match')));
    %folders{idx};
    if length(loc) == 1
        shortlist(b).record = 'Recorded';
    elseif length(loc) == 0
        shortlist(b).record = 'Not recorded yet';
    else 
        shortlist(b).record = '[ERROR]';
    end
    
    disp(sprintf('Candidate%u------------%s---------%udph---------Cage%u----------%s',...
    b,shortlist(b).id,shortlist(b).age,shortlist(b).cage,shortlist(b).record ));
    newline;
end


%%% show which bird are not recorded yet

newline;
notidx = find( strcmp({shortlist(:).record}.','Not recorded yet'));
notyet = shortlist(notidx);

for no = 1: length(notyet)
    

    disp(sprintf('Candidate%u------------%s---------%udph---------Cage%u----------%s',...
    no,notyet(no).id,notyet(no).age,notyet(no).cage,notyet(no).record ));
    newline;
end


