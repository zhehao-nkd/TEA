
path_txt = "D:\R705\P18\P07F002.txt"
path_plx = "D:\R705\P07\P07F002.plx"
path_folder = "D:\R705\P07\P07F2_SibDegs_R705"

ba= Batch(path_txt,path_plx,path_folder);
ba.select;
singleunits = ba.getn;


collect = {};
nametxt = {};
parfor p = 1: length(singleunits) % for each single unit
    % generate struct
    singlestruct = singleunits{p}.getImageInfo;
    collect{p} = singlestruct;
    %nametxt{p} = singleunits{p}.neuronname;
end


cated = horzcat(ccollect{:});
catedcated = horzcat(cated{:});

rows = unique({catedcated.channelunit}.');
[columns,cindex] = unique(cellstr({catedcated.soundpath}.'));

specrow = {};
for ind = 1: length(cindex)
    specrow{1,ind} = catedcated(cindex(ind)).specimg;
end

size3 = size(catedcated(1).rasterimg);
finalimg = repmat({uint8(255*ones(size3(1),size3(2),size3(3)))}, length(rows), length(columns));



for dd = 1: length(catedcated)
    
    rownum = find(~cellfun(@isempty,regexp(rows,catedcated(dd).channelunit)));
    [columnnum,~] = find(ismember(columns,convertStringsToChars(catedcated(dd).soundpath)));
    finalimg{rownum,columnnum} = catedcated(dd).rasterimg;
end

finalimg = vertcat(specrow,finalimg);
FINALIMG = cell2mat(finalimg);

[~,plxname,~] = fileparts(path_plx);
imwrite(FINALIMG,sprintf('SimuRecorded_Neurons_%s.tiff',plxname));