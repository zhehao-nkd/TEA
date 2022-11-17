addpath(genpath("C:\Users\Zhehao\Dropbox (OIST)\My_Matlab\TEA"))
path_txt = "D:\New_Stimuli_after220901\R707_Ephys_New\P14\P14F001.txt"
path_pl2 = "D:\New_Stimuli_after220901\R707_Ephys_New\P14\P14F001.pl2"
path_folder = "D:\New_Stimuli_after220901\R707_Ephys_New\P14\P14F1_ConsSibs"

cho = Chorus(path_txt,path_pl2,path_folder);

singleunits = cho.getn;

% 
% collect = {};
% nametxt = {};
% for p = 1: length(singleunits) % for each single unit
%     % generate struct
%     singlestruct = singleunits{p}.collectImages;
%     collect{p} = singlestruct;
%     %nametxt{p} = singleunits{p}.neuronname;
% end
% 
% 
% cated = horzcat(ccollect{:});
% catedcated = horzcat(cated{:});

% 
catedcated = cho.collectImages

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

[~,pl2name,~] = fileparts(path_pl2);
imwrite(FINALIMG,sprintf('SimuRecorded_Neurons_%s.tiff',pl2name));