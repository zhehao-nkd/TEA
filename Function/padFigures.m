%place matrices into a cell 
%matrices = {rand(5,13),rand(4,7),rand(5,10)};
function outimg = padFigures(matrices)
% 补齐不同的figure使它们拥有同样的size

sz1=cellfun('size',matrices,1);
mx = max(sz1);
%pad each matrix
padded = cellfun(...
    @(M)[M(1:end,:,:);256*ones(mx-size(M,1),size(M,2),size(M,3))],...
    matrices,...
    'UniformOutput', false...
    );

%find maxmum of number of columns of matrices
sz2=cellfun('size',padded,2);
mx = max(sz2);
%pad each matrix
paddedpadded = cellfun(...
    @(M)[M(:,1:end,:),256*ones(size(M,1),mx-size(M,2),size(M,3))],...
    padded,...
    'UniformOutput', false...
    );
            
%concatenate matrices
outimg = cell2mat(paddedpadded);


end