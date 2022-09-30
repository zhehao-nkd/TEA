addpath(genpath("C:\Users\Zhehao\Dropbox (OIST)\My_Matlab\TEA"));

T = readtable("C:\Users\Zhehao\Dropbox (OIST)\My_EphysInput\Y675.xlsx")


for k = 56
%for k = 30: height(T)
    
    temp = T.path_txt(k);
    txt = temp{1};
    temp = T.path_plx(k);
    plx = temp{1};
    temp = T.path_folder(k);
    folder = temp{1};
    
    try
        
    b = Chorus(txt,plx,folder);
    b.select;
    neuronlist = b.getn;
    
    for m = 1: length(neuronlist)
        thisn = neuronlist{m};
        thisn.rawthree;
        %thisn.Three;
    end
    catch error
    end


end