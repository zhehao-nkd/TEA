% R693 -Z01

path_txt = "D:\Zhehao\R695_0531\P01\R695_P01F001.txt";
path_plx = "D:\Zhehao\R695_0531\P01\R695_P01F001.plx";
path_folder = "D:\Zhehao\R695_0531\P01\R695_Cons4s";


path_txt = "D:\Zhehao\R695_0531\P01\R695_P01F002.txt";
path_plx = "D:\Zhehao\R695_0531\P01\R695_P01F002.plx";
path_folder = "D:\Zhehao\R695_0531\P01\R695_P01F2_Degs"


%P2 D:\Zhehao\R695_0531\P02\R695_P02F001.txt

path_txt = "D:\Zhehao\R695_0531\P02\R695_P02F001.txt";
path_plx = "D:\Zhehao\R695_0531\P02\R695_P02F001.plx";
path_folder = "D:\Zhehao\R695_0531\P01\R695_P01F1_Cons4s";


% addpath genpath "D:\ZhehaoScripts\Dropbox (OIST)\My_Matlab\TEA"
addpath(genpath("D:\ZhehaoScripts\Dropbox (OIST)\My_Matlab\TEA"))
% 核心分析代码
b = Batch(path_txt,path_plx,path_folder);
b.select;
neuronlist = b.getn;

for k = 1: length(neuronlist)
    thisn = neuronlist{k};
    %thisn.three;
    thisn.rawthree;
    %thisn.ResponseBasedOrderedThreePlots;
end

A = Analysis(neuronlist{1});
A.drawAlignedNormFragTwoPlots;

