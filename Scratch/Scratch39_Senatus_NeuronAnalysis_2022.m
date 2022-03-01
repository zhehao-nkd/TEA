% R693 -Z01

path_txt = "E:\2022Recordings\R693@02220224\R693_Z01_Neuron1_Cons_2d5s_sub.txt";
path_plx = "E:\2022Recordings\R693@02220224\R693_Z01_Neuron1_Cons_2d5s_sub.plx";
path_folder = "Z:\Yazaki-SugiyamaU\Zhehao\My_Stimuli\R693@02232022\R693-1-Cons-2d5s";


% R693 -Z02

path_txt = "E:\2022Recordings\R693@02220224\R693_Z02_Neuron1_Degs_Fcall_2d5s_sub.txt";
path_plx = "E:\2022Recordings\R693@02220224\R693_Z02_Neuron1_Degs_Fcall_2d5s_sub.plx";
path_folder = "Z:\Yazaki-SugiyamaU\Zhehao\My_Stimuli\R693@02232022\R693-2-Degs-Fcall-2d5s";

% R693 -Z03


path_txt = "E:\2022Recordings\R693@02220224\R693_Z03_Neuron1_Details_Fcall1st_2d5s_sub.txt";
path_plx = "E:\2022Recordings\R693@02220224\R693_Z03_Neuron1_Details_Fcall1st_2d5s_sub.plx";
path_folder = "Z:\Yazaki-SugiyamaU\Zhehao\My_Stimuli\R693@02232022\R693-3-Details-Fcall-1st-2d5s";



% R693 -Z04 E:\2022Recordings\R693@02220224\R693_Z04_Neuron2_Norms_4s_sub.txt

path_txt = "E:\2022Recordings\R693@02220224\R693_Z04_Neuron2_Cons_2d5s.txt";
path_plx = "E:\2022Recordings\R693@02220224\R693_Z04_Neuron2_Cons_2d5s_sub.plx";
path_folder = "Z:\Yazaki-SugiyamaU\Zhehao\My_Stimuli\R693@02232022\R693-1-Cons-2d5s";

%R693 -Z05 "E:\2022Recordings\R693@02220224\R693_Z05_Neuron2_Degs_2d5s_sub.plx"


path_txt = "E:\2022Recordings\R693@02220224\R693_Z05_Neuron2_Degs_Y616_2d5s.txt";
path_plx = "E:\2022Recordings\R693@02220224\R693_Z05_Neuron2_Degs_Y616_2d5s_sub.plx";
path_folder = "Z:\Yazaki-SugiyamaU\Zhehao\My_Stimuli\R693@02232022\R693-5-Degs-Y616-2d5s";

% Test 1 "E:\2022Recordings\R693@02220224\Test_R693_Z04_Neuron2_Cons_2d5s_sub.plx"
path_txt = "E:\2022Recordings\R693@02220224\Test_R693_Z04_Neuron2_Cons_2d5s.txt";
path_plx = "E:\2022Recordings\R693@02220224\Test_R693_Z04_Neuron2_Cons_2d5s_sub.plx";
path_folder = "Z:\Yazaki-SugiyamaU\Zhehao\My_Stimuli\R693@02232022\R693-1-Cons-2d5s";

% Test1 Deg "E:\2022Recordings\R693@02220224\Test1_R693_Degs_2d5s_sub.plx"
path_txt = "E:\2022Recordings\R693@02220224\Test1_R693_Degs_2d5s_sub.txt";
path_plx = "E:\2022Recordings\R693@02220224\Test1_R693_Degs_2d5s_sub.plx";
path_folder = "Z:\Yazaki-SugiyamaU\Zhehao\My_Stimuli\R693@02232022\Test1-R693-Degs-2d5s";

% addpath genpath
addpath(genpath("C:\Users\Zhehao\Dropbox (OIST)\My_Matlab\TEA"))
% 核心分析代码
b = Batch(path_txt,path_plx,path_folder);
b.select;
neuronlist = b.getn;

for k = 1: length(neuronlist)
    thisn = neuronlist{k};
    thisn.three;
    thisn.rawthree;
    %thisn.ResponseBasedOrderedThreePlots;
end

A = Analysis(neuronlist{1});
A.drawAlignedNormFragTwoPlots;

