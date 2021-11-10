 % Y687
 path_txt = "E:\2021Recordings\Z1027_Y687\Y687_Z31.txt";
 path_plx = "E:\2021Recordings\Z1027_Y687\Y687_Z31.plx";
 path_folder = "Z:\Yazaki-SugiyamaU\Zhehao\My_Stimuli\O686@11022021\O686-1";

b = Batch(path_txt,path_plx,path_folder);

b.select;

neuronlist = b.getn;
 
 for k = 1: length(neuronlist)
     thisn = neuronlist{k};
     thisn.rawthree;
     thisn.threesingle;
 end
 
 
 
 %%%%%%
 path_txt = "E:\2021Recordings\Z1102_O686\O686_33.txt";
 path_plx = "E:\2021Recordings\Z1102_O686\O686_33.plx";
 path_folder = "Z:\Yazaki-SugiyamaU\Zhehao\My_Stimuli\O686@11022021\O686-20";

b = Batch(path_txt,path_plx,path_folder);

b.select;

neuronlist = b.getn;
 
 for k = 1: length(neuronlist)
     thisn = neuronlist{k};
     %thisn.rawthree;
     thisn.threesingle;
 end