% scratch 6



%%%new batch for complex analysis


addpath(genpath("C:\Users\Zhehao\Dropbox (OIST)\My_Matlab\TEA"));
T = readtable("C:\Users\Zhehao\Dropbox (OIST)\My_EphysInput\selected.xlsx");
source = table2struct(T);

b = Batch(source(4).path_txt,source(4).path_plx,source(4).path_folder);
b.select;
neuronlist = b.getn;
 
 for k = 1: length(neuronlist)
     thisn = neuronlist{k};
     %thisn.rawthree;
     thisn.threesingle;
 end
 
 


 
 

% txt = 'E:\2021Recordings\Z0630_Y661\Y661_Z05.txt';
% plx = 'E:\2021Recordings\Z0630_Y661\Y661_Z05.plx';
% folder = 'Z:\Yazaki-SugiyamaU\Zhehao\My_Stimuli\Y661@06292021\Y661-1';

b1 = Batch(source(3).path_txt,source(3).path_plx,source(3).path_folder);
b1.select(1);
neuronlist = b1.getn;
list1 = neuronlist{1}.todisplay;




% txt = 'E:\2021Recordings\Z0630_Y661\Y661_Z06.txt';
% plx = 'E:\2021Recordings\Z0630_Y661\Y661_Z06.plx';
% folder = 'Z:\Yazaki-SugiyamaU\Zhehao\My_Stimuli\Y661@06292021\Y661-2';

b2 = Batch(source(4).path_txt,source(4).path_plx,source(4).path_folder);
b2.select(1);
neuronlist = b2.getn;
list2 = neuronlist{1}.todisplay;



listall = horzcat(list1,list2);


d = Display(listall);
d.showsyl('B606',mergedeleinf)
d.showdeg('B606')
d.showincre('B606')
d.showcatego('Y606')



% Display(listall).showdeg('Y606')
% 
% Display(listall).showsyl('G548')
% % figure;
% draw.spec( listall(5).y,32000)
% 
% figure;
% 
% draw.spec( listall(92).y,32000)

