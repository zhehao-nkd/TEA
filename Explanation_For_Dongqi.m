% Explanation

% Ephys相关数据被我储存在三种数据文件中
% .plx 文件是原始文件，External中的 readPLXFileC 是用来读取 .plx文件的,可以在 matlab file exchange上找到 
% .txt文件储存spikes的时间和波形 
% .wav文件其中一个通道是声刺激，另一个通道是用来标定声刺激播放时间的trigger signal


%文件储存结构
plxpath = "D:\Ephys-O706-W\P01\O706_P01F1_5sCons.plx"
txtpath = "D:\Ephys-O706-W\P01\O706_P01F1_5sCons.txt"
stimulidir = "D:\Ephys-O706-W\P01\Stimuli\O706-P01F1-5sCons"
% Bird---O706   P01----Plexon system Recording01 P01F1---对这批神经元进行的第一种刺激


% 数据用两种采集系统进行采集，P指的是plexon系统,sample frequency=40000 Z指的是zeus系统,sample rate=30000
% Trigger signal用方波间距来编码stimuli代号，从右向左读，长为1，短为0，二转十进制后得到stimuli代号
% 两个系统摄入Trigger signal的方式也不同,P系统通过Analog Input通道接收，Z系统通过digital Input接收
 t = Trigger(path_plx) % 可用来读取trigger信息
 
