addpath(genpath("Z:\Zhehao\Dropbox (OIST)\My_Matlab\TEA"))

dbstop if error;
subdirs = extract.folder("Y:\Yazaki-SugiyamaU\Bird-song").';


wb = waitbar(0,'Start processing');
num_files = length(subdirs);
utl.UpdateParforWaitbar(num_files, wb);
D = parallel.pool.DataQueue;
afterEach(D, @utl.UpdateParforWaitbar);

for k = 1: length(subdirs)
    getInf.moveFromBucket(subdirs{k});
    send(D,1);
end