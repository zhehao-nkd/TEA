path_txt = "D:\R707\P01\R707_P01F001_mrg.txt"
path_plx = "D:\R707\P01\R707_P01F001_mrg.plx"
path_wav = "D:\R707\Sib_Merge"

ba = Batch(path_txt,path_plx,path_wav);

ba.select;
Ns = ba.getn;

for k = 1:length(Ns)
   A = Analysis(Ns{k});
   save(A.neuronname,A,'-v7.3');
end