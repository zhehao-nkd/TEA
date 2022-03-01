filenames = extract.filename("C:\Users\Zhehao\Dropbox (OIST)\My_Stimuli\R677@11232021\senatusOneMotif\SegData",'*.mat')


for w = 1: length(filenames)
    reviseSeg(filenames{w});
end