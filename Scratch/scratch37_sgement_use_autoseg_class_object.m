dbstop if error

wav_dir = "E:\WavsCollection";
subdirs = extract.folder(wav_dir);
tic
for r = 1:length(subdirs)
    
    
    birdid  = split(subdirs{r},'\');
    birdid = regexp(birdid{end},'\d+','match');
    birdid = str2double(birdid{1});
    if birdid > 300
          autoseg(subdirs{r}).autosegmenter_restrict_bout(2);
    
    fprintf('Now__%u__of %u dirs are processed',r,length(subdirs));
%     else
%         
%         if exist(sprintf('%s\\SegData',subdirs{r}), 'dir')
%             rmdir(sprintf('%s\\SegData',subdirs{r}));
%         end
        
    end
    
%    if ~isempty(extract.filename(sprintf('%s\\SegData',subdirs{r}),'*.mat'))
%        continue
%    end
  
end
toc

send_mail_message('379380788@qq.com','First section finished','Oops')

tic
parfor r = 1:length(subdirs)
    
    %    if ~isempty(extract.filename(sprintf('%s\\SegDataIndi',subdirs{r}),'*.mat'))
    %        continue
    %    end
     birdid  = split(subdirs{r},'\');
    birdid = regexp(birdid{end},'\d+','match');
    birdid = str2double(birdid{1});
    if birdid > 300
          autoseg(subdirs{r}).autoSegmenter;
    
    fprintf('Now__%u__of %u dirs are processed',r,length(subdirs));
    end
   
    
    
  
end

send_mail_message('379380788@qq.com','Second section finished','Oops')
toc
