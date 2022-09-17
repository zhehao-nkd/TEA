classdef convert
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    
    methods(Static)
        
        function  c2mex()
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            dir = uigetdir("C:\Users\Zhehao\Dropbox (OIST)\My_Matlab\NewMatlab\SAT")
            cfiles = getFilenames(dir,'*.c');
            eval(['cd ','''',dir,'''']);

            for n = 1: length(cfiles)
                [~,name,ext] = fileparts(cfiles{n});
                filename = [name,ext];
                eval(['mex ', filename])
            end
            
        end
        
        function newid = bid(rawid,which_to_keep) % used for normalize the bird id to B345-like str
            
            captured = regexp(rawid,'(?<letter>[OGBYR][a-zA-Z]*)\s*(?<number>\d{3})?(?<spe>TUT|BOS|V2|Fcall|Mcall|WNS|Het)?','names');
            if isempty(captured) % 如果没有发现B521这种ID，那么继续找有没有特殊标签
                captured = regexp(rawid,'(?<letter>[OGBYR][a-zA-Z]*)\s*(?<number>\d{3})?|(?<spe>TUT|BOS|V2|Fcall|Mcall|WNS|Het)?','names');
            end
            % 之所以用上面letter的regexp方式是因为大小不同，枚举过难
            % 而 \s*是为了防止空格
            if isempty(captured)% If nothing is detected ruturn the warning message
                newid = rawid; fprintf('Message@ convert.bid: No change, keep raw input--%s\n',rawid); return;
            end
            
            if length(captured) > 1 % if More than 1, choose which to keep
                captured = captured(which_to_keep);
            end
            
            % Switch to 1
            newid = upper(strcat(captured.letter(1),captured.number,captured.spe)); % 取首字母
            
        end
            
        function two2one(dir) % convert 2 channel to 1 channel wav
            [parentdir,rawdir,~] = fileparts(dir);
            rawoutdir = sprintf('SingleChannel_%s',rawdir);
            outdir = fullfile(parentdir,rawoutdir);
            mkdir(outdir);
            
            files = extract.filename(dir,'*.wav');
            for idx = 1: length(files)
                [y,fs] = audioread(files{idx});
                ys = y(:,2);
                [~,name,ext] = fileparts(files{idx});
                newfilepath = sprintf('%s/%s%s',outdir,name,ext);
                audiowrite(newfilepath,ys,fs);
            end
        end
        
        function two2oneAllLevel(dir) % convert 2 channel to 1 channel wav, no matter how many levels of folders there are
            dbstop if error
            f = waitbar(0,'Please wait...');
            
            % generat the new parent dir
            [~,old,~] = fileparts(dir);
            new = sprintf('SingleChannel_%s',old);
            
            %  create the folders
            allfolders = extract.foldersAllLevel(dir);
            for k = 1: length(allfolders)
                New_allfolders{k} = replace(allfolders{k},old,new);
                mkdir(New_allfolders{k});
                waitbar(k/length(allfolders),f,'Creating folders...');
            end
            
            % re write the files
            waitbar(0,f,'Rewriting files...');
            
            allfiles = extract.filesAllLevel(dir,'*.wav');
            
            for k = 1: length(allfiles)
                New_allfiles{k} = replace(allfiles{k},old,new);
                [y,fs] = audioread(allfiles{k});
                ys = y(:,2); % y-signal
                audiowrite(New_allfiles{k},ys,fs);
                waitbar(k/length(allfiles),f,'Creating folders...');
            end
            
            close(f);
            
        end
          
        function new_sptimes = sptimesOnset2Zero(sptimes, onset_time)
            % convert the onset of sptimes to zeros
            for k = 1: length(sptimes)
                new_sptimes{k} = sptimes{k} - onset_time;
            end
            
        end
        
        function rename(dirpath,ext,old,new)
            
            
            files = dir(sprintf('%s%s%s',dirpath,'\',ext));
            % Loop through each
            for id = 1:length(files)
                % Get the file name (minus the extension)
                [~, f,ext] = fileparts(files(id).name);
                % replace _ with -
                
                targets = regexp(f,old);
                
                f(targets) = new;
                
                newname = sprintf('%s%s',f,ext);
                
                
                
                % If not the same, rename
                if ~strcmp(fullfile(files(id).folder,files(id).name), fullfile(files(id).folder,newname))
                    movefile(fullfile(files(id).folder,files(id).name), fullfile(files(id).folder,newname));
                end
                
            end
            
            
        end
        
        function cutwav(dir,initial,terminal)  % in order to cut off silence in wavs
            
            % read wav and remove the emptys  
         %dir = "C:\Users\Zhehao\3D Objects\wavs";initial = 2.99; terminal = 3.09;
          
            [~,rawdir,~] = fileparts(dir);
            outdir = sprintf('cutwav_%s',rawdir);
            mkdir(outdir);
  
            files = extract.filename(dir,'*.wav');
            
            for idx = 1: length(files)
                [yraw,fs] = audioread(files{idx});
                ynew = yraw(initial*fs:terminal*fs);
                 [~,name,ext] = fileparts(files{idx});
                newfilepath = sprintf('%s/%s%s',outdir,name,ext);
                audiowrite(newfilepath,ynew,fs);
            end
            
            
        end
        
        function destineydir = mergeSubfolders(complexfolder,ext)
            
            % "C:\Users\Zhehao\Desktop\ComplexFolder"
            
            FS = filesep;
            parts = split(complexfolder,FS);
            parts(length(parts)) = strcat('Merged_',parts(length(parts)));
            destineydir = fullfile(parts{:});
            mkdir(destineydir);
            files = extract.filesAllLevel(complexfolder,ext);
            for k = 1: length(files)
                
              [~,name,ext] = fileparts(files{k});
              
              if strcmp(ext, '.') || strcmp(ext, '..')
                  continue
              end
              
              nameext = strcat(name,ext);
              SourceFile = files{k};
              DestinyFile = fullfile(destineydir,nameext);
              copyfile(SourceFile, DestinyFile, 'f');
              
            end
        end
        
        function renameAllLevel(targetdir,ext,old,new)
            dbstop if error
            files = extract.filesAllLevel(targetdir, ext);
            
            wb = waitbar(0,'Processing...');
            
            for id = 1:length(files)
                
                waitbar(id/length(files),wb,sprintf('当前是所有%u个文件中的%u',length(files),id))
                % Get the file name
                oldname = files{id};
                
                [folder,oldfilename,ext] = fileparts(oldname);
                newfilename = strrep(oldfilename,old,new);
                
                newname = fullfile(folder,strcat(newfilename,ext));
                if ~strcmp(oldname,newname)
                movefile(oldname, newname,'f');
                end
            end
            
            close(wb);

        end   
        
        function newpath = path(oldpath,handle)
            switch handle
                case 'win' % 非常简陋，不一定是对的
                    newpath = strrep(oldpath,'/','\');
                case 'unix'
                    %convert to a unix version:
                    newpath = strrep(oldpath,'\','/');
            end
        end
        
        function colorifyImage(originalImage)
            
            figure
            subplot(2,1,1);
            imshow(originalImage);
            % loops are unnecessary
            % your mask does not depend on color selection
            % and your color selection does not select what you think it selects
            % these masks (very) roughly select the chips in the image
            maskR = originalImage(:,:,1) > 200 & originalImage(:,:,2) < 100 & originalImage(:,:,3) < 100;
            maskG = originalImage(:,:,1) < 50 & originalImage(:,:,2) > 100 & originalImage(:,:,3) < 150;
            maskB = originalImage(:,:,1) < 10 & originalImage(:,:,2) < 100 & originalImage(:,:,3) > 220;
            maskW = originalImage(:,:,1) == 255 & originalImage(:,:,2) == 255 & originalImage(:,:,3) == 255;
            % i'm not dealing with a tedious prompt.  feel free to change
            colormode = 'white';
            % you can also specify a color other than black
            newcolor = [220 0 0];
            outpict = originalImage;
            switch lower(colormode)
                case 'red'
                    selectedmask = repmat(maskR,[1 1 3]);
                case 'green'
                    selectedmask = repmat(maskG,[1 1 3]);
                case 'blue'
                    selectedmask = repmat(maskB,[1 1 3]);
                case 'white'
                    selectedmask = repmat(maskW,[1 1 3]);
                otherwise
                    error('ha ha you typed the wrong thing and got an error')
            end
            outpict(selectedmask) = 0;
            outpict = outpict + uint8(selectedmask.*permute(newcolor,[1 3 2]));
            % of course, that depends on the class of the input image
            subplot(2,1,2);
            imshow(outpict);
            
        end
        
        function originalImage = colorEdge(originalImage,color)
            
            switch color
                case 'r'
                    targetcolor = uint8([225,0,0]);
                case 'b'
                    targetcolor = uint8([0,0,225]);
                case 'g'
                    targetcolor = uint8([0,225,0]);
            end
            
            range = 9;
            
            for k = 1: size(originalImage,2)
                for kk = 1: range
                    originalImage(kk,k,:) = targetcolor;
                    originalImage(size(originalImage,1)-kk+1,k,:) = targetcolor;
                end
            end
            
            for k = 1: size(originalImage,1)
                for kk = 1: range
                    originalImage(k,kk,:) = targetcolor;
                    originalImage(k,size(originalImage,2)-kk+1,:) = targetcolor;
                end
            end
         
        end
        
        function Iall = mergeImage(Icollect)
            % Merge images in one, Icollect is {5x1 cell,4x1 cell, 3x1 cell}
            I_eachColumn = {};
            for k = 1:length(Icollect)
                temp = Icollect{k};
                I_eachColumn{k} = vertcat(temp{:});
            end
            
            size1 = [];
            for oo = 1: length(I_eachColumn)
                size1(oo) = size(I_eachColumn{oo},1);
            end
            
            [max_size1,max_oo] = max(size1);
            
            Ipad = {};
            for oo = 1: length(I_eachColumn)
                localI = I_eachColumn{oo};
                Ibase= uint8(256*ones( size(I_eachColumn{max_oo} ) ));
                Ibase(1:size(localI,1),1:size(localI,2),1:size(localI,3)) = localI;
                Ipad{oo} = Ibase;
            end
            
            Iall = horzcat(Ipad{:});
            
        end
        
        function breakdownEphysTxt(txtfile)
            
            [pathstr, oldname, ext] = fileparts(txtfile);
            
            bid = regexp(oldname,'[OGBYR]\d{3}','match');
            zpid = regexp(oldname,'[ZP]\d+F\d+','match');
            
            
            rawT = readtable(txtfile);% sorted plus unsorted
            fid = fopen(txtfile);
            frewind(fid);
            first_line = utl.deblankl(fgetl(fid)); % deblankl for removal of blanks, as VariableNames can't include blanks
            titlecell = split(first_line,',');
            lentitle = length(titlecell);
            rawT.Properties.VariableNames(1:lentitle) = titlecell;
            
            temp = unique(rawT.channelname);
            spkc_ids = find(~cellfun(@isempty,regexp(temp,'SPKC')));
            spkc_names = {temp{spkc_ids}}.';
            
            for k = 1:length(spkc_names)
                spkc_num = regexp(spkc_names{k},'(?<=SPKC)\d{2}','match');
                num_ids = find(~cellfun(@isempty,regexp(rawT.channelname,spkc_num{1})));
                subsetT = rawT(num_ids,:);
                unitnames = unique(subsetT.unit);
                unitnames( unitnames==0) = []; % Remove the unit zero,cause everyone has it
                %if length(unique(subsetT.unit)) > 1 % for this channel, there are at least some sorted spikes
                for kk = 1: length(unitnames)
                    
                    if ~isempty(bid)&&~isempty(zpid)
                    newname = strjoin({convertStringsToChars(bid),convertStringsToChars(zpid),spkc_names{k},...
                        sprintf('U%u',unitnames(kk))},'_'); % U means unit
                    else
                        newname = strjoin({convertStringsToChars(oldname),spkc_names{k},...
                        sprintf('U%u',unitnames(kk))},'_'); % if filename is not in the correct format, then keep using the oldname
                    end
                    subsetname = fullfile(pathstr,strcat(newname,ext));
                    writetable(subsetT,subsetname);
                end
               % end
            end
              
        end
        
    end
    
end

