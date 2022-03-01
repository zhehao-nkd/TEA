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
        
        function newid = bid(rawid) % used for normalize the bird id to B345-like str
            
            captured = regexp(rawid,'(?<letter>[a-zA-Z]+)(?<number>\d+)','names');
            special_characters  = regexp(rawid,'(?<letter>TUT|BOS|Tut|Bos)','names');
            if length(special_characters)~= 0
                disp('It is a TUT or BOS!');
                temp = convertStringsToChars(captured.letter); 
                newletter = upper(temp(1));
                newid = strcat(newletter,captured.number,special_characters.letter);
                return
            end
            % 3 to 6 letter, red to orange{3,6}
            if length(captured) ~= 1
                captured2 = regexp(rawid,'(?<letter>Fcall|Mcall|WNS|Het)','names');
                
                if length(captured2) == 1
                    fprintf('hey! current id is a %s',captured2.letter);
                    newid = captured2.letter;
                    return
                else
                    warning('%s may not be correctly converted!',rawid);
                    newid = rawid;
                    return
                end
            end
           
            temp = convertStringsToChars(captured.letter);
            
            newletter = upper(temp(1));
            newid = strcat(newletter,captured.number);
                    
        end
            
        function two2one(dir) % convert 2 channel to 1 channel wav
            [~,rawdir,~] = fileparts(dir);
            outdir = sprintf('SingleChannel_%s',rawdir);
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
        
    end
end

