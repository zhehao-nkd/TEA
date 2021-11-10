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
            % 3 to 6 letter, red to orange{3,6}
            if length(captured) ~= 1
                warning('%s may not be correctly converted!',rawid);
                newid = rawid;
                return
            end
           
            temp = convertStringsToChars(captured.letter);
            
            newletter = upper(temp(1));
            newid = [newletter,captured.number];
                    
        end
            
        function two2one(dir) % convert 2 channel to 1 channel wav
            [~,rawdir,~] = fileparts(dir);
            outdir = sprintf('1channel_%s',rawdir);
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

