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
            outdir = sprintf('1ch_%s',rawdir);
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
            
    end
end

