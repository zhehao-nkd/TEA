classdef Extract
    %EXTRACT Summary of this class goes here
    %   Detailed explanation goes here
    
    methods(Static)
        
        function features = feature(y,fs)
            
            features = SAT_hijack(y,fs).features;
            % fix the pitch
            temp = features.pitch;
            newpitch = features.pitch;
            newpitch(newpitch > 2000) = nan;
            newpitch(newpitch < 400) = nan;
            features.pitch = newpitch;
            features.rawpitch = temp;
            
            % fix the fm
            temp2 = features.FM;
            newfm = [diff(features.pitch),0];
            
            features.fm = newfm;
            features.rawfm = temp2;
            
            
        end
        
        function filenames = filename(dirpath,extension) % Extract filename only in dirpath folder, but not in target folder
            
            switch class(dirpath)
                
                case {'char','string'}
                    
                    fileExt = extension;
                    files = dir(fullfile(dirpath,fileExt));
                    filenames = {};
                    coder.varsize(filenames);
                    for n=1:size(files,1)
                        filenames{n} = fullfile(dirpath,files(n,1).name);
                    end
                    
                case 'cell'
                    
                    filenames = {};
                    for k = 1:length(dirpath)
                        fileExt = extension;
                        files = dir(fullfile(dirpath{k},fileExt));
                        for n=1:size(files,1)
                            filenames_sub{n} = fullfile(dirpath{k},files(n,1).name);
                        end
                        filenames{k} = filenames_sub;
                    end
                    filenames = horzcat(filenames{:});
            end
            
            filenames = filenames.'; % Trasnpose only for convinient checking variables
            
            fclose('all') ;
         end
        
        function filenames = Classical_filename(dirpath,extension)
            fileExt = extension;
            dirpath = convertCharsToStrings(dirpath);
            files = dir(fullfile(dirpath,fileExt));
            len = size(files,1);
            
            filenames = {};
            coder.varsize(filenames);
            
            for n=1:len
                filenames{n} = fullfile(dirpath,files(n,1).name);
            end
            
            filenames = filenames.'; % Trasnpose only for convinient checking variables
        end
        
        function subDirsNames = folder(parentDir)
            %GetFolders
            % 函数功能为获取父文件夹下所有子文件夹的路径
            % 函数的输入为ParentFolder：父文件夹路径。eg: 'D:\Program Files'
            % 函数的输出为SubFolders：子文件夹路径。为一个元胞数组，eg: {'D:\Program Files\FileZilla FTP Client\docs'}
            if ~strcmp(parentDir(end),'\')
                parentDir = strcat(parentDir,'\');
            end
            files = dir(parentDir);
            % Get a logical vector that tells which is a directory.
            dirFlags = [files.isdir];
            % Extract only those that are directories.
            subDirs = files(dirFlags); % A structure with extra info.
            % Get only the folder names into a cell array.
            subDirsNames = {subDirs(3:end).name}.';
            subDirsNames(ismember(subDirsNames ,'..')) = [];
            subDirsNames(ismember(subDirsNames ,'.')) = [];
            for k = 1:length(subDirsNames)
                subDirsNames{k} = strcat(parentDir,subDirsNames{k});
            end

            subDirsNames = cellstr(subDirsNames);

        end
        
        function newspt = cutspt(sptimes,zpt,ylen)
            sptimes = cellfun(@(x) x - zpt, sptimes, 'un', 0);
            sptimes = cellfun(@(v) v( v > 0), sptimes, 'UniformOutput', false);
            newspt = cellfun(@(v) v( v < ylen), sptimes, 'UniformOutput', false);
        end
        
        function birdid = birdid(dirpath)
            fileExt = '*.wav';
            files = dir(fullfile(dirpath,fileExt));
            len = size(files,1);
            
            filenames = {};
            coder.varsize(filenames);
            
            for n=1:len
                filenames = fullfile(dirpath,files(n,1).name);
                [~,rawid,~] = fileparts(filenames);
                parts = strsplit(rawid,'_');
                birdid{n} = parts{1};
            end
            
            birdid = birdid.'; % Trasnpose only for convinient checking variables
            
        end
        
        function value = json(json_path)
            % read json file
            fname = json_path;
            fid = fopen(fname);
            raw = fread(fid,inf);
            str = char(raw');
            fclose(fid);
            value = jsondecode(str);
            
        end
        
        function list = all(path) % Extract everything files/folders in the path folder
            list = dir(fullfile(path,'**/*.*'));
        end
        
        function folders = foldersAllLevel(path) % Extract everything files/folders in the path folder
            list = dir(fullfile(path,'**/*.*'));
            
            dirids = find([list.isdir].');
            
            temp = {list.name}.';
            dot_ids = find(~cellfun(@isempty,regexp({list.name}.','\.\.')));% Too bad code
            
            %dot_list = list(dot_ids);
            
            screened_ids = intersect(dirids,dot_ids);
            
            screened_list = list(screened_ids);
            
            folders = {};
            for k = 1: length( screened_list)
                
                folders{k} = screened_list(k).folder;
                
            end
            
            folders = folders.';
            
            
        end
        
        function files = filesAllLevel(path, extension)
            
            folders = Extract.foldersAllLevel(path);
            
            for k = 1:length(folders)
                
                subdir_files{k} = Extract.filename(folders{k},extension);
            end
            
            files = vertcat(subdir_files{:});
            
        end
        
        function info = fileDirRelation(path,extension) % Available for two layers of folders
            
            subfolders = Extract.folder(path).';
            
            info = struct;
            ids = 0;
            for k = 1:length(subfolders)
               
                wavfiles = Extract.filename(subfolders{k},extension);
                
                for p = 1: length(wavfiles)
                     ids = ids + 1;
                     info(ids).fileid = ids;
                     info(ids).filename = wavfiles{p};
                     info(ids).folder = subfolders{k};
                    
                end
                
            end

        end
        
        function new_sptimes = sptimes(spike_times, initial, terminal)
            
            new_sptimes = {};
            for k = 1: length(spike_times)
                local_spt = spike_times{k};
                if isempty(local_spt)
                    new_sptimes{k} = [];
                    continue
                end
                ids = find((local_spt>initial )& (local_spt<terminal));
                if isempty(ids)
                    new_sptimes{k} = [];
                    continue
                end
                new_sptimes{k} = local_spt(ids);
            end
            
        end
        
        function new_sptimes = sptimes_resetSP(spike_times, initial, terminal)
            new_sptimes = {};
            for k = 1: length(spike_times)
                local_spt = spike_times{k};
                if isempty(local_spt)
                    new_sptimes{k} = [];
                    continue
                end
                ids = find((local_spt>initial )& (local_spt<terminal));
                if isempty(ids)
                    new_sptimes{k} = [];
                    continue
                end
                new_sptimes{k} = local_spt(ids);
            end
            
            new_sptimes = Convert.sptimesOnset2Zero(new_sptimes, initial);
            
        end
        
        function parent_dir = parentDir(input_dir) % 找到一个dir 的 parent dir
            
           % input_dir = "C:\Users\Zhehao\Desktop\ComplexFolder\P01"
           f = filesep;
           parts = split(input_dir,f);
           parentparts = parts(1:length(parts)-1);
           parent_dir = fullfile(parentparts{:});

        end
              
    end  
   
end

