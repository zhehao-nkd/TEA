classdef extract
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
        
        function filenames = filename(dirpath,extension)
            fileExt = extension;
            files = dir(fullfile(dirpath,fileExt));
            len = size(files,1);
            
            filenames = {};
            coder.varsize(filenames);
            
            for n=1:len
                filenames{n} = fullfile(dirpath,files(n,1).name);
            end
            
            filenames = filenames.'; % Trasnpose only for convinient checking variables
        end
        
        function [SubFolders] = folder(ParentFolder)
            %GetFolders
            % 函数功能为获取父文件夹下所有子文件夹的路径
            % 函数的输入为ParentFolder：父文件夹路径。eg: 'D:\Program Files'
            % 函数的输出为SubFolders：子文件夹路径。为一个元胞数组，eg: {'D:\Program Files\FileZilla FTP Client\docs'}
            
            SubFolderNames = dir(ParentFolder);
            for i=1:length(SubFolderNames)
                if( isequal( SubFolderNames( i ).name, '.' )||...
                        isequal( SubFolderNames( i ).name, '..')||...
                        ~SubFolderNames( i ).isdir) % 如果不是目录则跳过
                    continue;
                end
                SubFolder(i).SubFolderName = fullfile( ParentFolder, SubFolderNames( i ).name );
            end
            
            temp = {SubFolder.SubFolderName};
            idx = cellfun(@(x)~isempty(x),temp,'UniformOutput',true); % 利用cellfun函数得到元胞数组中所有非空元素的下标
            SubFolders = temp(idx);
            
        end
        
        function newspt = cutspt(sptimes,zpt,ylen)
           sptimes = cellfun(@(x) x - zpt, sptimes, 'un', 0);
           sptimes = cellfun(@(v) v( v > 0), sptimes, 'UniformOutput', false);
           newspt = cellfun(@(v) v( v < ylen), sptimes, 'UniformOutput', false);
        end
    end
end

