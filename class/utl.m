classdef Utl
    % store some utilis functions that I have not yet decided where to put
    
  
    methods(Static)
        
        
        function padded = pad(y,targetlen)
            % pad zero to data in both side
            % here length(y) must < targetlen
            y = y(:); % force y to be column vector
            if mod(targetlen - length(y),2) == 0 % if even number
                frontside = (targetlen - length(y))/2;
                backside = (targetlen - length(y))/2; % length in one side
            else
                frontside = (targetlen - length(y)-1)/2;
                backside = targetlen - length(y)-frontside;
            end
            
            padded = [zeros(frontside,1);y;zeros(backside,1)];
            
            
        end
        
        function merged = mergestruct(dir,rows) % put the need to merge mat files into the same folder
            % rows specificy how many rows to merge
            files = Extract.filename(dir,'*.mat');
            summer = {};
            for idx= 1: length(files)
                eval(['load ',files{idx}]);
                if exist('rows','var')
                    summer{idx} = syllables(1:rows);
                else
                    summer{idx} = syllables;
                end
                merged = horzcat(summer{:});
                
            end
        end
        
        
        function coreThree(path_txt,path_pl2,path_folder)
            b = Chorus(path_txt,path_pl2,path_folder);
            b.select;
            neuronlist = b.getn;
            
            for k = 1: length(neuronlist)
                thisn = neuronlist{k};
                thisn.three;
                thisn.rawthree;
                %thisn.ResponseBasedOrderedThreePlots;
            end
        end
        
        function s = deblankl(x)
            if ~isempty(x)
                s = lower(x);
                s = s(find(s~=32));
            else
                s = [];
                
            end
        end
        
        function p = UpdateParforWaitbar(data, h)
%             D = parallel.pool.DataQueue;
%             h = waitbar(0, '开始生成 Neuron objects');
%             Utl.UpdateParforWaitbar(num_files, h);
%             afterEach(D, @Utl.UpdateParforWaitbar);
%             send(D, 1);
            
            persistent TOTAL COUNT H
            if nargin == 2
                % initialisation mode
                H = h;
                TOTAL = data;
                COUNT = 0;
            else
                % afterEach call, increment COUNT
                COUNT = 1 + COUNT;
                p = COUNT / TOTAL;
                waitbar(p, H,sprintf('此为总共%u个神经元中的%u',TOTAL,COUNT));
            end
        end
        
        function filename = fileparts(absolute_path)
            % to use the fileparts function in cellfun
            [~,filename,~] = fileparts(absolute_path);
        end
        
        function renamed = load(path)
            % load variable from path and rename
            loaded = load(path);
            renamed = loaded.(subsref(fieldnames(loaded),substruct('{}',{1})));
        end
        
        function bucketDriveletter = bucketletter(~)
            % get the drive letter of bucket server
            driveletters = cellfun( @(x) char(regexp(x,'[A-Za-z]','match')), getdrives );
            disknames = arrayfun(@DriveName,driveletters,'Uni',0);
            bucketDriveletter = driveletters(find( ~cellfun(@isempty, regexp(disknames,'bucket'))));
        end
         
    end
    
    
    
    
    
end

