classdef Utl
    % store some utilis functions that I have not yet decided where to put
    
  
    methods(Static)

        function c = intersect_allowrepeats(a,b)
            % a = [1,2,3,3],b = [3,3,4,5], the intersection c =  [3,3]
            if isnumeric(a) && isnumeric(b)
                c = [];
                for i = 1:length(a)
                    for j = 1:length(b)
                        if a(i) == b(j)
                            c(end+1) = a(i);
                            b(j) = NaN; % Mark the element as already counted
                            break;
                        end
                    end
                end
            elseif iscellstr(a) && iscellstr(b)
                c = {};
                for i = 1:length(a)
                    for j = 1:length(b)
                        if strcmp(a{i}, b{j})
                            c{end+1} = a{i};
                            b{j} = ''; % Mark the element as already counted
                            break;
                        end
                    end
                end
            else
                error('Input arrays must be both numeric or both cell arrays of character strings');
            end

        end
        
        
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
        
   
        function s = deblankl(x)
            % 去除字符串里的空白字符
            if ~isempty(x)
                s = lower(x);
                s = s(find(s~=32));
            else
                s = [];
                
            end
        end
        

        function filename = fileparts(absolute_path)
            % to use the fileparts function in cellfun
            [~,filename,~] = fileparts(absolute_path);
        end
        
        
        function bucketDriveletter = bucketletter(~)
            % get the drive letter of bucket server
            driveletters = cellfun( @(x) char(regexp(x,'[A-Za-z]','match')), getdrives );
            disknames = arrayfun(@DriveName,driveletters,'Uni',0);
            bucketDriveletter = driveletters(find( ~cellfun(@isempty, regexp(disknames,'bucket'))));
        end
         
    end
    

    methods(Static) % 已经冻结的方法

        function merged = Frozen_mergestruct(dir,rows) % put the need to merge mat files into the same folder
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

        function Frozen_coreThree(path_txt,path_pl2,path_folder)
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

        function p = Frozen_UpdateParforWaitbar(data, h)
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
        
        function renamed = Frozen_load(path)
            % load variable from path and rename
            loaded = load(path);
            renamed = loaded.(subsref(fieldnames(loaded),substruct('{}',{1})));
        end
    
    end
    
    
    
    
end

