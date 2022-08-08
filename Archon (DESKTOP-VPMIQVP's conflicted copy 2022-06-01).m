classdef Archon
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        bird_folders
        zpfolders
        all_data_info
    end
    
    methods
        function a = Archon(dirpath)
            
            dbstop if error
            a.bird_folders = [extract.folder(dirpath)].'; %#ok<NBRAK>
            
            pure_folname = {}; % 初始化
            for k = 1: length(a.bird_folders)
                temp1 = split(a.bird_folders{k},'\');
                pure_folname{k} = cellstr(temp1{end});
            end
            
            contain_ephys_ids = find(~cellfun(@isempty, regexp([pure_folname{:}].','Ephys')));
            if ~isempty(contain_ephys_ids)
                a.bird_folders =  a.bird_folders(contain_ephys_ids);
            else
                a.bird_folders =  dirpath;
            end
            
            infocollect = {};
            for k = 1:length(a.bird_folders)
                plxfiles = extract.filename(a.bird_folders{k},'*.plx');
                txtfiles = extract.filename(a.bird_folders{k},'*.txt');
                
                if exist(sprintf('%s\\Stimuli',a.bird_folders{k}), 'dir')
                    stidirs = extract.folder(sprintf('%s\\Stimuli',a.bird_folders{k})).';
                else % 如果未把所有stimuli放进stimuli folder里面
                    stidirs = extract.folder(sprintf('%s',a.bird_folders{k})).'
                end
                
                plxname = {}; % extract .plx filenames
                for m = 1: length(plxfiles)
                    [~,plxname{m},~] = fileparts(plxfiles{m});
                end
                
                txtname = {};% extract .txt filenames
                for m = 1: length(txtfiles)
                    [~,txtname{m},~] = fileparts(txtfiles{m});
                end
                
                stiname = {}; % extract stimuli folder names
                for m = 1: length( stidirs)
                    temps = split( stidirs{m},'\');
                    stiname{m} = temps{end};
                end
                
                % Based on .plx filenames, find the corresponding .txt and
                % stimuli folder names
                datainfo = struct;
                for m = 1: length(plxname)
                    
                    txt_hit = find(~cellfun(@isempty, regexp([cellstr(txtname)].',plxname{m})));
                    
                    sti_hit = find(~cellfun(@isempty, regexp([cellstr(stiname)].',plxname{m})));
                    
                    % fill the struct
                    sadtemp = regexp(plxname{m}, '[_-]', 'split');
                    datainfo(m).animal = sadtemp{end};
                    datainfo(m).plxname = plxname{m};
                    datainfo(m).keyword = regexp(plxname{m},'Con|Deg|Detail|Frag|Repla|Mir|Rev','match');
                    datainfo(m).sectionname = regexp(plxname{m},'Z\d{2}|P\d{2}','match');
                    fuck = regexp(cellstr(datainfo(m).sectionname),'\d{2}','match');
                    datainfo(m).sectionid = str2num(fuck{1}{1});
                    datainfo(m).subfilename = regexp(plxname{m},'F\d{1}','match');
                    if ~isempty(datainfo(m).subfilename) %如果文件名里存在F标识subfile的id
                        fuck2 = regexp(datainfo(m).subfilename,'\d','match');
                        datainfo(m).subfileid = str2num(fuck{1}{1});
                    elseif ~isempty(regexp(datainfo(m).plxname,'Mrg|mrg'))% 如果不存在F标识却存在‘Mrg|mrg’标识
                        datainfo(m).subfilename = "Mrg";
                    else % 如果不存在任何标识
                        datainfo(m).subfilename = [];
                    end
                    
                    datainfo(m).plxpath = plxfiles{m};
                    if ~isempty(txt_hit)
                        datainfo(m).txtpath = txtfiles{txt_hit};
                    else
                        datainfo(m).txtpath = [];
                    end
                    
                    if ~isempty(sti_hit)
                        datainfo(m).stidirpath = stidirs{sti_hit};
                    else
                        datainfo(m).stidirpath  = [];
                    end
                    
                end
                
                infocollect{k} = datainfo;
                
            end
            
            
            a.all_data_info = horzcat(infocollect{:});
            
            
            
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
    
    methods(Static)
        
        function A = Generate_Analysis_Object(birdname,ZorP,channelname,unit)
            arch = Archon('D:/');
            hitid = find(~cellfun(@isempty,regexp(cellstr(arch.bird_folders),birdname)));
            hitbird_folder = arch.bird_folders{hitid};
            
            
            % 找到hit folder 内部的 subfolders
            subdirs = extract.folder(hitbird_folder);
            hitZorP = find(~cellfun(@isempty,regexp(cellstr(subdirs.'),ZorP)));
            hit_subdir = subdirs {hitZorP};
            
            plxfiles = extract.filename(hit_subdir,'*.plx');
            zfids_plx = find(~cellfun(@isempty,regexp( cellfun(@(x)convertStringsToChars(x),{plxfiles{:}}.','Uni',0),'[ZP]\d+F\d+')));
            plxfiles = plxfiles(zfids_plx); %只考虑那些拥有ZP_F_ 格式的plx文件
                     
            txtfiles = extract.filename(hit_subdir,'*.txt');
            zfids = find(~cellfun(@isempty,regexp( cellfun(@(x)convertStringsToChars(x),{txtfiles{:}}.','Uni',0),'[ZP]\d+F\d+'))); %#ok<CCAT1>
            txtfiles = txtfiles(zfids); %只考虑那些拥有ZP_F_ 格式的txt文件
            
            stimuli_collect_dir = append(hit_subdir,'\Stimuli');
            stimuli_dirs = extract.folder(stimuli_collect_dir );
            
            Ns = {};
            
            for k = 1: length(plxfiles)
                
                fid = regexp(plxfiles{k},'(?<=F)\d*','match');
                
                % find corresponding txt file
                followF = regexp(cellstr(txtfiles),'(?<=F)\d*','match');
                duiying_txt_id = find(~cellfun(@isempty, regexp([followF{:}].', convertStringsToChars(fid))));
                
                followF = regexp(cellstr(stimuli_dirs),'(?<=F)\d*','match');
                followF = cellfun(@str2num, [followF{:}].','Uni',0);
                followF = [followF{:}].';
                duiying_stidir_id = find(followF == str2num(fid)) ;
                
                
                b = Batch(txtfiles{duiying_txt_id},plxfiles{k},stimuli_dirs{duiying_stidir_id});
                
                
                this_channel = channelname;
                this_unit = unit;
                neu_list = {b.nlist.neuronname}.';
                
                channel_ids = find(~cellfun(@isempty,regexp(neu_list,this_channel)));
                unit_ids = find(~cellfun(@isempty,regexp(neu_list,sprintf('_%u',this_unit))));
                neuron_ids = intersect(channel_ids,unit_ids);
                
                b.select(neuron_ids);
                N = b.getn{1};
                % allocate sap-based feature information to each neuron
                json_dir = append(hit_subdir,'\Json');
                
                json_files = extract.filename(json_dir,'*.json');
                
                datafile_id = find(~cellfun(@isempty, regexp(cellstr(json_files),'Data')));
                infofile_id = find(~cellfun(@isempty, regexp(cellstr(json_files),'Info')));
                % exist(json_dir,'dir')
                sorted_data = Neuron.extractFeaturesFromSapRawDataJson( json_files{datafile_id} , json_files{infofile_id} );
                N.setEachStimuliSapFeatures(sorted_data);
                N.calMeanFeatures;
            
                %duiying_song_folder
                Ns{k} = N;
                 
            end
            
            
           
            
            
            A = Analysis(Ns);
            A.birdname = birdname;
            A.zorp = ZorP;
            A.channelname = channelname;
            A.unitname = unit;
           
            A.image_naming_format = sprintf('%s_%s_%s_%u',A.birdname,A.zorp,A.channelname,A.unitname);
            save(sprintf('%s_%s_%s_%u',A.birdname,A.zorp,A.channelname,A.unitname),'A','-v7.3');
            
            
%             try
%             A.drawMeanFeaturesVsRespAsLineChart;
%             A.sort_frags_by_response_strength_and_then_draw;
%             A.drawAlignedNormDegsTwoPlots;
%             A.AlignReplasWithNormsThenDraw;
%             catch Error
%             end
            
        end
        
        function [neuroster,malfunctions] = BatchMode_Generate_Analysis_Object(datadir)
            % 分析每一个数据盘里的Z P
            % folder,生成Analysis对象，并且对于那些不可分析的（某些信息不够），生成记录文件以便补充这些信息
            %datadir = 'D:';
            dbstop if error
            
            birddir = extract.folder(datadir);
            ephysid = find(~cellfun(@isempty, regexp(cellstr(birddir),'Ephys')));
            birddir = birddir(ephysid);
            
            neuroster = struct; % Neuron 花名册
            malfunctions = struct; % malfunctions
            count = 0;
            malcount = 0;
            
            % waitbars
            wb1 = waitbar(0,'Bird folders');
            wb2 = waitbar(0,'Neuron of this bird');
            pos_w1=get(wb1,'position');
            pos_w2=[pos_w1(1) pos_w1(2)+pos_w1(4)+1 pos_w1(3) pos_w1(4)];
            set(wb2,'position',pos_w2,'doublebuffer','on');
                    
            for k = 1: length(birddir)
                
                waitbar(k/length(birddir),wb1,sprintf('当前Bird Folder是所有 %u 中的%u',length(birddir),k))
                
                %提取包含Z或者P的文件夹
                neurondir = extract.folder(birddir{k});
                neuronid = find(~cellfun(@isempty, regexp(cellstr(neurondir),'[ZP]\d+')));
                neurondir = neurondir(neuronid);
                
                for kk = 1: length(neurondir) 
                    
                    
                    waitbar(kk/length(birddir),wb2,sprintf('正处理 %u个神经元文件夹中的%u',length(neurondir) ,kk))
                    
                    
                     txtfile = extract.filename(neurondir{kk},'*.txt');
                     zfids = find(~cellfun(@isempty,regexp( cellfun(@(x)convertStringsToChars(x),{txtfile{:}}.','Uni',0),'[ZP]\d+F\d+'))); %#ok<CCAT1>
                     txtfile = txtfile(zfids); %只考虑那些拥有ZP_F_ 格式的txt文件
                     
                     plxfile = extract.filename(neurondir{kk},'*.plx');
                     zfids_plx = find(~cellfun(@isempty,regexp( cellfun(@(x)convertStringsToChars(x),{plxfile{:}}.','Uni',0),'[ZP]\d+F\d+')));
                     plxfile = plxfile(zfids_plx); %只考虑那些拥有ZP_F_ 格式的plx文件
                     
                     stidir = extract.folder(fullfile(neurondir{kk},'Stimuli'));
                     if isempty(stidir)
                         malcount = malcount + 1;
                         malfunctions(malcount).birdname = regexp(birddir{k},'[BRGYOX]\d{3}','match');
                         malfunctions(malcount).neuronid = regexp(neurondir{kk},'[ZP]\d+','match');
                         malfunctions(malcount).reason = 'Stimuli文件夹为空';
                         continue
                     end
                     zfids_sti = find(~cellfun(@isempty,regexp( cellfun(@(x)convertStringsToChars(x),{stidir{:}}.','Uni',0),'[ZP]\d+F\d+')));
                     stidir = stidir(zfids_sti);
                     
                     %判断.plx都已进行spike sorting并导出结果为.txt文件
                     if length(txtfile) < length(plxfile)
                         malcount = malcount + 1;
                         malfunctions(malcount).birdname = regexp(birddir{k},'[BRGYOX]\d{3}','match');
                         malfunctions(malcount).neuronid = regexp(neurondir{kk},'[ZP]\d+','match');
                         malfunctions(malcount).reason = 'txt文件缺失'; 
                         continue
                     end
                     
                     if isempty(txtfile) 
                         malcount = malcount + 1;
                         malfunctions(malcount).birdname = regexp(birddir{k},'[BRGYOX]\d{3}','match');
                         malfunctions(malcount).neuronid = regexp(neurondir{kk},'[ZP]\d+','match');
                         malfunctions(malcount).reason = 'txt文件不存在';
                         continue
                     end
                     
                     if length(stidir) < length(plxfile)
                         malcount = malcount + 1;
                         malfunctions(malcount).birdname = regexp(birddir{k},'[BRGYOX]\d{3}','match');
                         malfunctions(malcount).neuronid = regexp(neurondir{kk},'[ZP]\d+','match');
                         malfunctions(malcount).reason = 'stimuli folder缺失';
                         continue
                     end
                     
                     
                     
                    consistentunits = Archon.findConsistentNeurons(neurondir{kk}); 
                    % 判断是否存在每个文件都有的神经元
                    if isempty(consistentunits)
                        malcount = malcount + 1;
                        malfunctions(malcount).birdname = regexp(birddir{k},'[BRGYOX]\d{3}','match');
                        malfunctions(malcount).neuronid = regexp(neurondir{kk},'[ZP]\d+','match');
                        malfunctions(malcount).reason = '不存在每个文件都有的神经元';
                        continue
                    end
                    
                    
                    
                    
                    
                    for kkk = 1: length(consistentunits)
                    count = count + 1;
                    
                    neuroster(count).birdname = regexp(birddir{k},'[BRGYOX]\d{3}','match');
                    neuroster(count).neuronid = regexp(neurondir{kk},'[ZP]\d+','match');
                    neuroster(count).neuronpath = neurondir{kk};
                    neuroster(count).howManyFile = length(txtfile);
                    temp = regexp(consistentunits{kkk},'SPKC\d+','match');neuroster(count).channelname = temp{1};
                    temp = regexp(consistentunits{kkk},'(?<=_)\d+','match');neuroster(count).unit = str2num(temp{1});
                    
                    end
                    
                    
                end
                
                
                
            end
        end
        
        function channel_unit_info = findSortedNeurons(path_txt) % from one file
            
           % path_txt = "D:\Ephys-O709-W\P07\O709_P07F5.txt"
           channel_unit_info = struct;
           sorted_tables = Spike.split(path_txt);
           
           for k = 1: length(sorted_tables)
               
               channel_unit_info(k).channelname = sorted_tables{k}.('channelname'){1}
               channel_unit_info(k).unit =sorted_tables{k}.('unit')(1)
               channel_unit_info(k).channel_unit = sprintf('%s_%u',channel_unit_info(k).channelname,channel_unit_info(k).unit)
           end
            
        end
        
        function shared = findConsistentNeurons(dir_txt)  % 找到在每一个文件里都有的神经元
              
              %dir_txt = "D:\Ephys-O709-W\P07"
              txtfile = extract.filename(dir_txt,'*.txt');
              zfids = find(~cellfun(@isempty,regexp( cellfun(@(x)convertStringsToChars(x),{txtfile{:}}.','Uni',0),'[ZP]\d+F\d+'))); %#ok<CCAT1>
              txtfile = txtfile(zfids); %只考虑那些拥有ZP_F_ 格式的txt文件
                     
              for k = 1: length(txtfile)
                  info = Archon.findSortedNeurons(txtfile{k});
                  
                  if length(fieldnames(info))== 0
                     shared = [];
                     return
                  end
                      
                  consiss{k} = cellfun(@(x)convertCharsToStrings(x),{info.channel_unit}.','UniformOutput',1);
                  
              end
              
              
              shared = consiss{1}; %找到所有的txt文件都拥有的一个channel+unit
              for m = 1: length(consiss)
                  shared = intersect(shared,consiss{m});
              end
              
          end
              
    end
end

