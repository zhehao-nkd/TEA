classdef Fenxi < handle
    %FENXI Summary of this class goes here
    %   Detailed explanation goes here
    
    % The fenxi should be further edited
    
    properties
        group
       
        source    % read xlxs into struct
        list
        b1
        n1
        n2
        n3
        sylnames % a cell of unique syl names
        degnames
        increnames
        categonames
        transnames
        mergedeleinf
        
        
     
    end
    
    methods
        
        function f = Fenxi(group)
            %FENXI Construct an instance of this class
            %   Detailed explanation goes here
            f.group = group;
            f.source = table2struct(readtable(pa.selected));
            
            targets = find([f.source.group].' == group);
            prims =  find([f.source.type].' == 1);
            seconds =  find([f.source.type].' == 2);
            primpath = f.source(intersect(targets,prims));
            secondpath = f.source(intersect(targets,seconds));
            
            b1 = Chorus(primpath.path_txt,primpath.path_pl2,primpath.path_folder);
            f.b1 = b1;
            n1idx = find(~cellfun(@isempty,regexp({b1.nlist.neuronname}.',primpath.neuron)));
            
            
            
            b2 = Chorus(secondpath.path_txt,secondpath.path_pl2,secondpath.path_folder);
             % b2is the batch object for 2nd-trail stimuli
            n2idx = find(~cellfun(@isempty,regexp({b2.nlist.neuronname}.',secondpath.neuron)));
            
            b1.select(n1idx); b2.select(n2idx);
            
           temp = b1.getn;
           n1 = temp{1};
           f.n1 = n1;
          
           temp = b2.getn;
           n2 = temp{1};
          % f.n2 = n2; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           
           list1 = n1.toList;
           list2 = n2.toList;
           f.list = horzcat(list1,list2);
           
           % sort
           f.sort  % !!!!!!!!!!!!
            
        end
        
        function f = sort(f)
            % syl
            sylidx = find(~cellfun(@isempty, regexp({f.list.stimuliname}.','syl')));
            syllist = f.list(sylidx);
            for k = 1: length(syllist)
                
                temp = strsplit(syllist(k).stimuliname,'-');
                
                sylnames{k} = temp{2};
            end
            if exist('sylnames','var')
                f.sylnames = unique(sylnames);
            end
            
            
            % deg
            degidx = find(~cellfun(@isempty, regexp({f.list.stimuliname}.','deg')));
            deglist = f.list(degidx);
            for k = 1: length(deglist)
                
                temp = strsplit(deglist(k).stimuliname,'-');
                
                degnames{k} = temp{2};
            end
            if exist('degnames','var')
                f.degnames = unique(degnames);
            end
           
            
            % incre
            increidx = find(~cellfun(@isempty, regexp({f.list.stimuliname}.','incre')));
            increlist = f.list(increidx);
            for k = 1: length(increlist)
                
                temp = strsplit(increlist(k).stimuliname,'-');
                
                increnames{k} = temp{2};
            end
            if exist('increnames','var')
                f.increnames = unique(increnames);
            end
          
            
            % catego
            categoidx = find(~cellfun(@isempty, regexp({f.list.stimuliname}.','catego')));
            categolist = f.list(categoidx);
            for k = 1: length(categolist)
                
                temp = strsplit(categolist(k).stimuliname,'-');
                
                categonames{k} = ['-',temp{2},'-'];
            end
            if exist('categonames','var')
                f.categonames = unique(categonames);
            end
            
            % trans/transform
            transidx = find(~cellfun(@isempty, regexp({f.list.stimuliname}.','trans')));
            translist = f.list(transidx);
            for k = 1: length(translist)
                
                temp = strsplit(translist(k).stimuliname,'-');
                
                transnames{k} = ['-',temp{2},'-'];
            end
            if exist('transnames','var')
                f.transnames = unique(transnames);
            end
            
            
            
        end
        
        function info = toVAE(f)% this is a function to generate the input for VAE analysis
            %info.sylnames = f.sylnames
        end
        function primthree(f)
            
            f.b1.select
            selected = f.b1.getn;
            
            for k = 1: length(selected)
                selected{k}.three;
            end
         
        end
        
        function drawselect(f,ids,range)
            d = Display(f.list);
            d.showselect(ids,range);
        end
        function drawsyl(f,keyword,rangeratio,ids)
            
            load("C:\Users\Zhehao\Dropbox (OIST)\My_Stimuli\Y661@06282021\mergedeleinf.mat");
            d = Display(f.list);
            
            if exist('keyword','var')
                
                if exist('rangeratio','var')
                    if exist('ids','var')
                        d.showsyl(keyword,mergedeleinf,rangeratio,ids);
                    else
                        d.showsyl(keyword,mergedeleinf,rangeratio);
                    end
                    
                else
                    d.showsyl(keyword,mergedeleinf);
                end  
            else
                for k = 1: length(f.sylnames)
                    d.showsyl(f.sylnames{k},mergedeleinf);
                end
            end
            
        end
        
        function drawdeg(f)
            %load("C:\Users\Zhehao\Dropbox (OIST)\My_Stimuli\Y661@06282021\mergedeleinf.mat");
            d = Display(f.list);
            for k = 1: length(f.degnames)
                d.showdeg(f.degnames{k});
            end
        end
        
        function drawcatego(f)
       
            d = Display(f.list);
            
            for k = 1: length(f.categonames)
                d.showcatego(f.categonames{k});
            end
            
        end
        
        function drawpartcatego(f,initial ,terminal) % show part of all
       
            d = Display(f.list);
            catego = d.findcatego;
            
            catego = catego(initial: terminal);
            
            d = Display(catego);
            d.showcatego;
            
        end  
        
        function drawallcatego(f)
            d = Display(f.list);
            d.showcatego;
            
        end
        
        function drawincre(f)
            
            d = Display(f.list);
            for k = 1: length(f.increnames)
                d.showincre(f.increnames{k});
            end
            
        end
        
        function drawtrans(f) % draw the stimuli-response figure of transformed stimuli
            d = Display(f.list);
            for k = 1: length(f.transnames)
                d.showtrans(f.transnames{k});
            end
        end 
        
        function drawsinglenorm(f,keyword)
            d = Display(f.list);
            d.shownorm(keyword);
            
        end
        
        function drawsimplesyl(f,keyword)
            load("C:\Users\Zhehao\Dropbox (OIST)\My_Stimuli\Y661@06282021\mergedeleinf.mat");
            d = Display(f.list);
            
            if exist('keyword','var')
                d.showsimplesyl(mergedeleinf,keyword);
            else
                d.showsimplesyl(mergedeleinf);
            end
        end
        
        function drawtransform(f,keyword)
        end
        
        function syllist = toAcousticSpace(f)
            
           ids = find(~cellfun(@isempty, regexp({f.list.stimuliname}.','syl'))); % find all syls
           
           syllist = f.list(ids);
           
           for n = 1: length(syllist)
            tempsum = Cal.psth_syl(syllist(n).rawy,syllist(n).fs,syllist(n).rawsptimes);
            halfsum = sum(tempsum(end/2:end));
            fullsum = sum(tempsum);
            maxvalue = max(Cal.psth_syl(syllist(n).rawy,syllist(n).fs,syllist(n).rawsptimes))
            syllist(n).maxvalue = maxvalue;
            syllist(n).halfsum = halfsum;
            syllist(n).fullsum = fullsum;
            if maxvalue > 8 % here the threshold is very important
                syllist(n).label = 1;
            else
                syllist(n).label = 0;
            end    
           end
           
        end
        
        function syllist = to2ndAcousticSpace(f)
            ids = find(~cellfun(@isempty, regexp({f.list.stimuliname}.','catego'))); % find all syls
           
           syllist = f.list(ids);
           
           for n = 1: length(syllist)
            tempsum = Cal.psth_syl(syllist(n).rawy,syllist(n).fs,syllist(n).rawsptimes);
            range = 1 % the very fisrt 1 second
            beginmax = max(tempsum(1: ceil(length(tempsum)*range/(length(syllist(n).rawy)/syllist(n).fs)) ));% the maximum value of the begining 0.5 second
            halfsum = sum(tempsum(end/2:end));
            fullsum = sum(tempsum);
            maxvalue = max(Cal.psth_syl(syllist(n).rawy,syllist(n).fs,syllist(n).rawsptimes));
            syllist(n).maxvalue = maxvalue;
            syllist(n).halfsum = halfsum;
            syllist(n).fullsum = fullsum;
             syllist(n).beginmax = beginmax;
            if maxvalue > 8 % here the threshold is very important
                syllist(n).label = 1;
            else
                syllist(n).label = 0;
            end    
           end
           
        end
        
        
        function back_analysis_in_spectral_space(f)
            
            % to read the segmentation info
            
            
            % to read the coordinates in VRAE
            
            % to read the significance info
            
            % to plot it in 2-d spectral space
            
            
        end
        
        

    end
end

