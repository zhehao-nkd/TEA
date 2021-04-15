classdef bucketfolder
    
    methods(Static)
        function folders = rmbad(iso, folders)
            
            [~,fnames,~] = cellfun(@fileparts,folders,'UniformOutput',false);
            
            
            
            for n = 1: length(fnames)
                
                alphabets = regexp(fnames{n},'[A-Za-z]');
                initial = fnames{n}(alphabets(1));
                
                number = fnames{n}(regexp(fnames{n},'\d'));
                
                new_fnames{n} = [initial,number];
                
            end
            
            % remove iso
            for m = 1:length(iso)
                
                idxes = find (~cellfun(@isempty, regexpi(new_fnames, iso(m))));
                
                if ~isempty(idxes)
                    for bad = 1:length (idxes) % this is a bad code
                        folders{idxes(bad)} = NaN;
                    end
                end
                
            end
            
            % remove G123B234 like mixed marks
            not1Idx = find(  cellfun(@length, regexp(new_fnames, '[A-Za-z]+\d{3}'))~=1 );
            
            if ~isempty(not1Idx)
                
                for worst = 1:length(not1Idx)
                    folders{not1Idx(worst)} = NaN;
                end
                
            end
            
            
            % remove wierd folders, e.g. channel, test B345346 etc.
            
            wierdIdx = find (cellfun(@isempty, regexp(new_fnames, '[A-Za-z]+\d{3}')));
            
            if ~isempty(wierdIdx)
                
                for worse = 1:length(wierdIdx)
                    folders{wierdIdx(worse)} = NaN;
                end
                
            end
            
            folders(cellfun(@(x) any(isnan(x)),folders)) = [];
            
        end
        
    end
end