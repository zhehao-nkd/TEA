% figure sorter
classdef SFC %sylFigClassifier
    
    % Instruction: In order to classifier syllables based on teh writted figures
    
    methods(Static)
        
        
        function deleteFig(birdname)
            % When the segmentation of figures are wrongs, use this to delete all figures with this birdname
            files = Extract.filesAllLevel('./','*.png');
            ids = find(~cellfun(@isempty, regexp(files,birdname)));
            
            files_to_be_deleted = files(ids);
            for k = 1:length(files_to_be_deleted)
                delete(files_to_be_deleted{k});
            end
        end
        
    end
    
    
end