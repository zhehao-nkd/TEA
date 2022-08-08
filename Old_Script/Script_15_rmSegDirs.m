% remove segdata or syldata or eledata

subdirs = extract.folder("E:\Stimuli_Source\allBirdsSong")
for k = 1:length(subdirs)
    segdir = sprintf('%s\\SegData',subdirs{k});
    syldir = sprintf('%s\\SylData',subdirs{k});
    eledir = sprintf('%s\\EleData',subdirs{k});
    
%     if isfolder(segdir)
%     rmdir(  segdir );
%     end
%     
    if isfolder( syldir)
        
    cmd_rmdir(  syldir );
    
    end
    
    
    if isfolder( eledir)
    rmdir(  eledir );
    end
    
    disp('日暮途远，人间何世；将军一去，大树飘零')
 
    
end



function    [ st, msg ] = cmd_rmdir( folderspec )       
%   cmd_rmdir removes a directory and its contents 
%   
%   Removes all directories and files in the specified directory in
%   addition to the directory itself.  Used to remove a directory tree.
%   See also: xtests\generic_utilies_test.m
           
    narginchk( 1, 1 )
    
    dos_cmd = sprintf( 'rmdir /S /Q "%s"', folderspec );
    
    [ st, msg ] = system( dos_cmd );
end