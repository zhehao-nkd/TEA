% scratch to replace filenames


function rename(dirpath,ext,old,new)


    files = dir(sprintf('%s\%s',dirpath,ext));
    % Loop through each
    for id = 1:length(files)
        % Get the file name (minus the extension)
        [~, f,ext] = fileparts(files(id).name);
        % replace _ with -

        targets = regexp(f,old);

        f(targets) = new;

        newname = sprintf('%s%s',f,ext);



        % If numeric, rename
        movefile(fullfile(files(id).folder,files(id).name), fullfile(files(id).folder,newname));

    end


end