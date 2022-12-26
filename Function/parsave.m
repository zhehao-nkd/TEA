function parsave(filename,variablenames,variables,version)

if iscell(variablenames) % 如果input是包含变量名string的cell

    for k =1:length(variablenames)
        eval([variablenames{k},'=variables{k};']);
        try
            save('-mat',filename,variablenames{k},'-append'); % 每次循环保存一个变量
        catch
            if exist('version','var')
                save('-mat',filename,variablenames{k},version); % version = '-v7.3'
            else
                save('-mat',filename,variablenames{k}); % version = '-v7.3'
            end
        end
    end

elseif isstring(variablenames)||ischar(variables) % 当只有一个变量时，可以直接以string/char 变量名作为输入

    eval([variablenames,'=variables;']);

    if exist('version','var')
        save('-mat',filename,variablenames,version);
    else
        save('-mat',filename,variablenames);
    end

end

end