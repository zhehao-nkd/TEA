

sourcelist =  listDataSources();

datasource = 'sap';
username = 'root';
password = 'sap2011';
conn = database(datasource,username,password);

tablelist = table2struct(sqlfind(conn,'','Catalog','sap'));
alltablenames = {tablelist.Table}.';

infonames = alltablenames(~cellfun(@isempty,regexp(alltablenames, 'file_table_[ogbyr]\d{3}[zp]\d+')));
datanames = alltablenames(~cellfun(@isempty,regexp(alltablenames, 'raw_[ogbyr]\d{3}[zp]\d+')));

% featuredir = "C:\Users\Zhehao\Music\Forced_FeatureDir";
% mkdir(featuredir)


wb1 = waitbar(0,'开始写入Info文件'); % loop the query to export txt
for k = 1:length(infonames) % info 文件


    bid = upper(regexp(convertCharsToStrings(infonames{k}),'[ogbyr]\d{3}','match'));
    zp = upper(regexp(convertCharsToStrings(infonames{k}),'[zp]\d+','match'));
    [~,zpdir] = Archon.getDirpath(bid,zp);
    %featuredir = fullfile(zpdir,'Feature');
    %mkdir(featuredir); % 建立Feature文件夹
    infofilename = sprintf('%s-%s-Info.txt',bid,zp);
    outpath = convert.path(fullfile(featuredir,infofilename),'unix');
    sqlquery = sprintf(['(SELECT ''bird_ID'',''file_index'',''file_name'',''file_age'',''bird_age'')\n',...
        'UNION\n',...
        '(SELECT bird_ID,file_index,file_name,file_age,bird_age\n',...
        'FROM sap.%s\n',...
        'INTO OUTFILE ''%s''\n',...
        'FIELDS ENCLOSED BY ''"'' TERMINATED BY '','' ESCAPED BY ''"''\n',...
        'LINES TERMINATED BY '',\\r'');']...
        ,infonames{k},outpath);

    if ~isfile(outpath)
        execute(conn,sqlquery);
    else
        disp('此文件已经生成')
    end
    waitbar(k/length(infonames),wb1,sprintf('已完成全部%u中的第%u个Info文件',length(infonames),k));
end

wb2 = waitbar(0,'开始写入Data文件');
for k = 1:length(datanames)% data 文件
    bid = upper(regexp(convertCharsToStrings(datanames{k}),'[ogbyr]\d{3}','match'));
    zp = upper(regexp(convertCharsToStrings(datanames{k}),'[zp]\d+','match'));
    [~,zpdir] = Archon.getDirpath(bid,zp);
    
    featuredir = fullfile(zpdir,'Feature');
    mkdir(featuredir); % 建立Feature文件夹
    
    infofilename = sprintf('%s-%s-Data.txt',bid,zp);
    outpath = convert.path(fullfile(featuredir,infofilename),'unix');
    sqlquery = sprintf(...
        ['(SELECT ''time'',''file_index'',''amplitude'',''mean_frequency_amp'',''pitch'',''mean_frequency'',',...
        '''FM'',''am'',''goodness'',''entropy'',''peak_frequency'',''DAS'',''continuity_t'',''continuity_f'')\n',...
        'UNION\n',...
        '(SELECT time,file_index,amplitude,mean_frequency_amp,pitch,mean_frequency,',...
        'FM,am,goodness,entropy,peak_frequency,DAS,continuity_t,continuity_f\n',...
        'FROM sap.%s',...
        ' INTO OUTFILE ''%s''\n',...
        'FIELDS ENCLOSED BY ''"'' TERMINATED BY '','' ESCAPED BY ''"''\n',...
        'LINES TERMINATED BY '',\\r'');']...
        ,datanames{k},outpath);
    
    if ~isfile(outpath)
        execute(conn,sqlquery);
    else
        disp('此文件已经生成')
    end
    waitbar(k/length(datanames),wb2,sprintf('已完成全部%u中的第%u个Data文件',length(datanames),k));
end
