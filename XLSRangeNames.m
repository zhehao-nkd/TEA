function [ RangeNames , RangeStatus ]= XLSRangeNames(filename)
%XLSRANGENAMES Returns a list of named ranges in an excel workbook
%   RangeNames = XLSRangeNames(filename)
%   RangeNames is a (n x 2) cell array with the name of the named range in
%   the first column and the value (aka location) in the second column
%
%   The value is an excel formula similar to: '=Sheet1!$A$1'
%
%   RangeStatus returns a 1 if at least one named range to return, 0
%   otherwise
%
%   Author: Stephen Lienhard
%
%   Example:
%     filename = 'C:\mytest.xls'
%     [ RangeNames , RangeStats ] = XLSRangeNames(filename);

%% Validate filename data type
if nargin < 1
    error('ExcelStuff:XLSRangeNames:Nargin',...
        'Filename must be specified.');
end
if ~isstr(filename)
    error('ExcelStuff:XLSRangeNames:InputClass','Filename must be a Char aray.');
end

% Validate filename is not empty
if isempty(filename)
    error('ExcelStuff:XLSRangeNames:FileName',...
        'Filename must not be empty.');
end

%% handle requested Excel workbook filename
[ pathstr , name , ext ] = fileparts( filename );
if isempty( ext )
	ext = '.xls';
	filename = [ filename , '.xls' ];
end
if isempty( pathstr )
	filename = which( filename , '-all' );
	if size( filename , 1 ) ~= 1
		error([ 'File was either not located, or multiple locations ' ...
            'were found. Please reissue readfromexcel command, ' ... 
            'providing absolute path to the file of interest.']);
	end
end

%% Attempt to start Excel as ActiveX server process
try
    Excel = actxserver( 'excel.application' );
    Excel.Visible = 0;
catch
    if ispc
        warning('ExcelStuff:XLSRangeNames:ActiveX', ...
            ['Could not start Excel server. ' ...
            'See documentation for resulting limitations.'])
    end
    return;
end

%% Attempt to open workbook
try
    workbook = Excel.workbooks.Open( filename );
catch
    workbook.Close( false ); % close workbook without saving any changes.
    delete( Excel ); % delete COM server
    RangeNames{ 1 , 1 } = '';
    RangeNames{ 1 , 2 } = '';
    RangeStatus         = 0;
    return;
end

%% Get named range names
NamedRangeCount = Excel.ActiveWorkBook.names.count;
if NamedRangeCount > 0
    RangeStatus = 1;
    for ii=1:NamedRangeCount
        RName=get( Excel.ActiveWorkbook.names.Item( ii ) );
        RangeNames{ ii , 1 } = RName.Name;
        RangeNames{ ii , 2 } = RName.Value;
    end
else
    RangeNames{ 1 , 1 } = '';
    RangeNames{ 1 , 2 } = '';
    RangeStatus         = 0;
end

%% Close everything down
try
    workbook.Close( false ); % close workbook without saving any changes.
    delete( Excel ); % delete COM server
end

return;