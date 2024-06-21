function info = mmil_csv2struct(fname,parseflag,separator,quote)
%function info = mmil_csv2struct(fname,[parseflag],[separator],[quote])
%
% Required Input:
%  fname: input file name
%
% Optional Input:
%   parseflag: [0|1] whether to parse strings for arrays and matrices
%     {default = 1}
%   separator: separator between columns
%     {default =  ','}
%   quote: opening/closing quote character
%     {default = '"'}
%
% Created:  02/02/09 by Don Hagler
% Prev Mod: 04/25/19 by Feng Xue
% Prev Mod: 05/22/19 by Don Hagler
% Last Mod: 05/12/20 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist(fname,'file'), error('file %s not found',fname); end;
if ~exist('parseflag','var') || isempty(parseflag), parseflag = 1; end;
if ~exist('separator','var') || isempty(separator), separator = ','; end;
if ~exist('quote','var'), quote = '"'; end;

% load comma separated value text file
raw_info = mmil_readtext(fname,separator,'',quote);

% exclude columns with empty column headers
ind_keep = find(~cellfun(@isempty,raw_info(1,:)));
raw_info = raw_info(:,ind_keep);

% recognize arrays and matrices from strings
if parseflag, raw_info = mmil_parsecells(raw_info); end;

% replace illegal characters in column headers with underscore
tags = regexprep(raw_info(1,:),'[\.\s();\?-]','_');

% replace multiple underscores with single underscore
tags = regexprep(tags,'_+','_');

% convert cell array to struct array
info = cell2struct(raw_info(2:end,:),tags,2);

