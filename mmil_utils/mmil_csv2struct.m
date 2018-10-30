function info = mmil_csv2struct(fname)
%function info = mmil_csv2struct(fname)
%
% Created:  02/02/09 by Don Hagler
% Last Mod: 05/03/16 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist(fname,'file'), error('file %s not found',fname); end;

% load comma separated value text file
raw_info = mmil_readtext(fname);

% exclude columns with empty column headers
ind_keep = find(~cellfun(@isempty,raw_info(1,:)));
raw_info = raw_info(:,ind_keep);

% recognize arrays and matricies from strings
raw_info = mmil_parsecells(raw_info);

% replace illegal characters in column headers with underscore
tags = regexprep(raw_info(1,:),'[\.\s();\?-]','_');

% replace multiple underscores with single underscore
tags = regexprep(tags,'_+','_');

% convert cell array to struct array
info = cell2struct(raw_info(2:end,:),tags,2);

